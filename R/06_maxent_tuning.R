# ============================================================================
# 06_maxent_tuning.R
# Tunes MaxEnt hyperparameters (feature classes × regularisation multipliers)
# via manual grid search with spatial block CV, using year-matched covariate
# extraction. Fits the final model on all data and saves it for prediction.
#
# Year-matched extraction uses year-specific rasters for dynamic covariates
# (LST night, rainfall) so that each occurrence carries the environmental
# conditions from its observation year. This preserves 98 presences vs. 83
# under long-term mean extraction. The manual grid search bypasses ENMeval's
# hardcoded cell deduplication, which cannot accommodate year-matched values.
#
# Inputs:  data/processed/occurrences_thinned.csv
#          data/processed/background_points.csv
#          data/raw/ (covariate rasters — static + year-specific)
#          outputs/models/retained_vars.rds
#          outputs/models/spatial_cv_folds.rds
# Outputs: outputs/tables/enmeval_results.csv
#          outputs/models/selected_tuning.rds
#          outputs/models/maxent_final.rds
#          outputs/models/training_data.rds
#          outputs/figures/maxent_tuning_cbi.png
#          outputs/figures/maxent_tuning_cbi.pdf
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(maxnet)
  library(ecospat)
})

set.seed(SEED)

cat("Grid:", length(ENM_FC), "feature classes ×",
    length(ENM_RM), "regularisation multipliers =",
    length(ENM_FC) * length(ENM_RM), "combinations\n")
cat("Folds:", K_FOLDS, "\n")
cat("Total model fits:", length(ENM_FC) * length(ENM_RM) * K_FOLDS, "\n")

# ------------------------------ Load data -----------------------------------

occ <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
bg  <- read.csv(here::here("data", "processed", "background_points.csv"))

cat("Presences:", nrow(occ), "| Background:", nrow(bg), "\n")

retained_vars <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))

covs <- rast(file.path(DIR_COVARIATES, COV_FILES[retained_vars]))
names(covs) <- retained_vars

mask_r <- rast(file.path(DIR_COVARIATES, "ecological_mask_150mm.tif"))
covs_masked <- mask(covs, mask_r, maskvalues = 0)

cat("Covariates:", paste(retained_vars, collapse = ", "), "\n")

# ---------------------- Year-matched extraction -----------------------------
# Static covariates use a single raster; dynamic covariates use the
# year-specific annual raster for each point's observation year.

static_files <- c(
  slope      = "slope_1km.tif",
  river_dist = "river_distance_1km.tif",
  vertisols  = "vertisols_1km.tif"
)

dynamic_patterns <- c(
  lst_night = "lst_night_annual_{year}_1km.tif",
  rainfall  = "rainfall_{year}_1km.tif"
)

extract_year_matched <- function(pts, retained_vars) {
  static_vars  <- intersect(retained_vars, names(static_files))
  dynamic_vars <- intersect(retained_vars, names(dynamic_patterns))

  static_r <- rast(file.path(DIR_COVARIATES, static_files[static_vars]))
  names(static_r) <- static_vars
  static_r <- mask(static_r, mask_r, maskvalues = 0)

  static_vals <- terra::extract(static_r,
    as.matrix(pts[, c("longitude", "latitude")]))

  unique_years <- sort(unique(pts$year))
  dynamic_vals <- as.data.frame(matrix(NA, nrow = nrow(pts),
                                       ncol = length(dynamic_vars),
                                       dimnames = list(NULL, dynamic_vars)))

  for (yr in unique_years) {
    idx <- which(pts$year == yr)
    yr_files <- sapply(dynamic_vars, function(v) {
      file.path(DIR_COVARIATES, gsub("\\{year\\}", yr, dynamic_patterns[v]))
    })
    yr_r <- rast(yr_files)
    names(yr_r) <- dynamic_vars
    yr_r <- mask(yr_r, mask_r, maskvalues = 0)

    yr_vals <- terra::extract(yr_r,
      as.matrix(pts[idx, c("longitude", "latitude")]))
    dynamic_vals[idx, ] <- yr_vals
  }

  cbind(static_vals, dynamic_vals)
}

occ_env <- extract_year_matched(occ, retained_vars)
bg_env  <- extract_year_matched(bg, retained_vars)

cat("Occurrence extraction:", nrow(occ_env), "rows,",
    sum(complete.cases(occ_env)), "complete\n")
cat("Background extraction:", nrow(bg_env), "rows,",
    sum(complete.cases(bg_env)), "complete\n")

# ----------------------- Prepare partitions ---------------------------------

cv_folds <- readRDS(file.path(DIR_MODELS, "spatial_cv_folds.rds"))

n_occ <- nrow(occ)
occs_grp <- cv_folds$folds_ids[1:n_occ]
bg_grp   <- cv_folds$folds_ids[(n_occ + 1):length(cv_folds$folds_ids)]

# Remove NA rows
occ_na <- !complete.cases(occ_env)
bg_na  <- !complete.cases(bg_env)

occ_clean <- occ[!occ_na, ]
occ_env   <- occ_env[!occ_na, ]
occs_grp  <- occs_grp[!occ_na]

bg_clean <- bg[!bg_na, ]
bg_env   <- bg_env[!bg_na, ]
bg_grp   <- bg_grp[!bg_na]

cat("After NA removal — Presences:", nrow(occ_clean),
    "| Background:", nrow(bg_clean), "\n")

# Cell × year deduplication — only drop points sharing both cell and year
occ_cells <- cellFromXY(covs_masked, as.matrix(occ_clean[, c("longitude", "latitude")]))
cell_year <- paste(occ_cells, occ_clean$year, sep = "_")
dup_mask  <- duplicated(cell_year)

occ_clean <- occ_clean[!dup_mask, ]
occ_env   <- occ_env[!dup_mask, ]
occs_grp  <- occs_grp[!dup_mask]

cat("Cell×year deduplicated presences:", nrow(occ_clean),
    "(removed", sum(dup_mask), ")\n")

stopifnot(length(occs_grp) == nrow(occ_clean))
stopifnot(nrow(occ_env) == nrow(occ_clean))
stopifnot(length(bg_grp) == nrow(bg_clean))
stopifnot(nrow(bg_env) == nrow(bg_clean))

cat("Presence folds:", table(occs_grp), "\n")
cat("Background folds:", table(bg_grp), "\n")

# ----------------------------- Grid search ----------------------------------

pa <- c(rep(1, nrow(occ_clean)), rep(0, nrow(bg_clean)))
env_all <- rbind(occ_env, bg_env)
grp_all <- c(occs_grp, bg_grp)

results <- list()

for (fc in ENM_FC) {
  classes <- list(
    l = grepl("L", fc), q = grepl("Q", fc),
    h = grepl("H", fc), p = grepl("P", fc),
    t = grepl("T", fc)
  )

  for (rm in ENM_RM) {
    fold_metrics <- list()

    for (k in 1:K_FOLDS) {
      train_idx <- grp_all != k
      test_occ  <- which(grp_all == k & pa == 1)
      test_bg   <- which(grp_all == k & pa == 0)

      mod_k <- tryCatch(
        maxnet(p = pa[train_idx],
               data = env_all[train_idx, ],
               f = maxnet.formula(p = pa[train_idx],
                                  data = env_all[train_idx, ],
                                  classes = paste(names(classes)[unlist(classes)],
                                                  collapse = "")),
               regmult = rm),
        error = function(e) NULL
      )

      if (is.null(mod_k)) next

      pred_test_occ <- predict(mod_k, env_all[test_occ, ],
                               type = "cloglog") |> as.numeric()
      pred_test_bg  <- predict(mod_k, env_all[test_bg, ],
                               type = "cloglog") |> as.numeric()

      # AUC
      auc_k <- mean(sapply(pred_test_occ,
                           function(p) mean(p > pred_test_bg)))

      # CBI
      cbi_k <- tryCatch(
        ecospat::ecospat.boyce(
          fit = c(pred_test_occ, pred_test_bg),
          obs = pred_test_occ,
          nclass = 0, window.w = "default", res = 100,
          PEplot = FALSE
        )$cor,
        error = function(e) NA_real_
      )

      # Omission rate at 10th percentile
      train_occ_idx <- which(train_idx & pa == 1)
      pred_train_occ <- predict(mod_k, env_all[train_occ_idx, ],
                                type = "cloglog") |> as.numeric()
      thresh_10p <- quantile(pred_train_occ, probs = 0.1, na.rm = TRUE)
      or_10p <- mean(pred_test_occ < thresh_10p)

      fold_metrics[[k]] <- data.frame(auc = auc_k, cbi = cbi_k,
                                      or_10p = or_10p)
    }

    if (length(fold_metrics) == 0) next

    fm <- bind_rows(fold_metrics)
    results[[paste0(fc, "_", rm)]] <- data.frame(
      fc = fc, rm = rm,
      auc.val.avg = mean(fm$auc, na.rm = TRUE),
      auc.val.sd  = sd(fm$auc, na.rm = TRUE),
      cbi.val.avg = mean(fm$cbi, na.rm = TRUE),
      cbi.val.sd  = sd(fm$cbi, na.rm = TRUE),
      or.10p.avg  = mean(fm$or_10p, na.rm = TRUE),
      or.10p.sd   = sd(fm$or_10p, na.rm = TRUE)
    )

    cat(fc, "rm =", rm, "| CBI:",
        round(mean(fm$cbi, na.rm = TRUE), 3), "\n")
  }
}

res <- bind_rows(results)
cat("\nModels evaluated:", nrow(res), "\n")

# ----------------------------- Select model ---------------------------------
# LQH CBI is flat across regularisation (0.80–0.89 over 0.5–4.0 rm), so
# configurations are not meaningfully distinguishable on CV performance alone.
# Selection among equivalent-performing configurations uses response curve
# plausibility as a secondary criterion:
#   - rm=1.5 produces a sharp LST night thermal threshold and humped rainfall
#     response consistent with P. orientalis biology, both of which appear
#     independently in RF and GBT comparators.
#   - Higher rm (e.g. 2.0) penalises the hinge features that capture these
#     ecological shapes and compensates by shifting weight to river distance —
#     the one covariate all three algorithms disagree on and the accessibility
#     diagnostic flagged as a geographic proxy.
# Selecting rm=1.5 from within the range of equivalent CBI configurations
# prioritises ecological interpretability over a small, non-significant CBI
# difference.

best <- res |>
  filter(fc == "LQH", rm == 1.5)

cat("Selected model:\n")
cat("  Feature classes:", best$fc, "\n")
cat("  Regularization:", best$rm, "\n")
cat("  CBI (mean ± sd):", round(best$cbi.val.avg, 3), "±",
    round(best$cbi.val.sd, 3), "\n")
cat("  AUC (mean ± sd):", round(best$auc.val.avg, 3), "±",
    round(best$auc.val.sd, 3), "\n")
cat("  Omission (10p):", round(best$or.10p.avg, 3), "±",
    round(best$or.10p.sd, 3), "\n")

# ----------------------------- Tuning figure --------------------------------

p_tune <- res |>
  mutate(rm = as.numeric(rm)) |>
  ggplot(aes(x = rm, y = cbi.val.avg, colour = fc)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = cbi.val.avg - cbi.val.sd,
                    ymax = cbi.val.avg + cbi.val.sd),
                width = 0.15, alpha = 0.4) +
  geom_point(data = . %>% filter(fc == best$fc, rm == as.numeric(best$rm)),
             size = 4, shape = 1, stroke = 1.2, colour = "black") +
  scale_colour_brewer(palette = "Set1", name = "Feature\nclasses") +
  labs(x = "Regularization multiplier",
       y = "CBI (mean ± sd across folds)",
       title = "MaxEnt tuning: spatial block CV (year-matched extraction)",
       subtitle = paste0("Selected: ", best$fc, " rm=", best$rm,
                         " | CBI=", round(best$cbi.val.avg, 3))) +
  theme_minimal()

ggsave(file.path(DIR_FIGS, "maxent_tuning_cbi.pdf"), p_tune,
       width = 8, height = 5)
ggsave(file.path(DIR_FIGS, "maxent_tuning_cbi.png"), p_tune,
       width = 8, height = 5, dpi = 300)
cat("Saved tuning figure\n")

# ----------------------------- Fit final model ------------------------------
# Fit selected model on ALL presences + background (no hold-out).

best_classes <- tolower(best$fc)

final_mod <- maxnet(
  p    = pa,
  data = env_all,
  f    = maxnet.formula(p = pa, data = env_all, classes = best_classes),
  regmult = as.numeric(best$rm)
)

cat("Final model fitted on full dataset\n")
cat("  Presences:", sum(pa == 1), "| Background:", sum(pa == 0), "\n")
cat("  Coefficients:", length(final_mod$betas), "\n")

# --------------------------------- Save -------------------------------------

write.csv(res, file.path(DIR_TABLES, "enmeval_results.csv"),
          row.names = FALSE)

saveRDS(list(fc = best$fc, rm = as.numeric(best$rm)),
        file.path(DIR_MODELS, "selected_tuning.rds"))

saveRDS(final_mod, file.path(DIR_MODELS, "maxent_final.rds"))

saveRDS(list(occ_env = occ_env, bg_env = bg_env,
             occ_clean = occ_clean, bg_clean = bg_clean),
        file.path(DIR_MODELS, "training_data.rds"))

cat("\nSaved:\n")
cat("  Results:       ", file.path(DIR_TABLES, "enmeval_results.csv"), "\n")
cat("  Tuning params: ", file.path(DIR_MODELS, "selected_tuning.rds"), "\n")
cat("  Final model:   ", file.path(DIR_MODELS, "maxent_final.rds"), "\n")
cat("  Training data: ", file.path(DIR_MODELS, "training_data.rds"), "\n")

cat("06_maxent_tuning.R complete\n")