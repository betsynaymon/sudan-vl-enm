# 06_maxent_tuning.R
# Tunes MaxEnt (feature classes x regularisation multipliers) by grid search
# with spatial block CV and year-matched extraction, selects a configuration
# by the rule in params.R, and fits the final model on all data.
#
# Year-matched extraction gives each point the conditions of its observation
# year for dynamic covariates; prediction uses the long-term means. The manual
# grid search bypasses ENMeval's cell deduplication, which cannot accommodate
# year-matched values.
#
# Inputs:  OCC_FILE, BG_FILE, FOLD_TABLE_FILE, retained_vars.rds, covariates
# Outputs: TUNING_FILE, MODEL_FILE, TRAIN_FILE
#          outputs/tables/enmeval_results.csv
#          outputs/figures/maxent_tuning_cbi.png / .pdf
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))

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

occ <- read.csv(OCC_FILE)
bg  <- read.csv(BG_FILE)
retained_vars <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
cat("Presences:", nrow(occ), "| Background:", nrow(bg), "\n")
cat("Covariates:", paste(retained_vars, collapse = ", "), "\n")

# ---------------------- Year-matched extraction -----------------------------
# helpers.R: every point must be in the domain with complete values.

occ_env <- extract_year_matched(occ, retained_vars)
bg_env  <- extract_year_matched(bg,  retained_vars)

# Annual (fitting) and long-term (prediction) rasters must be on one scale
dyn   <- intersect(retained_vars, names(COV_ANNUAL))
ltm_r <- rast(file.path(DIR_COVARIATES, COV_FILES[dyn])); names(ltm_r) <- dyn
ltm   <- terra::extract(ltm_r, as.matrix(occ[, c("longitude", "latitude")]))
cat("\nAt presences, median [min, max]:\n")
for (v in dyn) cat(sprintf(
  "  %-10s annual %8.1f [%7.1f, %7.1f] | long-term %8.1f [%7.1f, %7.1f]\n", v,
  median(occ_env[[v]]), min(occ_env[[v]]), max(occ_env[[v]]),
  median(ltm[[v]]),     min(ltm[[v]]),     max(ltm[[v]])))

# ------------------------------- CV folds -----------------------------------
# Joined by point ID; each row carries its own fold from here on.

fold_table <- readRDS(FOLD_TABLE_FILE)
occ$fold <- attach_fold(occ$coordinate_id, 1, fold_table)
bg$fold  <- attach_fold(bg$bg_id,          0, fold_table)
cat("\nPresence folds:", table(occ$fold), "| Background folds:", table(bg$fold), "\n")

# Within-year thinning (02) means no two presences share a cell and year
occ_cy <- paste(cellFromXY(rast(DOMAIN_FILE),
                           as.matrix(occ[, c("longitude", "latitude")])), occ$year)
stopifnot("Presences share a cell and year" = !anyDuplicated(occ_cy))

# >= 150 mm region of the domain, for the within-belt underfitting check
sens_r  <- rast(SENS_MASK_FILE)
occ$wet <- terra::extract(sens_r, as.matrix(occ[, c("longitude", "latitude")]))[, 1] %in% 1
bg$wet  <- terra::extract(sens_r, as.matrix(bg[,  c("longitude", "latitude")]))[, 1] %in% 1
cat("In the >= 150 mm region: presences", sum(occ$wet), "| background", sum(bg$wet), "\n")
# ----------------------------- Grid search ----------------------------------

pa      <- c(rep(1, nrow(occ)), rep(0, nrow(bg)))
env_all <- rbind(occ_env, bg_env)
grp_all <- c(occ$fold, bg$fold)
wet_all <- c(occ$wet, bg$wet)

results  <- list()
fold_res <- list()

for (fc in ENM_FC) {
  for (rm in ENM_RM) {
    fold_metrics <- list()
    for (k in seq_len(K_FOLDS)) {
      train_idx <- grp_all != k
      mod_k <- tryCatch(fit_maxnet(pa[train_idx], env_all[train_idx, ], fc, rm),
                        error = function(e) NULL)
      if (is.null(mod_k)) next

      metrics_k <- eval_fold(mod_k,
        test_occ  = env_all[grp_all == k & pa == 1, ],
        test_bg   = env_all[grp_all == k & pa == 0, ],
        train_occ = env_all[train_idx & pa == 1, ])
        metrics_wet <- eval_fold(mod_k,
          test_occ  = env_all[grp_all == k & pa == 1 & wet_all, ],
          test_bg   = env_all[grp_all == k & pa == 0 & wet_all, ],
          train_occ = env_all[train_idx & pa == 1, ])
      names(metrics_wet) <- paste0(names(metrics_wet), "_wet")
      fold_metrics[[k]] <- cbind(fold = k, metrics_k, metrics_wet)
    }
    fm <- bind_rows(fold_metrics)
    if (nrow(fm) == 0) next

    fold_res[[paste0(fc, "_", rm)]] <- cbind(fc = fc, rm = rm, fm)

    results[[paste0(fc, "_", rm)]] <- data.frame(
      fc = fc, rm = rm, n_folds = nrow(fm), n_cbi = sum(!is.na(fm$cbi)),
      auc.val.avg = mean(fm$auc),              auc.val.sd = sd(fm$auc),
      cbi.val.avg = mean(fm$cbi, na.rm = TRUE), cbi.val.sd = sd(fm$cbi, na.rm = TRUE),
      or.10p.avg  = mean(fm$or_10p),           or.10p.sd  = sd(fm$or_10p),
      cbi_wet.avg = mean(fm$cbi_wet, na.rm = TRUE), cbi_wet.sd = sd(fm$cbi_wet, na.rm = TRUE))
    cat(fc, "rm =", rm, "| folds:", nrow(fm), "| CBI:",
        round(mean(fm$cbi, na.rm = TRUE), 3), "\n")
  }
}

res <- bind_rows(results)

fold_res <- bind_rows(fold_res)
write.csv(fold_res, file.path(DIR_TABLES, "tuning_by_fold.csv"), row.names = FALSE)
cat("\nConfigurations evaluated:", nrow(res), "of", length(ENM_FC) * length(ENM_RM), "\n")

# ----------------------------- Select model ---------------------------------
# Rule (params.R): highest mean validation CBI among configurations scored on
# every fold. Configurations within one standard error of the best (one-SE
# rule) are statistically indistinguishable: they are saved as the plateau and
# compared on response curves in 07, which can override via params.R.
# Ties go to the simpler configuration (fewer feature classes, then higher rm).

complete <- res |> filter(n_folds == K_FOLDS, n_cbi == K_FOLDS)
cat("Scored on all folds:", nrow(complete), "of", nrow(res), "\n")
stopifnot("No configuration was scored on every fold" = nrow(complete) > 0)

cat("\nTop configurations by CBI:\n")
complete |> arrange(desc(cbi.val.avg)) |> head(8) |>
  mutate(across(where(is.numeric), ~ round(., 3))) |> print()

top1    <- complete |> arrange(desc(cbi.val.avg)) |> slice(1)
se_top  <- top1$cbi.val.sd / sqrt(K_FOLDS)
plateau <- complete |> filter(cbi.val.avg >= top1$cbi.val.avg - se_top) |>
  arrange(desc(cbi.val.avg)) |> select(fc, rm, cbi.val.avg, cbi.val.sd)
cat("\nWithin one SE of the best (CBI >=", round(top1$cbi.val.avg - se_top, 3), "):\n")
print(plateau |> mutate(across(where(is.numeric), ~ round(., 3))))

if (is.null(SELECTED_FC) && is.null(SELECTED_RM)) {
  best <- complete |> arrange(desc(cbi.val.avg), nchar(fc), desc(rm)) |> slice(1)
  rule <- "highest CBI"
} else {
  stopifnot("Set both SELECTED_FC and SELECTED_RM, or neither" =
              !is.null(SELECTED_FC) && !is.null(SELECTED_RM))
  best <- complete |> filter(fc == SELECTED_FC, rm == SELECTED_RM)
  stopifnot("Override is not a configuration scored on every fold" = nrow(best) == 1)
  rule <- "override (params.R)"
}

# The grid must contain the optimum: at the top of ENM_RM, CBI for the
# selected feature classes must not still be rising by more than one SE.
if (best$rm == max(ENM_RM)) {
  lower <- complete |> filter(fc == best$fc, rm < best$rm)
  if (best$cbi.val.avg - max(lower$cbi.val.avg) > se_top)
    stop("CBI still rising at the top of ENM_RM: extend the grid")
}

# Underfitting check (params.R): does any configuration beat the selected one
# on within-belt CBI by more than one SE?
se_wet   <- best$cbi_wet.sd / sqrt(K_FOLDS)
best_wet <- complete |> arrange(desc(cbi_wet.avg), nchar(fc), desc(rm)) |> slice(1)
cat("\nSelected within-belt CBI:", round(best$cbi_wet.avg, 3),
    "| best within-belt:", best_wet$fc, best_wet$rm, round(best_wet$cbi_wet.avg, 3), "\n")
if (rule == "highest CBI" && best_wet$cbi_wet.avg - best$cbi_wet.avg > se_wet) {
  best <- best_wet
  rule <- "highest within-belt CBI (underfitting check)"
  cat("Selection moved to within-belt CBI\n")
}

cat("\nCBI by fold, plateau configurations:\n")
fold_res |> semi_join(plateau, by = c("fc", "rm")) |> select(fc, rm, fold, cbi) |>
  tidyr::pivot_wider(names_from = fold, values_from = cbi, names_prefix = "fold_") |>
  mutate(across(where(is.numeric), ~ round(., 3))) |> print()

cat("\nWithin-belt CBI (>= 150 mm), all configurations:\n")
complete |> select(fc, rm, cbi_wet.avg) |> mutate(cbi_wet.avg = round(cbi_wet.avg, 3)) |>
  tidyr::pivot_wider(names_from = rm, values_from = cbi_wet.avg) |> print()

cat("\nSelected configuration, all metrics by fold:\n")
fold_res |> filter(fc == best$fc, rm == best$rm) |>
  mutate(across(where(is.numeric), ~ round(., 3))) |> print()

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

final_mod <- fit_maxnet(pa, env_all, best$fc, best$rm)

cat("Final model fitted on full dataset\n")
cat("  Presences:", sum(pa == 1), "| Background:", sum(pa == 0), "\n")
cat("  Coefficients:", length(final_mod$betas), "\n")
cat("\nFinal model features (non-zero coefficients):\n")
print(signif(sort(final_mod$betas), 3))

# ------------------- Held-out presence predictions --------------------------
# Selected configuration: each presence's prediction from the model trained
# without its fold, and whether it falls below that model's OMISSION_Q
# threshold. Shows which records the model fails to transfer to.

held_out <- lapply(seq_len(K_FOLDS), function(k) {
  m_k <- fit_maxnet(pa[grp_all != k], env_all[grp_all != k, ], best$fc, best$rm)
  thr <- quantile(as.numeric(predict(m_k, env_all[grp_all != k & pa == 1, ],
                                     type = "cloglog")), OMISSION_Q)
  d <- occ[occ$fold == k, c("coordinate_id", "source", "year",
                            "longitude", "latitude", "fold")]
  d$pred <- as.numeric(predict(m_k, occ_env[occ$fold == k, ], type = "cloglog"))
  d$below_threshold <- d$pred < thr
  d
}) |> bind_rows()

write.csv(held_out, file.path(DIR_TABLES, "heldout_presence_predictions.csv"),
          row.names = FALSE)
cat("\nHeld-out presences below threshold (rows = fold):\n")
print(table(held_out$fold, held_out$below_threshold))

# --------------------------------- Save -------------------------------------

write.csv(res, file.path(DIR_TABLES, "enmeval_results.csv"),
          row.names = FALSE)

saveRDS(list(fc = best$fc, rm = best$rm, rule = rule,
             cbi = best$cbi.val.avg, cbi_sd = best$cbi.val.sd,
             auc = best$auc.val.avg, or_10p = best$or.10p.avg,
             plateau = plateau), TUNING_FILE)
saveRDS(final_mod, MODEL_FILE)
saveRDS(list(occ_env = occ_env, bg_env = bg_env,
             occ_clean = occ, bg_clean = bg), TRAIN_FILE)

cat("\nSaved:\n")
cat("  Results:       ", file.path(DIR_TABLES, "enmeval_results.csv"), "\n")
cat("  Tuning params: ", file.path(DIR_MODELS, "selected_tuning.rds"), "\n")
cat("  Final model:   ", file.path(DIR_MODELS, "maxent_final.rds"), "\n")
cat("  Training data: ", file.path(DIR_MODELS, "training_data.rds"), "\n")

cat("06_maxent_tuning.R complete\n")