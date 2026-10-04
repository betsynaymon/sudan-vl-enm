# ============================================================================
# 15_variable_sensitivity.R
# Is the estimate sensitive to which covariates represent each environmental
# axis? Refits the model with seven alternative covariate sets:
#   Seasonal   lst_dry, lst_wet        dry- / wet-season LST night for annual
#   Swap       ndvi_for_rain           NDVI for rainfall (excluded on collinearity)
#   Addition   plus_lst_day, plus_treecover, plus_elevation
#              (excluded in 03 on ecological and measurement grounds, not
#              collinearity, so each is tested by adding it)
#   Drop       minus_river             river distance removed
# Each variant is fitted at the primary settings (the main comparison: any
# change reflects the covariate set alone) and at its own highest-CBI
# configuration in the 06 grid (descriptive). The selected configuration's
# rank in each variant's grid shows whether the primary-settings comparison
# uses a reasonable model (24). The primary set runs through the same code and
# must reproduce 06-08.
#
# Inputs:  TRAIN_FILE, TUNING_FILE, MODEL_FILE, SUIT_FILE, POP_ALIGNED_FILE,
#          DOMAIN_FILE, SENS_MASK_FILE, retained_vars.rds, correlation_matrix.rds
#          (03), enmeval_results.csv (06), arp_summary.csv and
#          arp_plateau_candidates.csv (08), covariate rasters (COV_FILES,
#          COV_ANNUAL, SEASONAL_FILES, SEASONAL_ANNUAL)
# Outputs: VARIANT_SURFACES_FILE (layers {variant}_fixed, {variant}_tuned)
#          outputs/tables/variant_tuning.csv
#          outputs/tables/variable_sensitivity_comparison.csv
#          outputs/figures/sensitivity_response_shared.png
#          outputs/figures/sensitivity_response_unique.png
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(maxnet); library(ggplot2)
})

# ------------------------------ Load inputs ---------------------------------

train   <- readRDS(TRAIN_FILE)
tuning  <- readRDS(TUNING_FILE)
mod     <- readRDS(MODEL_FILE)
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
cor_03  <- readRDS(file.path(DIR_MODELS, "correlation_matrix.rds"))
enm_06  <- read.csv(file.path(DIR_TABLES, "enmeval_results.csv"))
arp_08  <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))
cand_08 <- read.csv(file.path(DIR_TABLES, "arp_plateau_candidates.csv"))
suit_r  <- rast(SUIT_FILE)
pop     <- rast(POP_ALIGNED_FILE)
occ     <- train$occ_clean
bg      <- train$bg_clean
files_all  <- c(COV_FILES, SEASONAL_FILES)
annual_all <- c(COV_ANNUAL, SEASONAL_ANNUAL)
cat("Model:", tuning$fc, "rm =", tuning$rm, "| presences:", nrow(occ),
    "| background:", nrow(bg), "\n")

# ------------------------------- Variants -----------------------------------

swap <- function(from, to) replace(vars, vars == from, to)
variants <- list(
  primary        = vars,
  lst_dry        = swap("lst_night", "lst_night_dry"),
  lst_wet        = swap("lst_night", "lst_night_wet"),
  ndvi_for_rain  = swap("rainfall", "ndvi"),
  plus_lst_day   = c(vars, "lst_day"),
  plus_treecover = c(vars, "treecover"),
  plus_elevation = c(vars, "elevation"),
  minus_river    = setdiff(vars, "river_dist")
)
stopifnot("A variant covariate has no file" = all(unlist(variants) %in% names(files_all)),
          "A variant has no colour" = setequal(names(variants), names(pal_variants)))

# Labels from 03's correlation matrix: a swap shows r with the covariate it
# replaces; an addition shows its strongest r with the five.
r_txt <- function(a, b) if (all(c(a, b) %in% rownames(cor_03))) sprintf(", r = %.2f", cor_03[a, b]) else ""
describe <- function(v) {
  add <- setdiff(v, vars); drop <- setdiff(vars, v)
  if (length(add) && length(drop)) return(paste0(add, " for ", drop, r_txt(add, drop)))
  if (length(add)) {
    j <- vars[which.max(abs(cor_03[add, vars]))]
    return(paste0("+ ", add, " (strongest with ", j, r_txt(add, j), ")"))
  }
  if (length(drop)) return(paste0("- ", drop))
  "primary"
}
for (nm in names(variants)) cat(sprintf("  %-15s %s\n", nm, describe(variants[[nm]])))

# ------------------------- Rows and covariate values ------------------------
# Same presences, background and folds as 06 (TRAIN_FILE), with every variant
# covariate extracted year-matched. Rows missing a variant's value are dropped
# for that variant only, and counted.

pts     <- function(d) d[, c("longitude", "latitude", "year")]
all_v   <- unique(unlist(variants))
occ_all <- extract_year_matched(pts(occ), all_v, annual = annual_all, require_complete = FALSE)
bg_all  <- extract_year_matched(pts(bg),  all_v, annual = annual_all, require_complete = FALSE)
stopifnot("Primary covariate values differ from 06" =
  isTRUE(all.equal(occ_all[, vars], train$occ_env[, vars], check.attributes = FALSE)) &&
  isTRUE(all.equal(bg_all[, vars],  train$bg_env[, vars],  check.attributes = FALSE)))
cat("\nRows missing each covariate:\n")
print(rbind(presences = colSums(is.na(occ_all)), background = colSums(is.na(bg_all))))

# -------------------------------- Refits ------------------------------------
# Per variant: the 06 grid (cv_maxnet), then refit_maxnet at the primary
# settings and at the variant's highest-CBI configuration, predicted over the
# domain with long-term means.

grid <- expand.grid(fc = ENM_FC, rm = ENM_RM, stringsAsFactors = FALSE)

run_variant <- function(nm) {
  t0 <- Sys.time()
  v  <- variants[[nm]]
  ko <- complete.cases(occ_all[, v]); kb <- complete.cases(bg_all[, v])
  oe <- occ_all[ko, v]; be <- bg_all[kb, v]
  pa <- c(rep(1, nrow(oe)), rep(0, nrow(be)))
  env  <- rbind(oe, be)
  fold <- c(occ$fold[ko], bg$fold[kb])
  tune <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    cv <- tryCatch(cv_maxnet(pa, env, fold, grid$fc[i], grid$rm[i]), error = function(e) NULL)
    data.frame(variant = nm, grid[i, ],
               cbi    = if (is.null(cv)) NA_real_ else mean(cv$cbi, na.rm = TRUE),
               cbi_sd = if (is.null(cv)) NA_real_ else sd(cv$cbi, na.rm = TRUE),
               auc    = if (is.null(cv)) NA_real_ else mean(cv$auc))
  }))
  best <- tune[which.max(tune$cbi), ]
  sel  <- tune[tune$fc == tuning$fc & tune$rm == tuning$rm, ]
  covs <- domain_covs(v, files = files_all)
  fits <- lapply(list(fixed = sel, tuned = best), function(cfg) {
    r <- refit_maxnet(oe, be, occ$fold[ko], bg$fold[kb], cfg$fc, cfg$rm)
    r$surf <- terra::predict(covs, r$mod, type = "cloglog", na.rm = TRUE)
    r$arp  <- arp_estimates(r$surf, pop, r$thr)
    r$cfg  <- cfg
    r
  })
  cat(sprintf("  %-15s %d / %d rows | selected config rank %d of %d (CBI %.3f) | best %s rm %s (%.3f) | %.1f min\n",
              nm, nrow(oe), nrow(be), sum(tune$cbi > sel$cbi, na.rm = TRUE) + 1,
              sum(!is.na(tune$cbi)), sel$cbi, best$fc, best$rm, best$cbi,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(tune = tune, fits = fits, n_occ = nrow(oe), n_bg = nrow(be),
       within_se = sel$cbi >= best$cbi - best$cbi_sd / sqrt(K_FOLDS))
}

cat("\nVariants (presence / background rows):\n")
res <- setNames(lapply(names(variants), run_variant), names(variants))

# ------------------------ Primary reproduces 06-08 --------------------------

p   <- res$primary$fits$fixed
chk <- inner_join(res$primary$tune, enm_06, by = c("fc", "rm")) |>
  filter(!is.na(cbi), !is.na(cbi.val.avg))
arp_08v <- setNames(arp_08$arp, arp_08$metric)
stopifnot(
  "Primary grid differs from 06"             = isTRUE(all.equal(chk$cbi, chk$cbi.val.avg)),
  "Coefficients differ from MODEL_FILE (06)" = isTRUE(all.equal(p$mod$betas, mod$betas)),
  "CV CBI differs from 06"                   = abs(mean(p$cv$cbi) - tuning$cbi) < 1e-8,
  "Surface differs from SUIT_FILE (07)" =
    global(abs(p$surf - suit_r), "max", na.rm = TRUE)[[1]] < 1e-6,
  "Estimates differ from 08" =
    all(round(arp_estimates(suit_r, pop, p$thr)) == arp_08v[names(p$arp)])
)
cat("Primary reproduces 06 (grid, coefficients, CV CBI), 07 (surface) and 08 (estimates)\n")

# ------------------------------ Comparison ----------------------------------
# Reading fixed before the run: at the primary settings, a national
# risk-weighted change within the plateau spread from 08 means the estimate is
# robust to that covariate choice; beyond it, the choice matters. Tuned rows
# are descriptive. A selected configuration outside one SE of a variant's best
# means that variant's primary-settings comparison uses a weak model (24).

plateau <- max(abs(cand_08$change_vs_selected[startsWith(cand_08$candidate, paste0(tuning$fc, " "))]))
belt    <- values(rast(SENS_MASK_FILE), mat = FALSE) %in% 1
v_prim  <- values(p$surf, mat = FALSE)

comparison <- bind_rows(lapply(names(res), function(nm)
  bind_rows(lapply(names(res[[nm]]$fits), function(rule) {
    r <- res[[nm]]$fits[[rule]]; v <- values(r$surf, mat = FALSE); ok <- !is.na(v) & !is.na(v_prim)
    data.frame(variant = nm, rule = rule, change = describe(variants[[nm]]),
               fc = r$cfg$fc, rm = r$cfg$rm,
               n_presences = res[[nm]]$n_occ, n_background = res[[nm]]$n_bg,
               n_coef = length(r$mod$betas),
               cbi = mean(r$cv$cbi, na.rm = TRUE), auc = mean(r$cv$auc),
               sel_within_se = res[[nm]]$within_se,
               rho_domain = cor(v[ok], v_prim[ok], method = "spearman"),
               rho_belt   = cor(v[ok & belt], v_prim[ok & belt], method = "spearman"),
               arp_weighted = r$arp[["risk_weighted"]],
               arp_maxsss = r$arp[["maxsss"]], arp_p10 = r$arp[["p10"]])
  })))) |>
  mutate(rw_change_pct  = 100 * (arp_weighted / p$arp[["risk_weighted"]] - 1),
         beyond_plateau = abs(rw_change_pct) > plateau)

show <- function(rule_nm) {
  comparison |> filter(rule == rule_nm) |>
    select(variant, change, fc, rm, cbi, auc, sel_within_se, rho_domain, rho_belt,
           arp_weighted, rw_change_pct, beyond_plateau) |>
    mutate(across(c(cbi, auc, rho_domain, rho_belt), ~ round(., 3)),
           arp_weighted = fmt(arp_weighted), rw_change_pct = round(rw_change_pct, 1)) |>
    print(row.names = FALSE, right = FALSE)
}
cat("\nPrimary settings (", tuning$fc, " rm ", tuning$rm, ") | plateau spread: +/-",
    round(plateau, 1), "%\n", sep = "")
show("fixed")
cat("\nEach variant's highest-CBI configuration (descriptive):\n")
show("tuned")

# ---------------------------- Response curves -------------------------------
# Primary-settings models; others at each variant's presence median (as 07),
# over that variant's data range. Seasonal LST-night curves share the LST
# night panel.

curves <- bind_rows(lapply(names(res), function(nm) {
  v  <- variants[[nm]]
  ko <- complete.cases(occ_all[, v]); kb <- complete.cases(bg_all[, v])
  data.frame(variant = nm,
             response_curves(res[[nm]]$fits$fixed$mod, sapply(occ_all[ko, v], median),
                             rbind(occ_all[ko, v], bg_all[kb, v]), v))
})) |>
  mutate(panel   = sub("_(dry|wet)$", "", variable),
         label   = factor(cov_labels[panel], levels = unique(cov_labels)),
         variant = factor(variant, levels = names(pal_variants)))

curve_plot <- function(d, nrow) {
  ggplot(d, aes(value, suit, colour = variant)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ label, scales = "free_x", nrow = nrow) +
    scale_colour_manual(values = pal_variants, name = NULL) +
    labs(x = NULL, y = "Suitability (cloglog)") +
    theme(legend.position = "bottom")
}
save_fig(file.path(DIR_FIGS, "sensitivity_response_shared.png"),
         curve_plot(filter(curves, panel %in% vars, variable != "vertisols"), 2),
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT)
save_fig(file.path(DIR_FIGS, "sensitivity_response_unique.png"),
         curve_plot(filter(curves, !panel %in% vars), 1),
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT * 0.6)

# --------------------------------- Save -------------------------------------

surf_out <- do.call(c, unname(unlist(lapply(res, function(x) lapply(x$fits, `[[`, "surf")),
                                     recursive = FALSE)))
names(surf_out) <- paste0(rep(names(res), each = 2), "_", c("fixed", "tuned"))
writeRaster(surf_out, VARIANT_SURFACES_FILE, overwrite = TRUE)
write.csv(bind_rows(lapply(res, `[[`, "tune")),
          file.path(DIR_TABLES, "variant_tuning.csv"), row.names = FALSE)
write.csv(comparison, file.path(DIR_TABLES, "variable_sensitivity_comparison.csv"), row.names = FALSE)
cat("15_variable_sensitivity.R complete\n")