# ============================================================================
# 15_variable_sensitivity.R
# Tests how sensitive the suitability surface and ARP are to covariate
# selection by refitting MaxEnt with nine alternative covariate sets. Set A
# (annual; slope, river distance, vertisols, LST night, rainfall) is the
# primary model; its results are loaded for comparison.
#
# Variants, by type:
#   Seasonal composite           A_dry, A_wet   dry- / wet-season LST night
#   Alternative representative   B_annual       NDVI for rainfall
#     (one covariate swapped     E_elevation    elevation for LST night
#      for another from the      F_lstday       LST day for rainfall
#      same correlated cluster)  G_treecover    tree cover for rainfall
#   Adds correlated covariate    C_annual       LST day alongside rainfall
#     (violates the |r| >= 0.7   D_annual       tree cover alongside rainfall
#      screening)
#   Drops covariate              H_noriver      river distance removed
#
# Each variant is fitted under two selection rules:
#   Highest CBI       — feature classes and regularisation chosen by the best
#                       spatial-CV CBI across the full tuning grid (the
#                       approach reported in the dissertation).
#   Primary settings  — fixed at the primary model's settings, so any change
#                       reflects the covariate set alone. CV scores come
#                       from the same tuning grid.
#
# All variants share the same occurrence points, background points and
# spatial CV folds, with year-matched covariate extraction. Predictions use
# long-term mean surfaces. ARP is reported on the full surface and within the
# ecological mask. Per-variant results are cached; delete
# variant_summary_{variant}.rds to refit a variant.
#
# Inputs:  data/processed/occurrences_thinned.csv
#          data/processed/background_points.csv
#          outputs/models/spatial_cv_folds.rds
#          outputs/models/maxent_final.rds
#          outputs/models/selected_tuning.rds
#          outputs/tables/enmeval_results.csv
#          outputs/tables/arp_summary.csv
#          outputs/surfaces/worldpop_2025_aligned.tif
#          data/raw/ (all covariate rasters)
# Outputs: outputs/models/maxent_final_{variant}.rds        (highest CBI)
#          outputs/models/maxent_final_{variant}_fixed.rds  (primary settings)
#          outputs/surfaces/maxent_suitability_{variant}.tif
#          outputs/surfaces/maxent_suitability_{variant}_fixed.tif
#          outputs/tables/enmeval_results_{variant}.csv     (tuning grids)
#          outputs/tables/variant_summary_{variant}.rds     (cached results)
#          outputs/tables/variable_sensitivity_comparison.csv
#          outputs/figures/sensitivity_response_shared.png
#          outputs/figures/sensitivity_response_unique.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(maxnet)
  library(ecospat)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(rnaturalearth)
  library(geodata)
})

set.seed(SEED)

# ------------------------------ Load shared data ----------------------------

occ      <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
bg       <- read.csv(here::here("data", "processed", "background_points.csv"))
cv_folds <- readRDS(file.path(DIR_MODELS, "spatial_cv_folds.rds"))
mask_r   <- rast(file.path(DIR_COVARIATES, "ecological_mask_150mm.tif"))

adm0        <- gadm(country = "SDN", level = 0, path = here::here("data", "raw"))
sudan       <- st_as_sf(adm0)
pop_aligned <- rast(file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"))
total_pop   <- global(pop_aligned, "sum", na.rm = TRUE)[[1]]

cat("Presences:", nrow(occ), "| Background:", nrow(bg), "\n")

# ----------------------------- Define variants ------------------------------

variants <- list(

  B_annual = list(
    type = "Alternative representative (NDVI for rainfall)",
    vars = c("slope", "river_dist", "vertisols", "lst_night", "ndvi"),
    static = c(slope = "slope_1km.tif", river_dist = "river_distance_1km.tif",
               vertisols = "vertisols_1km.tif"),
    dynamic_annual = c(lst_night = "lst_night_annual_{year}_1km.tif",
                       ndvi = "ndvi_annual_{year}_1km.tif"),
    dynamic_mean = c(lst_night = "lst_night_annual_mean_2000_2024_1km.tif",
                     ndvi = "ndvi_annual_mean_2000_2024_1km.tif")
  ),

  C_annual = list(
    type = "Adds correlated covariate (LST day with rainfall, |r| = 0.72)",
    vars = c("slope", "river_dist", "vertisols", "lst_night", "lst_day", "rainfall"),
    static = c(slope = "slope_1km.tif", river_dist = "river_distance_1km.tif",
               vertisols = "vertisols_1km.tif"),
    dynamic_annual = c(lst_night = "lst_night_annual_{year}_1km.tif",
                       lst_day = "lst_day_annual_{year}_1km.tif",
                       rainfall = "rainfall_{year}_1km.tif"),
    dynamic_mean = c(lst_night = "lst_night_annual_mean_2000_2024_1km.tif",
                     lst_day = "lst_day_annual_mean_2000_2024_1km.tif",
                     rainfall = "rainfall_mean_2000_2024_1km.tif")
  ),

  D_annual = list(
    type = "Adds correlated covariate (tree cover with rainfall, |r| = 0.70)",
    vars = c("slope", "river_dist", "vertisols", "lst_night", "rainfall", "treecover"),
    static = c(slope = "slope_1km.tif", river_dist = "river_distance_1km.tif",
               vertisols = "vertisols_1km.tif"),
    dynamic_annual = c(lst_night = "lst_night_annual_{year}_1km.tif",
                       rainfall = "rainfall_{year}_1km.tif",
                       treecover = "treecover_{year}_1km.tif"),
    dynamic_mean = c(lst_night = "lst_night_annual_mean_2000_2024_1km.tif",
                     rainfall = "rainfall_mean_2000_2024_1km.tif",
                     treecover = "treecover_mean_2000_2024_1km.tif")
  ),

  A_dry = list(
    type = "Seasonal composite (dry-season LST night)",
    vars = c("slope", "river_dist", "vertisols", "lst_night", "rainfall"),
    static = c(slope = "slope_1km.tif", river_dist = "river_distance_1km.tif",
               vertisols = "vertisols_1km.tif"),
    dynamic_annual = c(lst_night = "lst_night_dry_{year}_1km.tif",
                       rainfall = "rainfall_{year}_1km.tif"),
    dynamic_mean = c(lst_night = "lst_night_dry_mean_2000_2024_1km.tif",
                     rainfall = "rainfall_mean_2000_2024_1km.tif")
  ),

  A_wet = list(
    type = "Seasonal composite (wet-season LST night)",
    vars = c("slope", "river_dist", "vertisols", "lst_night", "rainfall"),
    static = c(slope = "slope_1km.tif", river_dist = "river_distance_1km.tif",
               vertisols = "vertisols_1km.tif"),
    dynamic_annual = c(lst_night = "lst_night_wet_{year}_1km.tif",
                       rainfall = "rainfall_{year}_1km.tif"),
    dynamic_mean = c(lst_night = "lst_night_wet_mean_2000_2024_1km.tif",
                     rainfall = "rainfall_mean_2000_2024_1km.tif")
  ),

  # ---- Clean within-cluster swaps and a drop test ----

  E_elevation = list(
    type = "Alternative representative (elevation for LST night)",
    vars = c("slope", "river_dist", "vertisols", "elevation", "rainfall"),
    static = c(slope = "slope_1km.tif", river_dist = "river_distance_1km.tif",
               vertisols = "vertisols_1km.tif", elevation = "elevation_1km.tif"),
    dynamic_annual = c(rainfall = "rainfall_{year}_1km.tif"),
    dynamic_mean = c(rainfall = "rainfall_mean_2000_2024_1km.tif")
  ),

  F_lstday = list(
    type = "Alternative representative (LST day for rainfall)",
    vars = c("slope", "river_dist", "vertisols", "lst_night", "lst_day"),
    static = c(slope = "slope_1km.tif", river_dist = "river_distance_1km.tif",
               vertisols = "vertisols_1km.tif"),
    dynamic_annual = c(lst_night = "lst_night_annual_{year}_1km.tif",
                       lst_day = "lst_day_annual_{year}_1km.tif"),
    dynamic_mean = c(lst_night = "lst_night_annual_mean_2000_2024_1km.tif",
                     lst_day = "lst_day_annual_mean_2000_2024_1km.tif")
  ),

  G_treecover = list(
    type = "Alternative representative (tree cover for rainfall)",
    vars = c("slope", "river_dist", "vertisols", "lst_night", "treecover"),
    static = c(slope = "slope_1km.tif", river_dist = "river_distance_1km.tif",
               vertisols = "vertisols_1km.tif"),
    dynamic_annual = c(lst_night = "lst_night_annual_{year}_1km.tif",
                       treecover = "treecover_{year}_1km.tif"),
    dynamic_mean = c(lst_night = "lst_night_annual_mean_2000_2024_1km.tif",
                     treecover = "treecover_mean_2000_2024_1km.tif")
  ),

  H_noriver = list(
    type = "Drops covariate (river distance)",
    vars = c("slope", "vertisols", "lst_night", "rainfall"),
    static = c(slope = "slope_1km.tif", vertisols = "vertisols_1km.tif"),
    dynamic_annual = c(lst_night = "lst_night_annual_{year}_1km.tif",
                       rainfall = "rainfall_{year}_1km.tif"),
    dynamic_mean = c(lst_night = "lst_night_annual_mean_2000_2024_1km.tif",
                     rainfall = "rainfall_mean_2000_2024_1km.tif")
  )
)

# ----------------------------- Helper functions -----------------------------

extract_year_matched <- function(pts, variant, mask_r) {
  static_vars  <- names(variant$static)
  dynamic_vars <- names(variant$dynamic_annual)

  static_r <- rast(file.path(DIR_COVARIATES, variant$static))
  names(static_r) <- static_vars
  static_r <- mask(static_r, mask_r, maskvalues = 0)

  static_vals <- terra::extract(static_r,
    as.matrix(pts[, c("longitude", "latitude")]))[, static_vars, drop = FALSE]

  dynamic_vals <- as.data.frame(matrix(
    NA, nrow = nrow(pts), ncol = length(dynamic_vars),
    dimnames = list(NULL, dynamic_vars)
  ))

  for (yr in sort(unique(pts$year))) {
    idx <- which(pts$year == yr)
    yr_files <- sapply(dynamic_vars, function(v) {
      file.path(DIR_COVARIATES, gsub("\\{year\\}", yr, variant$dynamic_annual[v]))
    })
    yr_r <- rast(yr_files)
    names(yr_r) <- dynamic_vars
    yr_r <- mask(yr_r, mask_r, maskvalues = 0)

    yr_vals <- terra::extract(yr_r,
      as.matrix(pts[idx, c("longitude", "latitude")]))[, dynamic_vars, drop = FALSE]
    dynamic_vals[idx, ] <- yr_vals
  }

  cbind(static_vals, dynamic_vals)
}

prepare_partitions <- function(occ, bg, occ_env, bg_env, cv_folds, mask_r) {
  # Fold vector: first nrow(occ) entries = presences, rest = background
  occ_folds <- cv_folds$folds_ids[1:nrow(occ)]
  bg_folds  <- cv_folds$folds_ids[(nrow(occ) + 1):length(cv_folds$folds_ids)]

  # Remove NA rows (tracking which survive)
  occ_na <- !complete.cases(occ_env)
  bg_na  <- !complete.cases(bg_env)

  occ_clean <- occ[!occ_na, ]
  occ_env   <- occ_env[!occ_na, ]
  occ_folds <- occ_folds[!occ_na]
  bg_clean  <- bg[!bg_na, ]
  bg_env    <- bg_env[!bg_na, ]
  bg_folds  <- bg_folds[!bg_na]

  cat("  After NA removal — Presences:", nrow(occ_clean),
      "| Background:", nrow(bg_clean), "\n")

  # Cell × year deduplication
  occ_cells <- cellFromXY(mask_r, as.matrix(occ_clean[, c("longitude", "latitude")]))
  cell_year <- paste(occ_cells, occ_clean$year, sep = "_")
  dup_mask  <- duplicated(cell_year)

  occ_clean <- occ_clean[!dup_mask, ]
  occ_env   <- occ_env[!dup_mask, ]
  occ_folds <- occ_folds[!dup_mask]

  cat("  Cell×year dedup — Presences:", nrow(occ_clean),
      "(removed", sum(dup_mask), ")\n")
  cat("  Fold sizes (presences):", paste(table(occ_folds), collapse = ", "), "\n")

  list(
    occ_clean = occ_clean, occ_env = occ_env, occs_grp = occ_folds,
    bg_clean = bg_clean, bg_env = bg_env, bg_grp = bg_folds
  )
}

tune_maxent <- function(prep, variant) {
  pa <- c(rep(1, nrow(prep$occ_env)), rep(0, nrow(prep$bg_env)))
  env_all <- rbind(prep$occ_env, prep$bg_env)
  grp_all <- c(prep$occs_grp, prep$bg_grp)

  results <- list()

  for (fc in ENM_FC) {
    fc_map <- c(L = "l", Q = "q", H = "h", P = "p", T = "t")
    classes <- paste(fc_map[strsplit(fc, "")[[1]]], collapse = "")

    for (rm in ENM_RM) {
      fold_metrics <- list()

      for (k in 1:K_FOLDS) {
        train_idx <- grp_all != k
        test_occ <- which(grp_all == k & pa == 1)
        test_bg  <- which(grp_all == k & pa == 0)

        mod_k <- tryCatch(
          maxnet(p = pa[train_idx], data = env_all[train_idx, ],
                 f = maxnet.formula(p = pa[train_idx],
                                    data = env_all[train_idx, ],
                                    classes = classes),
                 regmult = rm),
          error = function(e) NULL
        )

        if (is.null(mod_k)) next

        pred_test_occ <- predict(mod_k, env_all[test_occ, ],
                                 type = "cloglog") |> as.numeric()
        pred_test_bg  <- predict(mod_k, env_all[test_bg, ],
                                 type = "cloglog") |> as.numeric()

        auc_k <- mean(sapply(pred_test_occ,
                             function(p) mean(p > pred_test_bg)))

        cbi_k <- tryCatch(
          ecospat::ecospat.boyce(
            fit = c(pred_test_occ, pred_test_bg),
            obs = pred_test_occ,
            nclass = 0, window.w = "default", res = 100,
            PEplot = FALSE
          )$cor,
          error = function(e) NA_real_
        )

        train_occ_idx <- which(train_idx & pa == 1)
        pred_train_occ <- predict(mod_k, env_all[train_occ_idx, ],
                                  type = "cloglog") |> as.numeric()
        or_10p <- mean(pred_test_occ <
                       quantile(pred_train_occ, probs = 0.1, na.rm = TRUE))

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
    }
  }

  bind_rows(results)
}

fit_and_predict <- function(prep, tuning_res, variant, variant_name, fixed = NULL) {
  # Select settings: highest CBI (default), or the primary model's settings
  if (is.null(fixed)) {
    best <- tuning_res[which.max(tuning_res$cbi.val.avg), ]
    rule <- "Highest CBI"
  } else {
    best <- tuning_res[tuning_res$fc == fixed$fc & tuning_res$rm == fixed$rm, ]
    stopifnot("Primary settings missing from tuning grid" = nrow(best) == 1)
    rule <- "Primary settings"
    variant_name <- paste0(variant_name, "_fixed")
  }
  cat("  [", rule, "] Selected:", best$fc, "rm =", best$rm,
      "| CBI:", round(best$cbi.val.avg, 3),
      "| AUC:", round(best$auc.val.avg, 3), "\n")

  # Fit final model on all data
  pa      <- c(rep(1, nrow(prep$occ_env)), rep(0, nrow(prep$bg_env)))
  env_all <- rbind(prep$occ_env, prep$bg_env)

  fc_map <- c(L = "l", Q = "q", H = "h", P = "p", T = "t")
  best_classes <- paste(fc_map[strsplit(best$fc, "")[[1]]], collapse = "")

  final_mod <- maxnet(
    p = pa, data = env_all,
    f = maxnet.formula(p = pa, data = env_all, classes = best_classes),
    regmult = as.numeric(best$rm)
  )

  cat("  Coefficients:", length(final_mod$betas), "\n")

  # Prediction surface (long-term means)
  static_r  <- rast(file.path(DIR_COVARIATES, variant$static))
  names(static_r) <- names(variant$static)
  dynamic_r <- rast(file.path(DIR_COVARIATES, variant$dynamic_mean))
  names(dynamic_r) <- names(variant$dynamic_mean)

  pred_stack <- mask(c(static_r, dynamic_r), vect(sudan))
  suit_r <- terra::predict(pred_stack, final_mod, type = "cloglog", na.rm = TRUE)

  cat("  Suitability range:", round(minmax(suit_r)[1], 4), "\u2013",
      round(minmax(suit_r)[2], 4), "\n")

  # Thresholds and ARP
  pred_occ <- predict(final_mod, prep$occ_env, type = "cloglog")[, 1]
  pred_bg  <- predict(final_mod, prep$bg_env, type = "cloglog")[, 1]

  p10 <- unname(quantile(pred_occ, 0.10))
  candidates <- sort(unique(c(pred_occ, pred_bg)))
  sens <- sapply(candidates, function(t) mean(pred_occ >= t))
  spec <- sapply(candidates, function(t) mean(pred_bg < t))
  maxsss <- candidates[which.max(sens + spec)]

  arp_p10      <- global(pop_aligned * (suit_r >= p10), "sum", na.rm = TRUE)[[1]]
  arp_maxsss   <- global(pop_aligned * (suit_r >= maxsss), "sum", na.rm = TRUE)[[1]]
  arp_weighted <- global(pop_aligned * suit_r, "sum", na.rm = TRUE)[[1]]
  arp_weighted_masked <- global(pop_aligned * mask(suit_r, mask_r, maskvalues = 0),
                                "sum", na.rm = TRUE)[[1]]

  cat("  ARP (risk-weighted):", format(round(arp_weighted), big.mark = ","), "\n")

  # Save
  saveRDS(final_mod, file.path(DIR_MODELS,
          paste0("maxent_final_", variant_name, ".rds")))
  writeRaster(suit_r, file.path(DIR_SURFACES,
          paste0("maxent_suitability_", variant_name, ".tif")),
          overwrite = TRUE)
  write.csv(tuning_res, file.path(DIR_TABLES,
          paste0("enmeval_results_", variant_name, ".csv")),
          row.names = FALSE)

  data.frame(
    variant = variant_name, type = variant$type, fc = best$fc, rm = best$rm,
    cbi = round(best$cbi.val.avg, 3), cbi_sd = round(best$cbi.val.sd, 3),
    auc = round(best$auc.val.avg, 3), auc_sd = round(best$auc.val.sd, 3),
    p10_thresh = round(p10, 4), maxsss_thresh = round(maxsss, 4),
    arp_p10 = round(arp_p10), arp_maxsss = round(arp_maxsss),
    arp_weighted = round(arp_weighted),
    arp_weighted_masked = round(arp_weighted_masked),
    n_coefs = length(final_mod$betas)
  )
}

# ================================ RUN VARIANTS ==============================
primary_tuning <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))
primary_fixed  <- list(fc = primary_tuning$fc, rm = primary_tuning$rm)
cat("Fixed-settings refits use:", primary_fixed$fc, "rm =", primary_fixed$rm, "\n")

results_list <- list()

for (vname in names(variants)) {
  summary_path <- file.path(DIR_TABLES, paste0("variant_summary_", vname, ".rds"))

  if (file.exists(summary_path)) {
    cat("\nLoading cached results for", vname, "\n")
    results_list[[vname]] <- readRDS(summary_path)
    next
  }

  cat("\n", strrep("=", 50), "\n")
  cat("Running variant:", vname, "\n")
  cat("Type:", variants[[vname]]$type, "\n")
  cat("Variables:", paste(variants[[vname]]$vars, collapse = ", "), "\n")
  cat(strrep("=", 50), "\n")

  set.seed(SEED)
  occ_env <- extract_year_matched(occ, variants[[vname]], mask_r)
  bg_env  <- extract_year_matched(bg, variants[[vname]], mask_r)

  prep <- prepare_partitions(occ, bg, occ_env, bg_env, cv_folds, mask_r)

  set.seed(SEED)
  tuning_res <- tune_maxent(prep, variants[[vname]])

  res_auto  <- fit_and_predict(prep, tuning_res, variants[[vname]], vname)
  res_fixed <- fit_and_predict(prep, tuning_res, variants[[vname]], vname,
                               fixed = primary_fixed)
  results_list[[vname]] <- bind_rows(res_auto, res_fixed)
  saveRDS(results_list[[vname]], summary_path)

  cat("  Done.\n")
}

# ========================= COMPARISON TABLE =================================

# Load Set A baseline
baseline_tuning <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))
baseline_mod    <- readRDS(file.path(DIR_MODELS, "maxent_final.rds"))
baseline_res    <- read.csv(file.path(DIR_TABLES, "enmeval_results.csv"))
baseline_arp    <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))

baseline_best <- baseline_res |>
  filter(fc == baseline_tuning$fc, rm == baseline_tuning$rm)

baseline_row <- data.frame(
  variant = "A_annual", fc = baseline_tuning$fc, rm = baseline_tuning$rm,
  type = "Primary",
  rule = "Ecological plausibility",
  cbi = round(baseline_best$cbi.val.avg, 3),
  cbi_sd = round(baseline_best$cbi.val.sd, 3),
  auc = round(baseline_best$auc.val.avg, 3),
  auc_sd = round(baseline_best$auc.val.sd, 3),
  p10_thresh = round(baseline_arp$threshold[baseline_arp$metric == "p10"], 4),
  maxsss_thresh = round(baseline_arp$threshold[baseline_arp$metric == "maxSSS"], 4),
  arp_p10 = baseline_arp$arp[baseline_arp$metric == "p10"],
  arp_maxsss = baseline_arp$arp[baseline_arp$metric == "maxSSS"],
  arp_weighted = baseline_arp$arp[baseline_arp$metric == "risk_weighted"],
  arp_weighted_masked = global(pop_aligned * mask(rast(file.path(DIR_SURFACES,
                          "maxent_suitability.tif")), mask_r, maskvalues = 0),
                          "sum", na.rm = TRUE)[[1]],
  n_coefs = length(baseline_mod$betas)
)

comparison <- bind_rows(baseline_row, bind_rows(results_list))

comparison <- comparison |>
  mutate(
    rule = case_when(
      variant == "A_annual"     ~ "Ecological plausibility",
      grepl("_fixed$", variant) ~ "Primary settings",
      TRUE                      ~ "Highest CBI"
    ),
    arp_weighted_masked = round(arp_weighted_masked)
  )

cat("\n--- Variable sensitivity comparison ---\n")
comparison |>
  mutate(arp_weighted_m = round(arp_weighted / 1e6, 1)) |>
  select(variant, type, rule, fc, rm, cbi, auc, arp_weighted_m, arp_weighted_masked) |>
  arrange(variant, rule) |>
  as_tibble() |>
  print(n = Inf, width = Inf)

write.csv(comparison, file.path(DIR_TABLES, "variable_sensitivity_comparison.csv"),
          row.names = FALSE)

# ========================= RESPONSE CURVES ==================================

all_variants <- c(
  list(A_annual = list(
    vars = c("slope", "river_dist", "vertisols", "lst_night", "rainfall"),
    static = c(slope = "slope_1km.tif", river_dist = "river_distance_1km.tif",
               vertisols = "vertisols_1km.tif"),
    dynamic_mean = c(lst_night = "lst_night_annual_mean_2000_2024_1km.tif",
                     rainfall = "rainfall_mean_2000_2024_1km.tif")
  )),
  variants
)

model_files <- c(
  A_annual = "maxent_final.rds",
  setNames(paste0("maxent_final_", names(variants), ".rds"), names(variants))
)

variant_order <- c("A_annual", names(variants))

variant_colours <- c(
  A_annual = "grey30", A_dry = "#B15928", A_wet = "#5E3C99",
  B_annual = "#1B9E77", C_annual = "#D95F02", D_annual = "#7570B3",
  E_elevation = "#E7298A", F_lstday = "#66A61E",
  G_treecover = "#E6AB02", H_noriver = "#A6761D"
)

n_pts <- 200
set.seed(SEED)

response_all <- list()

for (vname in names(all_variants)) {
  v   <- all_variants[[vname]]
  mod <- readRDS(file.path(DIR_MODELS, model_files[vname]))

  s_r <- rast(file.path(DIR_COVARIATES, v$static))
  names(s_r) <- names(v$static)
  d_r <- rast(file.path(DIR_COVARIATES, v$dynamic_mean))
  names(d_r) <- names(v$dynamic_mean)
  stack_r <- mask(c(s_r, d_r), mask_r, maskvalues = 0)

  bg_sample <- spatSample(stack_r, size = 5000, method = "random",
                          na.rm = TRUE, values = TRUE)
  bg_medians <- apply(bg_sample, 2, median)
  bg_ranges  <- apply(bg_sample, 2, range)

  for (var in v$vars) {
    if (var == "vertisols") {
      newdata <- as.data.frame(t(replicate(2, bg_medians)))
      newdata[[var]] <- c(0, 1)
    } else {
      newdata <- as.data.frame(t(replicate(n_pts, bg_medians)))
      newdata[[var]] <- seq(bg_ranges[1, var], bg_ranges[2, var],
                            length.out = n_pts)
    }

    newdata$suit <- as.numeric(
      predict(mod, newdata, clamp = TRUE, type = "cloglog"))

    response_all[[paste0(vname, "_", var)]] <- tibble(
      variant = vname, variable = var,
      value = newdata[[var]], suit = newdata$suit
    )
  }
}

response_df <- bind_rows(response_all)

var_labels <- c(
  slope = "Slope (degrees)", river_dist = "Distance to river (m)",
  vertisols = "Vertisols (0/1)", lst_night = "LST night (\u00b0C)",
  rainfall = "Rainfall (mm/yr)", ndvi = "NDVI",
  lst_day = "LST day (\u00b0C)", treecover = "Tree cover (%)",
  elevation = "Elevation (m)"
)

response_df <- response_df |>
  mutate(var_label = var_labels[variable],
         variant = factor(variant, levels = variant_order))

# Shared variables
shared_vars <- c("slope", "river_dist", "lst_night", "rainfall")

p_shared <- ggplot(response_df |> filter(variable %in% shared_vars),
       aes(x = value, y = suit, colour = variant)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ var_label, scales = "free_x", nrow = 2) +
  scale_colour_manual(values = variant_colours, name = "Variant") +
  labs(x = NULL, y = "Habitat suitability (cloglog)") +
       #title = "Marginal response curves \u2014 shared covariates across variants",
       #subtitle = "Each covariate varied across background range; others at median") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(DIR_FIGS, "sensitivity_response_shared.png"), p_shared,
       width = 10, height = 7, dpi = 300)
cat("Saved sensitivity_response_shared.png\n")

# Unique variables
unique_vars <- c("ndvi", "lst_day", "treecover", "elevation")

p_unique <- ggplot(response_df |> filter(variable %in% unique_vars),
       aes(x = value, y = suit, colour = variant)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ var_label, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = variant_colours, name = "Variant") +
  labs(x = NULL, y = "Habitat suitability (cloglog)",
       title = "Marginal response curves \u2014 variant-specific covariates") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(DIR_FIGS, "sensitivity_response_unique.png"), p_unique,
       width = 10, height = 4, dpi = 300)
cat("Saved sensitivity_response_unique.png\n")

cat("\n15_variable_sensitivity.R complete\n")