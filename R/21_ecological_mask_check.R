# ============================================================================
# 21_mask_sensitivity.R
# "Maximally permissive" sensitivity check: full-country background (no
# ecological mask) with all occurrence points retained (including Khartoum
# referral records). Tests whether the ecological mask threshold and the
# point-filtering decisions drive the reported suitability surface and ARP.
#
# Design: same MaxEnt configuration (LQH, rm = 1.5), same spatial CV
# structure (100 km blocks, k = 4), same evaluation metrics. Only the
# background extent and presence set differ from the primary model.
#
# Inputs:  data/raw/compiled_vl_presences.csv
#          data/raw/ (covariate rasters — static + year-specific)
#          data/raw/gadm/ (cached GADM boundary via geodata)
#          outputs/models/retained_vars.rds
#          outputs/surfaces/maxent_suitability.tif  (primary, for comparison)
#          outputs/tables/arp_summary.csv            (primary, for comparison)
#          outputs/surfaces/worldpop_2025_aligned.tif
# Outputs: outputs/sensitivity/mask/suitability_no_mask.tif
#          outputs/sensitivity/mask/arp_summary_no_mask.csv
#          outputs/sensitivity/mask/arp_by_state_no_mask.csv
#          outputs/sensitivity/mask/comparison_summary.csv
#          outputs/sensitivity/mask/suitability_comparison.png
#          outputs/sensitivity/mask/occurrences_thinned_no_mask.csv
#          outputs/sensitivity/mask/maxent_no_mask.rds
#          outputs/sensitivity/mask/spatial_cv_folds_no_mask.rds
#          outputs/sensitivity/mask/training_data_no_mask.rds
#          outputs/sensitivity/mask/background_points_no_mask.csv
#          outputs/sensitivity/mask/cv_metrics_no_mask.csv
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(maxnet)
  library(ecospat)
  library(spThin)
  library(blockCV)
  library(geodata)
  library(patchwork)
})

set.seed(SEED)

# Output directory
DIR_SENS <- file.path(DIR_OUTPUTS, "sensitivity", "mask")
dir.create(DIR_SENS, showWarnings = FALSE, recursive = TRUE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SENSITIVITY CHECK: full-country background, all presences\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# =============================== SETUP =======================================

# Sudan boundary (GADM) — used as mask instead of ecological mask
adm0  <- gadm(country = "SDN", level = 0, path = here::here("data", "raw"))
adm1  <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))
sudan <- st_as_sf(adm0)

# Covariates — unmasked (no ecological mask applied)
retained_vars <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
covs <- rast(file.path(DIR_COVARIATES, COV_FILES[retained_vars]))
names(covs) <- retained_vars

# Create full-country mask from Sudan boundary:
# rasterize GADM to the covariate grid, 1 = inside Sudan, NA = outside
sudan_mask <- rasterize(adm0, covs[[1]], field = 1, background = NA)

# Mask covariates to Sudan boundary (not ecological mask)
covs_sudan <- mask(covs, sudan_mask)

cat("Retained covariates:", paste(retained_vars, collapse = ", "), "\n")
cat("Full-country mask cells:",
    format(global(sudan_mask, "notNA")[[1]], big.mark = ","), "\n\n")

# ============================ 1. THINNING ====================================
# Load ALL raw occurrences — no referral filtering

occ_raw <- read.csv(here::here("data", "raw", "compiled_vl_presences.csv"))
cat("Raw occurrence records:", nrow(occ_raw), "\n")

# Within-year spatial thinning at 5 km (identical to primary pipeline)
set.seed(SEED)

years <- unique(occ_raw$year)

thinned_by_year <- lapply(years, function(y) {
  yr_data <- occ_raw %>% filter(year == y)

  if (nrow(yr_data) < 2) return(yr_data)

  thin_input <- yr_data %>%
    transmute(LAT = latitude, LONG = longitude, SPEC = "VL")

  res <- thin(
    loc.data   = thin_input,
    lat.col    = "LAT",
    long.col   = "LONG",
    spec.col   = "SPEC",
    thin.par   = THIN_KM,
    reps       = 100,
    locs.thinned.list.return = TRUE,
    write.files = FALSE, write.log.file = FALSE, verbose = FALSE
  )

  best <- res[[which.max(sapply(res, nrow))]]

  best_keys <- paste(round(best$Longitude, 5), round(best$Latitude, 5))
  yr_data %>%
    filter(paste(round(longitude, 5), round(latitude, 5)) %in% best_keys)
})

occ_all <- bind_rows(thinned_by_year) %>%
  mutate(year = as.integer(year))

cat("After within-year thinning:", nrow(occ_all), "records\n")
cat("Unique locations:", occ_all %>% distinct(longitude, latitude) %>% nrow(), "\n")

# Compare with primary (which excluded referral IDs)
occ_primary <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
cat("Primary model presences:", nrow(occ_primary), "\n")
cat("Sensitivity presences:  ", nrow(occ_all), "\n")
cat("Net new points:         ", nrow(occ_all) - nrow(occ_primary), "\n\n")

# Save for reference
write.csv(occ_all, file.path(DIR_SENS, "occurrences_thinned_no_mask.csv"),
          row.names = FALSE)

# ========================= 2. BACKGROUND SAMPLING ============================
# 10k uniform-random points from all of Sudan (not restricted to eco mask)

set.seed(SEED)

bg_pts <- spatSample(sudan_mask, size = N_BACKGROUND, method = "random",
                     na.rm = TRUE, as.points = TRUE)

bg <- as.data.frame(bg_pts, geom = "XY") %>%
  rename(longitude = x, latitude = y) %>%
  select(longitude, latitude)

# Assign years from occurrence-year distribution
year_weights <- occ_all %>% count(year, name = "weight")
bg$year <- sample(year_weights$year, size = nrow(bg), replace = TRUE,
                  prob = year_weights$weight)

cat("Background points sampled:", nrow(bg), "\n\n")

# ======================== 3. SPATIAL CV FOLDS ================================
# Same structure as primary: 100 km blocks, k = 4, random assignment

pts <- bind_rows(
  occ_all %>% select(longitude, latitude) %>% mutate(pa = 1),
  bg      %>% select(longitude, latitude) %>% mutate(pa = 0)
)

pts_sf <- st_as_sf(pts, coords = c("longitude", "latitude"), crs = 4326)

folds <- cv_spatial(
  x         = pts_sf,
  column    = "pa",
  size      = BLOCK_SIZE_M,
  k         = K_FOLDS,
  hexagon   = FALSE,
  selection = "random",
  iteration = 200,
  seed      = 1238L
)

cat("Fold balance:\n")
print(folds$records)
cat("\n")

# ================= 4. YEAR-MATCHED COVARIATE EXTRACTION ======================
# Same logic as 06, but masking by Sudan boundary instead of ecological mask

static_files <- c(
  slope      = "slope_1km.tif",
  river_dist = "river_distance_1km.tif",
  vertisols  = "vertisols_1km.tif"
)

dynamic_patterns <- c(
  lst_night = "lst_night_annual_{year}_1km.tif",
  rainfall  = "rainfall_{year}_1km.tif"
)

extract_year_matched_full <- function(pts_df, retained_vars) {
  static_vars  <- intersect(retained_vars, names(static_files))
  dynamic_vars <- intersect(retained_vars, names(dynamic_patterns))

  # Static covariates — masked by Sudan boundary, not eco mask
  static_r <- rast(file.path(DIR_COVARIATES, static_files[static_vars]))
  names(static_r) <- static_vars
  static_r <- mask(static_r, sudan_mask)

  static_vals <- terra::extract(static_r,
    as.matrix(pts_df[, c("longitude", "latitude")]))

  # Dynamic covariates — year-matched, masked by Sudan boundary
  unique_years <- sort(unique(pts_df$year))
  dynamic_vals <- as.data.frame(matrix(NA, nrow = nrow(pts_df),
                                       ncol = length(dynamic_vars),
                                       dimnames = list(NULL, dynamic_vars)))

  for (yr in unique_years) {
    idx <- which(pts_df$year == yr)
    yr_files <- sapply(dynamic_vars, function(v) {
      file.path(DIR_COVARIATES, gsub("\\{year\\}", yr, dynamic_patterns[v]))
    })
    yr_r <- rast(yr_files)
    names(yr_r) <- dynamic_vars
    yr_r <- mask(yr_r, sudan_mask)

    yr_vals <- terra::extract(yr_r,
      as.matrix(pts_df[idx, c("longitude", "latitude")]))
    dynamic_vals[idx, ] <- yr_vals
  }

  cbind(static_vals, dynamic_vals)
}

occ_env <- extract_year_matched_full(occ_all, retained_vars)
bg_env  <- extract_year_matched_full(bg, retained_vars)

cat("Occurrence extraction:", nrow(occ_env), "rows,",
    sum(complete.cases(occ_env)), "complete\n")
cat("Background extraction:", nrow(bg_env), "rows,",
    sum(complete.cases(bg_env)), "complete\n")

# ======================== 5. PREPARE PARTITIONS ==============================

n_occ <- nrow(occ_all)
occs_grp <- folds$folds_ids[1:n_occ]
bg_grp   <- folds$folds_ids[(n_occ + 1):length(folds$folds_ids)]

# Remove NA rows
occ_na <- !complete.cases(occ_env)
bg_na  <- !complete.cases(bg_env)

if (any(occ_na)) cat("Dropping", sum(occ_na), "NA occurrence rows\n")
if (any(bg_na))  cat("Dropping", sum(bg_na), "NA background rows\n")

occ_clean <- occ_all[!occ_na, ]
occ_env   <- occ_env[!occ_na, ]
occs_grp  <- occs_grp[!occ_na]

bg_clean <- bg[!bg_na, ]
bg_env   <- bg_env[!bg_na, ]
bg_grp   <- bg_grp[!bg_na]

# Cell × year deduplication
occ_cells <- cellFromXY(covs_sudan, as.matrix(occ_clean[, c("longitude", "latitude")]))
cell_year <- paste(occ_cells, occ_clean$year, sep = "_")
dup_mask  <- duplicated(cell_year)

occ_clean <- occ_clean[!dup_mask, ]
occ_env   <- occ_env[!dup_mask, ]
occs_grp  <- occs_grp[!dup_mask]

cat("After cleanup — Presences:", nrow(occ_clean),
    "| Background:", nrow(bg_clean), "\n\n")

stopifnot(length(occs_grp) == nrow(occ_clean))
stopifnot(nrow(occ_env)    == nrow(occ_clean))
stopifnot(length(bg_grp)   == nrow(bg_clean))
stopifnot(nrow(bg_env)     == nrow(bg_clean))

# ========================== 6. SPATIAL CV ====================================
# Fixed at LQH rm = 1.5 — this is a sensitivity check, not a re-tuning

pa      <- c(rep(1, nrow(occ_clean)), rep(0, nrow(bg_clean)))
env_all <- rbind(occ_env, bg_env)
grp_all <- c(occs_grp, bg_grp)

cat("Running 4-fold spatial CV (LQH, rm = 1.5)...\n")

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
                              classes = "lqh"),
           regmult = 1.5),
    error = function(e) { cat("  Fold", k, "failed:", e$message, "\n"); NULL }
  )

  if (is.null(mod_k)) next

  pred_test_occ <- predict(mod_k, env_all[test_occ, ],
                           type = "cloglog") |> as.numeric()
  pred_test_bg  <- predict(mod_k, env_all[test_bg, ],
                           type = "cloglog") |> as.numeric()

  # AUC
  auc_k <- mean(sapply(pred_test_occ, function(p) mean(p > pred_test_bg)))

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

  fold_metrics[[k]] <- data.frame(fold = k, auc = auc_k, cbi = cbi_k)
  cat("  Fold", k, "| CBI:", round(cbi_k, 3), "| AUC:", round(auc_k, 3), "\n")
}

fm <- bind_rows(fold_metrics)
cat("\nCV summary:\n")
cat("  CBI:", round(mean(fm$cbi, na.rm = TRUE), 3), "\u00b1",
    round(sd(fm$cbi, na.rm = TRUE), 3), "\n")
cat("  AUC:", round(mean(fm$auc, na.rm = TRUE), 3), "\u00b1",
    round(sd(fm$auc, na.rm = TRUE), 3), "\n\n")

# ========================= 7. FIT FINAL MODEL ================================

final_mod <- maxnet(
  p    = pa,
  data = env_all,
  f    = maxnet.formula(p = pa, data = env_all, classes = "lqh"),
  regmult = 1.5
)

cat("Final model fitted\n")
cat("  Presences:", sum(pa == 1), "| Background:", sum(pa == 0), "\n")
cat("  Coefficients:", length(final_mod$betas), "\n\n")

# ======================= 8. PREDICTION SURFACE ===============================
# Predict to all of Sudan (same extent as primary model)

cat("Predicting suitability surface...\n")
covs_pred <- mask(covs, vect(sudan))
sens_suit <- terra::predict(covs_pred, final_mod, type = "cloglog", na.rm = TRUE)

cat("Prediction range:", round(minmax(sens_suit)[1], 4), "\u2013",
    round(minmax(sens_suit)[2], 4), "\n\n")

writeRaster(sens_suit, file.path(DIR_SENS, "suitability_no_mask.tif"),
            overwrite = TRUE)

# ======================== 9. ARP ESTIMATES ===================================

# Thresholds from THIS model's training predictions
pred_occ_final <- predict(final_mod, occ_env, type = "cloglog")[, 1]
pred_bg_final  <- predict(final_mod, bg_env, type = "cloglog")[, 1]

p10_sens    <- unname(quantile(pred_occ_final, 0.10))
candidates  <- sort(unique(c(pred_occ_final, pred_bg_final)))
sens_spec   <- sapply(candidates, function(t) mean(pred_occ_final >= t)) +
               sapply(candidates, function(t) mean(pred_bg_final < t))
maxsss_sens <- candidates[which.max(sens_spec)]

cat("Thresholds: p10 =", round(p10_sens, 4),
    "| maxSSS =", round(maxsss_sens, 4), "\n")

# Population overlay — reuse aligned WorldPop from primary pipeline
pop_aligned <- rast(file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"))

arp_p10_s      <- global(pop_aligned * (sens_suit >= p10_sens), "sum", na.rm = TRUE)[[1]]
arp_maxsss_s   <- global(pop_aligned * (sens_suit >= maxsss_sens), "sum", na.rm = TRUE)[[1]]
arp_weighted_s <- global(pop_aligned * sens_suit, "sum", na.rm = TRUE)[[1]]
total_pop      <- global(pop_aligned, "sum", na.rm = TRUE)[[1]]

cat("\n--- Sensitivity ARP ---\n")
cat("p10:          ", format(round(arp_p10_s), big.mark = ","), "\n")
cat("maxSSS:       ", format(round(arp_maxsss_s), big.mark = ","), "\n")
cat("Risk-weighted:", format(round(arp_weighted_s), big.mark = ","), "\n\n")

# State-level breakdown
risk_weighted_r <- pop_aligned * sens_suit

state_arp_sens <- data.frame(
  state        = adm1$NAME_1,
  total_pop    = terra::extract(pop_aligned, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]],
  arp_p10      = terra::extract(pop_aligned * (sens_suit >= p10_sens), adm1,
                                fun = "sum", na.rm = TRUE, ID = FALSE)[[1]],
  arp_maxsss   = terra::extract(pop_aligned * (sens_suit >= maxsss_sens), adm1,
                                fun = "sum", na.rm = TRUE, ID = FALSE)[[1]],
  arp_weighted = terra::extract(risk_weighted_r, adm1,
                                fun = "sum", na.rm = TRUE, ID = FALSE)[[1]]
) %>%
  mutate(
    pct_p10      = round(100 * arp_p10 / total_pop, 1),
    pct_maxsss   = round(100 * arp_maxsss / total_pop, 1),
    pct_weighted = round(100 * arp_weighted / total_pop, 1)
  ) %>%
  arrange(desc(arp_weighted))

# Save ARP tables
arp_summary_sens <- data.frame(
  metric     = c("p10", "maxSSS", "risk_weighted"),
  threshold  = c(round(p10_sens, 4), round(maxsss_sens, 4), NA),
  arp        = c(round(arp_p10_s), round(arp_maxsss_s), round(arp_weighted_s)),
  pct_of_pop = c(round(100 * arp_p10_s / total_pop, 1),
                 round(100 * arp_maxsss_s / total_pop, 1),
                 round(100 * arp_weighted_s / total_pop, 1)),
  total_pop  = round(total_pop)
)

write.csv(arp_summary_sens, file.path(DIR_SENS, "arp_summary_no_mask.csv"),
          row.names = FALSE)
write.csv(state_arp_sens, file.path(DIR_SENS, "arp_by_state_no_mask.csv"),
          row.names = FALSE)

# ==================== 10. COMPARISON WITH PRIMARY ============================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPARISON WITH PRIMARY MODEL\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Load primary surface and ARP
primary_suit <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))
primary_arp  <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))

# Surface correlation (50k random cells)
set.seed(SEED)
sample_cells <- spatSample(c(primary_suit, sens_suit), size = 50000,
                           method = "random", na.rm = TRUE, as.df = TRUE)
names(sample_cells) <- c("primary", "sensitivity")

pearson_r  <- cor(sample_cells$primary, sample_cells$sensitivity,
                  method = "pearson")
spearman_r <- cor(sample_cells$primary, sample_cells$sensitivity,
                  method = "spearman")

cat("Surface correlation:\n")
cat("  Pearson r:  ", round(pearson_r, 3), "\n")
cat("  Spearman \u03c1:", round(spearman_r, 3), "\n\n")

# ARP comparison
primary_rw <- primary_arp$arp[primary_arp$metric == "risk_weighted"]
sens_rw    <- round(arp_weighted_s)
pct_change <- round(100 * (sens_rw - primary_rw) / primary_rw, 1)

cat("Risk-weighted ARP:\n")
cat("  Primary:     ", format(primary_rw, big.mark = ","), "\n")
cat("  Sensitivity: ", format(sens_rw, big.mark = ","), "\n")
cat("  Change:      ", pct_change, "%\n\n")

# Summary table
comparison <- data.frame(
  metric = c("presences", "background_extent",
             "cbi_mean", "auc_mean",
             "pearson_r", "spearman_rho",
             "arp_rw_primary", "arp_rw_sensitivity", "arp_rw_pct_change",
             "p10_threshold_primary", "p10_threshold_sensitivity",
             "maxsss_threshold_primary", "maxsss_threshold_sensitivity"),
  primary = c(
    nrow(occ_primary), "ecological mask (150mm)",
    NA, NA,
    NA, NA,
    primary_rw,
    NA,
    NA,
    primary_arp$threshold[primary_arp$metric == "p10"],
    NA,
    primary_arp$threshold[primary_arp$metric == "maxSSS"],
    NA
  ),
  sensitivity = c(
    nrow(occ_clean), "full country (Sudan boundary)",
    round(mean(fm$cbi, na.rm = TRUE), 3),
    round(mean(fm$auc, na.rm = TRUE), 3),
    round(pearson_r, 3),
    round(spearman_r, 3),
    NA,
    sens_rw,
    pct_change,
    NA,
    round(p10_sens, 4),
    NA,
    round(maxsss_sens, 4)
  )
)

write.csv(comparison, file.path(DIR_SENS, "comparison_summary.csv"),
          row.names = FALSE)

cat("Comparison table saved\n\n")

# ====================== 11. COMPARISON FIGURE ================================
# Side-by-side suitability maps + difference

states_sf <- st_as_sf(adm1)

# Convert to data frames for plotting
primary_df <- as.data.frame(primary_suit, xy = TRUE)
names(primary_df) <- c("x", "y", "suitability")
primary_df <- primary_df[!is.na(primary_df$suitability), ]

sens_df <- as.data.frame(sens_suit, xy = TRUE)
names(sens_df) <- c("x", "y", "suitability")
sens_df <- sens_df[!is.na(sens_df$suitability), ]

# Difference surface
diff_r <- sens_suit - primary_suit
diff_df <- as.data.frame(diff_r, xy = TRUE)
names(diff_df) <- c("x", "y", "diff")
diff_df <- diff_df[!is.na(diff_df$diff), ]

suit_colours <- c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7",
                  "#EF8A62", "#B2182B")

map_xlim <- c(21.5, 39)
map_ylim <- c(8, 24.5)

# Panel A: Primary
p_a <- ggplot() +
  geom_sf(data = sudan, fill = "grey95", colour = NA) +
  geom_raster(data = primary_df, aes(x = x, y = y, fill = suitability)) +
  scale_fill_gradientn(colours = suit_colours, limits = c(0, 1),
                       name = "Suitability", na.value = "transparent") +
  geom_sf(data = states_sf, fill = NA, colour = "grey40", linewidth = 0.15) +
  geom_sf(data = sudan, fill = NA, colour = "black", linewidth = 0.3) +
  geom_point(data = occ_primary, aes(x = longitude, y = latitude),
             shape = 21, size = 1, stroke = 0.3,
             fill = "white", colour = "black") +
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
  labs(title = paste0("(a) Primary model (n = ", nrow(occ_primary), ")"),
       subtitle = "Ecological mask, referral points excluded") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.3, "cm"))

# Panel B: Sensitivity
p_b <- ggplot() +
  geom_sf(data = sudan, fill = "grey95", colour = NA) +
  geom_raster(data = sens_df, aes(x = x, y = y, fill = suitability)) +
  scale_fill_gradientn(colours = suit_colours, limits = c(0, 1),
                       name = "Suitability", na.value = "transparent") +
  geom_sf(data = states_sf, fill = NA, colour = "grey40", linewidth = 0.15) +
  geom_sf(data = sudan, fill = NA, colour = "black", linewidth = 0.3) +
  geom_point(data = occ_all, aes(x = longitude, y = latitude),
             shape = 21, size = 1, stroke = 0.3,
             fill = "white", colour = "black") +
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
  labs(title = paste0("(b) Full-country variant (n = ", nrow(occ_clean), ")"),
       subtitle = "No mask, all points retained") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.3, "cm"))

# Panel C: Difference (sensitivity - primary)
max_abs_diff <- max(abs(range(diff_df$diff, na.rm = TRUE)))

p_c <- ggplot() +
  geom_sf(data = sudan, fill = "grey95", colour = NA) +
  geom_raster(data = diff_df, aes(x = x, y = y, fill = diff)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0,
                       limits = c(-max_abs_diff, max_abs_diff),
                       name = "\u0394 Suitability",
                       na.value = "transparent") +
  geom_sf(data = states_sf, fill = NA, colour = "grey40", linewidth = 0.15) +
  geom_sf(data = sudan, fill = NA, colour = "black", linewidth = 0.3) +
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
  labs(title = "(c) Difference (sensitivity \u2212 primary)",
       subtitle = paste0("Pearson r = ", round(pearson_r, 3),
                         " | \u0394 ARP = ", pct_change, "%")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.3, "cm"))

# Combine
fig_comparison <- p_a + p_b + p_c +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Mask sensitivity: ecological mask vs. full-country background",
    subtitle = paste0("Risk-weighted ARP: primary ",
                      format(primary_rw, big.mark = ","),
                      " vs. sensitivity ",
                      format(sens_rw, big.mark = ","),
                      " (", pct_change, "% change)"),
    theme = theme(plot.title = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 10))
  )

ggsave(file.path(DIR_SENS, "suitability_comparison.png"), fig_comparison,
       width = 18, height = 8, dpi = 300)
cat("Saved suitability_comparison.png\n")

# ====================== 12. SAVE MODEL OBJECTS ===============================
 
saveRDS(final_mod, file.path(DIR_SENS, "maxent_no_mask.rds"))
saveRDS(folds, file.path(DIR_SENS, "spatial_cv_folds_no_mask.rds"))
 
saveRDS(list(occ_env = occ_env, bg_env = bg_env,
             occ_clean = occ_clean, bg_clean = bg_clean),
        file.path(DIR_SENS, "training_data_no_mask.rds"))
 
write.csv(bg, file.path(DIR_SENS, "background_points_no_mask.csv"),
          row.names = FALSE)
 
write.csv(fm, file.path(DIR_SENS, "cv_metrics_no_mask.csv"),
          row.names = FALSE)
 
cat("\nSaved model objects:\n")
cat("  Model:          ", file.path(DIR_SENS, "maxent_no_mask.rds"), "\n")
cat("  CV folds:       ", file.path(DIR_SENS, "spatial_cv_folds_no_mask.rds"), "\n")
cat("  Training data:  ", file.path(DIR_SENS, "training_data_no_mask.rds"), "\n")
cat("  Background pts: ", file.path(DIR_SENS, "background_points_no_mask.csv"), "\n")
cat("  CV metrics:     ", file.path(DIR_SENS, "cv_metrics_no_mask.csv"), "\n")

# ========================= FINAL SUMMARY =====================================

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SENSITIVITY CHECK COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("Presences:        primary", nrow(occ_primary),
    "| sensitivity", nrow(occ_clean), "\n")
cat("Background:       ecological mask → full country\n")
cat("Surface Pearson:  r =", round(pearson_r, 3), "\n")
cat("Surface Spearman: \u03c1 =", round(spearman_r, 3), "\n")
cat("ARP (risk-weighted): primary",
    format(primary_rw, big.mark = ","),
    "→ sensitivity",
    format(sens_rw, big.mark = ","),
    "(", pct_change, "%)\n")
cat("CBI: ", round(mean(fm$cbi, na.rm = TRUE), 3), "\n")
cat("AUC: ", round(mean(fm$auc, na.rm = TRUE), 3), "\n")
cat("\nOutputs in:", DIR_SENS, "\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")

cat("\n21_mask_sensitivity.R complete\n")