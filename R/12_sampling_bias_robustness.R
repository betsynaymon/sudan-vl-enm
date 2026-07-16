# ============================================================================
# 12_sampling_bias_robustness.R
# Resamples background proportional to accessibility (Weiss et al. 2018
# travel-time surface), refits MaxEnt, and compares to the uniform-background
# model. Tests both sqrt and log bias transforms for robustness. The 13M-to-
# 17.7M divergence quantifies the population rendered invisible by geographic
# bias in occurrence records.
#
# Inputs:  outputs/models/training_data.rds
#          outputs/models/spatial_cv_folds.rds
#          outputs/models/maxent_final.rds
#          outputs/models/retained_vars.rds
#          outputs/surfaces/maxent_suitability.tif
#          outputs/surfaces/worldpop_2025_aligned.tif
#          outputs/tables/arp_summary.csv
#          data/raw/weiss_travel_time.tif
#          data/raw/ecological_mask_150mm.tif
#          data/processed/occurrences_thinned.csv
# Outputs: outputs/models/maxent_biased_background.rds
#          outputs/models/maxent_log_bias_background.rds
#          outputs/surfaces/maxent_suitability_biased_bg.tif
#          outputs/tables/bias_correction_summary.csv
#          outputs/tables/arp_by_state_bias_comparison.csv
#          outputs/figures/suitability_biased_bg.png
#          outputs/figures/response_curves_bias_comparison.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(maxnet)
  library(ecospat)
  library(ggplot2)
  library(ggspatial)
  library(rnaturalearth)
  library(geodata)
})

# ------------------------------ Load inputs ---------------------------------

train   <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
folds   <- readRDS(file.path(DIR_MODELS, "spatial_cv_folds.rds"))
suit_mx <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
mod     <- readRDS(file.path(DIR_MODELS, "maxent_final.rds"))
tuning  <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))
tt_raw  <- rast(here::here("data", "raw", "weiss_travel_time.tif"))
mask_r  <- rast(file.path(DIR_COVARIATES, "ecological_mask_150mm.tif"))

best_classes <- tolower(tuning$fc)

cat("Presences:", nrow(train$occ_env), "\n")
cat("Original background:", nrow(train$bg_env), "\n")
cat("Tuning:", tuning$fc, "rm =", tuning$rm, "\n")

# -------------------- Bias surface and resample -----------------------------

tt_aligned <- resample(tt_raw, mask_r, method = "bilinear")

# Bias weight: inverse sqrt of travel time (accessible = high weight)
bias_r <- 1 / sqrt(1 + tt_aligned)
mask_r <- subst(mask_r, 0, NA)
bias_r <- mask(bias_r, mask_r)

cat("Bias weight range:", round(global(bias_r, "min", na.rm = TRUE)[[1]], 4),
    "\u2013", round(global(bias_r, "max", na.rm = TRUE)[[1]], 4), "\n")

set.seed(SEED)
bg_biased <- spatSample(bias_r, size = N_BACKGROUND, method = "weights",
                        na.rm = TRUE, as.points = TRUE)

bg_biased_df <- as.data.frame(bg_biased, geom = "XY") |>
  rename(longitude = x, latitude = y) |>
  select(longitude, latitude)

cat("Biased background sampled:", nrow(bg_biased_df), "\n")

# Assign years from occurrence-year distribution
year_weights <- train$occ_clean |> count(year, name = "weight")
bg_biased_df$year <- sample(year_weights$year, size = nrow(bg_biased_df),
                            replace = TRUE, prob = year_weights$weight)

# Accessibility comparison
tt_orig   <- terra::extract(tt_raw, vect(train$bg_clean, geom = c("longitude", "latitude"),
                                         crs = "EPSG:4326"))[, 2]
tt_biased <- terra::extract(tt_raw, vect(bg_biased_df, geom = c("longitude", "latitude"),
                                         crs = "EPSG:4326"))[, 2]

cat("\nBackground accessibility comparison:\n")
cat("Original  \u2014 median:", round(median(tt_orig, na.rm = TRUE)), "min\n")
cat("Biased    \u2014 median:", round(median(tt_biased, na.rm = TRUE)), "min\n")

# ------------------- Extract covariates and refit ---------------------------

cov_stack_mean <- rast(file.path(DIR_COVARIATES, COV_FILES[vars]))
names(cov_stack_mean) <- vars

bg_biased_pts <- vect(bg_biased_df, geom = c("longitude", "latitude"), crs = "EPSG:4326")
bg_biased_env <- terra::extract(cov_stack_mean, bg_biased_pts, ID = FALSE)

complete <- complete.cases(bg_biased_env)
bg_biased_env <- bg_biased_env[complete, ]
bg_biased_df  <- bg_biased_df[complete, ]
cat("Biased background after NA drop:", nrow(bg_biased_env), "\n")

p_mat <- as.matrix(train$occ_env[, vars])
b_mat <- as.matrix(bg_biased_env[, vars])

mod_biased <- maxnet(
  p    = c(rep(1, nrow(p_mat)), rep(0, nrow(b_mat))),
  data = as.data.frame(rbind(p_mat, b_mat)),
  f    = maxnet.formula(
    p    = c(rep(1, nrow(p_mat)), rep(0, nrow(b_mat))),
    data = as.data.frame(rbind(p_mat, b_mat)),
    classes = best_classes
  ),
  regmult = tuning$rm
)

cat("Biased-bg model:", sum(mod_biased$betas != 0), "non-zero /",
    length(mod_biased$betas), "total coefficients\n")

saveRDS(mod_biased, file.path(DIR_MODELS, "maxent_biased_background.rds"))

# ----------------------- Spatial CV evaluation -------------------------------
# Presences use the dedup-fixed fold vector; biased background gets fold
# assignments via nearest-block matching (new points not in original fold set).

occ_orig   <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
occ_keys   <- paste(round(occ_orig$longitude, 5), round(occ_orig$latitude, 5), occ_orig$year)
train_keys <- paste(round(train$occ_clean$longitude, 5), round(train$occ_clean$latitude, 5), train$occ_clean$year)
dedup_idx  <- which(!occ_keys %in% train_keys)

pres_folds <- folds$folds_ids[-(dedup_idx)][1:nrow(train$occ_env)]

# Assign biased background to nearest block
bg_sf <- st_as_sf(bg_biased_df, coords = c("longitude", "latitude"), crs = 4326)
nearest_idx <- st_nearest_feature(bg_sf, folds$blocks)
bg_folds <- folds$blocks$folds[nearest_idx]

df_biased <- bind_rows(
  as.data.frame(train$occ_env[, vars]) |> mutate(pa = 1),
  as.data.frame(bg_biased_env[, vars]) |> mutate(pa = 0)
)
df_biased$fold <- c(pres_folds, bg_folds)

cat("NA folds:", sum(is.na(df_biased$fold)), "\n")
cat("Presences per fold:", table(df_biased$fold[df_biased$pa == 1]), "\n")

fold_cbi <- fold_auc <- numeric(4)

for (k in 1:4) {
  idx_train <- df_biased$fold != k
  idx_test  <- df_biased$fold == k

  p_tr <- as.matrix(df_biased[idx_train & df_biased$pa == 1, vars])
  b_tr <- as.matrix(df_biased[idx_train & df_biased$pa == 0, vars])

  mod_k <- maxnet(
    p    = c(rep(1, nrow(p_tr)), rep(0, nrow(b_tr))),
    data = as.data.frame(rbind(p_tr, b_tr)),
    f    = maxnet.formula(
      p    = c(rep(1, nrow(p_tr)), rep(0, nrow(b_tr))),
      data = as.data.frame(rbind(p_tr, b_tr)),
      classes = best_classes
    ),
    regmult = tuning$rm
  )

  test_data <- df_biased[idx_test, vars]
  pred_test <- predict(mod_k, newdata = test_data, type = "cloglog")[, 1]

  pres_pred <- pred_test[df_biased$pa[idx_test] == 1]
  bg_pred   <- pred_test[df_biased$pa[idx_test] == 0]

  boyce <- ecospat.boyce(fit = pred_test, obs = pres_pred,
                         nclass = 0, PEplot = FALSE)
  fold_cbi[k] <- boyce$cor

  n1 <- length(pres_pred)
  n0 <- length(bg_pred)
  fold_auc[k] <- (sum(rank(c(pres_pred, bg_pred))[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cat("\nSpatial CV \u2014 biased background:\n")
cat("CBI:", round(mean(fold_cbi, na.rm = TRUE), 3),
    "\u00b1", round(sd(fold_cbi, na.rm = TRUE), 3), "\n")
cat("AUC:", round(mean(fold_auc), 3),
    "\u00b1", round(sd(fold_auc), 3), "\n")

# ------------------- Predict and compare surfaces ---------------------------

suit_biased_path <- file.path(DIR_SURFACES, "maxent_suitability_biased_bg.tif")

if (file.exists(suit_biased_path)) {
  cat("Loading cached biased-bg surface\n")
  suit_biased <- rast(suit_biased_path)
} else {
  cov_stack <- rast(file.path(DIR_COVARIATES, COV_FILES[vars]))
  names(cov_stack) <- vars

  suit_biased <- predict(cov_stack, mod_biased, type = "cloglog",
                         clamp = TRUE, na.rm = TRUE)

  writeRaster(suit_biased, suit_biased_path, overwrite = TRUE)
  cat("Computed and saved biased-bg surface\n")
}

# Mask both to ecological mask for comparison
suit_mx_masked     <- mask(suit_mx, mask_r)
suit_biased_masked <- mask(suit_biased, mask_r)

# Surface correlation
set.seed(SEED)
valid_cells <- which(!is.na(values(suit_mx_masked)) &
                     !is.na(values(suit_biased_masked)))
samp_idx <- sample(valid_cells, min(50000, length(valid_cells)))

r_pearson  <- round(cor(values(suit_mx_masked)[samp_idx],
                        values(suit_biased_masked)[samp_idx]), 3)
r_spearman <- round(cor(values(suit_mx_masked)[samp_idx],
                        values(suit_biased_masked)[samp_idx],
                        method = "spearman"), 3)

cat("\nSurface correlation (50k sample, within ecological mask):\n")
cat("Pearson: ", r_pearson, "\n")
cat("Spearman:", r_spearman, "\n")

# ----------------------------- ARP comparison -------------------------------

pop_aligned <- rast(file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"))
total_pop   <- global(pop_aligned, "sum", na.rm = TRUE)[[1]]

# Biased-background thresholds
pred_occ_biased <- predict(mod_biased, newdata = train$occ_env[, vars],
                           type = "cloglog")[, 1]
pred_bg_biased  <- predict(mod_biased, newdata = bg_biased_env[, vars],
                           type = "cloglog")[, 1]

biased_p10 <- unname(quantile(pred_occ_biased, 0.10))

candidates <- sort(unique(c(pred_occ_biased, pred_bg_biased)))
sens <- sapply(candidates, function(t) mean(pred_occ_biased >= t))
spec <- sapply(candidates, function(t) mean(pred_bg_biased < t))
biased_maxsss <- candidates[which.max(sens + spec)]

biased_arp_weighted <- global(pop_aligned * suit_biased_masked, "sum", na.rm = TRUE)[[1]]
biased_arp_maxsss   <- global(pop_aligned * (suit_biased_masked >= biased_maxsss), "sum", na.rm = TRUE)[[1]]
biased_arp_p10      <- global(pop_aligned * (suit_biased_masked >= biased_p10), "sum", na.rm = TRUE)[[1]]

# Original ARP (recomputed on masked surface for consistency)
orig_pred_occ <- predict(mod, newdata = train$occ_env[, vars], type = "cloglog")[, 1]
orig_pred_bg  <- predict(mod, newdata = train$bg_env[, vars], type = "cloglog")[, 1]
orig_p10 <- unname(quantile(orig_pred_occ, 0.10))

candidates_orig <- sort(unique(c(orig_pred_occ, orig_pred_bg)))
sens_orig <- sapply(candidates_orig, function(t) mean(orig_pred_occ >= t))
spec_orig <- sapply(candidates_orig, function(t) mean(orig_pred_bg < t))
orig_maxsss <- candidates_orig[which.max(sens_orig + spec_orig)]

orig_arp_weighted <- global(pop_aligned * suit_mx_masked, "sum", na.rm = TRUE)[[1]]
orig_arp_maxsss   <- global(pop_aligned * (suit_mx_masked >= orig_maxsss), "sum", na.rm = TRUE)[[1]]
orig_arp_p10      <- global(pop_aligned * (suit_mx_masked >= orig_p10), "sum", na.rm = TRUE)[[1]]

cat("\n--- ARP comparison (ecological mask applied) ---\n")
cat(sprintf("%-20s %14s %14s\n", "", "Original", "Biased-bg"))
cat(sprintf("%-20s %14s %14s\n", "Risk-weighted",
    format(round(orig_arp_weighted), big.mark = ","),
    format(round(biased_arp_weighted), big.mark = ",")))
cat(sprintf("%-20s %14s %14s\n", "maxSSS binary",
    format(round(orig_arp_maxsss), big.mark = ","),
    format(round(biased_arp_maxsss), big.mark = ",")))
cat(sprintf("%-20s %14s %14s\n", "p10 binary",
    format(round(orig_arp_p10), big.mark = ","),
    format(round(biased_arp_p10), big.mark = ",")))

pct_shift <- round(100 * (biased_arp_weighted - orig_arp_weighted) / orig_arp_weighted, 1)
cat("Risk-weighted shift (masked):", pct_shift, "%\n")

# Unmasked ARP — full prediction surface across Sudan
# This produces the 13M-to-17.7M bracket: the original model predicts across
# the full covariate extent (effectively self-masking since desert gets ~0),
# while the bias-corrected model assigns suitability to Nile-corridor areas
# that share fragments of the Gedaref environmental profile.
orig_arp_unmasked   <- global(pop_aligned * suit_mx, "sum", na.rm = TRUE)[[1]]
biased_arp_unmasked <- global(pop_aligned * suit_biased, "sum", na.rm = TRUE)[[1]]

pct_shift_unmasked <- round(100 * (biased_arp_unmasked - orig_arp_unmasked) / orig_arp_unmasked, 1)

cat("\n--- ARP comparison (full surface, unmasked) ---\n")
cat(sprintf("%-20s %14s %14s\n", "", "Original", "Biased-bg"))
cat(sprintf("%-20s %14s %14s\n", "Risk-weighted",
    format(round(orig_arp_unmasked), big.mark = ","),
    format(round(biased_arp_unmasked), big.mark = ",")))
cat("Risk-weighted shift (unmasked):", pct_shift_unmasked, "%\n")

# ----------------------------- Save summary ---------------------------------

bias_summary <- data.frame(
  metric     = c("Presences", "Background", "Coefficients",
                  "CBI (spatial CV)", "AUC (spatial CV)",
                  "Mean suitability", "ARP risk-weighted",
                  "ARP maxSSS binary", "ARP p10 binary",
                  "p10 threshold", "maxSSS threshold",
                  "Surface correlation (Pearson)", "Bias transform"),
  original   = c(98, 10000, 13,
                  0.857, 0.776,
                  round(global(suit_mx, "mean", na.rm = TRUE)[[1]], 4),
                  round(orig_arp_weighted), round(orig_arp_maxsss),
                  round(orig_arp_p10),
                  round(orig_p10, 4), round(orig_maxsss, 4),
                  NA, "none"),
  biased_bg  = c(98, nrow(bg_biased_env),
                  sum(mod_biased$betas != 0),
                  round(mean(fold_cbi, na.rm = TRUE), 3),
                  round(mean(fold_auc), 3),
                  round(global(suit_biased, "mean", na.rm = TRUE)[[1]], 4),
                  round(biased_arp_weighted), round(biased_arp_maxsss),
                  round(biased_arp_p10),
                  round(biased_p10, 4), round(biased_maxsss, 4),
                  r_pearson, "1/sqrt(1+tt)")
)

write.csv(bias_summary, file.path(DIR_TABLES, "bias_correction_summary.csv"),
          row.names = FALSE)

# -------------------- Log-transform sensitivity -----------------------------

bias_log <- 1 / log1p(1 + tt_aligned)
bias_log <- mask(bias_log, mask_r)

set.seed(SEED)
bg_log <- spatSample(bias_log, size = N_BACKGROUND, method = "weights",
                     na.rm = TRUE, as.points = TRUE)
bg_log_df <- as.data.frame(bg_log, geom = "XY") |>
  rename(longitude = x, latitude = y) |>
  select(longitude, latitude)

bg_log_df$year <- sample(year_weights$year, size = nrow(bg_log_df),
                         replace = TRUE, prob = year_weights$weight)
bg_log_pts <- vect(bg_log_df, geom = c("longitude", "latitude"), crs = "EPSG:4326")
bg_log_env <- terra::extract(cov_stack_mean, bg_log_pts, ID = FALSE)

complete <- complete.cases(bg_log_env)
bg_log_env <- bg_log_env[complete, ]

p_mat <- as.matrix(train$occ_env[, vars])
b_mat <- as.matrix(bg_log_env[, vars])

mod_log <- maxnet(
  p    = c(rep(1, nrow(p_mat)), rep(0, nrow(b_mat))),
  data = as.data.frame(rbind(p_mat, b_mat)),
  f    = maxnet.formula(
    p    = c(rep(1, nrow(p_mat)), rep(0, nrow(b_mat))),
    data = as.data.frame(rbind(p_mat, b_mat)),
    classes = best_classes
  ),
  regmult = tuning$rm
)

suit_log <- predict(cov_stack_mean, mod_log, type = "cloglog",
                    clamp = TRUE, na.rm = TRUE)

log_arp_weighted <- global(pop_aligned * mask(suit_log, mask_r), "sum", na.rm = TRUE)[[1]]

tt_log <- terra::extract(tt_raw, vect(bg_log_df[complete, ], geom = c("longitude", "latitude"),
                                      crs = "EPSG:4326"))[, 2]

cat("\n--- Bias transform sensitivity ---\n")
cat("Sqrt \u2014 bg median:", round(median(tt_biased, na.rm = TRUE)), "min, ARP:",
    format(round(biased_arp_weighted), big.mark = ","), "\n")
cat("Log  \u2014 bg median:", round(median(tt_log, na.rm = TRUE)), "min, ARP:",
    format(round(log_arp_weighted), big.mark = ","), "\n")

saveRDS(mod_log, file.path(DIR_MODELS, "maxent_log_bias_background.rds"))

# -------------------- State-level ARP comparison ----------------------------

adm1 <- geodata::gadm(country = "SDN", level = 1, path = here::here("data", "raw"))

rw_orig   <- pop_aligned * suit_mx_masked
rw_biased <- pop_aligned * suit_biased_masked

state_comparison <- data.frame(
  state      = adm1$NAME_1,
  total_pop  = terra::extract(pop_aligned, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]],
  arp_orig   = terra::extract(rw_orig, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]],
  arp_biased = terra::extract(rw_biased, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]]
) |>
  mutate(
    delta      = arp_biased - arp_orig,
    pct_change = round(100 * delta / arp_orig, 1),
    pct_orig   = round(100 * arp_orig / total_pop, 1),
    pct_biased = round(100 * arp_biased / total_pop, 1)
  ) |>
  arrange(desc(delta))

cat("\nState-level ARP: original vs accessibility-corrected (risk-weighted):\n")
state_comparison |>
  mutate(across(c(total_pop, arp_orig, arp_biased, delta),
                ~ format(round(.), big.mark = ","))) |>
  print(right = FALSE)

write.csv(state_comparison, file.path(DIR_TABLES, "arp_by_state_bias_comparison.csv"),
          row.names = FALSE)

# -------------------- Suitability map figure --------------------------------

adm0  <- gadm(country = "SDN", level = 0, path = here::here("data", "raw"))
sudan <- st_as_sf(adm0)
occ   <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))

suit_biased_sudan <- mask(suit_biased, vect(sudan))

pred_df <- as.data.frame(suit_biased_sudan, xy = TRUE)
names(pred_df) <- c("x", "y", "suitability")
pred_df <- pred_df[!is.na(pred_df$suitability), ]

p_suit <- ggplot() +
  geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
  geom_raster(data = pred_df, aes(x = x, y = y, fill = suitability)) +
  scale_fill_gradientn(
    colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7",
                "#EF8A62", "#B2182B"),
    na.value = "transparent",
    name = "Habitat\nsuitability",
    limits = c(0, 1)
  ) +
  geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
  geom_point(data = occ, aes(x = longitude, y = latitude),
             colour = "black", fill = "white",
             shape = 21, size = 1.5, stroke = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_minimal()) +
  labs(title = "VL habitat suitability \u2014 accessibility-corrected background",
       subtitle = "MaxEnt (LQH, rm = 1.0); background sampled \u221d 1/\u221a(1 + travel time)") +
  coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 23), crs = 4326) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank())

ggsave(file.path(DIR_FIGS, "suitability_biased_bg.png"), p_suit,
       width = 10, height = 8, dpi = 300)
cat("Saved suitability_biased_bg.png\n")

# ------------------- Response curves comparison -----------------------------

bg_medians_orig   <- apply(train$bg_env[, vars], 2, median)
bg_medians_biased <- apply(bg_biased_env[, vars], 2, median)

bg_all <- rbind(train$bg_env[, vars], bg_biased_env[, vars])
n_pts  <- 200

var_labels <- c(
  slope      = "Slope (degrees)",
  river_dist = "Distance to river (m)",
  vertisols  = "Vertisols (0/1)",
  lst_night  = "LST night (\u00b0C)",
  rainfall   = "Rainfall (mm/yr)"
)

build_curves <- function(model, bg_medians, label) {
  bind_rows(lapply(vars, function(var) {
    if (var == "vertisols") {
      newdata <- as.data.frame(t(replicate(2, bg_medians)))
      newdata[[var]] <- c(0, 1)
    } else {
      newdata <- as.data.frame(t(replicate(n_pts, bg_medians)))
      newdata[[var]] <- seq(min(bg_all[[var]]), max(bg_all[[var]]),
                            length.out = n_pts)
    }
    newdata$suitability <- predict(model, newdata, clamp = TRUE,
                                   type = "cloglog")
    tibble(
      variable = var,
      value    = newdata[[var]],
      suit     = as.numeric(newdata$suitability),
      model    = label
    )
  }))
}

curves <- bind_rows(
  build_curves(mod, bg_medians_orig, "Uniform background"),
  build_curves(mod_biased, bg_medians_biased, "Accessibility-corrected")
) |>
  mutate(var_label = var_labels[variable],
         var_label = factor(var_label, levels = var_labels))

p_curves <- ggplot(curves, aes(x = value, y = suit, colour = model)) +
  geom_line(data = curves |> filter(variable != "vertisols"),
            linewidth = 0.9) +
  geom_point(data = curves |> filter(variable == "vertisols"),
             size = 3) +
  facet_wrap(~ var_label, scales = "free_x", nrow = 2) +
  scale_colour_manual(
    values = c("Uniform background" = "#2166AC",
               "Accessibility-corrected" = "#B2182B"),
    name = NULL
  ) +
  labs(x = NULL,
       y = "Habitat suitability (cloglog)",
       title = "Marginal response curves \u2014 uniform vs accessibility-corrected background",
       subtitle = "Each covariate varied across shared range; others held at respective background median") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(DIR_FIGS, "response_curves_bias_comparison.png"), p_curves,
       width = 10, height = 6, dpi = 300)
cat("Saved response_curves_bias_comparison.png\n")

cat("12_sampling_bias_robustness.R complete\n")