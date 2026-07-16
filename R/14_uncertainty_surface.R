# ============================================================================
# 14_uncertainty_surface.R
# Refits MaxEnt on each spatial CV fold's training set, predicts each to the
# full raster, and computes pixelwise mean and SD across the 4 surfaces.
# Also computes fold-excluded ARP under all three estimation methods to
# quantify how much the estimate depends on any single geographic cluster.
#
# Inputs:  outputs/models/training_data.rds
#          outputs/models/spatial_cv_folds.rds
#          outputs/models/selected_tuning.rds
#          outputs/models/retained_vars.rds
#          outputs/surfaces/worldpop_2025_aligned.tif
#          outputs/tables/arp_summary.csv
#          data/raw/ (covariate rasters)
#          data/processed/occurrences_thinned.csv
# Outputs: outputs/surfaces/maxent_cv_mean.tif
#          outputs/surfaces/maxent_cv_sd.tif
#          outputs/models/cv_fold_predictions.rds  (cached)
#          outputs/tables/arp_fold_uncertainty.csv
#          outputs/figures/prediction_uncertainty_sd.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(maxnet)
  library(ggplot2)
  library(geodata)
  library(ggspatial)
})

# ------------------------------ Load inputs ---------------------------------

train   <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
folds   <- readRDS(file.path(DIR_MODELS, "spatial_cv_folds.rds"))
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
tuning  <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))

cov_stack <- rast(file.path(DIR_COVARIATES, COV_FILES))
names(cov_stack) <- names(COV_FILES)

cat("Presences:", nrow(train$occ_env), "| Background:", nrow(train$bg_env), "\n")
cat("Refitting with:", tuning$fc, "rm =", tuning$rm, "\n")

# ---------------------- Fold alignment --------------------------------------

df <- bind_rows(
  bind_cols(train$occ_env, train$occ_clean[, c("longitude", "latitude")]) |>
    mutate(pa = 1),
  bind_cols(train$bg_env, train$bg_clean[, c("longitude", "latitude")]) |>
    mutate(pa = 0)
)

occ_orig   <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
occ_keys   <- paste(round(occ_orig$longitude, 5), round(occ_orig$latitude, 5), occ_orig$year)
train_keys <- paste(round(train$occ_clean$longitude, 5), round(train$occ_clean$latitude, 5), train$occ_clean$year)
dedup_idx  <- which(!occ_keys %in% train_keys)

fold_ids <- folds$folds_ids[-(dedup_idx)]

stopifnot(
  "Fold vector length doesn't match training data" =
    length(fold_ids) == nrow(df)
)

df$fold <- fold_ids

cat("Presences per fold:", table(df$fold[df$pa == 1]), "\n")

# -------------------- CV refit and predict ----------------------------------

best_classes <- tolower(tuning$fc)

cv_pred_path <- file.path(DIR_MODELS, "cv_fold_predictions.rds")

if (file.exists(cv_pred_path)) {
  cat("Loading cached fold predictions\n")
  pred_stack <- readRDS(cv_pred_path)
} else {
  pred_stack <- list()

  for (k in 1:K_FOLDS) {
    cat("Fold", k, "\u2014 training on folds",
        paste(setdiff(1:K_FOLDS, k), collapse = ","), "...")

    train_idx <- df$fold != k

    p_vec   <- as.numeric(as.character(df$pa[train_idx]))
    env_mat <- df[train_idx, vars]

    mod_k <- maxnet(
      p       = p_vec,
      data    = env_mat,
      f       = maxnet.formula(p = p_vec, data = env_mat, classes = best_classes),
      regmult = tuning$rm
    )

    cat(" fitted (", length(mod_k$betas), "coefs) ...")

    pred_k <- predict(cov_stack, mod_k, type = "cloglog",
                      clamp = TRUE, na.rm = TRUE)
    pred_stack[[k]] <- pred_k

    cat(" predicted.\n")
  }

  saveRDS(pred_stack, cv_pred_path)
  cat("Computed and saved fold predictions\n")
}

# Compute mean and SD
pred_all  <- rast(pred_stack)
suit_mean <- app(pred_all, mean, na.rm = TRUE)
suit_sd   <- app(pred_all, sd, na.rm = TRUE)

cat("\nMean suitability \u2014 mean:", round(global(suit_mean, "mean", na.rm = TRUE)[[1]], 4),
    " range:", round(global(suit_mean, "min", na.rm = TRUE)[[1]], 4), "\u2013",
    round(global(suit_mean, "max", na.rm = TRUE)[[1]], 4), "\n")
cat("SD surface      \u2014 mean:", round(global(suit_sd, "mean", na.rm = TRUE)[[1]], 4),
    " range:", round(global(suit_sd, "min", na.rm = TRUE)[[1]], 4), "\u2013",
    round(global(suit_sd, "max", na.rm = TRUE)[[1]], 4), "\n")

writeRaster(suit_mean, file.path(DIR_SURFACES, "maxent_cv_mean.tif"), overwrite = TRUE)
writeRaster(suit_sd, file.path(DIR_SURFACES, "maxent_cv_sd.tif"), overwrite = TRUE)

# ---------------------- Uncertainty map -------------------------------------

adm0  <- gadm(country = "SDN", level = 0, path = here::here("data", "raw"))
sudan <- st_as_sf(adm0)

sd_masked <- mask(suit_sd, vect(sudan))
sd_df <- as.data.frame(sd_masked, xy = TRUE, na.rm = TRUE)
names(sd_df)[3] <- "sd"

occ_pts <- train$occ_clean[, c("longitude", "latitude")]

p_sd <- ggplot() +
  geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
  geom_raster(data = sd_df, aes(x, y, fill = sd)) +
  scale_fill_gradientn(
    colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7",
                "#EF8A62", "#B2182B"),
    na.value = "transparent",
    name = "Prediction\nSD",
    limits = c(0, max(sd_df$sd))
  ) +
  geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
  geom_point(data = occ_pts, aes(longitude, latitude),
             colour = "black", fill = "white",
             shape = 21, size = 1.5, stroke = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 22.5)) +
  labs(title = "Prediction uncertainty across spatial CV folds",
       subtitle = "SD of suitability across 4 fold-excluded refits") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.title = element_blank())

ggsave(file.path(DIR_FIGS, "prediction_uncertainty_sd.png"), p_sd,
       width = 7, height = 8, dpi = 300, bg = "white")
cat("Saved prediction_uncertainty_sd.png\n")

# -------------------- Fold-excluded ARP ------------------------------------

pop_aligned <- rast(file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"))
occ_xy <- as.matrix(train$occ_clean[, c("longitude", "latitude")])
bg_xy  <- as.matrix(train$bg_clean[, c("longitude", "latitude")])

fold_arps <- data.frame(
  fold_excluded = 1:K_FOLDS,
  arp_weighted  = numeric(K_FOLDS),
  arp_p10       = numeric(K_FOLDS),
  arp_maxsss    = numeric(K_FOLDS)
)

for (k in 1:K_FOLDS) {
  pred_k <- pred_stack[[k]]

  fold_arps$arp_weighted[k] <- global(pop_aligned * pred_k,
                                      "sum", na.rm = TRUE)[[1]]

  pred_at_occ <- terra::extract(pred_k, occ_xy)[[1]]
  pred_at_bg  <- terra::extract(pred_k, bg_xy)[[1]]

  p10_k <- quantile(pred_at_occ, 0.10, na.rm = TRUE)
  fold_arps$arp_p10[k] <- global(pop_aligned * (pred_k >= p10_k),
                                 "sum", na.rm = TRUE)[[1]]

  candidates <- sort(unique(c(pred_at_occ, pred_at_bg)))
  sens <- sapply(candidates, function(t) mean(pred_at_occ >= t, na.rm = TRUE))
  spec <- sapply(candidates, function(t) mean(pred_at_bg < t, na.rm = TRUE))
  maxsss_k <- candidates[which.max(sens + spec)]
  fold_arps$arp_maxsss[k] <- global(pop_aligned * (pred_k >= maxsss_k),
                                    "sum", na.rm = TRUE)[[1]]
}

mx_arp <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))

cat("\n--- ARP across fold-excluded models ---\n\n")
cat(sprintf("%-15s %12s %12s %12s\n",
            "Fold excluded", "Risk-weighted", "p10", "maxSSS"))
for (k in 1:K_FOLDS) {
  cat(sprintf("%-15s %12s %12s %12s\n",
              paste("Fold", k),
              format(round(fold_arps$arp_weighted[k]), big.mark = ","),
              format(round(fold_arps$arp_p10[k]), big.mark = ","),
              format(round(fold_arps$arp_maxsss[k]), big.mark = ",")))
}

cat(sprintf("\n%-15s %12s %12s %12s\n", "Mean",
    format(round(mean(fold_arps$arp_weighted)), big.mark = ","),
    format(round(mean(fold_arps$arp_p10)), big.mark = ","),
    format(round(mean(fold_arps$arp_maxsss)), big.mark = ",")))
cat(sprintf("%-15s %12s %12s %12s\n", "SD",
    format(round(sd(fold_arps$arp_weighted)), big.mark = ","),
    format(round(sd(fold_arps$arp_p10)), big.mark = ","),
    format(round(sd(fold_arps$arp_maxsss)), big.mark = ",")))

cat("\nFull-data model: weighted =",
    format(round(mx_arp$arp[mx_arp$metric == "risk_weighted"]), big.mark = ","),
    "| p10 =",
    format(round(mx_arp$arp[mx_arp$metric == "p10"]), big.mark = ","),
    "| maxSSS =",
    format(round(mx_arp$arp[mx_arp$metric == "maxSSS"]), big.mark = ","), "\n")

write.csv(fold_arps, file.path(DIR_TABLES, "arp_fold_uncertainty.csv"),
          row.names = FALSE)

cat("\n14_uncertainty_surface.R complete\n")