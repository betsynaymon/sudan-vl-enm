# ============================================================================
# 09_rf_comparator.R
# Fits a random forest comparator using the same data, covariates, and spatial
# CV folds as MaxEnt. Tests whether the suitability surface and ARP are robust
# to algorithmic assumptions. Produces partial dependence plots comparing
# learned ecological relationships across both algorithms.
#
# Inputs:  outputs/models/training_data.rds
#          outputs/models/spatial_cv_folds.rds
#          outputs/models/maxent_final.rds
#          outputs/surfaces/maxent_suitability.tif
#          outputs/surfaces/worldpop_2025_aligned.tif
#          outputs/models/retained_vars.rds
#          outputs/tables/arp_summary.csv
# Outputs: outputs/models/rf_final.rds
#          outputs/surfaces/rf_suitability.tif
#          outputs/tables/arp_full_comparison.csv
#          outputs/figures/pdp_maxent_vs_rf.png
#          outputs/figures/suitability_maxent_vs_rf.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ranger)
  library(ggplot2)
  library(patchwork)
  library(ggspatial)
  library(ecospat)
  library(geodata)
  library(maxnet)
})

# ------------------------------ Load inputs ---------------------------------

train   <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
folds   <- readRDS(file.path(DIR_MODELS, "spatial_cv_folds.rds"))
suit_mx <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
mod     <- readRDS(file.path(DIR_MODELS, "maxent_final.rds"))

cat("Presences:", nrow(train$occ_env),
    "| Background:", nrow(train$bg_env), "\n")

# ---------------------- Training / fold alignment ---------------------------
# The fold vector (length 10,099) is one longer than the training data
# (10,098) because MaxEnt's cell×year deduplication dropped one row.
# Identify the dropped row by coordinates and remove it from the fold
# vector so both algorithms use identical fold assignments.

df <- bind_rows(
  bind_cols(train$occ_env, train$occ_clean[, c("longitude", "latitude")]) |>
    mutate(pa = 1),
  bind_cols(train$bg_env, train$bg_clean[, c("longitude", "latitude")]) |>
    mutate(pa = 0)
)

# Find which presence was dropped during cell×year dedup
occ_orig <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
occ_keys   <- paste(round(occ_orig$longitude, 5), round(occ_orig$latitude, 5), occ_orig$year)
train_keys <- paste(round(train$occ_clean$longitude, 5), round(train$occ_clean$latitude, 5), train$occ_clean$year)
dropped_idx <- which(!occ_keys %in% train_keys)

cat("Dropped row:", dropped_idx,
    "| Point:", occ_orig$longitude[dropped_idx], occ_orig$latitude[dropped_idx],
    "year:", occ_orig$year[dropped_idx], "\n")

# Remove dropped row from fold vector
fold_ids <- folds$folds_ids[-(dropped_idx)]

stopifnot(
  "Fold vector length doesn't match training data" =
    length(fold_ids) == nrow(df)
)

df$fold <- fold_ids

cat("Combined:", nrow(df), "rows (",
    sum(df$pa == 1), "pres,", sum(df$pa == 0), "bg)\n")
cat("NA fold assignments:", sum(is.na(df$fold)), "\n")
cat("Presences per fold:", table(df$fold[df$pa == 1]), "\n")

# ---------------------- Spatial CV tuning ------------------------------------
# Tune mtry with CBI under the same spatial block CV folds.
# Case weights balance classes: each presence counts as ~102 background points.

n_pres <- sum(df$pa == 1)
n_bg   <- sum(df$pa == 0)
df$weight <- ifelse(df$pa == 1, n_bg / n_pres, 1)

df$pa <- factor(df$pa, levels = c("0", "1"))

mtry_grid <- 1:5
n_trees   <- 1000

results <- data.frame()

for (m in mtry_grid) {
  fold_cbi <- fold_auc <- numeric(4)

  for (k in 1:4) {
    idx_train <- df$fold != k
    idx_test  <- df$fold == k

    rf_k <- ranger(
      pa ~ slope + river_dist + vertisols + lst_night + rainfall,
      data         = df[idx_train, ],
      case.weights = df$weight[idx_train],
      num.trees    = n_trees,
      mtry         = m,
      probability  = TRUE,
      seed         = SEED,
      verbose      = FALSE
    )

    pred_test <- predict(rf_k, data = df[idx_test, ])$predictions[, 2]
    pres_pred <- pred_test[df$pa[idx_test] == "1"]
    bg_pred   <- pred_test[df$pa[idx_test] == "0"]

    boyce <- ecospat.boyce(fit = pred_test, obs = pres_pred,
                           nclass = 0, PEplot = FALSE)
    fold_cbi[k] <- boyce$cor

    n1 <- length(pres_pred)
    n0 <- length(bg_pred)
    fold_auc[k] <- (sum(rank(c(pres_pred, bg_pred))[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }

  results <- rbind(results, data.frame(
    mtry     = m,
    mean_cbi = mean(fold_cbi, na.rm = TRUE),
    sd_cbi   = sd(fold_cbi, na.rm = TRUE),
    mean_auc = mean(fold_auc),
    sd_auc   = sd(fold_auc)
  ))
  cat("mtry =", m, " | CBI =", round(mean(fold_cbi, na.rm = TRUE), 3),
      "\u00b1", round(sd(fold_cbi, na.rm = TRUE), 3),
      " | AUC =", round(mean(fold_auc), 3), "\n")
}

best <- results[which.max(results$mean_cbi), ]
cat("\nBest: mtry =", best$mtry,
    "| CBI =", round(best$mean_cbi, 3),
    "| AUC =", round(best$mean_auc, 3), "\n")

# ----------------------- Fit final RF model ---------------------------------

best_mtry <- best$mtry

rf_final <- ranger(
  pa ~ slope + river_dist + vertisols + lst_night + rainfall,
  data         = df,
  case.weights = df$weight,
  num.trees    = 1000,
  mtry         = best_mtry,
  probability  = TRUE,
  importance   = "permutation",
  seed         = SEED
)

cat("RF final model: mtry =", best_mtry,
    "| OOB error:", round(rf_final$prediction.error, 4), "\n")

saveRDS(rf_final, file.path(DIR_MODELS, "rf_final.rds"))

# ----------------------- Predict RF surface ---------------------------------

suit_rf_path <- file.path(DIR_SURFACES, "rf_suitability.tif")

if (file.exists(suit_rf_path)) {
  cat("Loading cached RF suitability surface\n")
  suit_rf <- rast(suit_rf_path)
} else {
  cov_stack <- rast(file.path(DIR_COVARIATES, COV_FILES[vars]))
  names(cov_stack) <- vars

  suit_rf <- predict(cov_stack, rf_final, fun = function(model, ...) {
    predict(model, ...)$predictions[, 2]
  }, na.rm = TRUE)

  writeRaster(suit_rf, suit_rf_path, overwrite = TRUE)
  cat("Computed and saved RF suitability surface\n")
}

cat("RF suitability range:", round(global(suit_rf, "min", na.rm = TRUE)[[1]], 4),
    "\u2013", round(global(suit_rf, "max", na.rm = TRUE)[[1]], 4), "\n")

# ---------------------- Surface comparison ----------------------------------

set.seed(SEED)
valid_cells <- which(!is.na(values(suit_mx)) & !is.na(values(suit_rf)))
samp_idx <- sample(valid_cells, min(50000, length(valid_cells)))

mx_vals <- values(suit_mx)[samp_idx]
rf_vals <- values(suit_rf)[samp_idx]

cat("\nSurface correlation (50k sample):\n")
cat("  Pearson: ", round(cor(mx_vals, rf_vals, method = "pearson"), 3), "\n")
cat("  Spearman:", round(cor(mx_vals, rf_vals, method = "spearman"), 3), "\n")

# ----------------------- ARP comparison -------------------------------------

pop_aligned <- rast(file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"))
total_pop   <- global(pop_aligned, "sum", na.rm = TRUE)[[1]]

# RF thresholds
rf_pred_occ <- predict(rf_final, data = df[df$pa == "1", ])$predictions[, 2]
rf_pred_bg  <- predict(rf_final, data = df[df$pa == "0", ])$predictions[, 2]

rf_p10 <- unname(quantile(rf_pred_occ, 0.10))

candidates <- sort(unique(c(rf_pred_occ, rf_pred_bg)))
sens <- sapply(candidates, function(t) mean(rf_pred_occ >= t))
spec <- sapply(candidates, function(t) mean(rf_pred_bg  <  t))
rf_maxsss <- candidates[which.max(sens + spec)]

cat("\nRF p10:   ", round(rf_p10, 4), "\n")
cat("RF maxSSS:", round(rf_maxsss, 4), "\n")

rf_arp_p10     <- global(pop_aligned * (suit_rf >= rf_p10), "sum", na.rm = TRUE)[[1]]
rf_arp_maxsss  <- global(pop_aligned * (suit_rf >= rf_maxsss), "sum", na.rm = TRUE)[[1]]
rf_arp_weighted <- global(pop_aligned * suit_rf, "sum", na.rm = TRUE)[[1]]

# Load MaxEnt ARP for comparison
mx_arp <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))

cat("\n--- ARP comparison ---\n")
cat(sprintf("%-20s %12s %12s\n", "", "MaxEnt", "RF"))
cat(sprintf("%-20s %12s %12s\n", "p10 binary",
    format(round(mx_arp$arp[mx_arp$metric == "p10"]), big.mark = ","),
    format(round(rf_arp_p10), big.mark = ",")))
cat(sprintf("%-20s %12s %12s\n", "maxSSS binary",
    format(round(mx_arp$arp[mx_arp$metric == "maxSSS"]), big.mark = ","),
    format(round(rf_arp_maxsss), big.mark = ",")))
cat(sprintf("%-20s %12s %12s\n", "Risk-weighted",
    format(round(mx_arp$arp[mx_arp$metric == "risk_weighted"]), big.mark = ","),
    format(round(rf_arp_weighted), big.mark = ",")))

# Save comparison table
arp_comp <- data.frame(
  metric = c("p10", "maxSSS", "risk_weighted"),
  maxent = mx_arp$arp,
  rf     = c(round(rf_arp_p10), round(rf_arp_maxsss), round(rf_arp_weighted))
)
write.csv(arp_comp, file.path(DIR_TABLES, "arp_full_comparison.csv"), row.names = FALSE)

# -------------------- Partial dependence plots ------------------------------

set.seed(SEED)
samp_rows <- sample(nrow(df), 500)
base_data <- df[samp_rows, vars]

var_labels <- c(
  slope      = "Slope (\u00b0)",
  river_dist = "River distance (m)",
  vertisols  = "Vertisols (binary)",
  lst_night  = "LST night (\u00b0C)",
  rainfall   = "Rainfall (mm)"
)

compute_pdp <- function(var_name, grid_n = 50) {
  var_seq <- seq(min(df[[var_name]]), max(df[[var_name]]), length.out = grid_n)

  mx_means <- rf_means <- numeric(grid_n)

  for (i in seq_along(var_seq)) {
    modified <- base_data
    modified[[var_name]] <- var_seq[i]

    mx_means[i] <- mean(predict(mod, modified, type = "cloglog")[, 1])
    rf_means[i] <- mean(predict(rf_final, data = modified)$predictions[, 2])
  }

  data.frame(
    variable    = var_name,
    value       = rep(var_seq, 2),
    suitability = c(mx_means, rf_means),
    model       = rep(c("MaxEnt", "RF"), each = grid_n)
  )
}

cat("Computing PDPs...\n")
pdp_all <- bind_rows(lapply(vars, function(v) {
  cat("  ", v, "\n")
  compute_pdp(v)
}))

pdp_all$facet_label <- var_labels[pdp_all$variable]

p_pdp <- ggplot(pdp_all, aes(x = value, y = suitability, colour = model)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_label, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = c("MaxEnt" = "steelblue", "RF" = "firebrick")) +
  labs(x = NULL, y = "Predicted suitability", colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(DIR_FIGS, "pdp_maxent_vs_rf.png"), p_pdp,
       width = 10, height = 7, dpi = 300, bg = "white")
cat("Saved pdp_maxent_vs_rf.png\n")

# -------------------- Side-by-side suitability maps -------------------------

adm0 <- gadm(country = "SDN", level = 0, path = here::here("data", "raw"))
sudan <- st_as_sf(adm0)

suit_mx_masked <- mask(suit_mx, vect(sudan))
suit_rf_masked <- mask(suit_rf, vect(sudan))

mx_df <- as.data.frame(suit_mx_masked, xy = TRUE, na.rm = TRUE)
rf_df <- as.data.frame(suit_rf_masked, xy = TRUE, na.rm = TRUE)
names(mx_df)[3] <- names(rf_df)[3] <- "suitability"

suit_colours <- c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7",
                  "#EF8A62", "#B2182B")

p_mx <- ggplot() +
  geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
  geom_raster(data = mx_df, aes(x, y, fill = suitability)) +
  scale_fill_gradientn(colours = suit_colours, limits = c(0, 1),
                       na.value = "transparent", name = "Habitat\nsuitability") +
  geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
  coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 22.5), crs = 4326) +
  labs(title = "MaxEnt (LQH, rm = 1.0)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "none")

p_rf <- ggplot() +
  geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
  geom_raster(data = rf_df, aes(x, y, fill = suitability)) +
  scale_fill_gradientn(colours = suit_colours, limits = c(0, 1),
                       na.value = "transparent", name = "Habitat\nsuitability") +
  geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 22.5), crs = 4326) +
  labs(title = "Random Forest (mtry = 1)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank())

p_maps <- p_mx + p_rf +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(DIR_FIGS, "suitability_maxent_vs_rf.png"), p_maps,
       width = 12, height = 6, dpi = 300, bg = "white")
cat("Saved suitability_maxent_vs_rf.png\n")

cat("09_rf_comparator.R complete\n")