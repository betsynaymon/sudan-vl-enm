# ============================================================================
# 17_gbt_comparator.R
# Fits a gradient boosted tree comparator as a third independent algorithm
# alongside MaxEnt and RF. GBT is a sequential boosted ensemble — 
# mathematically distinct from both. Produces three-way surface comparison,
# ARP estimates, and partial dependence plots.
#
# Inputs:  outputs/models/training_data.rds
#          outputs/models/spatial_cv_folds.rds
#          outputs/models/maxent_final.rds
#          outputs/models/rf_final.rds
#          outputs/models/retained_vars.rds
#          outputs/surfaces/maxent_suitability.tif
#          outputs/surfaces/rf_suitability.tif
#          outputs/surfaces/worldpop_2025_aligned.tif
#          outputs/tables/arp_full_comparison.csv
#          data/processed/occurrences_thinned.csv
# Outputs: outputs/models/gbt_final.rds
#          outputs/surfaces/gbt_suitability.tif
#          outputs/figures/suitability_three_models.png
#          outputs/figures/pdp_three_models.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(gbm)
  library(ggplot2)
  library(patchwork)
  library(ggspatial)
  library(ecospat)
  library(maxnet)
  library(ranger)
  library(geodata)
})

# ------------------------------ Load inputs ---------------------------------

train    <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
folds    <- readRDS(file.path(DIR_MODELS, "spatial_cv_folds.rds"))
suit_mx  <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))
suit_rf  <- rast(file.path(DIR_SURFACES, "rf_suitability.tif"))
vars     <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))

cat("Presences:", nrow(train$occ_env),
    "| Background:", nrow(train$bg_env), "\n")

# ---------------------- Training / fold alignment ---------------------------

df <- bind_rows(
  bind_cols(train$occ_env, train$occ_clean[, c("longitude", "latitude")]) |>
    mutate(pa = 1),
  bind_cols(train$bg_env, train$bg_clean[, c("longitude", "latitude")]) |>
    mutate(pa = 0)
)

occ_orig   <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
occ_keys   <- paste(round(occ_orig$longitude, 5), round(occ_orig$latitude, 5), occ_orig$year)
train_keys <- paste(round(train$occ_clean$longitude, 5), round(train$occ_clean$latitude, 5), train$occ_clean$year)
dropped_idx <- which(!occ_keys %in% train_keys)

fold_ids <- folds$folds_ids[-(dropped_idx)]

stopifnot(length(fold_ids) == nrow(df))
df$fold <- fold_ids

n_pres <- sum(df$pa == 1)
n_bg   <- sum(df$pa == 0)
df$weight <- ifelse(df$pa == 1, n_bg / n_pres, 1)

cat("Presences per fold:", table(df$fold[df$pa == 1]), "\n")

# ---------------------- Spatial CV tuning ------------------------------------

lr_grid    <- c(0.001, 0.005, 0.01)
depth_grid <- c(1, 3, 5)
n_trees    <- 5000
tree_steps <- seq(500, n_trees, by = 500)
bag_frac   <- 0.75

results <- data.frame()

for (lr in lr_grid) {
  for (depth in depth_grid) {
    fold_cbi <- matrix(NA, nrow = K_FOLDS, ncol = length(tree_steps))
    fold_auc <- matrix(NA, nrow = K_FOLDS, ncol = length(tree_steps))

    for (k in 1:K_FOLDS) {
      idx_train <- df$fold != k
      idx_test  <- df$fold == k

      set.seed(SEED)
      gbm_k <- gbm(
        pa ~ slope + river_dist + vertisols + lst_night + rainfall,
        data              = df[idx_train, ],
        distribution      = "bernoulli",
        weights           = df$weight[idx_train],
        n.trees           = n_trees,
        interaction.depth = depth,
        shrinkage         = lr,
        bag.fraction      = bag_frac,
        n.minobsinnode    = 10,
        verbose           = FALSE
      )

      for (t in seq_along(tree_steps)) {
        pred_test <- predict(gbm_k, newdata = df[idx_test, ],
                             n.trees = tree_steps[t], type = "response")

        pres_pred <- pred_test[df$pa[idx_test] == 1]
        bg_pred   <- pred_test[df$pa[idx_test] == 0]

        boyce <- tryCatch(
          ecospat::ecospat.boyce(fit = pred_test, obs = pres_pred,
                                nclass = 0, PEplot = FALSE),
          error = function(e) list(cor = NA_real_)
        )
        fold_cbi[k, t] <- boyce$cor

        n1 <- length(pres_pred)
        n0 <- length(bg_pred)
        fold_auc[k, t] <- (sum(rank(c(pres_pred, bg_pred))[1:n1]) -
                            n1 * (n1 + 1) / 2) / (n1 * n0)
      }
    }

    for (t in seq_along(tree_steps)) {
      results <- rbind(results, data.frame(
        lr = lr, depth = depth, n.trees = tree_steps[t],
        mean_cbi = mean(fold_cbi[, t], na.rm = TRUE),
        sd_cbi   = sd(fold_cbi[, t], na.rm = TRUE),
        mean_auc = mean(fold_auc[, t], na.rm = TRUE),
        sd_auc   = sd(fold_auc[, t], na.rm = TRUE)
      ))
    }

    best_t <- tree_steps[which.max(colMeans(fold_cbi, na.rm = TRUE))]
    best_cbi <- max(colMeans(fold_cbi, na.rm = TRUE), na.rm = TRUE)
    cat(sprintf("lr=%.3f depth=%d | best n.trees=%d | CBI=%.3f\n",
                lr, depth, best_t, best_cbi))
  }
}

best <- results[which.max(results$mean_cbi), ]
cat("\nSelected: lr =", best$lr, "depth =", best$depth,
    "n.trees =", best$n.trees,
    "| CBI =", round(best$mean_cbi, 3),
    "| AUC =", round(best$mean_auc, 3), "\n")

# ----------------------- Fit final GBT model --------------------------------

set.seed(SEED)
gbt_final <- gbm(
  pa ~ slope + river_dist + vertisols + lst_night + rainfall,
  data              = df,
  distribution      = "bernoulli",
  weights           = df$weight,
  n.trees           = best$n.trees,
  interaction.depth = best$depth,
  shrinkage         = best$lr,
  bag.fraction      = bag_frac,
  n.minobsinnode    = 10,
  verbose           = FALSE
)

imp <- summary(gbt_final, plotit = FALSE)
cat("\nRelative influence:\n")
for (i in 1:nrow(imp)) {
  cat(sprintf("  %-12s %.1f%%\n", imp$var[i], imp$rel.inf[i]))
}

saveRDS(gbt_final, file.path(DIR_MODELS, "gbt_final.rds"))

# ----------------------- Predict GBT surface --------------------------------

suit_gbt_path <- file.path(DIR_SURFACES, "gbt_suitability.tif")

if (file.exists(suit_gbt_path)) {
  cat("Loading cached GBT suitability surface\n")
  suit_gbt <- rast(suit_gbt_path)
} else {
  cov_stack <- rast(file.path(DIR_COVARIATES, COV_FILES[vars]))
  names(cov_stack) <- vars

  suit_gbt <- predict(cov_stack, gbt_final, n.trees = best$n.trees,
                      type = "response", na.rm = TRUE)

  writeRaster(suit_gbt, suit_gbt_path, overwrite = TRUE)
  cat("Computed and saved GBT suitability surface\n")
}

cat("GBT suitability range:", round(global(suit_gbt, "min", na.rm = TRUE)[[1]], 4),
    "\u2013", round(global(suit_gbt, "max", na.rm = TRUE)[[1]], 4), "\n")

# ---------------------- Surface comparison ----------------------------------

set.seed(SEED)
valid_cells <- which(!is.na(values(suit_mx)) & !is.na(values(suit_gbt)))
samp_idx <- sample(valid_cells, min(50000, length(valid_cells)))

cat("\nSurface correlation (GBT vs MaxEnt, 50k sample):\n")
cat("  Pearson: ", round(cor(values(suit_mx)[samp_idx],
                             values(suit_gbt)[samp_idx]), 3), "\n")
cat("  Spearman:", round(cor(values(suit_mx)[samp_idx],
                             values(suit_gbt)[samp_idx],
                             method = "spearman"), 3), "\n")

# ----------------------- ARP comparison -------------------------------------

pop_aligned <- rast(file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"))

# GBT thresholds
gbt_pred_occ <- predict(gbt_final, newdata = df[df$pa == 1, ],
                        n.trees = best$n.trees, type = "response")
gbt_pred_bg  <- predict(gbt_final, newdata = df[df$pa == 0, ],
                        n.trees = best$n.trees, type = "response")

gbt_p10 <- unname(quantile(gbt_pred_occ, 0.10))
candidates <- sort(unique(c(gbt_pred_occ, gbt_pred_bg)))
sens <- sapply(candidates, function(t) mean(gbt_pred_occ >= t))
spec <- sapply(candidates, function(t) mean(gbt_pred_bg  <  t))
gbt_maxsss <- candidates[which.max(sens + spec)]

gbt_arp_p10      <- global(pop_aligned * (suit_gbt >= gbt_p10), "sum", na.rm = TRUE)[[1]]
gbt_arp_maxsss   <- global(pop_aligned * (suit_gbt >= gbt_maxsss), "sum", na.rm = TRUE)[[1]]
gbt_arp_weighted <- global(pop_aligned * suit_gbt, "sum", na.rm = TRUE)[[1]]

# Load MaxEnt and RF ARPs
arp_comp <- read.csv(file.path(DIR_TABLES, "arp_full_comparison.csv"))

cat("\n--- ARP comparison (three models) ---\n")
cat(sprintf("%-20s %12s %12s %12s\n", "", "MaxEnt", "RF", "GBT"))
cat(sprintf("%-20s %12s %12s %12s\n", "Risk-weighted",
    format(round(arp_comp$maxent[arp_comp$metric == "risk_weighted"]), big.mark = ","),
    format(round(arp_comp$rf[arp_comp$metric == "risk_weighted"]), big.mark = ","),
    format(round(gbt_arp_weighted), big.mark = ",")))
cat(sprintf("%-20s %12s %12s %12s\n", "maxSSS binary",
    format(round(arp_comp$maxent[arp_comp$metric == "maxSSS"]), big.mark = ","),
    format(round(arp_comp$rf[arp_comp$metric == "maxSSS"]), big.mark = ","),
    format(round(gbt_arp_maxsss), big.mark = ",")))
cat(sprintf("%-20s %12s %12s %12s\n", "p10 binary",
    format(round(arp_comp$maxent[arp_comp$metric == "p10"]), big.mark = ","),
    format(round(arp_comp$rf[arp_comp$metric == "p10"]), big.mark = ","),
    format(round(gbt_arp_p10), big.mark = ",")))

# -------------------- Three-way suitability maps ----------------------------

tuning <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))
adm0   <- gadm(country = "SDN", level = 0, path = here::here("data", "raw"))
sudan  <- st_as_sf(adm0)

suit_mx_m  <- mask(suit_mx, vect(sudan))
suit_rf_m  <- mask(suit_rf, vect(sudan))
suit_gbt_m <- mask(suit_gbt, vect(sudan))

mx_df  <- as.data.frame(suit_mx_m, xy = TRUE, na.rm = TRUE)
rf_df  <- as.data.frame(suit_rf_m, xy = TRUE, na.rm = TRUE)
gbt_df <- as.data.frame(suit_gbt_m, xy = TRUE, na.rm = TRUE)
names(mx_df)[3] <- names(rf_df)[3] <- names(gbt_df)[3] <- "suitability"

suit_colours <- c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7",
                  "#EF8A62", "#B2182B")

p_mx <- ggplot() +
  geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
  geom_raster(data = mx_df, aes(x, y, fill = suitability)) +
  scale_fill_gradientn(colours = suit_colours, limits = c(0, 1),
                       na.value = "transparent", name = "Habitat\nsuitability") +
  geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
  coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 22.5), crs = 4326) +
  labs(title = paste0("MaxEnt (LQH, rm = ", tuning$rm, ")")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "none")

p_rf <- ggplot() +
  geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
  geom_raster(data = rf_df, aes(x, y, fill = suitability)) +
  scale_fill_gradientn(colours = suit_colours, limits = c(0, 1),
                       na.value = "transparent", name = "Habitat\nsuitability") +
  geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
  coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 22.5), crs = 4326) +
  labs(title = "Random Forest (mtry = 1)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        legend.position = "none")

p_gbt <- ggplot() +
  geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
  geom_raster(data = gbt_df, aes(x, y, fill = suitability)) +
  scale_fill_gradientn(colours = suit_colours, limits = c(0, 1),
                       na.value = "transparent", name = "Habitat\nsuitability") +
  geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 22.5), crs = 4326) +
  labs(title = paste0("GBT (lr = ", best$lr, ", depth = ", best$depth, ")")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank())

p_maps <- p_mx + p_rf + p_gbt +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(DIR_FIGS, "suitability_three_models.png"), p_maps,
       width = 15, height = 6, dpi = 300, bg = "white")
cat("Saved suitability_three_models.png\n")

# -------------------- Partial dependence plots ------------------------------

mod      <- readRDS(file.path(DIR_MODELS, "maxent_final.rds"))
rf_final <- readRDS(file.path(DIR_MODELS, "rf_final.rds"))

set.seed(SEED)
samp_rows <- sample(nrow(df), 500)
base_data <- df[samp_rows, vars]

compute_pdp_3 <- function(var_name, grid_n = 50) {
  var_seq <- seq(min(df[[var_name]]), max(df[[var_name]]), length.out = grid_n)

  mx_means <- rf_means <- gbt_means <- numeric(grid_n)

  for (i in seq_along(var_seq)) {
    modified <- base_data
    modified[[var_name]] <- var_seq[i]

    mx_means[i]  <- mean(predict(mod, modified, type = "cloglog")[, 1])
    rf_means[i]  <- mean(predict(rf_final, data = modified)$predictions[, 2])
    gbt_means[i] <- mean(predict(gbt_final, newdata = modified,
                                 n.trees = best$n.trees, type = "response"))
  }

  data.frame(
    variable    = var_name,
    value       = rep(var_seq, 3),
    suitability = c(mx_means, rf_means, gbt_means),
    model       = rep(c("MaxEnt", "RF", "GBT"), each = grid_n)
  )
}

cat("Computing PDPs...\n")
pdp_all <- bind_rows(lapply(vars, function(v) {
  cat("  ", v, "\n")
  compute_pdp_3(v)
}))

var_labels <- c(
  slope      = "Slope (\u00b0)",
  river_dist = "River distance (m)",
  vertisols  = "Vertisols (binary)",
  lst_night  = "LST night (\u00b0C)",
  rainfall   = "Rainfall (mm)"
)
pdp_all$facet_label <- var_labels[pdp_all$variable]
pdp_all$model <- factor(pdp_all$model, levels = c("MaxEnt", "RF", "GBT"))

p_pdp <- ggplot(pdp_all, aes(x = value, y = suitability, colour = model)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ facet_label, scales = "free", ncol = 3) +
  scale_colour_manual(values = c("MaxEnt" = "steelblue",
                                 "RF" = "firebrick",
                                 "GBT" = "darkorange")) +
  labs(x = NULL, y = "Predicted suitability", colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(DIR_FIGS, "pdp_three_models.png"), p_pdp,
       width = 10, height = 7, dpi = 300, bg = "white")
cat("Saved pdp_three_models.png\n")

cat("\n17_gbt_comparator.R complete\n")