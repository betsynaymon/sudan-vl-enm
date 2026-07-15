# ============================================================================
# 03_collinearity.R
# Screens the full candidate covariate set for collinearity, selects the final
# predictor set on ecological grounds, and writes the retained variable names
# to disk for all downstream scripts. Screening uses the prediction-surface
# layers (static covariates + 2000–2024 long-term means).
#
# Inputs:  data/raw/ (9 candidate covariate rasters + ecological mask)
# Outputs: outputs/retained_vars.rds
#          outputs/correlation_matrix.rds  (cached, speeds up reruns)
#          outputs/collinearity_sample.rds (cached, speeds up reruns)
#          outputs/figures/collinearity_correlation_matrix.png
#          outputs/figures/collinearity_correlation_matrix.pdf
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(usdm)
  library(ggplot2)
})

set.seed(SEED)

# ----------------------- Load candidate covariates --------------------------

covariate_files <- c(
  elevation  = "elevation_1km.tif",
  slope      = "slope_1km.tif",
  river_dist = "river_distance_1km.tif",
  vertisols  = "vertisols_1km.tif",
  lst_day    = "lst_day_annual_mean_2000_2024_1km.tif",
  lst_night  = "lst_night_annual_mean_2000_2024_1km.tif",
  ndvi       = "ndvi_annual_mean_2000_2024_1km.tif",
  rainfall   = "rainfall_mean_2000_2024_1km.tif",
  treecover  = "treecover_mean_2000_2024_1km.tif"
)

covs <- rast(file.path(DIR_COVARIATES, covariate_files))
names(covs) <- names(covariate_files)

# Apply ecological mask so screening covers only the modelled area
mask_r <- rast(file.path(DIR_COVARIATES, "ecological_mask_150mm.tif"))
covs_masked <- mask(covs, mask_r, maskvalues = 0)

cat("Non-NA cells per layer:\n")
print(global(covs_masked, "notNA"))

# ----------------------- Sample and correlate -------------------------------
# This step is slow (~100k raster cells). Results are cached to disk so the
# figure and VIF sections can be rerun without resampling.

cor_mat_path <- file.path(DIR_MODELS, "correlation_matrix.rds")
sample_vals_path <- file.path(DIR_MODELS, "collinearity_sample.rds")

if (file.exists(sample_vals_path) && file.exists(cor_mat_path)) {
  cat("Loading cached collinearity sample and correlation matrix\n")
  sample_vals <- readRDS(sample_vals_path)
  cor_mat <- readRDS(cor_mat_path)
} else {
  set.seed(SEED)

  sample_vals <- spatSample(
    covs_masked,
    size = COLLIN_SAMPLE_N,
    method = "random",
    na.rm = TRUE,
    values = TRUE
  )

  cor_mat <- cor(sample_vals, use = "complete.obs", method = "pearson")

  saveRDS(sample_vals, sample_vals_path)
  saveRDS(cor_mat, cor_mat_path)
  cat("Computed and cached collinearity results\n")
}

cat("Sample:", nrow(sample_vals), "cells\n")
cat("\nPearson correlation matrix:\n")
print(round(cor_mat, 2))

# Flag pairs exceeding threshold
high_pairs <- which(abs(cor_mat) >= COR_THRESHOLD & upper.tri(cor_mat), arr.ind = TRUE)
if (nrow(high_pairs) > 0) {
  cat("\nPairs exceeding |r| >=", COR_THRESHOLD, ":\n")
  for (i in seq_len(nrow(high_pairs))) {
    r <- high_pairs[i, ]
    cat("  ", rownames(cor_mat)[r[1]], " — ", colnames(cor_mat)[r[2]],
        ": ", round(cor_mat[r[1], r[2]], 2), "\n")
  }
}

# ------------------------ VIF candidate sets --------------------------------
# Two collinearity clusters: greenness-wetness (NDVI, tree cover, rainfall)
# and temperature-elevation (LST day, LST night, elevation). Each candidate
# set resolves both clusters differently.

candidate_sets <- list(
  A_rain         = c("slope", "river_dist", "vertisols", "lst_night", "rainfall"),
  B_ndvi         = c("slope", "river_dist", "vertisols", "lst_night", "ndvi"),
  C_rain_daytemp = c("slope", "river_dist", "vertisols", "lst_night", "lst_day", "rainfall"),
  D_rain_tree    = c("slope", "river_dist", "vertisols", "lst_night", "rainfall", "treecover")
)

for (set_name in names(candidate_sets)) {
  cat("\n---", set_name, "---\n")
  print(usdm::vif(sample_vals[, candidate_sets[[set_name]]]))
}

# ---------------------- Finalise variable set -------------------------------
# Set A selected on ecological grounds:
#   slope       topographic steepness; independent of elevation (r = 0.56)
#   river_dist  distance to drainage incl. seasonal khors; independent (|r| <= 0.21)
#   vertisols   black-cotton soil; independent (|r| <= 0.25)
#   lst_night   night temperature — conditions the nocturnal vector experiences
#   rainfall    represents the moisture-vegetation axis (most stable single measure)
#
# Dropped (4), each redundant with a retained variable:
#   elevation   distal proxy; acts through temperature (r = -0.77 with lst_night)
#   lst_day     daytime thermal window; tracks lst_night (0.62) and rainfall (-0.72)
#   ndvi        near-duplicate of rainfall (0.86); conflates crop and woodland
#   treecover   redundant with rainfall (0.70); forest-built metric, near-zero in belt

retained_vars <- c("slope", "river_dist", "vertisols", "lst_night", "rainfall")

cat("\nFinal variable set VIF:\n")
print(usdm::vif(sample_vals[, retained_vars]))

saveRDS(retained_vars, file.path(DIR_MODELS, "retained_vars.rds"))
cat("Saved retained_vars.rds\n")

# -------------------- Correlation matrix figure -----------------------------

var_order  <- c("rainfall", "ndvi", "treecover",
                "lst_day", "lst_night", "elevation",
                "slope", "river_dist", "vertisols")
var_labels <- c(rainfall = "Rainfall", ndvi = "NDVI", treecover = "Tree cover",
                lst_day = "Day LST", lst_night = "Night LST", elevation = "Elevation",
                slope = "Slope", river_dist = "Distance to river", vertisols = "Vertisols")

cor_mat <- cor_mat[var_order, var_order]
cor_df <- as.data.frame(as.table(cor_mat))
names(cor_df) <- c("var1", "var2", "r")
cor_df$var1 <- factor(cor_df$var1, levels = var_order)
cor_df$var2 <- factor(cor_df$var2, levels = var_order)
cor_df <- cor_df[as.integer(cor_df$var1) >= as.integer(cor_df$var2), ]
cor_df$flag <- abs(cor_df$r) >= COR_THRESHOLD & cor_df$var1 != cor_df$var2

p <- ggplot(cor_df, aes(var2, var1, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_tile(data = cor_df[cor_df$flag, ], fill = NA,
            colour = "grey10", linewidth = 1.1) +
  geom_text(aes(label = sprintf("%.2f", r),
                colour = abs(r) > 0.5), size = 3.2) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"),
                      guide = "none") +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = seq(-1, 1, 0.2), name = "Pearson r") +
  scale_x_discrete(limits = var_order, labels = var_labels) +
  scale_y_discrete(limits = rev(var_order), labels = var_labels) +
  coord_fixed() + labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(1.6, "cm"))

ggsave(file.path(DIR_FIGS, "collinearity_correlation_matrix.pdf"),
       p, width = 7.5, height = 6.5)
ggsave(file.path(DIR_FIGS, "collinearity_correlation_matrix.png"),
       p, width = 7.5, height = 6.5, dpi = 300)
cat("Saved correlation matrix figure\n")

cat("03_collinearity.R complete\n")