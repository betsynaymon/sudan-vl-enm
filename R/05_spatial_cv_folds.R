# ============================================================================
# 05_spatial_cv_folds.R
# Estimates covariate autocorrelation range, builds spatial block CV folds for
# presences and background points, and saves a single fold object that all
# downstream scripts share.
#
# Block size (100 km) was selected after testing larger sizes (706, 350, 200 km)
# which all produced severely imbalanced folds due to Gedaref clustering.
# 100 km with k=4 and random assignment (200 iterations) gives viable
# per-fold presence balance while maintaining spatial independence.
#
# Inputs:  data/processed/occurrences_thinned.csv
#          data/processed/background_points.csv
#          data/raw/ (retained covariate rasters + ecological mask)
#          outputs/retained_vars.rds
# Outputs: outputs/models/spatial_cv_folds.rds  (committed to repo)
#          outputs/figures/spatial_cv_fold_map.png
#          outputs/figures/spatial_cv_fold_map.pdf
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(blockCV)
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
})

set.seed(SEED)

# ------------------------------ Load data -----------------------------------

occ <- read.csv(here::here("data", "processed", "occurrences_thinned.csv")) %>%
  mutate(pa = 1)

bg <- read.csv(here::here("data", "processed", "background_points.csv")) %>%
  mutate(pa = 0)

pts <- bind_rows(
  occ %>% select(longitude, latitude, pa),
  bg  %>% select(longitude, latitude, pa)
)

cat("Presences:", sum(pts$pa == 1),
    "| Background:", sum(pts$pa == 0),
    "| Total:", nrow(pts), "\n")

pts_sf <- st_as_sf(pts, coords = c("longitude", "latitude"), crs = 4326)

# --------------------------- Load covariates --------------------------------

retained_vars <- readRDS(file.path(DIR_OUTPUTS, "retained_vars.rds"))
cat("Retained:", paste(retained_vars, collapse = ", "), "\n")

covs <- rast(file.path(DIR_COVARIATES, COV_FILES[retained_vars]))
names(covs) <- retained_vars

mask_r <- rast(file.path(DIR_COVARIATES, "ecological_mask_150mm.tif"))
covs_masked <- mask(covs, mask_r, maskvalues = 0)

cat("Stack:", nlyr(covs_masked), "layers at", res(covs_masked)[1], "deg\n")

# ----------------------- Autocorrelation range ------------------------------
# Empirical variograms for each retained covariate. The median range informs
# the minimum block size for spatial independence.

sac <- cv_spatial_autocor(r = covs_masked, num_sample = 5000)

ranges_km <- sac$range_table$range / 1000
names(ranges_km) <- sac$range_table$layers

cat("Per-covariate autocorrelation ranges (km):\n")
print(round(sort(ranges_km)))
cat("\nMedian:", round(median(ranges_km)), "km\n")

# ------------------------------ Assign folds --------------------------------
# 100 km blocks, k=4, random assignment with 200 iterations for best balance.
# Median autocorrelation range is 706 km, but blocks at that scale produce
# severely imbalanced folds because most presences cluster in the Gedaref
# corridor. 100 km is a pragmatic compromise: enough blocks to distribute
# presences across folds while still grouping nearby points together.
#
# The fold object is committed to the repo so that cloning reproduces
# the exact CV metrics reported in the dissertation. Delete the file to
# regenerate from scratch — results will differ slightly because
# cv_spatial's random block-to-fold assignment is sensitive to R's
# global RNG state, which varies across sessions and platforms even
# with the same seed parameter. The spatial structure and balance
# will be comparable; only the specific fold membership changes.

folds_path <- file.path(DIR_MODELS, "spatial_cv_folds.rds")

if (file.exists(folds_path)) {
  cat("Loading existing fold assignments from", folds_path, "\n")
  folds <- readRDS(folds_path)
} else {
  folds <- cv_spatial(
    x         = pts_sf,
    column    = "pa",
    size      = BLOCK_SIZE_M,
    k         = K_FOLDS,
    hexagon   = FALSE,
    selection = "random",
    iteration = 200,
    seed      = SEED
  )

  saveRDS(folds, folds_path)
  cat("Generated and saved new fold assignments\n")
}

cat("\nFinal configuration:\n")
cat("  Block size:", BLOCK_SIZE_M / 1000, "km\n")
cat("  Folds:", folds$k, "\n")
cat("  Selection: random (200 iterations)\n\n")
cat("Per-fold balance:\n")
print(folds$records)

# ------------------------------- Fold map -----------------------------------

pts_sf$fold <- folds$folds_ids

pres_sf <- pts_sf[pts_sf$pa == 1, ]
bg_sf   <- pts_sf[pts_sf$pa == 0, ]

p <- ggplot() +
  geom_sf(data = folds$blocks, fill = NA, colour = "grey60", linewidth = 0.3) +
  geom_sf(data = bg_sf, aes(colour = factor(fold)),
          size = 0.3, alpha = 0.15) +
  geom_sf(data = pres_sf, aes(colour = factor(fold)),
          size = 2, shape = 17) +
  scale_colour_brewer(palette = "Set1", name = "Fold") +
  labs(title = "Spatial block CV fold assignment",
       subtitle = paste0("100 km blocks, k = 4 | ",
                         sum(pts_sf$pa == 1), " presences, ",
                         sum(pts_sf$pa == 0), " background")) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(DIR_FIGS, "spatial_cv_fold_map.pdf"), p,
       width = 8, height = 6)
ggsave(file.path(DIR_FIGS, "spatial_cv_fold_map.png"), p,
       width = 8, height = 6, dpi = 300)
cat("Saved fold map\n")

# --------------------------------- Summary ----------------------------------

cat("Fold object:", file.path(DIR_MODELS, "spatial_cv_folds.rds"), "\n")
cat("Contains:", length(folds$folds_list), "folds\n")
cat("Fold sizes (test presences):",
    sapply(folds$folds_list, function(f) sum(pts_sf$pa[f[[2]]] == 1)), "\n")

cat("05_spatial_cv_folds.R complete\n")