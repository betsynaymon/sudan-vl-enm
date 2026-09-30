# ============================================================================
# 05_spatial_cv_folds.R
# Estimates covariate autocorrelation range, shows presence balance across
# candidate block sizes, builds spatial block CV folds, and saves them keyed
# by point ID so downstream scripts join folds to data rather than relying
# on row order.
#
# Inputs:  OCC_FILE, BG_FILE, DOMAIN_FILE, ADM0_FILE, retained_vars.rds
# Outputs: FOLD_TABLE_FILE (pa, id, fold): what downstream scripts use
#          FOLDS_FILE (blockCV object; its blocks assign folds to new points)
#          outputs/figures/spatial_cv_fold_map.png / .pdf
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
occ <- read.csv(OCC_FILE)
bg  <- read.csv(BG_FILE)

stopifnot(
  "Presence IDs are not unique"   = !anyDuplicated(occ$coordinate_id),
  "Background IDs are not unique" = !anyDuplicated(bg$bg_id)
)

pts <- bind_rows(
  occ %>% transmute(pa = 1, id = coordinate_id, longitude, latitude),
  bg  %>% transmute(pa = 0, id = bg_id,         longitude, latitude)
)
cat("Presences:", sum(pts$pa == 1), "| Background:", sum(pts$pa == 0), "\n")

pts_sf <- st_as_sf(pts, coords = c("longitude", "latitude"), crs = 4326)

# --------------------------- Load covariates --------------------------------
retained_vars <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
cat("Retained:", paste(retained_vars, collapse = ", "), "\n")

covs <- rast(file.path(DIR_COVARIATES, COV_FILES[retained_vars]))
names(covs) <- retained_vars

covs_dom <- mask(covs, rast(DOMAIN_FILE), maskvalues = 0)

cat("Stack:", nlyr(covs_dom), "layers at", res(covs_dom)[1], "deg\n")

# ----------------------- Autocorrelation range ------------------------------
# Empirical variograms for each retained covariate, reported for context.
# Block size is set by fold balance (params.R): ranges for covariates that
# vary as regional gradients do not define a feasible independence distance.

set.seed(SEED)
sac <- cv_spatial_autocor(r = covs_dom, num_sample = SAC_SAMPLE_N)

ranges_km <- sac$range_table$range / 1000
names(ranges_km) <- sac$range_table$layers

cat("Per-covariate autocorrelation ranges (km):\n")
print(round(sort(ranges_km)))
cat("\nMedian:", round(median(ranges_km)), "km\n")

# ------------------------ Block size vs balance -----------------------------
# Presence balance across folds at each candidate block size.

for (km in BLOCK_TEST_KM) {
  f <- tryCatch(
    cv_spatial(x = pts_sf, column = "pa", size = km * 1000, k = K_FOLDS,
               hexagon = FALSE, selection = "random", iteration = FOLD_ITER,
               seed = SEED, plot = FALSE, report = FALSE, progress = FALSE),
    error = function(e) NULL)
  if (is.null(f)) {
    cat(sprintf("  %4d km: could not form %d folds\n", km, K_FOLDS)); next
  }
  n <- table(factor(f$folds_ids[pts$pa == 1], levels = seq_len(K_FOLDS)))
  cat(sprintf("  %4d km: test presences per fold %s\n", km, paste(n, collapse = " / ")))
}

# ------------------------------ Assign folds --------------------------------
# cv_spatial keeps the most balanced of FOLD_ITER random block-to-fold
# assignments. Regenerated every run from the current points.

folds <- cv_spatial(
  x = pts_sf, column = "pa", size = BLOCK_SIZE_M, k = K_FOLDS,
  hexagon = FALSE, selection = "random", iteration = FOLD_ITER, seed = SEED
)
stopifnot("Fold vector does not match the points" =
  length(folds$folds_ids) == nrow(pts) && !anyNA(folds$folds_ids))

fold_table <- pts %>% select(pa, id) %>% mutate(fold = folds$folds_ids)
saveRDS(fold_table, FOLD_TABLE_FILE)
saveRDS(folds, FOLDS_FILE)

cat("Block size:", BLOCK_SIZE_M / 1000, "km | folds:", K_FOLDS,
    "| iterations:", FOLD_ITER, "\n")
print(folds$records)

test_pres <- table(factor(fold_table$fold[fold_table$pa == 1], levels = seq_len(K_FOLDS)))
cat("Test presences per fold:", test_pres, "\n")
stopifnot("A fold has fewer than MIN_TEST_PRES test presences" =
  min(test_pres) >= MIN_TEST_PRES)

# ------------------------------- Fold map -----------------------------------
pts_sf$fold <- folds$folds_ids

pres_sf <- pts_sf[pts_sf$pa == 1, ]
bg_sf   <- pts_sf[pts_sf$pa == 0, ]

p <- ggplot() +
  geom_sf(data = folds$blocks, fill = NA, colour = "grey60", linewidth = 0.3) +
  geom_sf(data = sf::st_as_sf(vect(ADM0_FILE)), fill = NA, colour = "grey30") +
  geom_sf(data = bg_sf, aes(colour = factor(fold)),
          size = 0.3, alpha = 0.15) +
  geom_sf(data = pres_sf, aes(colour = factor(fold)),
          size = 2, shape = 17) +
  scale_colour_brewer(palette = "Set1", name = "Fold") +
  labs(title = "Spatial block CV fold assignment",
       subtitle = paste0(BLOCK_SIZE_M / 1000, " km blocks, k = ", K_FOLDS, " | ",
                  sum(pts_sf$pa == 1), " presences, ",
                  sum(pts_sf$pa == 0), " background")) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(DIR_FIGS, "spatial_cv_fold_map.pdf"), p,
       width = 8, height = 6)
ggsave(file.path(DIR_FIGS, "spatial_cv_fold_map.png"), p,
       width = 8, height = 6, dpi = 300)
cat("Saved fold map\n")

cat("05_spatial_cv_folds.R complete\n")