# ============================================================================
# 04_background_sampling.R
# Samples uniform-random background points within the ecological mask and
# assigns each a year drawn from the occurrence-year distribution for
# year-matched covariate extraction downstream.
#
# Inputs:  data/processed/occurrences_thinned.csv
#          data/raw/ecological_mask_150mm.tif
# Outputs: data/processed/background_points.csv
#          outputs/figures/background_vs_occurrences.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------ Load inputs ---------------------------------

occ <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
cat("Occurrence records:", nrow(occ), "\n")

stopifnot(
  "No occurrence records found — run 02_spatial_thinning.R first" =
    nrow(occ) > 0
)

# Ecological mask: 1 = inside study area, 0 = excluded
mask <- rast(file.path(DIR_COVARIATES, "ecological_mask_150mm.tif"))
mask <- subst(mask, 0, NA)
cat("Mask cells (valid):", global(mask, fun = "notNA") %>% pull(notNA), "\n")

# -------------------------- Sample background -------------------------------

set.seed(SEED)

bg_pts <- spatSample(mask, size = N_BACKGROUND, method = "random",
                     na.rm = TRUE, as.points = TRUE)

bg <- as.data.frame(bg_pts, geom = "XY") %>%
  rename(longitude = x, latitude = y) %>%
  select(longitude, latitude)

cat("Background points sampled:", nrow(bg), "\n")

stopifnot(
  "Background sample is empty — check mask" = nrow(bg) > 0
)

# Assign years weighted by occurrence-year distribution
year_weights <- occ %>% count(year, name = "weight")

bg$year <- sample(year_weights$year, size = nrow(bg), replace = TRUE,
                  prob = year_weights$weight)

cat("\nYear distribution comparison:\n")
occ_pct <- occ %>% count(year) %>% mutate(occ_pct = round(100 * n / sum(n), 1))
bg_pct  <- bg  %>% count(year) %>% mutate(bg_pct  = round(100 * n / sum(n), 1))
left_join(occ_pct, bg_pct, by = "year", suffix = c("_occ", "_bg")) %>%
  select(year, n_occ, occ_pct, n_bg, bg_pct) %>%
  as.data.frame() %>% print()

# ----------------------------- Save and map ---------------------------------

dir.create(here::here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
write.csv(bg, here::here("data", "processed", "background_points.csv"),
          row.names = FALSE)
cat("Saved", nrow(bg), "background points\n")

p <- ggplot() +
  geom_point(data = bg, aes(longitude, latitude),
             colour = "grey70", size = 0.3, alpha = 0.3) +
  geom_point(data = occ, aes(longitude, latitude),
             colour = "firebrick", size = 1.5) +
  coord_fixed() +
  labs(title = paste0("Background (n=", nrow(bg),
                      ") and occurrences (n=", nrow(occ), ")"),
       x = "Longitude", y = "Latitude") +
  theme_minimal()

ggsave(file.path(DIR_FIGS, "background_vs_occurrences.png"), p,
       width = 8, height = 7, dpi = 300)
cat("Saved background_vs_occurrences.png\n")

cat("04_background_sampling.R complete\n")