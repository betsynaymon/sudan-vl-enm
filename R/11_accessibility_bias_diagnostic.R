# ============================================================================
# 11_accessibility_bias_diagnostic.R
# Tests whether occurrence records are systematically biased toward accessible
# areas by comparing Weiss et al. (2018) travel-time-to-city values at
# occurrence vs. background locations. Establishes the empirical fact of the
# bias, motivating the sampling bias correction (12) and accessibility-weighted
# null model (13).
#
# Inputs:  outputs/models/training_data.rds
#          data/raw/weiss_travel_time.tif (downloaded via Malaria Atlas Project)
# Outputs: outputs/figures/accessibility_bias_diagnostic.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(httr)
})

# ------------------------------ Load inputs ---------------------------------

train <- readRDS(file.path(DIR_MODELS, "training_data.rds"))

cat("Presences:", nrow(train$occ_clean), "\n")
cat("Background:", nrow(train$bg_clean), "\n")

# ---------------------- Load travel-time surface ----------------------------

tt_path <- here::here("data", "raw", "weiss_travel_time.tif")

if (!file.exists(tt_path)) {
  wcs_url <- paste0(
    "https://data.malariaatlas.org/geoserver/Accessibility/ows?",
    "service=WCS&version=2.0.1&request=GetCoverage&format=image/geotiff&",
    "CoverageId=Accessibility__201501_Global_Travel_Time_to_Cities&",
    "subset=Long(21.5,39)&subset=Lat(8.5,22.5)"
  )

  response <- GET(
    wcs_url,
    user_agent("R - MSc dissertation, e.p.naymon@lse.ac.uk"),
    write_disk(tt_path, overwrite = TRUE),
    progress()
  )
  Sys.sleep(1)
  stopifnot("Travel-time download failed" = status_code(response) == 200)
  cat("Downloaded:", basename(tt_path), "\n")
} else {
  cat("Travel-time raster already present\n")
}

tt_raw <- rast(tt_path)

cat("Dimensions:", ncol(tt_raw), "x", nrow(tt_raw), "\n")
cat("Resolution:", res(tt_raw), "\n")
cat("Range (minutes):", round(global(tt_raw, "min", na.rm = TRUE)[[1]]),
    "\u2013", round(global(tt_raw, "max", na.rm = TRUE)[[1]]), "\n")

# --------------------- Accessibility comparison -----------------------------

occ_pts <- vect(train$occ_clean, geom = c("longitude", "latitude"), crs = "EPSG:4326")
bg_pts  <- vect(train$bg_clean,  geom = c("longitude", "latitude"), crs = "EPSG:4326")

tt_occ <- terra::extract(tt_raw, occ_pts)[, 2]
tt_bg  <- terra::extract(tt_raw, bg_pts)[, 2]

cat("\nTravel time to nearest city (minutes):\n")
cat("At occurrences (n =", sum(!is.na(tt_occ)), "):\n")
print(round(summary(tt_occ)))
cat("\nAt background (n =", sum(!is.na(tt_bg)), "):\n")
print(round(summary(tt_bg)))

wt <- wilcox.test(tt_occ, tt_bg, alternative = "two.sided")
cat("\nWilcoxon rank-sum p-value:", format.pval(wt$p.value, digits = 3), "\n")
cat("Median occ:", round(median(tt_occ, na.rm = TRUE)),
    "| Median bg:", round(median(tt_bg, na.rm = TRUE)), "minutes\n")

# ----------------------------- Density plot ---------------------------------

plot_df <- bind_rows(
  data.frame(travel_time = tt_occ, group = "Occurrences"),
  data.frame(travel_time = tt_bg,  group = "Background")
)

p_acc <- ggplot(plot_df, aes(x = travel_time, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0, quantile(tt_bg, 0.99, na.rm = TRUE))) +
  scale_fill_manual(values = c("Background" = "grey60", "Occurrences" = "firebrick")) +
  labs(x = "Travel time to nearest city (minutes)",
       y = "Density", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave(file.path(DIR_FIGS, "accessibility_bias_diagnostic.png"), p_acc,
       width = 7, height = 5, dpi = 300, bg = "white")
cat("Saved accessibility_bias_diagnostic.png\n")

cat("\n11_accessibility_bias_diagnostic.R complete\n")