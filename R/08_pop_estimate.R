# ============================================================================
# 08_pop_estimate.R
# Converts the continuous suitability surface to binary predictions at two
# thresholds (p10, maxSSS), downloads and aligns WorldPop 2025 population,
# and computes at-risk population estimates at national and state level.
#
# The risk-weighted estimate (suitability × population, no threshold) is
# the primary metric. Binary thresholds are reported for comparison and
# to demonstrate threshold sensitivity.
#
# Inputs:  outputs/models/maxent_final.rds
#          outputs/models/selected_tuning.rds
#          outputs/models/training_data.rds
#          outputs/surfaces/maxent_suitability.tif
# Outputs: outputs/surfaces/binary_p10.tif
#          outputs/surfaces/binary_maxsss.tif
#          outputs/surfaces/worldpop_2025_aligned.tif
#          outputs/tables/arp_summary.csv
#          outputs/tables/arp_by_state.csv
#          outputs/figures/threshold_sensitivity_curve.png
#          outputs/figures/threshold_sensitivity_curve.pdf
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(maxnet)
  library(ggplot2)
  library(httr)
  library(geodata)
})

# ------------------------------ Load inputs ---------------------------------

mod    <- readRDS(file.path(DIR_MODELS, "maxent_final.rds"))
tuning <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))
train  <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
suit_r <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))

cat("Tuning: fc =", tuning$fc, ", rm =", tuning$rm, "\n")
cat("Presences:", nrow(train$occ_env),
    "| Background:", nrow(train$bg_env), "\n")
cat("Suitability range:", round(global(suit_r, "min", na.rm = TRUE)[[1]], 4),
    "\u2013", round(global(suit_r, "max", na.rm = TRUE)[[1]], 4), "\n")

# ----------------------- Threshold selection --------------------------------

pred_occ <- predict(mod, train$occ_env, type = "cloglog")[, 1]
pred_bg  <- predict(mod, train$bg_env, type = "cloglog")[, 1]

# p10: 10th percentile of training presences
p10 <- unname(quantile(pred_occ, 0.10))

# maxSSS: maximum sensitivity + specificity
candidates <- sort(unique(c(pred_occ, pred_bg)))
sens <- sapply(candidates, function(t) mean(pred_occ >= t))
spec <- sapply(candidates, function(t) mean(pred_bg < t))
maxsss <- candidates[which.max(sens + spec)]

cat("p10 threshold:    ", round(p10, 4), "\n")
cat("maxSSS threshold: ", round(maxsss, 4), "\n")
cat("Presences >= p10:   ", sum(pred_occ >= p10), "/", length(pred_occ),
    "(", round(100 * mean(pred_occ >= p10), 1), "%)\n")
cat("Presences >= maxSSS:", sum(pred_occ >= maxsss), "/", length(pred_occ),
    "(", round(100 * mean(pred_occ >= maxsss), 1), "%)\n")

# ----------------------- Binary surfaces ------------------------------------

suit_p10    <- suit_r >= p10
suit_maxsss <- suit_r >= maxsss

cell_area_km2 <- cellSize(suit_r, unit = "km")
area_p10    <- global(suit_p10 * cell_area_km2, "sum", na.rm = TRUE)[[1]]
area_maxsss <- global(suit_maxsss * cell_area_km2, "sum", na.rm = TRUE)[[1]]
total_area  <- global(mask(cell_area_km2, !is.na(suit_r), maskvalues = 0),
                      "sum", na.rm = TRUE)[[1]]

cat("\np10: ", format(round(area_p10), big.mark = ","), "km\u00b2",
    "(", round(100 * area_p10 / total_area, 1), "% of Sudan)\n")
cat("maxSSS:", format(round(area_maxsss), big.mark = ","), "km\u00b2",
    "(", round(100 * area_maxsss / total_area, 1), "% of Sudan)\n")

# ----------------------- Population overlay ---------------------------------

pop_dir <- here::here("data", "raw", "population")
if (!dir.exists(pop_dir)) dir.create(pop_dir, recursive = TRUE)

pop_100m_path <- file.path(pop_dir, "sdn_pop_2025_100m_constrained.tif")

if (!file.exists(pop_100m_path)) {
  wp_url <- paste0(
    "https://data.worldpop.org/GIS/Population/Global_2015_2030/",
    "R2025A/2025/SDN/v1/100m/constrained/",
    "sdn_pop_2025_CN_100m_R2025A_v1.tif"
  )

  response <- GET(
    wp_url,
    user_agent("R - MSc dissertation, e.p.naymon@lse.ac.uk"),
    write_disk(pop_100m_path, overwrite = TRUE),
    progress()
  )
  Sys.sleep(1)
  stopifnot("WorldPop download failed" = status_code(response) == 200)
  cat("Downloaded:", basename(pop_100m_path), "\n")
} else {
  cat("Population raster already present\n")
}

pop_100m <- rast(pop_100m_path)

cat("Raw population total:",
    format(round(global(pop_100m, "sum", na.rm = TRUE)[[1]]), big.mark = ","), "\n")

# Aggregate 100m → 1km by summing, then align to suitability grid
agg_factor <- round(res(suit_r)[1] / res(pop_100m)[1])
pop_1km <- aggregate(pop_100m, fact = agg_factor, fun = "sum", na.rm = TRUE)
pop_aligned <- resample(pop_1km, suit_r, method = "sum")

cat("Aligned population total:",
    format(round(global(pop_aligned, "sum", na.rm = TRUE)[[1]]), big.mark = ","), "\n")

# -------------------- At-risk population estimates --------------------------

arp_p10      <- global(pop_aligned * suit_p10, "sum", na.rm = TRUE)[[1]]
arp_maxsss   <- global(pop_aligned * suit_maxsss, "sum", na.rm = TRUE)[[1]]
arp_weighted <- global(pop_aligned * suit_r, "sum", na.rm = TRUE)[[1]]
total_pop    <- global(pop_aligned, "sum", na.rm = TRUE)[[1]]

cat("\n--- At-Risk Population ---\n")
cat("Total Sudan population (WorldPop 2025):",
    format(round(total_pop), big.mark = ","), "\n\n")
cat("p10 ARP:          ", format(round(arp_p10), big.mark = ","),
    "(", round(100 * arp_p10 / total_pop, 1), "%)\n")
cat("maxSSS ARP:       ", format(round(arp_maxsss), big.mark = ","),
    "(", round(100 * arp_maxsss / total_pop, 1), "%)\n")
cat("Risk-weighted ARP:", format(round(arp_weighted), big.mark = ","),
    "(", round(100 * arp_weighted / total_pop, 1), "%)\n")

# -------------------- Threshold sensitivity curve ---------------------------

thresholds <- seq(0, 0.95, by = 0.01)

arp_by_thresh <- sapply(thresholds, function(t) {
  global(pop_aligned * (suit_r >= t), "sum", na.rm = TRUE)[[1]]
})

thresh_df <- data.frame(
  threshold = thresholds,
  arp       = arp_by_thresh,
  pct_pop   = 100 * arp_by_thresh / total_pop
)

p_thresh <- ggplot(thresh_df, aes(x = threshold, y = arp / 1e6)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = p10,    linetype = "dashed", colour = col_p10) +
  geom_vline(xintercept = maxsss, linetype = "dashed", colour = col_maxsss) +
  annotate("text", x = p10 + 0.02,    y = max(arp_by_thresh / 1e6) * 0.9,
           label = paste0("p10 (", round(p10, 3), ")"),
           hjust = 0, size = 3.2, colour = "steelblue") +
  annotate("text", x = maxsss + 0.02, y = max(arp_by_thresh / 1e6) * 0.8,
           label = paste0("maxSSS (", round(maxsss, 3), ")"),
           hjust = 0, size = 3.2, colour = "firebrick") +
  labs(title = "Threshold sensitivity of at-risk population estimate",
       x = "Suitability threshold",
       y = "At-risk population (millions)") +
  theme_dissertation(gridlines = "both") + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        plot.title = element_text(hjust = 0.5))

save_fig(file.path(DIR_FIGS, "threshold_sensitivity_curve.png"), p_thresh,
       width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT)
save_fig(file.path(DIR_FIGS, "threshold_sensitivity_curve.pdf"), p_thresh,
       width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT)
cat("Saved threshold_sensitivity_curve.png and threshold_sensitivity_curve.pdf\n")

# ----------------------- State-level breakdown ------------------------------

adm1 <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))

risk_weighted_r <- pop_aligned * suit_r

state_arp <- data.frame(
  state        = adm1$NAME_1,
  total_pop    = terra::extract(pop_aligned, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]],
  arp_p10      = terra::extract(pop_aligned * suit_p10, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]],
  arp_maxsss   = terra::extract(pop_aligned * suit_maxsss, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]],
  arp_weighted = terra::extract(risk_weighted_r, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]]
)

state_arp <- state_arp |>
  mutate(
    pct_p10      = round(100 * arp_p10 / total_pop, 1),
    pct_maxsss   = round(100 * arp_maxsss / total_pop, 1),
    pct_weighted = round(100 * arp_weighted / total_pop, 1)
  ) |>
  arrange(desc(arp_weighted))

cat("\nState-level ARP (sorted by risk-weighted):\n")
state_arp |>
  mutate(across(c(total_pop, arp_p10, arp_maxsss, arp_weighted),
                ~ format(round(.), big.mark = ","))) |>
  print(right = FALSE)

# --------------------------------- Save -------------------------------------

write.csv(state_arp, file.path(DIR_TABLES, "arp_by_state.csv"), row.names = FALSE)

writeRaster(suit_p10,    file.path(DIR_SURFACES, "binary_p10.tif"), overwrite = TRUE)
writeRaster(suit_maxsss, file.path(DIR_SURFACES, "binary_maxsss.tif"), overwrite = TRUE)
writeRaster(pop_aligned, file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"), overwrite = TRUE)

arp_summary <- data.frame(
  metric     = c("p10", "maxSSS", "risk_weighted"),
  threshold  = c(round(p10, 4), round(maxsss, 4), NA),
  arp        = c(round(arp_p10), round(arp_maxsss), round(arp_weighted)),
  pct_of_pop = c(round(100 * arp_p10 / total_pop, 1),
                 round(100 * arp_maxsss / total_pop, 1),
                 round(100 * arp_weighted / total_pop, 1)),
  total_pop  = round(total_pop)
)
write.csv(arp_summary, file.path(DIR_TABLES, "arp_summary.csv"), row.names = FALSE)

cat("\nSaved:\n")
cat("  ", file.path(DIR_TABLES, "arp_by_state.csv"), "\n")
cat("  ", file.path(DIR_TABLES, "arp_summary.csv"), "\n")
cat("  ", file.path(DIR_SURFACES, "binary_p10.tif"), "\n")
cat("  ", file.path(DIR_SURFACES, "binary_maxsss.tif"), "\n")
cat("  ", file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"), "\n")

cat("08_pop_estimate.R complete\n")