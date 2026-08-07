# ============================================================================
# 22_state_accessibility.R
# Calculates the percentage of each state's population living beyond the
# background median travel time (324 minutes) to the nearest city.
# Produces a summary table for inclusion in the state-level comparison table.
#
# Inputs:  data/raw/weiss_travel_time.tif
#          data/raw/population/sdn_pop_2025_100m_constrained.tif
#          GADM admin1 boundaries
#          outputs/surfaces/maxent_suitability.tif (as alignment template)
# Outputs: outputs/tables/state_accessibility.csv
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(geodata)
})

# ----- Threshold -----
# Background median travel time from 11_accessibility_bias_diagnostic.R
TT_THRESHOLD <- 324  # minutes

# ----- Load inputs -----

# Travel-time surface (Weiss et al. 2018)
tt <- rast(here::here("data", "raw", "weiss_travel_time.tif"))

# Population raster 
pop <- rast(here::here("data", "raw", "population", "sdn_pop_2025_100m_constrained.tif"))

# Suitability surface as alignment template
suit <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))

# GADM admin1 boundaries
adm1 <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))

cat("States:", length(adm1), "\n")
cat("Population total:", round(global(pop, "sum", na.rm = TRUE)[[1]]), "\n")

# ----- Crop and align travel time to population grid -----

tt <- rast(here::here("data", "raw", "weiss_travel_time.tif"))
tt <- crop(tt, pop)  # trim to same extent first
tt_aligned <- resample(tt, pop, method = "bilinear")

# Verify alignment
cat("Pop dims:", ncol(pop), "x", nrow(pop), "\n")
cat("TT dims:", ncol(tt_aligned), "x", nrow(tt_aligned), "\n")

# ----- Zonal calculation per state -----

# Total population per state
state_total <- terra::extract(pop, adm1, fun = sum, na.rm = TRUE)

# Population beyond threshold: multiply pop by binary indicator
beyond <- ifel(tt_aligned > TT_THRESHOLD, 1, 0)
pop_beyond <- pop * beyond
state_beyond <- terra::extract(pop_beyond, adm1, fun = sum, na.rm = TRUE)

# Population with NA travel time (no data — typically remote areas)
tt_na <- ifel(is.na(tt_aligned) & !is.na(pop), 1, 0)
pop_na <- pop * tt_na
state_na <- terra::extract(pop_na, adm1, fun = sum, na.rm = TRUE)

# ----- Build table -----

results <- data.frame(
  state       = adm1$NAME_1,
  total_pop   = round(state_total[, 2]),
  pop_beyond  = round(state_beyond[, 2]),
  pop_na_tt   = round(state_na[, 2])
) |>
  mutate(
    pct_beyond = round(pop_beyond / total_pop * 100, 1),
    pct_na     = round(pop_na_tt / total_pop * 100, 1)
  ) |>
  arrange(desc(pct_beyond))

cat("\n--- Population beyond", TT_THRESHOLD, "minutes per state ---\n")
print(results, row.names = FALSE)

cat("\n--- Summary ---\n")
cat("Total pop beyond threshold:", sum(results$pop_beyond), "\n")
cat("Total pop with NA travel time:", sum(results$pop_na_tt), "\n")
cat("National pct beyond:", 
    round(sum(results$pop_beyond) / sum(results$total_pop) * 100, 1), "%\n")

# ----- Save -----

write.csv(results, file.path(DIR_TABLES, "state_accessibility.csv"), 
          row.names = FALSE)
cat("\nSaved state_accessibility.csv\n")

cat("\n22_state_accessibility.R complete\n")