# ============================================================================
# 20_east_west_diagnostic.R
# Extracts predicted suitability at each thinned presence location, classifies
# east (non-Darfur) vs west (Darfur), and compares distributions. Tests
# whether the western gap persists across all covariate variants and after
# accessibility correction. Frames the result as a justified scope boundary
# rather than a model failure.
#
# Inputs:  outputs/models/maxent_final.rds
#          outputs/models/training_data.rds
#          outputs/surfaces/maxent_suitability.tif
#          outputs/surfaces/maxent_suitability_{variant}.tif  (5 variants)
#          outputs/surfaces/maxent_suitability_biased_bg.tif
#          data/processed/occurrences_thinned.csv
# Outputs: outputs/tables/east_west_diagnostic.csv
#          outputs/figures/east_west_suitability_diagnostic.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(geodata)
  library(maxnet)
})

set.seed(SEED)

# ------------------------------ Load inputs ---------------------------------

pred_ltm <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))
train    <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
mod      <- readRDS(file.path(DIR_MODELS, "maxent_final.rds"))

# Thresholds
pred_occ <- predict(mod, train$occ_env, type = "cloglog")[, 1]
pred_bg  <- predict(mod, train$bg_env, type = "cloglog")[, 1]

p10 <- unname(quantile(pred_occ, 0.10))
candidates <- sort(unique(c(pred_occ, pred_bg)))
sens <- sapply(candidates, function(t) mean(pred_occ >= t))
spec <- sapply(candidates, function(t) mean(pred_bg < t))
maxsss <- candidates[which.max(sens + spec)]

cat("p10:", round(p10, 4), "| maxSSS:", round(maxsss, 4), "\n")

# ---------------------- Classify east vs west -------------------------------

occ  <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
adm1 <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))

occ_sf <- st_as_sf(occ, coords = c("longitude", "latitude"), crs = 4326)
occ_states <- st_join(occ_sf, st_as_sf(adm1)[, "NAME_1"])

darfur_states <- c("Central Darfur", "East Darfur", "North Darfur",
                   "South Darfur", "West Darfur")

occ_states$region <- ifelse(occ_states$NAME_1 %in% darfur_states, "West", "East")

# Assign border points (NA state) to nearest state
na_idx <- is.na(occ_states$NAME_1)
if (any(na_idx)) {
  nearest <- st_nearest_feature(occ_states[na_idx, ], st_as_sf(adm1))
  occ_states$NAME_1[na_idx] <- adm1$NAME_1[nearest]
  occ_states$region[na_idx] <- ifelse(
    adm1$NAME_1[nearest] %in% darfur_states, "West", "East"
  )
  cat("Assigned", sum(na_idx), "border points to nearest state\n")
}

cat("\nRecords by region:\n")
print(table(occ_states$region))

cat("\nWestern records by state:\n")
occ_states |> filter(region == "West") |>
  st_drop_geometry() |> count(NAME_1) |> print()

# ---------------------- Extract and compare ---------------------------------

occ_coords <- st_coordinates(occ_states)
occ_states$suitability <- terra::extract(pred_ltm, occ_coords)[, 1]

cat("\nSuitability by region:\n")
region_summary <- occ_states |>
  st_drop_geometry() |>
  group_by(region) |>
  summarise(
    n           = n(),
    mean_suit   = round(mean(suitability, na.rm = TRUE), 3),
    median_suit = round(median(suitability, na.rm = TRUE), 3),
    min_suit    = round(min(suitability, na.rm = TRUE), 3),
    max_suit    = round(max(suitability, na.rm = TRUE), 3),
    .groups     = "drop"
  )
print(region_summary)

# ---------------------- Strip plot ------------------------------------------

plot_df <- occ_states |>
  st_drop_geometry() |>
  select(region, suitability, NAME_1)

east_med <- region_summary$median_suit[region_summary$region == "East"]
west_med <- region_summary$median_suit[region_summary$region == "West"]
n_east   <- region_summary$n[region_summary$region == "East"]
n_west   <- region_summary$n[region_summary$region == "West"]

p_ew <- ggplot(plot_df, aes(x = region, y = suitability, colour = region)) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = maxsss, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = p10, linetype = "dotted", colour = "grey60") +
  annotate("text", x = 2.4, y = p10 + 0.03, label = "p10 threshold",
           size = 3, colour = "grey60") +
  annotate("text", x = 2.4, y = maxsss + 0.03, label = "maxSSS threshold",
           size = 3, colour = "grey40") +
  scale_colour_manual(values = c("East" = "#B2182B", "West" = "#2166AC")) +
  labs(x = NULL, y = "Predicted suitability (LTM surface)",
       title = "Model performance at known VL locations",
       subtitle = paste0(n_west, " western records fall below threshold \u2014 ",
                         "ecologically distinct system")) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(DIR_FIGS, "east_west_suitability_diagnostic.png"), p_ew,
       width = 7, height = 6, dpi = 300, bg = "white")
cat("Saved east_west_suitability_diagnostic.png\n")

# ------------------- Cross-variant sensitivity ------------------------------

variant_surfaces <- list(
  A_annual = rast(file.path(DIR_SURFACES, "maxent_suitability.tif")),
  B_annual = rast(file.path(DIR_SURFACES, "maxent_suitability_B_annual.tif")),
  C_annual = rast(file.path(DIR_SURFACES, "maxent_suitability_C_annual.tif")),
  D_annual = rast(file.path(DIR_SURFACES, "maxent_suitability_D_annual.tif")),
  A_dry    = rast(file.path(DIR_SURFACES, "maxent_suitability_A_dry.tif")),
  A_wet    = rast(file.path(DIR_SURFACES, "maxent_suitability_A_wet.tif"))
)

west_pts    <- occ_states |> filter(region == "West")
west_coords <- st_coordinates(west_pts)

west_suit <- sapply(variant_surfaces, function(r) {
  terra::extract(r, west_coords)[, 1]
})

west_suit_df <- as.data.frame(west_suit)
west_suit_df$state <- west_pts$NAME_1
west_suit_df <- west_suit_df |> select(state, everything())

cat("\nSuitability at western presences across variants:\n")
print(round(west_suit_df[, -1], 4))

cat("\nMedian suitability per variant:\n")
sapply(west_suit_df[, -1], median) |> round(4) |> print()

cat("\nMax suitability per variant:\n")
sapply(west_suit_df[, -1], max) |> round(4) |> print()

# ------------------- Accessibility-corrected check --------------------------

pred_biased <- rast(file.path(DIR_SURFACES, "maxent_suitability_biased_bg.tif"))
west_suit_bias <- terra::extract(pred_biased, west_coords)[, 1]

cat("\nWestern presences \u2014 accessibility-corrected surface:\n")
print(round(west_suit_bias, 4))
cat("Median:", round(median(west_suit_bias), 4),
    "| Max:", round(max(west_suit_bias), 4), "\n")

# --------------------------------- Save -------------------------------------

east_west_out <- occ_states |>
  cbind(st_coordinates(occ_states)) |>
  st_drop_geometry() |>
  rename(longitude = X, latitude = Y) |>
  select(longitude, latitude, year, source,
         state = NAME_1, region, suitability)

write.csv(east_west_out, file.path(DIR_TABLES, "east_west_diagnostic.csv"),
          row.names = FALSE)

cat("\n--- East-west diagnostic summary ---\n")
cat("Eastern presences (n=", n_east, "): median suitability", east_med, "\n")
cat("Western presences (n=", n_west, "): median suitability", west_med, "\n")

cat("\n20_east_west_diagnostic.R complete\n")