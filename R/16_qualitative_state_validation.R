# ============================================================================
# 16_qualitative_state_validation.R
# Compares model-predicted suitability against independently documented VL
# endemic status per state. Tests whether the model recovers the broad
# geographic hierarchy of VL burden without any state-level information
# entering the training process.
#
# Inputs:  outputs/surfaces/maxent_suitability.tif
#          data/raw/gadm/ (state boundaries)
# Outputs: outputs/tables/state_validation.csv
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(geodata)
})

# -------------------- Literature-based classification -----------------------

lit <- tribble(
  ~state,            ~status,           ~evidence,
  "Al Qadarif",      "Core endemic",    "Principal endemic focus",
  "Blue Nile",       "Core endemic",    "Historical endemic belt",
  "Sennar",          "Core endemic",    "Historical endemic belt",
  "Kassala",         "Core endemic",    "Historical endemic belt",
  "White Nile",      "Reported",        "Evidence of shift to endemic status",
  "South Kurdufan",  "Reported",        "Scattered foci",
  "North Kurdufan",  "Reported",        "Scattered foci",
  "West Kurdufan",   "Reported",        "Cited under regional Kordofan grouping",
  "North Darfur",    "Reported",        "Secondary VL focus",
  "South Darfur",    "Reported",        "Scattered foci",
  "Central Darfur",  "Reported",        "Created 2012 from West and South Darfur; both documented",
  "East Darfur",     "Reported",        "Darfur foci",
  "West Darfur",     "Reported",        "Scattered historically-reported foci",
  "Al Jazirah",      "Reported",        "Historical outbreaks; current scattered endemicity",
  "Red Sea",         "Reported",        "New foci added to control strategy in 2017",
  "Khartoum",        "Uncertain",       "Sporadic cases in 1960s; potentially imported currently",
  "River Nile",      "No documented",   "No published VL transmission",
  "Northern",        "No documented",   "No published VL transmission"
)

cat("States classified:", nrow(lit), "\n")

# --------------------- Extract suitability per state ------------------------

suit <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))
adm1 <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))
states <- st_as_sf(adm1)

state_stats <- data.frame()

for (i in seq_len(nrow(states))) {
  state_poly <- vect(states[i, ])
  vals <- terra::extract(suit, state_poly, ID = FALSE)[[1]]
  vals <- vals[!is.na(vals)]

  state_stats <- rbind(state_stats, data.frame(
    state       = states$NAME_1[i],
    mean_suit   = round(mean(vals), 4),
    median_suit = round(median(vals), 4),
    max_suit    = round(max(vals), 4),
    n_cells     = length(vals)
  ))
}

state_stats <- state_stats |> arrange(desc(mean_suit))

cat("\nState-level suitability summary:\n")
cat(sprintf("%-18s %10s %10s %10s\n", "State", "Mean", "Median", "Max"))
for (i in seq_len(nrow(state_stats))) {
  cat(sprintf("%-18s %10.4f %10.4f %10.4f\n",
              state_stats$state[i],
              state_stats$mean_suit[i],
              state_stats$median_suit[i],
              state_stats$max_suit[i]))
}

# ----------------------- Join and summarise ---------------------------------

validation <- lit |>
  left_join(state_stats, by = "state") |>
  arrange(desc(mean_suit))

if (any(is.na(validation$mean_suit))) {
  cat("WARNING — unmatched states:\n")
  print(validation$state[is.na(validation$mean_suit)])
}

cat("\n--- Mean suitability by endemic status ---\n")
validation |>
  group_by(status) |>
  summarise(
    n          = n(),
    mean_suit  = round(mean(mean_suit), 3),
    range      = paste0(round(min(mean_suit), 3), "\u2013", round(max(mean_suit), 3)),
    .groups    = "drop"
  ) |>
  arrange(desc(mean_suit)) |>
  print()

cat("\n--- Full table ---\n")
validation |>
  select(state, status, mean_suit, median_suit, max_suit) |>
  print(n = 18)

# --------------------------------- Save -------------------------------------

write.csv(validation |> select(state, status, evidence, mean_suit, median_suit, max_suit, n_cells),
          file.path(DIR_TABLES, "state_validation.csv"), row.names = FALSE)

cat("\nSaved:", file.path(DIR_TABLES, "state_validation.csv"), "\n")
cat("16_qualitative_state_validation.R complete\n")