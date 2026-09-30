# ============================================================================
# 02_spatial_thinning.R
# Diagnoses spatial clustering in occurrence records, draws a retention curve
# across candidate distances (THIN_TEST_KM), and thins within each year at
# THIN_KM. Repeat observations at one location in different years (with
# different annual covariates) are retained.
#
# Inputs:  OCC_RAW_FILE
# Outputs: OCC_FILE
#          outputs/figures/nn_distances.png
#          outputs/figures/retention_curve.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(spThin)
})

# ----------------------------- Load and filter ------------------------------
# Drop referral records (REFERRAL_IDS): a national referral centre's location
# carries no information on where infection occurred. Khartoum case records
# (KHARTOUM_CASE_IDS) and ID 121 (vector-positive trap site) are retained;
# script 10 refits without the Khartoum case records.

occ <- read.csv(OCC_RAW_FILE)

stopifnot(
  "A referral ID is not in the data"      = all(REFERRAL_IDS %in% occ$coordinate_id),
  "A Khartoum case ID is not in the data" = all(KHARTOUM_CASE_IDS %in% occ$coordinate_id)
)

occ_fit <- occ %>% filter(!coordinate_id %in% REFERRAL_IDS)

cat("Total records:", nrow(occ),
    "| After referral drop:", nrow(occ_fit), "\n")
stopifnot(
  "No records remain after referral drop" = nrow(occ_fit) > 0
)

# ---------------------- Nearest-neighbor diagnostic ------------------------
# Clustering diagnosed on unique locations (collapsing year-duplicates) since
# the concern is spatial autocorrelation in geographic space. spThin uses
# great-circle distances.

occ_unique <- occ_fit %>% distinct(longitude, latitude, .keep_all = TRUE)
cat("Unique fitting locations:", nrow(occ_unique), "\n")

pts <- st_as_sf(occ_unique, coords = c("longitude", "latitude"), crs = 4326)

dmat <- st_distance(pts)
diag(dmat) <- NA
nn_dist <- apply(dmat, 1, min, na.rm = TRUE)
nn_km <- as.numeric(nn_dist) / 1000

cat("Nearest-neighbor summary:\n")
print(summary(nn_km))
within <- sapply(THIN_TEST_KM, function(d) sum(nn_km < d))
names(within) <- paste0("<", THIN_TEST_KM, "km")
cat("Locations with a neighbour within each test distance:\n")
print(within)

p_nn <- ggplot(data.frame(nn_km), aes(nn_km)) +
  geom_histogram(binwidth = 5, boundary = 0) +
  labs(x = "Distance to nearest neighbour (km)",
       y = "Number of locations",
       title = paste("Nearest-neighbour distances:",
                     nrow(occ_unique), "unique fitting locations")) +
  theme_minimal()

ggsave(file.path(DIR_FIGS, "nn_distances.png"), p_nn,
       width = 8, height = 5, dpi = 300)
cat("Saved nn_distances.png\n")

# ----------------------------- Retention curve ------------------------------
# Tests multiple thinning thresholds. The "knee" — where retained count drops
# sharply — marks the transition from dropping near-duplicates to losing real
# spatial coverage.

set.seed(SEED)

retention <- lapply(THIN_TEST_KM, function(d) {
  res <- thin(
    loc.data = occ_unique %>% transmute(LAT = latitude, LONG = longitude, SPEC = "VL"),
    lat.col = "LAT", long.col = "LONG", spec.col = "SPEC",
    thin.par = d, reps = THIN_REPS,
    locs.thinned.list.return = TRUE,
    write.files = FALSE, write.log.file = FALSE, verbose = FALSE
  )
  counts <- sapply(res, nrow)
  data.frame(
    thin_km   = d,
    max_kept  = max(counts),
    min_kept  = min(counts),
    mean_kept = round(mean(counts), 1)
  )
}) %>% bind_rows()

retention$pct_kept <- round(100 * retention$max_kept / nrow(occ_unique), 1)
print(retention)

p_ret <- ggplot(retention, aes(thin_km, max_kept)) +
  geom_line() + geom_point(size = 2) +
  geom_text(aes(label = paste0(max_kept, " (", pct_kept, "%)")),
            vjust = -0.8, size = 3) +
  scale_x_continuous(breaks = THIN_TEST_KM) +
  labs(x = "Thinning distance (km)",
       y = paste0("Locations retained (max over ", THIN_REPS, " reps)"),
       title = "Retention curve across thinning distances") +
  theme_minimal()

ggsave(file.path(DIR_FIGS, "retention_curve.png"), p_ret,
       width = 8, height = 5, dpi = 300)
cat("Saved retention_curve.png\n")

# -------------------------- Within-year thinning ----------------------------
# A location observed in 2005 and again in 2011 represents two ecologically
# distinct observations (different annual covariates), so both are retained.
# Only within-year spatial redundancy is removed.
years <- unique(occ_fit$year)

thinned_by_year <- lapply(years, function(y) {
  set.seed(SEED + y)   # per-year seed: one year's records can't change another's thinning
  yr_data <- occ_fit %>% filter(year == y)

  # spThin needs >1 point; single-point years pass through
  if (nrow(yr_data) < 2) return(yr_data)

  thin_input <- yr_data %>%
    transmute(LAT = latitude, LONG = longitude, SPEC = "VL")

  res <- thin(
    loc.data = thin_input,
    lat.col = "LAT", long.col = "LONG", spec.col = "SPEC",
    thin.par = THIN_KM, reps = THIN_REPS,
    locs.thinned.list.return = TRUE,
    write.files = FALSE, write.log.file = FALSE, verbose = FALSE
  )

  best <- res[[which.max(sapply(res, nrow))]]

  # Match back to full rows by rounded coordinates
  best_keys <- paste(round(best$Longitude, 5), round(best$Latitude, 5))
  yr_data %>%
    filter(paste(round(longitude, 5), round(latitude, 5)) %in% best_keys)
})

occ_thinned <- bind_rows(thinned_by_year)
# Rows are matched back by rounded coordinates, so records sharing a location
# in the same year would all survive. This also guarantees 06's cell x year
# dedup drops no presences.
stopifnot("Duplicate location-year after thinning" = !any(duplicated(
  data.frame(round(occ_thinned$longitude, 5), round(occ_thinned$latitude, 5),
             occ_thinned$year))))

cat("Before:", nrow(occ_fit), "rows | After within-year thin at",
    THIN_KM, "km:", nrow(occ_thinned), "rows\n")

# Per-year summary
occ_thinned %>%
  count(year, name = "kept") %>%
  left_join(occ_fit %>% count(year, name = "original"), by = "year") %>%
  mutate(dropped = original - kept) %>%
  print()

cat("Khartoum case records retained:",
    sum(occ_thinned$coordinate_id %in% KHARTOUM_CASE_IDS), "of",
    length(KHARTOUM_CASE_IDS), "\n")
# ----------------------------------- Save -----------------------------------

occ_thinned <- occ_thinned %>% mutate(year = as.integer(year))

write.csv(occ_thinned, OCC_FILE, row.names = FALSE)

cat("Saved", nrow(occ_thinned), "thinned occurrence records\n")
cat("Unique locations:", occ_thinned %>% distinct(longitude, latitude) %>% nrow(), "\n")
cat("Year range:", range(occ_thinned$year), "\n")
cat("\nYear distribution:\n")
occ_thinned %>% count(year, name = "n") %>% print()

cat("02_spatial_thinning.R complete\n")