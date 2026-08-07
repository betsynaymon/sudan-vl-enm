# ============================================================================
# 21_unverified_low_prediction.R
# ----------------------------------------------------------------------------
# Quantifies the population living where the model predicts LOW suitability AND
# where geographic accessibility is so poor that this low prediction has not
# been -- and cannot easily be -- checked against ground observation.
#
# Definition (per cell, within the >=150 mm ecological domain):
#   below-threshold : suitability <  {p10, maxSSS}    (classified not-at-risk)
#   unreachable     : travel time  >  cutoff          (beyond observation reach)
#   invisible pop   : WorldPop 2025 summed where both hold
#
# Extent = ecological mask (>=150 mm). 
#
# "Unreachable" is anchored to the evidence base: the travel-time distribution
# at the 98 occurrence records defines how far observation actually reached.
# Cells beyond its upper tail (75th / 90th / 95th percentile) lie beyond that
# frontier. The (below-p10, tt > occ-95th) corner is the conservative floor
# ("at least N"); lower cutoffs and maxSSS give larger, still-defensible counts.
#
# Inputs:  outputs/models/maxent_final.rds
#          outputs/models/training_data.rds
#          outputs/surfaces/maxent_suitability.tif
#          outputs/surfaces/worldpop_2025_aligned.tif      (from 08)
#          data/raw/ecological_mask_150mm.tif
#          data/raw/weiss_travel_time.tif
# Outputs: outputs/tables/unverified_low_prediction_summary.csv
#          outputs/tables/unverified_low_prediction_by_state.csv
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(maxnet)
  library(geodata)
})

# ------------------------------ Load inputs ---------------------------------

mod    <- readRDS(file.path(DIR_MODELS, "maxent_final.rds"))
train  <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
suit_r <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))
pop_r  <- rast(file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"))
mask_r <- rast(file.path(DIR_COVARIATES, "ecological_mask_150mm.tif"))
tt_raw <- rast(here::here("data", "raw", "weiss_travel_time.tif"))

# 1 inside the >=150 mm domain, NA outside (matches masking convention in 12)
mask_bin <- subst(mask_r, 0, NA)

# ----------------------- Thresholds (identical to 08) -----------------------

pred_occ <- predict(mod, train$occ_env, type = "cloglog")[, 1]
pred_bg  <- predict(mod, train$bg_env,  type = "cloglog")[, 1]

p10 <- unname(quantile(pred_occ, 0.10))

candidates <- sort(unique(c(pred_occ, pred_bg)))
sens <- sapply(candidates, function(t) mean(pred_occ >= t))
spec <- sapply(candidates, function(t) mean(pred_bg  <  t))
maxsss <- candidates[which.max(sens + spec)]

cat("Thresholds (verify against 08 -- expect ~0.1541 / ~0.3857):\n")
cat("  p10   :", round(p10, 4),    "\n")
cat("  maxSSS:", round(maxsss, 4), "\n")

# ----------------- Reachability envelope of the evidence base ---------------
# Travel time at the 98 occurrences = how far observation actually reached.

occ_pts <- vect(train$occ_clean, geom = c("longitude", "latitude"), crs = "EPSG:4326")
bg_pts  <- vect(train$bg_clean,  geom = c("longitude", "latitude"), crs = "EPSG:4326")

tt_occ <- terra::extract(tt_raw, occ_pts)[, 2]
tt_bg  <- terra::extract(tt_raw, bg_pts)[, 2]

occ_q  <- quantile(tt_occ, c(0.75, 0.90, 0.95), na.rm = TRUE)
bg_med <- median(tt_bg, na.rm = TRUE)

cat("\nTravel time at occurrences (min; verify median ~117 against 11):\n")
cat("  median:", round(median(tt_occ, na.rm = TRUE)),
    "| 75th:", round(occ_q[1]),
    "| 90th:", round(occ_q[2]),
    "| 95th:", round(occ_q[3]), "\n")
cat("  background median (verify ~324):", round(bg_med), "\n")

# Candidate "unreachable" cutoffs, in minutes
cutoffs <- c(
  occ_75th  = unname(occ_q[1]),
  occ_90th  = unname(occ_q[2]),
  occ_95th  = unname(occ_q[3]),
  bg_median = bg_med
)

# ----------------------- Align travel time to grid --------------------------

tt_aligned <- resample(tt_raw, suit_r, method = "bilinear")   # matches 11 / 12

# ------------------------ Below-threshold masks -----------------------------
# "Below" = classified not-at-risk. 08 defines at-risk as suitability >= t.

below_p10    <- mask(suit_r < p10,    mask_bin)
below_maxsss <- mask(suit_r < maxsss, mask_bin)

# ----------------------------- Denominators ---------------------------------

total_pop     <- global(pop_r, "sum", na.rm = TRUE)[[1]]
pop_in_mask   <- global(mask(pop_r, mask_bin), "sum", na.rm = TRUE)[[1]]
pop_below_p10 <- global(pop_r * below_p10,    "sum", na.rm = TRUE)[[1]]
pop_below_mss <- global(pop_r * below_maxsss, "sum", na.rm = TRUE)[[1]]

# ------------------- Invisible population across cutoffs ---------------------

invisible_pop <- function(below_layer, cutoff_min) {
  remote <- mask(tt_aligned > cutoff_min, mask_bin)
  global(pop_r * (below_layer & remote), "sum", na.rm = TRUE)[[1]]
}

results <- do.call(rbind, lapply(names(cutoffs), function(nm) {
  cut <- cutoffs[[nm]]
  data.frame(
    cutoff_name             = nm,
    cutoff_min              = round(cut),
    pop_below_p10_remote    = invisible_pop(below_p10,    cut),
    pop_below_maxsss_remote = invisible_pop(below_maxsss, cut)
  )
}))

results$pct_total_p10    <- round(100 * results$pop_below_p10_remote    / total_pop, 1)
results$pct_total_maxsss <- round(100 * results$pop_below_maxsss_remote / total_pop, 1)

cat("\n=== Population below threshold AND beyond observation reach ===\n")
cat("(within the >=150 mm ecological domain)\n\n")
cat("Denominators:\n")
cat("  Total Sudan (WorldPop 2025):", format(round(total_pop),     big.mark = ","), "\n")
cat("  Within ecological mask:     ", format(round(pop_in_mask),   big.mark = ","), "\n")
cat("  Below p10 (in mask):        ", format(round(pop_below_p10), big.mark = ","), "\n")
cat("  Below maxSSS (in mask):     ", format(round(pop_below_mss), big.mark = ","), "\n\n")

print(results |>
        mutate(across(c(pop_below_p10_remote, pop_below_maxsss_remote),
                      ~ format(round(.), big.mark = ","))),
      right = FALSE)

floor_val <- results$pop_below_p10_remote[results$cutoff_name == "occ_95th"]
cat("\nConservative floor (below-p10, tt > occ 95th pct):",
    format(round(floor_val), big.mark = ","), "\n")

# ------------- Transparency: full-Sudan (unmasked), primary cutoff ----------

primary_cut <- cutoffs[["occ_90th"]]
remote_full <- tt_aligned > primary_cut
u_p10 <- global(pop_r * ((suit_r < p10)    & remote_full), "sum", na.rm = TRUE)[[1]]
u_mss <- global(pop_r * ((suit_r < maxsss) & remote_full), "sum", na.rm = TRUE)[[1]]

cat("\n[Transparency] Full-Sudan, unmasked, cutoff = occ 90th pct (",
    round(primary_cut), "min):\n")
cat("  below-p10 & remote:   ", format(round(u_p10), big.mark = ","), "\n")
cat("  below-maxSSS & remote:", format(round(u_mss), big.mark = ","), "\n")

# ---------------- State-level breakdown at primary cutoff -------------------
# Primary cutoff = occurrence 90th percentile.

remote_p       <- mask(tt_aligned > primary_cut, mask_bin)
invis_p10_r    <- (below_p10    & remote_p) * pop_r
invis_maxsss_r <- (below_maxsss & remote_p) * pop_r

adm1 <- geodata::gadm(country = "SDN", level = 1, path = here::here("data", "raw"))

state_tbl <- data.frame(
  state                   = adm1$NAME_1,
  pop_below_p10_remote    = terra::extract(invis_p10_r,    adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]],
  pop_below_maxsss_remote = terra::extract(invis_maxsss_r, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]]
) |>
  arrange(desc(pop_below_maxsss_remote))

cat("\nState-level invisible population (cutoff = occ 90th pct =",
    round(primary_cut), "min):\n")
state_tbl |>
  mutate(across(c(pop_below_p10_remote, pop_below_maxsss_remote),
                ~ format(round(.), big.mark = ","))) |>
  print(right = FALSE)

# --------------------------------- Save -------------------------------------

summary_out <- results
summary_out$total_pop        <- round(total_pop)
summary_out$pop_in_mask      <- round(pop_in_mask)
summary_out$pop_below_p10    <- round(pop_below_p10)
summary_out$pop_below_maxsss <- round(pop_below_mss)

write.csv(summary_out, file.path(DIR_TABLES, "unverified_low_prediction_summary.csv"),
          row.names = FALSE)
write.csv(state_tbl,   file.path(DIR_TABLES, "unverified_low_prediction_by_state.csv"),
          row.names = FALSE)

cat("\nSaved:\n")
cat("  ", file.path(DIR_TABLES, "unverified_low_prediction_summary.csv"), "\n")
cat("  ", file.path(DIR_TABLES, "unverified_low_prediction_by_state.csv"), "\n")
cat("21_unverified_low_prediction.R complete\n")