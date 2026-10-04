# ============================================================================
# 20_east_west_diagnostic.R
# Does the model fail at the Darfur presences, and is the failure robust?
# Compares Darfur presences with the rest by (1) suitability on the published
# map, against 08's thresholds (below them, a presence lies outside the
# binary maps), and (2) held-out predictions from 06 (each presence predicted
# by the model trained without its fold, against that model's p10): whether
# the model could predict the record without its neighbours. Then checks the
# gap on every alternative surface from 10, 12, 15 and 24.
#
#
# Inputs:  TRAIN_FILE, SUIT_FILE, DOMAIN_FILE, ADM1_FILE, arp_summary.csv (08),
#          heldout_presence_predictions.csv (06), DQ_SURFACES_FILE (10),
#          BIAS_SURFACES_FILE (12), VARIANT_SURFACES_FILE (15),
#          TT_COV_SURFACES_FILE (24)
# Outputs: outputs/tables/east_west_diagnostic.csv
#          outputs/tables/east_west_robustness.csv
#          outputs/figures/fig_east_west.png / .pdf
#          outputs/figures/presence_suitability_by_longitude.png
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(ggplot2)
})

# ------------------------------ Load inputs ---------------------------------

train  <- readRDS(TRAIN_FILE)
suit_r <- rast(SUIT_FILE)
arp_08 <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))
ho     <- read.csv(file.path(DIR_TABLES, "heldout_presence_predictions.csv"))
thr    <- setNames(arp_08$threshold, arp_08$metric)

occ <- train$occ_clean
xy  <- as.matrix(occ[, c("longitude", "latitude")])
occ$state    <- as.character(terra::extract(state_zones(suit_r), xy)[, 1])
occ$region   <- factor(ifelse(grepl("Darfur", occ$state), "Darfur", "Rest of Sudan"),
                       levels = names(pal_region))
occ$cell     <- cellFromXY(suit_r, xy)
occ$map_suit <- terra::extract(suit_r, xy)[, 1]
occ <- left_join(occ, select(ho, coordinate_id, heldout_pred = pred,
                             heldout_below = below_threshold), by = "coordinate_id")
stopifnot(
  "A presence has no state" = !anyNA(occ$state),
  "Held-out predictions don't match the presences" =
    nrow(ho) == nrow(occ) && !anyNA(occ$heldout_pred)
)
cat("Presences:", nrow(occ), "| Darfur:", sum(occ$region == "Darfur"),
    "| thresholds (08): p10", round(thr[["p10"]], 3), "maxSSS", round(thr[["maxsss"]], 3), "\n")

# ------------------------- By region and state ------------------------------
# Map suitability against 08's thresholds; held-out against each fold model's p10.

summ <- function(d) summarise(d,
  n = n(), cells = n_distinct(cell),
  map_median = median(map_suit), map_min = min(map_suit), map_max = max(map_suit),
  outside_p10_map = sum(map_suit < thr[["p10"]]),
  outside_maxsss_map = sum(map_suit < thr[["maxsss"]]),
  heldout_median = median(heldout_pred), heldout_below_p10 = sum(heldout_below),
  .groups = "drop")
rnd <- function(d) mutate(d, across(where(is.double), ~ round(., 3)))

cat("\nBy region:\n");  occ |> group_by(region) |> summ() |> rnd() |> print(width = Inf)
cat("\nBy state:\n")
occ |> group_by(state, region) |> summ() |> arrange(desc(map_median)) |> rnd() |>
  print(n = Inf, width = Inf)

# --------------------------- Across surfaces --------------------------------
# Every alternative surface: 10's data-quality refits, 12's corrected
# backgrounds, 15's covariate variants (primary settings and tuned) and 24's
# averaged surface. Threshold-free: the Darfur median against the other
# presences' 10th percentile.

layers_of <- function(file, drop = NULL, keep = NULL, prefix) {
  r <- rast(file)
  keep <- if (is.null(keep)) setdiff(names(r), drop) else keep
  r <- r[[keep]]; names(r) <- paste0(prefix, keep); r
}
surfs <- c(setNames(suit_r, "primary"),
           layers_of(DQ_SURFACES_FILE,      drop = "primary",       prefix = "dq_"),
           layers_of(BIAS_SURFACES_FILE,    keep = c("half", "matched"), prefix = "bias_"),
           layers_of(VARIANT_SURFACES_FILE, drop = "primary_fixed", prefix = "var_"),
           layers_of(TT_COV_SURFACES_FILE,  keep = "averaged",      prefix = "ttcov_"))

vals <- terra::extract(surfs, xy)
dar  <- occ$region == "Darfur"
rob  <- bind_rows(lapply(names(surfs), function(s) {
  v <- vals[[s]]
  data.frame(surface = s, darfur_median = median(v[dar]), darfur_max = max(v[dar]),
             rest_median = median(v[!dar]), rest_p10 = quantile(v[!dar], 0.1, names = FALSE))
})) |> mutate(gap_holds = darfur_median < rest_p10)

cat("\nDarfur gap holds on", sum(rob$gap_holds), "of", nrow(rob), "surfaces\n")
rob |> rnd() |> print(row.names = FALSE)

# ----------------------------- Source types ---------------------------------
# Exploratory: the grouping was chosen after the dissertation inspected the
# records below p10, and source type is confounded with region, so counts
# only, no test.

occ$source_type <- case_when(
  grepl("radio_dabanga|IFRC|reliefweb|sudan_tribune", occ$source, ignore.case = TRUE) ~ "Humanitarian/news",
  grepl("Pigott", occ$source) ~ "Pigott database",
  TRUE                        ~ "Peer-reviewed")
cat("\nSource type by region:\n"); print(table(occ$source_type, occ$region))
occ |> group_by(source_type) |>
  summarise(n = n(), heldout_below_p10 = sum(heldout_below),
            outside_p10_map = sum(map_suit < thr[["p10"]]),
            map_median = round(median(map_suit), 3), .groups = "drop") |> print()

# -------------------------------- Figures -----------------------------------

p_ew <- ggplot(occ, aes(region, map_suit, colour = region)) +
  geom_point(position = position_jitter(width = 0.15, seed = SEED), size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = thr[["p10"]],    linetype = "dotted", colour = "grey50") +
  geom_hline(yintercept = thr[["maxsss"]], linetype = "dashed", colour = "grey40") +
  annotate("text", x = 2.45, y = thr[["p10"]] + 0.03, label = "p10 threshold",
           size = 3, colour = "grey50", hjust = 1) +
  annotate("text", x = 2.45, y = thr[["maxsss"]] + 0.03, label = "maxSSS threshold",
           size = 3, colour = "grey40", hjust = 1) +
  scale_colour_manual(values = pal_region) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Predicted suitability at the presence") +
  theme(legend.position = "none")
for (ext in c("png", "pdf"))
  save_fig(file.path(DIR_FIGS, paste0("fig_east_west.", ext)), p_ew,
           width = FIG_WIDTH_HALF * 1.5, height = FIG_HEIGHT_PLOT)

p_lon <- ggplot(occ, aes(longitude, map_suit, colour = region, shape = heldout_below)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = thr[["p10"]], linetype = "dotted", colour = "grey50") +
  scale_colour_manual(values = pal_region, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), name = NULL,
                     labels = c(`FALSE` = "Held out: at or above p10",
                                `TRUE`  = "Held out: below p10")) +
  labs(x = "Longitude (\u00b0E)", y = "Predicted suitability at the presence") +
  theme(legend.position = "bottom")
save_fig(file.path(DIR_FIGS, "presence_suitability_by_longitude.png"), p_lon,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT)

# --------------------------------- Save -------------------------------------

write.csv(select(occ, coordinate_id, source, source_type, presence_type, year, longitude,
                 latitude, state, region, fold, map_suit, heldout_pred, heldout_below),
          file.path(DIR_TABLES, "east_west_diagnostic.csv"), row.names = FALSE)
write.csv(rob, file.path(DIR_TABLES, "east_west_robustness.csv"), row.names = FALSE)
cat("20_east_west_diagnostic.R complete\n")