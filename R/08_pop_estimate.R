# ============================================================================
# 08_pop_estimate.R
# Population at risk from the selected model's surface (07): risk-weighted
# (population x suitability; the primary estimate, an index rather than a
# count) and binary at p10 and maxSSS, nationally and by state, with
# diagnostics on where the estimate rests and how much it depends on the
# plateau choice.
#
# Inputs:  MODEL_FILE, TUNING_FILE, TRAIN_FILE, SUIT_FILE, MESS_FILE,
#          MESS_VARS_FILE, CANDIDATES_FILE, POP_FILE, DOMAIN_FILE, ADM1_FILE,
#          SENS_MASK_FILE
# Outputs: POP_ALIGNED_FILE, outputs/surfaces/binary_p10.tif, binary_maxsss.tif
#          outputs/tables/arp_summary.csv, arp_by_state.csv,
#          arp_plateau_candidates.csv, mess_limiting_variable_by_state.csv
#          outputs/figures/threshold_sensitivity_curve.png / .pdf
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(maxnet); library(ggplot2); library(httr)
})

mod    <- readRDS(MODEL_FILE)
tuning <- readRDS(TUNING_FILE)
train  <- readRDS(TRAIN_FILE)
suit_r <- rast(SUIT_FILE)
cat("Model:", tuning$fc, "rm =", tuning$rm, "| rule:", tuning$rule, "\n")

# ----------------------------- Thresholds -----------------------------------
# From training predictions (year-matched values), as in 06. Two diagnostics:
# p10 from map values at the presence sites (the map uses long-term means),
# and maxSSS against within-belt background only (half the full background is
# desert, which makes specificity cheap).

pred_occ <- as.numeric(predict(mod, train$occ_env, type = "cloglog"))
pred_bg  <- as.numeric(predict(mod, train$bg_env,  type = "cloglog"))
thr <- thresholds_from(pred_occ, pred_bg)

occ_xy  <- as.matrix(train$occ_clean[, c("longitude", "latitude")])
map_p10 <- quantile(terra::extract(suit_r, occ_xy)[, 1], OMISSION_Q, names = FALSE)

bg_wet  <- terra::extract(rast(SENS_MASK_FILE),
                          as.matrix(train$bg_clean[, c("longitude", "latitude")]))[, 1] %in% 1
thr_wet <- thresholds_from(pred_occ, pred_bg[bg_wet])

cat("p10:   ", round(thr[["p10"]], 4), "| from map values at presences:", round(map_p10, 4), "\n")
cat("maxSSS:", round(thr[["maxsss"]], 4), "| against within-belt background:",
    round(thr_wet[["maxsss"]], 4), "\n")
cat("Presences at or above: p10", sum(pred_occ >= thr[["p10"]]), "| maxSSS",
    sum(pred_occ >= thr[["maxsss"]]), "of", length(pred_occ), "\n")

# ----------------------- Binary surfaces ------------------------------------

suit_p10    <- suit_r >= thr[["p10"]]
suit_maxsss <- suit_r >= thr[["maxsss"]]

cell_area_km2 <- cellSize(suit_r, unit = "km")
area_p10    <- global(suit_p10 * cell_area_km2, "sum", na.rm = TRUE)[[1]]
area_maxsss <- global(suit_maxsss * cell_area_km2, "sum", na.rm = TRUE)[[1]]
total_area  <- global(mask(cell_area_km2, !is.na(suit_r), maskvalues = 0),
                      "sum", na.rm = TRUE)[[1]]

cat("\np10: ", format(round(area_p10), big.mark = ","), "km\u00b2",
    "(", round(100 * area_p10 / total_area, 1), "% of Sudan)\n")
cat("maxSSS:", format(round(area_maxsss), big.mark = ","), "km\u00b2",
    "(", round(100 * area_maxsss / total_area, 1), "% of Sudan)\n")

# ----------------------------- Population -----------------------------------

if (!file.exists(POP_FILE)) {
  response <- GET(POP_URL, user_agent("R - MSc dissertation, e.p.naymon@lse.ac.uk"),
                  write_disk(POP_FILE, overwrite = TRUE), progress())
  stopifnot("WorldPop download failed" = status_code(response) == 200)
}
pop_100m  <- rast(POP_FILE)
raw_total <- global(pop_100m, "sum", na.rm = TRUE)[[1]]

# 100 m -> covariate grid by summing; the total must be conserved
agg_factor  <- round(res(suit_r)[1] / res(pop_100m)[1])
pop_aligned <- resample(aggregate(pop_100m, fact = agg_factor, fun = "sum", na.rm = TRUE),
                        suit_r, method = "sum")
aligned_total <- global(pop_aligned, "sum", na.rm = TRUE)[[1]]
stopifnot("Population total not conserved in alignment" =
  abs(aligned_total - raw_total) / raw_total < POP_TOL)
writeRaster(pop_aligned, POP_ALIGNED_FILE, overwrite = TRUE)

# People in cells with no prediction (outside the domain or a missing covariate)
pop_pred <- global(mask(pop_aligned, suit_r), "sum", na.rm = TRUE)[[1]]
cat("\nWorldPop 2025:", fmt(raw_total), "| in cells with a prediction:", fmt(pop_pred),
    "| without:", fmt(raw_total - pop_pred), "\n")

# ------------------------- National estimates -------------------------------

est <- arp_estimates(suit_r, pop_aligned, thr)
cat("\nPopulation at risk:\n")
for (k in names(est)) cat(sprintf("  %-14s %12s (%.1f%% of WorldPop)\n",
                                  k, fmt(est[[k]]), 100 * est[[k]] / raw_total))

# ---------------------- Where the estimate rests ----------------------------
# Share of the risk-weighted estimate in cells outside the presence range, in
# cells with nights hotter than any presence (suitability saturates there
# without support), and in the >= 150 mm region (supplement comparison with
# the dissertation's within-mask figure).

rw_r   <- pop_aligned * suit_r
mess_r <- rast(MESS_FILE)
hot_r  <- domain_covs("lst_night") > max(train$occ_env$lst_night)
wet_r  <- rast(SENS_MASK_FILE) == 1

part  <- function(cond) global(rw_r * cond, "sum", na.rm = TRUE)[[1]]
parts <- c(outside_presence_range = part(mess_r < 0),
           nights_hotter_than_presences = part(hot_r),
           within_150mm_region = part(wet_r))
cat("\nRisk-weighted estimate by support:\n")
for (k in names(parts)) cat(sprintf("  %-30s %12s (%.1f%% of headline)\n",
                                    k, fmt(parts[[k]]), 100 * parts[[k]] / est[["risk_weighted"]]))

# -------------------------- Plateau candidates ------------------------------
# Risk-weighted estimate from each 06b candidate. The selected model's 06b
# surface must match 07's exactly.

cand     <- rast(CANDIDATES_FILE)
sel_name <- paste0(tuning$fc, " rm ", tuning$rm)
stopifnot("06b surface for the selected model differs from SUIT_FILE" =
  global(abs(cand[[sel_name]] - suit_r), "max", na.rm = TRUE)[[1]] < 1e-6)

cand_rw <- sapply(names(cand), function(nm)
  global(pop_aligned * cand[[nm]], "sum", na.rm = TRUE)[[1]])
cand_df <- data.frame(candidate = names(cand), risk_weighted = round(cand_rw),
                      change_vs_selected = round(100 * (cand_rw / est[["risk_weighted"]] - 1), 1))
cat("\nPlateau candidates:\n"); print(cand_df, row.names = FALSE)

# -------------------------------- States ------------------------------------

zones  <- state_zones(suit_r)
layers <- c(pop_aligned, rw_r,
            pop_aligned * (suit_r >= thr[["p10"]]),
            pop_aligned * (suit_r >= thr[["maxsss"]]),
            rw_r * (mess_r < 0))
names(layers) <- c("total_pop", "arp_weighted", "arp_p10", "arp_maxsss", "arp_outside_range")
cand_r <- pop_aligned * cand
names(cand_r) <- paste0("rw_", gsub(" ", "_", names(cand)))

state_arp <- zonal(c(layers, cand_r), zones, fun = "sum", na.rm = TRUE)
names(state_arp)[1] <- "state"

stopifnot("State estimates do not sum to the national estimate" =
  abs(sum(state_arp$arp_weighted) - est[["risk_weighted"]]) / est[["risk_weighted"]] < POP_TOL)

pres_state <- as.character(terra::extract(zones, occ_xy)[, 1])
state_arp$n_presences <- as.integer(table(factor(pres_state, levels = state_arp$state)))
stopifnot("Not every presence was assigned to a state" =
  sum(state_arp$n_presences) == nrow(train$occ_clean))

alt <- setdiff(names(cand_r), paste0("rw_", gsub(" ", "_", sel_name)))
state_arp <- state_arp |>
  mutate(pct_weighted      = round(100 * arp_weighted / total_pop, 1),
         pct_outside_range = round(100 * arp_outside_range / arp_weighted, 1),
         across(all_of(alt), ~ round(100 * (. / arp_weighted - 1), 1),
                .names = "{.col}_pct_change")) |>
  arrange(desc(arp_weighted))

cat("\nStates (risk-weighted; support; plateau sensitivity):\n")
state_arp |>
  select(state, n_presences, total_pop, arp_weighted, pct_weighted, arp_maxsss,
         arp_p10, pct_outside_range, ends_with("_pct_change")) |>
  mutate(across(c(total_pop, arp_weighted, arp_maxsss, arp_p10), fmt)) |>
  print(right = FALSE)

# -------------- Outside the presence range: which covariate -----------------

mess_vars <- rast(MESS_VARS_FILE)
lim_r <- which.min(mess_vars)
lim_layers <- do.call(c, lapply(seq_len(nlyr(mess_vars)), function(i)
  rw_r * (mess_r < 0) * (lim_r == i)))
names(lim_layers) <- names(mess_vars)
lim_state <- zonal(lim_layers, zones, fun = "sum", na.rm = TRUE)
names(lim_state)[1] <- "state"
cat("\nRisk-weighted estimate outside the presence range, by limiting covariate:\n")
lim_state |> mutate(across(-state, fmt)) |> print(right = FALSE)

# -------------------- Threshold sensitivity curve ---------------------------
v <- data.frame(s = values(suit_r, mat = FALSE), p = values(pop_aligned, mat = FALSE))
v <- v[complete.cases(v), ]
thresh_df <- data.frame(threshold = seq(0, 0.95, by = 0.01))
thresh_df$arp <- sapply(thresh_df$threshold, function(t) sum(v$p[v$s >= t]))


p_thresh <- ggplot(thresh_df, aes(x = threshold, y = arp / 1e6)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = thr[['p10']],    linetype = "dashed", colour = col_p10) +
  geom_vline(xintercept = thr[['maxsss']], linetype = "dashed", colour = col_maxsss) +
  annotate("text", x = thr[["p10"]] + 0.02,    y = max(thresh_df$arp / 1e6) * 0.9,
           label = paste0("p10 (", round(thr[["p10"]], 3), ")"),
           hjust = 0, size = 3.2, colour = "steelblue") +
  annotate("text", x = thr[["maxsss"]] + 0.02, y = max(thresh_df$arp / 1e6) * 0.8,
           label = paste0("maxSSS (", round(thr[["maxsss"]], 3), ")"),
           hjust = 0, size = 3.2, colour = "firebrick") +
  labs(#title = "Threshold sensitivity of at-risk population estimate",
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

# --------------------------------- Save -------------------------------------

writeRaster(suit_r >= thr[["p10"]],    file.path(DIR_SURFACES, "binary_p10.tif"),    overwrite = TRUE)
writeRaster(suit_r >= thr[["maxsss"]], file.path(DIR_SURFACES, "binary_maxsss.tif"), overwrite = TRUE)

write.csv(data.frame(metric = names(est), arp = round(est),
                     pct_of_pop = round(100 * est / raw_total, 1),
                     threshold = c(NA, thr[["p10"]], thr[["maxsss"]])),
          file.path(DIR_TABLES, "arp_summary.csv"), row.names = FALSE)
write.csv(state_arp, file.path(DIR_TABLES, "arp_by_state.csv"), row.names = FALSE)
write.csv(cand_df,   file.path(DIR_TABLES, "arp_plateau_candidates.csv"), row.names = FALSE)
write.csv(lim_state, file.path(DIR_TABLES, "mess_limiting_variable_by_state.csv"), row.names = FALSE)
cat("08_pop_estimate.R complete\n")