# ============================================================================
# 12_sampling_bias_robustness.R
# Does handling sampling bias change the surface and the estimate? Refits the
# selected model with the background drawn in proportion to accessibility,
# (1 + travel time)^-b, so presences and background share the bias, at two
# strengths: b = 0.5 (the dissertation's weighting) and b matched so the
# background is as close to cities as the presences (the calibration of 13).
# Matching treats all of the presences' extra accessibility as sampling bias,
# so it is an upper bound on the correction: it also removes any ecological
# association with accessibility. The uniform background (06) runs through
# the same code and must reproduce 06 and 07.
#
# Inputs:  TRAIN_FILE, TUNING_FILE, MODEL_FILE, SUIT_FILE, POP_ALIGNED_FILE,
#          FOLDS_FILE, TT_FILE, DOMAIN_FILE, ADM1_FILE, DISPLAY_ADM0_FILE,
#          SENS_MASK_FILE, retained_vars.rds, arp_plateau_candidates.csv (08),
#          covariate rasters (COV_FILES, COV_ANNUAL)
# Outputs: BIAS_SURFACES_FILE (layers: uniform, half, matched)
#          outputs/tables/bias_correction_summary.csv
#          outputs/tables/bias_correction_by_state.csv
#          outputs/figures/response_curves_bias_comparison.png
#          outputs/figures/fig_bias_comparison.png / .pdf
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(maxnet)
  library(ggplot2); library(patchwork)
})

# ------------------------------ Load inputs ---------------------------------

train   <- readRDS(TRAIN_FILE)
tuning  <- readRDS(TUNING_FILE)
mod     <- readRDS(MODEL_FILE)
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
suit_r  <- rast(SUIT_FILE)
pop     <- rast(POP_ALIGNED_FILE)
cand_08 <- read.csv(file.path(DIR_TABLES, "arp_plateau_candidates.csv"))
occ     <- train$occ_clean
bg      <- train$bg_clean
occ_env <- train$occ_env[, vars]
bg_env  <- train$bg_env[, vars]
cat("Model:", tuning$fc, "rm =", tuning$rm, "| presences:", nrow(occ),
    "| background:", nrow(bg), "\n")

# ------------------------- Corrected backgrounds ----------------------------
# Frame: domain cells with every long-term covariate (04). Cells without travel
# time get weight 0 (no fill; 11). Drawn with draw_points(): with replacement,
# years from the presences' distribution, folds from 05's blocks.

covs    <- domain_covs(vars)
frame   <- !any(is.na(covs))
frame   <- mask(frame, frame, maskvalues = 0)
tt      <- travel_time(frame)
cells   <- which(!is.na(values(frame, mat = FALSE)))
tt_c    <- values(tt, mat = FALSE)[cells]
has_tt  <- !is.na(tt_c)
occ_tt  <- terra::extract(tt, as.matrix(occ[, c("longitude", "latitude")]))[, 1]
b_match <- tt_exponent(tt_c, median(occ_tt))
cat("Frame cells:", length(cells), "| matched b =", b_match, "\n")

year_w  <- table(occ$year)
draw_bg <- function(b, seed)
  draw_points(frame, cells, N_BACKGROUND, ifelse(has_tt, (1 + tt_c)^(-b), 0),
              year_w, vars, seed)

bgs <- list(
  uniform = list(env = bg_env, fold = bg$fold, dropped = 0,
                 cell = cellFromXY(frame, as.matrix(bg[, c("longitude", "latitude")]))),
  half    = draw_bg(0.5,     SEED + 2),
  matched = draw_bg(b_match, SEED + 3)
)
bias_b      <- c(uniform = NA, half = 0.5, matched = b_match)
bias_labels <- c(uniform = "Uniform background", half = "Corrected, b = 0.5",
                 matched = paste0("Corrected, as accessible as presences (b = ", b_match, ")"))

q3 <- function(x) paste(round(quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)), collapse = " / ")
cat("\nTravel time q25 / median / q75 (min). Presences:", q3(occ_tt), "\n")
for (nm in names(bgs)) {
  t_bg <- values(tt, mat = FALSE)[bgs[[nm]]$cell]
  cat(sprintf("%-8s %d points, %d distinct cells, %d dropped | %s | P[presence closer] %.3f\n",
              nm, nrow(bgs[[nm]]$env), length(unique(bgs[[nm]]$cell)), bgs[[nm]]$dropped,
              q3(t_bg), auc_ties(-occ_tt, -t_bg[!is.na(t_bg)])))
}

# -------------------------------- Refits ------------------------------------
# Selected configuration; presences and their folds unchanged. CV scores each
# model against its own background, so scores are not comparable across rows
# (under the matched background they describe separating presences from
# equally accessible places). Thresholds are taken against the uniform
# background, so binary estimates share one reference.

refits <- setNames(lapply(names(bgs), function(nm) {
  r <- refit_maxnet(occ_env, bgs[[nm]]$env, occ$fold, bgs[[nm]]$fold, tuning$fc, tuning$rm)
  r$thr  <- thresholds_from(as.numeric(predict(r$mod, occ_env, type = "cloglog")),
                            as.numeric(predict(r$mod, bg_env,  type = "cloglog")))
  r$surf <- terra::predict(covs, r$mod, type = "cloglog", na.rm = TRUE)
  r$arp  <- arp_estimates(r$surf, pop, r$thr)
  cat("  ", nm, ": ", length(r$mod$betas), " features\n", sep = "")
  r
}), names(bgs))

p <- refits$uniform
stopifnot(
  "Coefficients differ from MODEL_FILE (06)" = isTRUE(all.equal(p$mod$betas, mod$betas)),
  "CV CBI differs from 06"                   = abs(mean(p$cv$cbi) - tuning$cbi) < 1e-8,
  "Surface differs from SUIT_FILE (07)" =
    global(abs(p$surf - suit_r), "max", na.rm = TRUE)[[1]] < 1e-6
)
cat("Uniform background reproduces 06 (coefficients, CV CBI) and 07 (surface)\n")

# ------------------------------ Comparison ----------------------------------
# Reading fixed before the run: the backgrounds bracket the handling of
# sampling bias, from none (uniform) to full (matched, an upper bound). A
# national change within the plateau spread from 08 means accessibility
# handling matters less than model choice; beyond it, the estimate depends on
# how sampling bias is treated, and the uniform-to-matched range is reported
# beside the data-quality range (10).

plateau <- max(abs(cand_08$change_vs_selected[startsWith(cand_08$candidate, paste0(tuning$fc, " "))]))
belt_r  <- rast(SENS_MASK_FILE) == 1
belt    <- values(belt_r, mat = FALSE) %in% TRUE
v_prim  <- values(p$surf, mat = FALSE)
rw_part <- function(s, region) global(pop * s * region, "sum", na.rm = TRUE)[[1]]

summary_df <- bind_rows(lapply(names(refits), function(nm) {
  r <- refits[[nm]]; v <- values(r$surf, mat = FALSE); ok <- !is.na(v) & !is.na(v_prim)
  data.frame(
    background = nm, b = bias_b[[nm]], n_coef = length(r$mod$betas),
    cv_cbi_own_bg = mean(r$cv$cbi, na.rm = TRUE), cv_auc_own_bg = mean(r$cv$auc),
    rho_domain = cor(v[ok], v_prim[ok], method = "spearman"),
    rho_belt   = cor(v[ok & belt], v_prim[ok & belt], method = "spearman"),
    arp_weighted   = r$arp[["risk_weighted"]],
    rw_in_150mm    = rw_part(r$surf, belt_r),
    rw_below_150mm = rw_part(r$surf, !belt_r),
    arp_maxsss = r$arp[["maxsss"]], maxsss = r$thr[["maxsss"]],
    arp_p10    = r$arp[["p10"]],    p10    = r$thr[["p10"]])
})) |>
  mutate(rw_change_pct      = 100 * (arp_weighted / arp_weighted[background == "uniform"] - 1),
         beyond_plateau     = abs(rw_change_pct) > plateau,
         change_in_150mm    = rw_in_150mm    - rw_in_150mm[background == "uniform"],
         change_below_150mm = rw_below_150mm - rw_below_150mm[background == "uniform"])

cat("\nFit and surface (CV against each model's own background):\n")
summary_df |>
  select(background, b, n_coef, cv_cbi_own_bg, cv_auc_own_bg, rho_domain, rho_belt) |>
  mutate(across(where(is.double), ~ round(., 3))) |> print(row.names = FALSE)

cat("\nPopulation at risk | plateau spread: +/-", round(plateau, 1), "%\n")
summary_df |>
  select(background, arp_weighted, rw_change_pct, beyond_plateau, change_in_150mm,
         change_below_150mm, arp_maxsss, arp_p10) |>
  mutate(across(c(arp_weighted, change_in_150mm, change_below_150mm, arp_maxsss, arp_p10), fmt),
         rw_change_pct = round(rw_change_pct, 1)) |>
  print(row.names = FALSE)

# -------------------------------- States ------------------------------------

zones    <- state_zones(suit_r)
surf_all <- do.call(c, unname(lapply(refits, `[[`, "surf")))
names(surf_all) <- names(refits)
rw_all   <- pop * surf_all
names(rw_all) <- names(refits)

state_df <- zonal(rw_all, zones, fun = "sum", na.rm = TRUE)
names(state_df)[1] <- "state"
stopifnot("States do not sum to the national estimates" =
  all(abs(colSums(state_df[names(refits)]) / summary_df$arp_weighted - 1) < POP_TOL))
state_df <- state_df |>
  mutate(across(all_of(names(refits)[-1]), ~ round(100 * (. / uniform - 1), 1),
                .names = "{.col}_pct")) |>
  arrange(desc(uniform))

cat("\nRisk-weighted estimate by state (uniform; % change under each correction):\n")
state_df |> select(state, uniform, ends_with("_pct")) |>
  mutate(uniform = fmt(uniform)) |> print(right = FALSE, row.names = FALSE)

# ------------------------------- Figures ------------------------------------
# Response curves with the others at the presence median (as 07), over the
# uniform background's range, so the models share one reference.

ref     <- sapply(occ_env, median)
env_all <- rbind(occ_env, bg_env)
curves  <- bind_rows(lapply(names(refits), function(nm)
  data.frame(background = nm, response_curves(refits[[nm]]$mod, ref, env_all, vars)))) |>
  mutate(background = factor(background, levels = names(pal_bias)))

p_curves <- ggplot(curves, aes(value, suit, colour = background)) +
  geom_line(data = filter(curves, variable != "vertisols"), linewidth = 0.7) +
  geom_point(data = filter(curves, variable == "vertisols"), size = 2) +
  facet_wrap(~ variable, scales = "free_x", nrow = 2) +
  scale_colour_manual(values = pal_bias, labels = bias_labels, name = NULL) +
  labs(x = NULL, y = "Suitability (cloglog)") +
  theme(legend.position = "bottom")
save_fig(file.path(DIR_FIGS, "response_curves_bias_comparison.png"), p_curves,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT)

display  <- st_as_sf(vect(DISPLAY_ADM0_FILE))
excluded <- st_as_sf(vect(EXCLUDED_FILE))
states   <- st_as_sf(vect(ADM1_FILE))
lim      <- display_limits(display)

map_panel <- function(nm, title) {
  df <- as.data.frame(surf_all[[nm]], xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "suitability")
  ggplot() +
    geom_raster(data = df, aes(x, y, fill = suitability)) +
    scale_fill_suitability() +
    layer_excluded(excluded) +
    layer_admin1(data = states, colour = "black", linewidth = 0.15) +
    layer_country(data = display, colour = "black", linewidth = 0.3) +
    labs(title = title) +
    coord_display(lim) +
    theme_map() +
    theme(plot.title = element_text(size = 9, hjust = 0))
}
fig_maps <- map_panel("uniform", "(a) Uniform background") +
  map_panel("matched", "(b) Background matched to presence accessibility") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
for (ext in c("png", "pdf"))
  save_fig(file.path(DIR_FIGS, paste0("fig_bias_comparison.", ext)), fig_maps,
           width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT)

# --------------------------------- Save -------------------------------------

writeRaster(surf_all, BIAS_SURFACES_FILE, overwrite = TRUE)
write.csv(summary_df, file.path(DIR_TABLES, "bias_correction_summary.csv"),  row.names = FALSE)
write.csv(state_df,   file.path(DIR_TABLES, "bias_correction_by_state.csv"), row.names = FALSE)
cat("12_sampling_bias_robustness.R complete\n")



