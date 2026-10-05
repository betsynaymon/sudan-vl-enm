# ============================================================================
# 07_model_visualization.R
# Response curves, rainfall x night-temperature surface, environmental spread,
# permutation importance, the suitability surface, and extrapolation relative
# to the presences (MESS), for the model selected in 06.
#
# Inputs:  MODEL_FILE, TUNING_FILE, TRAIN_FILE, retained_vars.rds, OCC_FILE,
#          DOMAIN_FILE, ADM0_FILE, ADM1_FILE, SENS_MASK_FILE, COV_FILES
# Outputs: SUIT_FILE, MESS_FILE, MESS_VARS_FILE
#          outputs/tables/response_curve_features.csv, covariate_spread.csv,
#                         permutation_importance.csv
#          outputs/figures/ (curves, surfaces, maps)
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))


suppressPackageStartupMessages({
  library(terra)
  library(maxnet)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggspatial)
  library(sf)
})

set.seed(SEED)

# ------------------------------ Load inputs ---------------------------------
occ   <- read.csv(OCC_FILE)
mod   <- readRDS(MODEL_FILE)
sel   <- readRDS(TUNING_FILE)
train <- readRDS(TRAIN_FILE)
retained_vars <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))

display <- st_as_sf(vect(DISPLAY_ADM0_FILE))
excluded <- st_as_sf(vect(EXCLUDED_FILE))

occ_env <- train$occ_env[, retained_vars]
bg_env  <- train$bg_env[, retained_vars]
env_all <- rbind(occ_env, bg_env)

cat("Selected config:", sel$fc, "rm =", sel$rm, "| rule:", sel$rule, "\n")
cat("Model features:", length(mod$betas), "\n")
cat("Training presences:", nrow(occ_env), "| Background:", nrow(bg_env), "\n")

covs_dom <- domain_covs(retained_vars)

sudan  <- st_as_sf(vect(ADM0_FILE))
states <- st_as_sf(vect(ADM1_FILE))

# Map extent from the display outline, with a margin (plotting_theme.R)
lim <- display_limits(display)

# Reference for curves and surfaces: the presence median. With the desert in
# the background, the background median is a desert cell.
ref <- sapply(occ_env, median)

# ========================== RESPONSE CURVES =================================

response_data <- response_curves(mod, ref, env_all, retained_vars) |>
  mutate(cov_labels = factor(cov_labels[variable], levels = cov_labels))

pres_range <- occ_env |> pivot_longer(everything(), names_to = "variable") |>
  group_by(variable) |> summarise(lo = min(value), hi = max(value)) |>
  mutate(cov_labels = factor(cov_labels[variable], levels = cov_labels))

p_resp <- ggplot(response_data, aes(x = value, y = suit)) +
  geom_rect(data = filter(pres_range, variable != "vertisols"),
            aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey92") +
  geom_line(linewidth = 0.9, colour = "#2166AC") +
  geom_point(data = filter(response_data, variable == "vertisols"),
             size = 3, colour = "#2166AC") +
  facet_wrap(~ cov_labels, scales = "free_x", nrow = 2) +
  labs(x = NULL, y = "Habitat suitability (cloglog)",
       title = paste0("Response curves, selected MaxEnt (", sel$fc, ", rm = ", sel$rm, ")"),
       subtitle = "Others at the presence median; grey = range of presences") +
  theme_minimal() + theme(strip.text = element_text(face = "bold"))

ggsave(file.path(DIR_FIGS, "response_curves.png"), p_resp,
       width = 10, height = 6, dpi = 300)
cat("Saved response_curves.png\n")

# -------------- Response-curve features ----------
# The text quotes positions read off these marginal curves: the LST-night
# threshold, rainfall peak and decline, slope peak. Other covariates are held
# at the presence median, so x-positions are more robust than heights.
#   rise_XX = lowest value where suitability reaches XX% of the curve's peak
#   fall_XX = highest value where suitability is still at XX% of the peak

resp_features <- curve_table(response_data)
vert <- response_data |> filter(variable == "vertisols")

cat("\nResponse-curve features:\n")
resp_features |>
  mutate(across(where(is.numeric), ~ signif(., 3))) |>
  print(width = Inf)
cat("Vertisols: suitability", round(vert$suit[vert$value == 0], 3),
    "(absent) vs", round(vert$suit[vert$value == 1], 3), "(present)\n")

write.csv(resp_features, file.path(DIR_TABLES, "response_curve_features.csv"),
          row.names = FALSE)

# ------------- Rainfall x night temperature (interaction) -------------------

surf_rl <- pair_surface(mod, ref, env_all, "rainfall", "lst_night")
p_rl <- ggplot(surf_rl, aes(rainfall, lst_night, fill = suit)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Suitability", limits = c(0, 1)) +
  geom_point(data = bg_env, aes(rainfall, lst_night), inherit.aes = FALSE,
             colour = "grey80", size = 0.2, alpha = 0.15) +
  geom_point(data = occ_env, aes(rainfall, lst_night), inherit.aes = FALSE,
             shape = 21, fill = "white", size = 1.3) +
  labs(x = cov_labels[["rainfall"]], y = cov_labels[["lst_night"]],
       title = "Rainfall x night temperature, selected model",
       subtitle = "Others at the presence median; white = presences, grey = background") +
  theme_minimal()
ggsave(file.path(DIR_FIGS, "rain_lst_surface.png"), p_rl, width = 6, height = 5, dpi = 300)

# ======================= ENVIRONMENTAL SPREAD ===============================

occ_spread <- occ_env |> mutate(type = "Presence")
bg_spread  <- bg_env  |> mutate(type = "Background")

env_both <- bind_rows(occ_spread, bg_spread) |>
  pivot_longer(cols = all_of(retained_vars), names_to = "variable",
               values_to = "value") |>
  mutate(var_label = cov_labels[variable],
         var_label = factor(var_label, levels = cov_labels))

p_spread <- ggplot(env_both |> filter(variable != "vertisols"),
                   aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ var_label, scales = "free", nrow = 2) +
  scale_fill_manual(values = c(Background = "grey60", Presence = "#B2182B")) +
  labs(x = NULL, y = "Density", fill = NULL,
       title = "Environmental spread: presences vs. background",
       subtitle = "Narrow overlap = potential geographic proxy") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(DIR_FIGS, "env_spread_density.png"), p_spread,
       width = 10, height = 6, dpi = 300)
cat("Saved env_spread_density.png\n")

# --------------- Covariate spread ----------------
# Values at training presences vs background (year-matched extractions).

spread_q <- env_both |>
  filter(variable != "vertisols") |>
  group_by(variable, type) |>
  summarise(
    n = n(), min = min(value),
    q05 = quantile(value, 0.05), q25 = quantile(value, 0.25),
    median = median(value),
    q75 = quantile(value, 0.75), q95 = quantile(value, 0.95),
    max = max(value), .groups = "drop"
  ) |>
  arrange(variable, desc(type))

cat("\nCovariate spread, presences vs background:\n")
spread_q |>
  mutate(across(min:max, ~ signif(., 4))) |>
  print(n = Inf, width = Inf)

vert_share <- env_both |>
  filter(variable == "vertisols") |>
  group_by(type) |>
  summarise(pct_on_vertisols = round(100 * mean(value == 1), 1))
cat("\nShare of points on vertisols (%):\n")
print(vert_share)

write.csv(spread_q, file.path(DIR_TABLES, "covariate_spread.csv"),
          row.names = FALSE)

# ====================== PERMUTATION IMPORTANCE ==============================

# In-sample AUC drop when one covariate is shuffled across presences and
# background, over all of Sudan and within the >= 150 mm region. Over all of
# Sudan importance mostly reflects belt vs desert; within the region it shows
# what discriminates inside the belt.

sens_r  <- rast(SENS_MASK_FILE)
occ_wet <- terra::extract(sens_r, as.matrix(train$occ_clean[, c("longitude", "latitude")]))[, 1] %in% 1
bg_wet  <- terra::extract(sens_r, as.matrix(train$bg_clean[,  c("longitude", "latitude")]))[, 1] %in% 1

set.seed(SEED); imp_all <- perm_importance(mod, occ_env, bg_env, retained_vars)
set.seed(SEED); imp_wet <- perm_importance(mod, occ_env[occ_wet, ], bg_env[bg_wet, ], retained_vars)

perm_df <- bind_rows(
  data.frame(scope = "All Sudan",        variable = retained_vars, imp_all, row.names = NULL),
  data.frame(scope = "Within >= 150 mm", variable = retained_vars, imp_wet, row.names = NULL)) |>
  mutate(var_label = cov_labels[variable])
write.csv(perm_df, file.path(DIR_TABLES, "permutation_importance.csv"), row.names = FALSE)

cat("\nPermutation importance (AUC drop,", N_PERM, "permutations):\n")
perm_df |> mutate(across(c(mean, sd), ~ round(., 4))) |>
  select(scope, variable, mean, sd) |> arrange(scope, desc(mean)) |> print()

p_imp <- perm_df |>
  mutate(var_label = factor(var_label, levels = rev(cov_labels))) |>
  ggplot(aes(x = mean, y = var_label)) +
  geom_point(size = 3, colour = "#2166AC") +
  geom_errorbar(aes(xmin = mean - sd, xmax = mean + sd),
                width = 0.2, orientation = "y", colour = "#2166AC") +
  facet_wrap(~ scope) +
  labs(x = "AUC drop (permutation importance)",
       y = NULL,
       title = paste0("Variable importance \u2014 selected MaxEnt (",
                      sel$fc, ", rm = ", sel$rm, ")"),
       subtitle = paste0("Mean \u00b1 SD across ", N_PERM, " permutations")) +
  theme_minimal()

ggsave(file.path(DIR_FIGS, "variable_importance.png"), p_imp,
       width = 8, height = 5, dpi = 300)
cat("Saved variable_importance.png\n")

# ======================= PREDICTION SURFACE =================================
# Predicted over the domain: the same extent the model was calibrated on.
pred_r <- terra::predict(covs_dom, mod, type = "cloglog", na.rm = TRUE)
names(pred_r) <- "suitability"
writeRaster(pred_r, SUIT_FILE, overwrite = TRUE)

n_pred <- global(!is.na(pred_r), "sum")[[1]]
n_ok   <- global(!any(is.na(covs_dom)), "sum", na.rm = TRUE)[[1]]
stopifnot("A complete domain cell has no prediction" = n_pred == n_ok)
cat("Prediction surface:", n_pred, "cells | range",
    paste(round(minmax(pred_r)[, 1], 4), collapse = " to "), "\n")

pred_df <- as.data.frame(pred_r, xy = TRUE)
names(pred_df) <- c("x", "y", "suitability")
pred_df <- pred_df[!is.na(pred_df$suitability), ]

cat("Saved suitability map\n")

# ============================ DARFUR ZOOM ===================================
darfur <- states[grepl("Darfur", states$NAME_1), ]
p_darfur <- ggplot() +
  geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
  geom_raster(data = pred_df, aes(x = x, y = y, fill = suitability)) +
  scale_fill_gradientn(
    colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7",
                "#EF8A62", "#B2182B"),
    na.value = "transparent",
    name = "Habitat\nsuitability",
    limits = c(0, 1)
  ) +
  geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
  geom_point(data = occ, aes(x = longitude, y = latitude),
             colour = "black", fill = "white",
             shape = 21, size = 3, stroke = 0.7) +
  coord_sf(xlim = st_bbox(darfur)[c("xmin", "xmax")],
           ylim = st_bbox(darfur)[c("ymin", "ymax")], crs = 4326) +
  labs(title = "Darfur states, zoomed",
       subtitle = "Presence points over suitability surface") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank())

ggsave(file.path(DIR_FIGS, "suitability_darfur_zoom.png"), p_darfur,
       width = 8, height = 7, dpi = 300)
cat("Saved suitability_darfur_zoom.png\n")

# ============================ EXTRAPOLATION =================================
# MESS against the presences (year-matched training values). The background
# spans Sudan, so almost nothing is novel relative to presences plus
# background; the informative question is whether a cell lies within the
# range where presences support the fitted response.
mess_full <- mask(rast(dismo::mess(x = raster::stack(covs_dom), v = occ_env, full = TRUE)),
                  pred_r)
names(mess_full) <- c(retained_vars, "mess")
mess_r    <- mess_full[["mess"]]
mess_vars <- mess_full[[retained_vars]]
writeRaster(mess_r,    MESS_FILE,      overwrite = TRUE)
writeRaster(mess_vars, MESS_VARS_FILE, overwrite = TRUE)
stopifnot("MESS has values outside the complete domain cells" =
  global(!is.na(mess_r), "sum")[[1]] == n_ok)

mess_vals <- values(mess_r, na.rm = TRUE)
cat("\nMESS vs presences: outside the presence range in", sum(mess_vals < 0), "of",
    length(mess_vals), "cells (", round(100 * mean(mess_vals < 0), 1), "%)\n")

novel <- which(values(mess_r) < 0)
lim   <- apply(values(mess_vars)[novel, , drop = FALSE], 1, which.min)
cat("Limiting covariate among those cells:\n"); print(table(retained_vars[lim]))

mess_diff <- global(abs(min(mess_vars) - mess_r), "max", na.rm = TRUE)[[1]]
cat("Max |min(per-variable) - overall MESS|:", signif(mess_diff, 3), "\n")

# ================== DISSERTATION FIGURE ======================================
# Produces polished versions of the suitability and MESS maps using the
# shared dissertation theme. 
#
# Outputs:
#   outputs/figures/fig_suitability_mess.png  — combined panel, main text
#   outputs/figures/fig_suitability_mess.pdf
#   outputs/figures/fig_suitability.png       — standalone, full width
#   outputs/figures/fig_suitability.pdf
# =============================================================================

library(patchwork)
source(here::here("R", "plotting_theme.R"))
 
# ---- Occurrences as sf (for geom_sf consistency) ----
occ_sf <- st_as_sf(occ, coords = c("longitude", "latitude"), crs = 4326)
 
# ---- Clip MESS to Sudan boundary  ----
mess_df_clipped <- as.data.frame(mess_r, xy = TRUE)
names(mess_df_clipped) <- c("x", "y", "mess")
mess_df_clipped <- mess_df_clipped[!is.na(mess_df_clipped$mess), ]
mess_df_clipped$type <- ifelse(mess_df_clipped$mess < 0,
                               "Extrapolation", "Interpolation")
 
# ---------- Panel (a): Suitability surface -----------------------------------
p_suit_a <- ggplot() +
  geom_sf(data = sudan, fill = "grey95", colour = NA) +
  geom_raster(data = pred_df, aes(x = x, y = y, fill = suitability)) +
  layer_excluded(excluded) +
  scale_fill_suitability() +
  guides(fill = guide_colourbar(
    barheight = unit(0.4, "cm"),
    barwidth  = unit(3, "cm"),
    title.position = "left",
    title.theme = element_text(size = 7, face = "plain", vjust = 0.8),
    label.theme = element_text(size = 6)
  )) +
  layer_admin1(data = states, colour = "black", linewidth = 0.15) +
  layer_country(data = display, colour = "black", linewidth = 0.3) +
  geom_sf(data = occ_sf, shape = 21, size = 1.0, stroke = 0.3,
          fill = "white", colour = "black") +
  ggspatial::annotation_scale(
    location = "tl", width_hint = 0.15, text_cex = 0.5,
    line_width = 0.3, pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm")
  ) +
  labs(title = "(a) Predicted suitability") +
  coord_display(lim) +
  theme_map() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = 9, face = "plain", hjust = 0)
  )
 
 
# ---------- Panel (b): MESS (interpolation vs. extrapolation) ----------------
p_mess_b <- ggplot() +
  geom_sf(data = sudan, fill = "grey95", colour = NA) +
  geom_raster(data = mess_df_clipped, aes(x = x, y = y, fill = type)) +
  layer_excluded(excluded) +
  scale_fill_manual(values = pal_mess_binary, name = NULL) +
  guides(fill = guide_legend(
    keywidth  = unit(0.5, "cm"),
    keyheight = unit(0.4, "cm"),
    direction = "horizontal"
  )) +
  layer_admin1(data = states, colour = "black", linewidth = 0.15) +
  layer_country(data = display, colour = "black", linewidth = 0.3) +
  geom_sf(data = occ_sf, shape = 21, size = 1.0, stroke = 0.3,
          fill = "white", colour = "black") +
  labs(title = "(b) Outside the presence range") +
  coord_display(lim) +
  theme_map() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = 9, face = "plain", hjust = 0)
  )
 
 
# ---------- Combined figure --------------------------------------------------
 
fig_suit_mess <- p_suit_a + p_mess_b
 
save_fig(file.path(DIR_FIGS, "fig_suitability_mess.png"), fig_suit_mess,
         width = FIG_WIDTH_FULL, height = 12)
save_fig(file.path(DIR_FIGS, "fig_suitability_mess.pdf"), fig_suit_mess,
         width = FIG_WIDTH_FULL, height = 12)
cat("Saved fig_suitability_mess\n")
 
 
# ---------- Standalone suitability (full width) ------------------------------
 
p_suit_full <- ggplot() +
  geom_sf(data = sudan, fill = "grey95", colour = NA) +
  geom_raster(data = pred_df, aes(x = x, y = y, fill = suitability)) +
  layer_excluded(excluded) +
  scale_fill_suitability() +
  layer_admin1(data = states, colour = "black", linewidth = 0.15) +
  layer_country(data = display, colour = "black", linewidth = 0.3) +
  geom_sf(data = occ_sf, shape = 21, size = 1.5, stroke = 0.4,
          fill = "white", colour = "black") +
  add_scalebar() +
  add_north_arrow() +
  coord_display(lim) +
  theme_map()
 
save_fig(file.path(DIR_FIGS, "fig_suitability.png"), p_suit_full,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)
save_fig(file.path(DIR_FIGS, "fig_suitability.pdf"), p_suit_full,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)
cat("Saved fig_suitability (standalone)\n")

# ---------- Response curves + variable importance (2x3 grid) -----------------
# Five response curves fill positions 1-5; variable importance fills the 6th
# slot (bottom-right).
# Objects needed: response_data, perm_df, cov_labels (from earlier in script)

# Covariate color palette 
pal_covariates <- c(
  "Slope (degrees)"          = "#CC79A7",
  "Distance to river (m)"    = "#009E73",
  "Vertisols (0/1)"          = "#D55E00",
  "LST night (\u00b0C)"      = "#E69F00",
  "Rainfall (mm/yr)"         = "#0072B2"
)

# Helper: one response curve panel
make_response <- function(var, show_ylab = FALSE) {
  lab <- cov_labels[var]
  col <- pal_covariates[lab]
  d <- response_data |> filter(variable == var)

  p <- ggplot(d, aes(x = value, y = suit)) +
    geom_line(linewidth = 0.7, colour = col)

  if (var == "vertisols") p <- p + geom_point(size = 2, colour = col)

  p +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = lab, x = NULL,
         y = if (show_ylab) "Predicted suitability" else NULL) +
    theme_dissertation(gridlines = "both") +
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          panel.grid.major = element_line(colour = "grey92"),
          axis.text = element_text(size = 7))
}

# Build the five response panels 
p1 <- make_response("slope",      show_ylab = TRUE)
p2 <- make_response("river_dist")
p3 <- make_response("vertisols")
p4 <- make_response("lst_night",  show_ylab = TRUE)
p5 <- make_response("rainfall")

# Variable importance 
p_imp_grid <- perm_df |>
  filter(scope == "All Sudan") |>
  mutate(var_label = factor(var_label, levels = rev(cov_labels))) |>
  ggplot(aes(x = mean, y = var_label, colour = var_label)) +
  geom_segment(aes(x = mean - sd, xend = mean + sd, yend = var_label),
               linewidth = 0.5, show.legend = FALSE) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = pal_covariates, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3))) +
  labs(title = "Variable importance", x = "AUC drop", y = NULL) +
  theme_dissertation() +
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        panel.grid.major = element_line(colour = "grey92"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 5),
        axis.title.x = element_text(size = 7))

# 2x3 grid
fig_resp_imp <- (p1 + p2 + p3) / (p4 + p5 + p_imp_grid) &
  theme(legend.position = "none")

save_fig(file.path(DIR_FIGS, "fig_response_importance.png"), fig_resp_imp,
         width = FIG_WIDTH_FULL, height = 12)
save_fig(file.path(DIR_FIGS, "fig_response_importance.pdf"), fig_resp_imp,
         width = FIG_WIDTH_FULL, height = 12)
cat("Saved fig_response_importance\n")

cat("07_model_visualization.R complete\n")