# ============================================================================
# 07_model_visualization.R
# Produces response curves, variable importance, the continuous suitability
# surface, a Darfur zoom, and a MESS extrapolation check from the fitted
# MaxEnt model.
#
# Expensive computations (prediction surface, MESS, permutation importance)
# are cached to disk. To iterate on figures without recomputing, simply
# rerun — cached results load instantly.
#
# Inputs:  outputs/models/maxent_final.rds
#          outputs/models/selected_tuning.rds
#          outputs/models/training_data.rds
#          outputs/models/retained_vars.rds
#          data/raw/ (LTM covariate rasters + ecological mask)
#          data/processed/occurrences_thinned.csv
# Outputs: outputs/surfaces/maxent_suitability.tif
#          outputs/tables/permutation_importance.csv
#          outputs/figures/response_curves.png
#          outputs/figures/env_spread_density.png
#          outputs/figures/variable_importance.png
#          outputs/figures/maxent_suitability_map.pdf
#          outputs/figures/maxent_suitability_map.png
#          outputs/figures/suitability_darfur_zoom.png
#          outputs/figures/mess_extrapolation.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(maxnet)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggspatial)
  library(rnaturalearth)
  library(ecospat)
  library(geodata)
  library(sf)
})

set.seed(SEED)

# ------------------------------ Load inputs ---------------------------------

occ <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))

mod   <- readRDS(file.path(DIR_MODELS, "maxent_final.rds"))
sel   <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))
train <- readRDS(file.path(DIR_MODELS, "training_data.rds"))

occ_env   <- train$occ_env
bg_env    <- train$bg_env
occ_clean <- train$occ_clean

cat("Selected config:", sel$fc, "rm =", sel$rm, "\n")
cat("Model coefficients:", length(mod$betas), "\n")
cat("Training presences:", nrow(occ_env),
    "| Background:", nrow(bg_env), "\n")

retained_vars <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))

covs <- rast(file.path(DIR_COVARIATES, COV_FILES[retained_vars]))
names(covs) <- retained_vars

mask_r <- rast(file.path(DIR_COVARIATES, "ecological_mask_150mm.tif"))
covs_masked <- mask(covs, mask_r, maskvalues = 0)

adm0  <- gadm(country = "SDN", level = 0, path = here::here("data", "raw"))
sudan <- st_as_sf(adm0)

# Shared label mapping for figures
var_labels <- c(
  slope      = "Slope (degrees)",
  river_dist = "Distance to river (m)",
  vertisols  = "Vertisols (0/1)",
  lst_night  = "LST night (\u00b0C)",
  rainfall   = "Rainfall (mm/yr)"
)

# ========================== RESPONSE CURVES =================================

bg_medians <- apply(bg_env, 2, median)
n_pts <- 200

response_data <- bind_rows(lapply(retained_vars, function(var) {
  if (var == "vertisols") {
    newdata <- as.data.frame(t(replicate(2, bg_medians)))
    newdata[[var]] <- c(0, 1)
  } else {
    newdata <- as.data.frame(t(replicate(n_pts, bg_medians)))
    newdata[[var]] <- seq(min(bg_env[[var]]), max(bg_env[[var]]),
                         length.out = n_pts)
  }

  newdata$suitability <- predict(mod, newdata, clamp = TRUE, type = "cloglog")

  tibble(
    variable = var,
    value    = newdata[[var]],
    suit     = as.numeric(newdata$suitability)
  )
}))

response_data <- response_data |>
  mutate(var_label = var_labels[variable],
         var_label = factor(var_label, levels = var_labels))

p_resp <- ggplot(response_data, aes(x = value, y = suit)) +
  geom_line(linewidth = 0.9, colour = "#2166AC") +
  geom_point(data = response_data |> filter(variable == "vertisols"),
             size = 3, colour = "#2166AC") +
  facet_wrap(~ var_label, scales = "free_x", nrow = 2) +
  labs(x = NULL,
       y = "Habitat suitability (cloglog)",
       title = paste0("Marginal response curves \u2014 selected MaxEnt (",
                      sel$fc, ", rm = ", sel$rm, ")"),
       subtitle = "Each covariate varied across background range; others held at median") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(DIR_FIGS, "response_curves.png"), p_resp,
       width = 10, height = 6, dpi = 300)
cat("Saved response_curves.png\n")

# ======================= ENVIRONMENTAL SPREAD ===============================

occ_spread <- occ_env |> mutate(type = "Presence")
bg_spread  <- bg_env  |> mutate(type = "Background")

env_both <- bind_rows(occ_spread, bg_spread) |>
  pivot_longer(cols = all_of(retained_vars), names_to = "variable",
               values_to = "value") |>
  mutate(var_label = var_labels[variable],
         var_label = factor(var_label, levels = var_labels))

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

# ====================== PERMUTATION IMPORTANCE ==============================

perm_path <- file.path(DIR_TABLES, "permutation_importance.csv")

if (file.exists(perm_path)) {
  cat("Loading cached permutation importance\n")
  perm_df <- read.csv(perm_path)
} else {
  pred_occ <- predict(mod, occ_env, type = "cloglog") |> as.numeric()
  pred_bg  <- predict(mod, bg_env, type = "cloglog") |> as.numeric()
  baseline_auc <- mean(sapply(pred_occ, function(p) mean(p > pred_bg)))
  cat("Baseline AUC:", round(baseline_auc, 4), "\n")

  set.seed(SEED)
  n_perm <- 50

  perm_imp <- sapply(retained_vars, function(var) {
    drops <- replicate(n_perm, {
      occ_shuf <- occ_env
      bg_shuf  <- bg_env
      all_vals <- c(occ_env[[var]], bg_env[[var]])
      shuffled <- sample(all_vals)
      occ_shuf[[var]] <- shuffled[1:nrow(occ_env)]
      bg_shuf[[var]]  <- shuffled[(nrow(occ_env) + 1):length(shuffled)]

      p_occ <- predict(mod, occ_shuf, type = "cloglog") |> as.numeric()
      p_bg  <- predict(mod, bg_shuf, type = "cloglog") |> as.numeric()

      perm_auc <- mean(sapply(p_occ, function(p) mean(p > p_bg)))
      baseline_auc - perm_auc
    })
    c(mean = mean(drops), sd = sd(drops))
  })

  perm_df <- as.data.frame(t(perm_imp)) |>
    mutate(variable = retained_vars,
           var_label = var_labels[variable]) |>
    arrange(desc(mean))

  write.csv(perm_df, perm_path, row.names = FALSE)
  cat("Computed and saved permutation importance\n")
}

cat("Permutation importance (AUC drop, 50 reps):\n")
for (i in seq_len(nrow(perm_df))) {
  cat("  ", perm_df$var_label[i], ":",
      round(perm_df$mean[i], 4), "\u00b1",
      round(perm_df$sd[i], 4), "\n")
}

p_imp <- perm_df |>
  mutate(var_label = factor(var_label, levels = rev(var_label))) |>
  ggplot(aes(x = mean, y = var_label)) +
  geom_point(size = 3, colour = "#2166AC") +
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd),
                 height = 0.2, colour = "#2166AC") +
  labs(x = "AUC drop (permutation importance)",
       y = NULL,
       title = paste0("Variable importance \u2014 selected MaxEnt (",
                      sel$fc, ", rm = ", sel$rm, ")"),
       subtitle = "Mean \u00b1 SD across 50 permutations") +
  theme_minimal()

ggsave(file.path(DIR_FIGS, "variable_importance.png"), p_imp,
       width = 8, height = 5, dpi = 300)
cat("Saved variable_importance.png\n")

# ======================= PREDICTION SURFACE =================================

suit_path <- file.path(DIR_SURFACES, "maxent_suitability.tif")

if (file.exists(suit_path)) {
  cat("Loading cached suitability surface\n")
  pred_r <- rast(suit_path)
} else {
  covs_sudan <- mask(covs, vect(sudan))
  pred_r <- terra::predict(covs_sudan, mod, type = "cloglog", na.rm = TRUE)

  writeRaster(pred_r, suit_path, overwrite = TRUE)
  cat("Computed and saved suitability surface\n")
}

cat("Prediction surface:\n")
cat("  Non-NA:", sum(!is.na(values(pred_r))), "\n")
cat("  Range:", round(minmax(pred_r)[1], 4), "\u2013",
    round(minmax(pred_r)[2], 4), "\n")

pred_df <- as.data.frame(pred_r, xy = TRUE)
names(pred_df) <- c("x", "y", "suitability")
pred_df <- pred_df[!is.na(pred_df$suitability), ]

p_suit <- ggplot() +
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
             shape = 21, size = 1.5, stroke = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_minimal()) +
  #labs(title = paste0("VL habitat suitability \u2014 MaxEnt (",
  #                   sel$fc, ", rm = ", sel$rm, ")"),
  #     subtitle = "Continuous prediction across Sudan") +
  coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 23), crs = 4326) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank())

ggsave(file.path(DIR_FIGS, "maxent_suitability_map.pdf"), p_suit,
       width = 10, height = 8)
ggsave(file.path(DIR_FIGS, "maxent_suitability_map.png"), p_suit,
       width = 10, height = 8, dpi = 300)
cat("Saved suitability map\n")

# ============================ DARFUR ZOOM ===================================

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
  coord_sf(xlim = c(24, 28), ylim = c(11, 16), crs = 4326) +
  labs(title = "Darfur / West Kordofan \u2014 zoomed",
       subtitle = "Presence points over suitability surface") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank())

ggsave(file.path(DIR_FIGS, "suitability_darfur_zoom.png"), p_darfur,
       width = 8, height = 7, dpi = 300)
cat("Saved suitability_darfur_zoom.png\n")

# =============================== MESS ======================================

mess_path <- file.path(DIR_SURFACES, "mess_surface.tif")

if (file.exists(mess_path)) {
  cat("Loading cached MESS surface\n")
  mess_r <- rast(mess_path)
} else {
  ref_data <- rbind(occ_env[, retained_vars], bg_env[, retained_vars])
  covs_sudan <- mask(covs, vect(sudan))
  mess_r <- dismo::mess(x = raster::stack(covs_sudan), v = ref_data, full = FALSE)
  mess_r <- rast(mess_r)

  writeRaster(mess_r, mess_path, overwrite = TRUE)
  cat("Computed and saved MESS surface\n")
}

mess_vals <- values(mess_r, na.rm = TRUE)
cat("MESS summary:\n")
cat("  Range:", round(min(mess_vals), 1), "to", round(max(mess_vals), 1), "\n")
cat("  Novel environments (MESS < 0):",
    sum(mess_vals < 0), "of", length(mess_vals),
    "(", round(100 * sum(mess_vals < 0) / length(mess_vals), 1), "%)\n")

mess_df <- as.data.frame(mess_r, xy = TRUE)
names(mess_df) <- c("x", "y", "mess")
mess_df <- mess_df[!is.na(mess_df$mess), ]
mess_df$type <- ifelse(mess_df$mess < 0, "Extrapolation", "Interpolation")

p_mess <- ggplot() +
  geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
  geom_raster(data = mess_df, aes(x = x, y = y, fill = type)) +
  scale_fill_manual(values = c(Extrapolation = "#D73027",
                               Interpolation = "#4575B4"),
                    name = NULL) +
  geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
  geom_point(data = occ, aes(x = longitude, y = latitude),
             colour = "white", fill = "black",
             shape = 21, size = 1.5, stroke = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 23), crs = 4326) +
  labs(title = "MESS: extrapolation risk",
       subtitle = "Red = at least one covariate outside training range") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank())

ggsave(file.path(DIR_FIGS, "mess_extrapolation.png"), p_mess,
       width = 10, height = 8, dpi = 300)
cat("Saved mess_extrapolation.png\n")

cat("07_model_visualization.R complete\n")

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
library(sf)
source(here::here("R", "plotting_theme.R"))
 
# ---- State boundaries (GADM level 1) ----
adm1   <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))
states <- st_as_sf(adm1)
 
# ---- Occurrences as sf (for geom_sf consistency) ----
occ_sf <- st_as_sf(occ, coords = c("longitude", "latitude"), crs = 4326)
 
# ---- Clip MESS to Sudan boundary  ----
mess_r_clipped <- mask(mess_r, vect(sudan))
mess_df_clipped <- as.data.frame(mess_r_clipped, xy = TRUE)
names(mess_df_clipped) <- c("x", "y", "mess")
mess_df_clipped <- mess_df_clipped[!is.na(mess_df_clipped$mess), ]
mess_df_clipped$type <- ifelse(mess_df_clipped$mess < 0,
                               "Extrapolation", "Interpolation")
 
# ---- Shared extent ----
map_xlim <- c(21.5, 39)
map_ylim <- c(8, 24.5)
 
 
# ---------- Panel (a): Suitability surface -----------------------------------
 
p_suit_a <- ggplot() +
  geom_sf(data = sudan, fill = "grey95", colour = NA) +
  geom_raster(data = pred_df, aes(x = x, y = y, fill = suitability)) +
  scale_fill_suitability() +
  guides(fill = guide_colourbar(
    barheight = unit(0.4, "cm"),
    barwidth  = unit(3, "cm"),
    title.position = "left",
    title.theme = element_text(size = 7, face = "plain", vjust = 0.8),
    label.theme = element_text(size = 6)
  )) +
  layer_admin1(data = states, colour = "black", linewidth = 0.15) +
  layer_country(data = sudan, colour = "black", linewidth = 0.3) +
  geom_sf(data = occ_sf, shape = 21, size = 1.0, stroke = 0.3,
          fill = "white", colour = "black") +
  ggspatial::annotation_scale(
    location = "tl", width_hint = 0.15, text_cex = 0.5,
    line_width = 0.3, pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm")
  ) +
  labs(title = "(a) Predicted suitability") +
  coord_sf(xlim = map_xlim, ylim = map_ylim, crs = 4326, expand = FALSE) +
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
  scale_fill_manual(values = pal_mess_binary, name = NULL) +
  guides(fill = guide_legend(
    keywidth  = unit(0.5, "cm"),
    keyheight = unit(0.4, "cm"),
    direction = "horizontal"
  )) +
  layer_admin1(data = states, colour = "black", linewidth = 0.15) +
  layer_country(data = sudan, colour = "black", linewidth = 0.3) +
  geom_sf(data = occ_sf, shape = 21, size = 1.0, stroke = 0.3,
          fill = "white", colour = "black") +
  labs(title = "(b) MESS analysis") +
  coord_sf(xlim = map_xlim, ylim = map_ylim, crs = 4326, expand = FALSE) +
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
  scale_fill_suitability() +
  layer_admin1(data = states, colour = "black", linewidth = 0.15) +
  layer_country(data = sudan, colour = "black", linewidth = 0.3) +
  geom_sf(data = occ_sf, shape = 21, size = 1.5, stroke = 0.4,
          fill = "white", colour = "black") +
  add_scalebar() +
  add_north_arrow() +
  coord_sf(xlim = map_xlim, ylim = map_ylim, crs = 4326, expand = FALSE) +
  theme_map()
 
save_fig(file.path(DIR_FIGS, "fig_suitability.png"), p_suit_full,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)
save_fig(file.path(DIR_FIGS, "fig_suitability.pdf"), p_suit_full,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)
cat("Saved fig_suitability (standalone)\n")

# ---------- Response curves + variable importance (2x3 grid) -----------------
# Five response curves fill positions 1-5; variable importance fills the 6th
# slot (bottom-right).
# Objects needed: response_data, perm_df, var_labels (from earlier in script)

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
  lab <- var_labels[var]
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
  mutate(var_label = factor(var_label, levels = rev(perm_df$var_label))) |>
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