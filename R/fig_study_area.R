# ============================================================================
# fig_study_area.R
# Study area map: Sudan with ecological mask shown as rainfall texture,
# occurrence points by type, admin boundaries, and Africa context inset.
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ggplot2)
  library(geodata)
  library(rnaturalearth)
  library(patchwork)
  library(ggspatial)
  library(dplyr)
})

# ------------------------------ Load data ------------------------------------

adm0   <- gadm(country = "SDN", level = 0, path = here::here("data", "raw"))
sudan  <- st_as_sf(adm0)
adm1   <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))
states <- st_as_sf(adm1)

# Ecological mask -> polygon for dashed outline
mask_r <- rast(here::here("data", "raw", "ecological_mask_150mm.tif"))
mask_r <- mask(mask_r, vect(sudan))
mask_poly <- st_as_sf(as.polygons(mask_r, dissolve = TRUE))
mask_col <- names(mask_poly)[1]
mask_poly <- mask_poly[mask_poly[[mask_col]] == 1, ]

# Occurrences
occ <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
occ$type_label <- case_when(
  occ$presence_type == "vl_case"          ~ "Human VL Case",
  occ$presence_type == "facility_report"  ~ "Facility Report",
  occ$presence_type == "vector_pos"       ~ "Vector"
)
occ$type_label <- factor(occ$type_label,
                         levels = c("Human VL Case", "Facility Report", "Vector"))
occ_sf <- st_as_sf(occ, coords = c("longitude", "latitude"), crs = 4326)

# Africa countries for inset
africa <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")

cat("Data loaded | Occurrences:", nrow(occ_sf), "\n")

# ------------------------------ Main map -------------------------------------

map_xlim <- c(21.5, 38.5)
map_ylim <- c(8, 24.5)

p_main <- ggplot() +
  # Base country fill (grey outside mask)
  geom_sf(data = sudan, fill = "grey92", colour = NA) +
  # Ecological mask (plain fill + dashed outline)
  geom_sf(data = mask_poly, fill = "#C8D6A0", colour = "grey30",
          linetype = "dashed", linewidth = 0.25, alpha = 0.6) +
  # Admin boundaries
  geom_sf(data = states, fill = NA, colour = "black", linewidth = 0.15) +
  geom_sf(data = sudan, fill = NA, colour = "black", linewidth = 0.3) +
  # Occurrence points (only these use fill -> clean legend)
  geom_sf(data = occ_sf, aes(fill = type_label),
          shape = 21, size = 2, stroke = 0.3, colour = "black") +
  scale_fill_manual(values = pal_occurrences, name = NULL) +
  guides(fill = guide_legend(
    direction = "horizontal", nrow = 1,
    override.aes = list(size = 3)
  )) +
  # Scale bar
  ggspatial::annotation_scale(
    location = "tl", width_hint = 0.15, text_cex = 0.5,
    line_width = 0.3, pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm")
  ) +
  coord_sf(xlim = map_xlim, ylim = map_ylim, crs = 4326, expand = FALSE) +
  theme_map() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0)
  )

# ------------------------------ Africa inset ---------------------------------

p_inset <- ggplot() +
  geom_sf(data = africa, fill = "grey88", colour = "white", linewidth = 0.2) +
  geom_sf(data = sudan, fill = "#B2182B", colour = "black", linewidth = 0.3) +
  coord_sf(xlim = c(-18, 55), ylim = c(-5, 38), crs = 4326) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4)
  )

# ------------------------------ Combine --------------------------------------

fig_study <- p_main +
  inset_element(p_inset,
                left = 0.01, bottom = 0.58,
                right = 0.28, top = 0.95)

save_fig(file.path(DIR_FIGS, "fig_study_area.png"), p_main,
          width = FIG_WIDTH_FULL, height = 14)
save_fig(file.path(DIR_FIGS, "fig_study_area.pdf"), p_main,
          width = FIG_WIDTH_FULL, height = 14)

save_fig(file.path(DIR_FIGS, "fig_study_area_inset.png"), fig_study,
         width = FIG_WIDTH_FULL, height = 14)
save_fig(file.path(DIR_FIGS, "fig_study_area_insent.pdf"), fig_study,
         width = FIG_WIDTH_FULL, height = 14)
cat("Saved fig_study_area\n")

############# SUDAN MAP ###############
# Literature-classified VL endemic status by state, with Africa context inset.
#   Core endemic             -> solid fill
#   Reported cases/foci      -> hatched (45 deg stripe) over light fill
#   No documented cases/foci -> plain grey
# ============================================================================
 
source(here::here("R", "params.R"))
source(here::here("R", "plotting_theme.R"))
 
suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(geodata)
  library(rnaturalearth)
  library(ggpattern)      # hatching / dot fills for sf polygons
  library(ggrepel)        # non-overlapping state labels
  library(patchwork)
  library(dplyr)
  library(tibble)
})
 
# ---------------------------- Figure constants -------------------------------
 
# Endemic-status palette. Same hue family as pal_suitability so the figure reads
# alongside the suitability map. Move into plotting_theme.R if reused.
pal_endemic <- c(
  "Core endemic"             = "#4393C3",   
  "Reported cases/foci"      = "grey88",   # light salmon base under hatching
  "No documented cases/foci" = "grey88"
)
 
HATCH_COL     <- "#4393C3"   # stripe colour for the "Reported" class
HATCH_ANGLE   <- 45
HATCH_SPACING <- 0.005       # npc units: smaller = denser hatching
HATCH_DENSITY <- 0.28        # fraction of area covered by stripes
HATCH_SIZE    <- 0.01        # stripe line width
KEY_SCALE     <- 0.7         # pattern scaling inside the legend key
 
map_xlim <- c(21.5, 38.5)
map_ylim <- c(8, 24.5)
 
# --------------------------- State classification ----------------------------
# gadm_name must match GADM NAME_1 exactly; label is what appears on the map.
status_levels <- c("Core endemic", "Reported cases/foci",
                   "No documented cases/foci")
 
state_status <- tribble(
  ~gadm_name,       ~label,           ~status,
  "Al Qadarif",     "Gedaref",        "Core endemic",
  "Sennar",         "Sennar",         "Core endemic",
  "Blue Nile",      "Blue Nile",      "Core endemic",
  "White Nile",     "White Nile",     "Core endemic",
  "Kassala",        "Kassala",        "Core endemic",
  "South Kurdufan", "South Kordofan", "Reported cases/foci",
  "West Kurdufan",  "West Kordofan",  "Reported cases/foci",
  "North Kurdufan", "North Kordofan", "Reported cases/foci",
  "East Darfur",    "East Darfur",    "Reported cases/foci",
  "South Darfur",   "South Darfur",   "Reported cases/foci",
  "West Darfur",    "West Darfur",    "Reported cases/foci",
  "North Darfur",   "North Darfur",   "Reported cases/foci",
  "Khartoum",       "Khartoum",       "Reported cases/foci",
  "Al Jazirah",     "Gezira",         "Reported cases/foci",
  "Red Sea",        "Red Sea",        "Reported cases/foci",
  "Northern",       "Northern",       "No documented cases/foci",
  "River Nile",     "River Nile",     "No documented cases/foci",
  "Central Darfur", "Central Darfur", "No documented cases/foci"
) |>
  mutate(status = factor(status, levels = status_levels))
 
# Per-state label offsets in decimal degrees. 
label_nudge <- tribble(
  ~gadm_name,   ~nudge_x, ~nudge_y,
  "Khartoum",        0.0,      0.0
)
 
# ------------------------------ Load data ------------------------------------
 
adm0  <- gadm(country = "SDN", level = 0, path = here::here("data", "raw"))
sudan <- st_as_sf(adm0)
adm1  <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))
states <- st_as_sf(adm1)
 
# Verify the classification covers GADM exactly before plotting.
missing_in_gadm <- setdiff(state_status$gadm_name, states$NAME_1)
unclassified    <- setdiff(states$NAME_1, state_status$gadm_name)
 
if (length(missing_in_gadm) > 0) {
  stop("Classified names absent from GADM NAME_1: ",
       paste(missing_in_gadm, collapse = ", "),
       "\nGADM returned: ", paste(sort(states$NAME_1), collapse = ", "))
}
if (length(unclassified) > 0) {
  stop("GADM states with no classification: ",
       paste(unclassified, collapse = ", "),
       "\nAdd them to state_status (e.g. a disputed-territory polygon).")
}
 
states_cls <- states |>
  left_join(state_status, by = c("NAME_1" = "gadm_name")) |>
  left_join(label_nudge, by = c("NAME_1" = "gadm_name")) |>
  mutate(nudge_x = coalesce(nudge_x, 0), nudge_y = coalesce(nudge_y, 0))
 
cat("States classified:",
    paste(names(table(states_cls$status)), table(states_cls$status),
          sep = " = ", collapse = " | "), "\n")
 
# Label anchors: point_on_surface sits inside concave polygons, unlike centroid.
# Warning about lon/lat is expected and irrelevant for label placement.
lab_xy <- suppressWarnings(
  st_coordinates(st_point_on_surface(st_geometry(states_cls)))
)
lab_df <- data.frame(
  label = states_cls$label,
  X     = lab_xy[, "X"] + states_cls$nudge_x,
  Y     = lab_xy[, "Y"] + states_cls$nudge_y
)
 
# Africa countries for the inset
africa <- ne_countries(scale = "medium", continent = "Africa",
                       returnclass = "sf")
 
# ------------------------------ Main map -------------------------------------
 
p_main <- ggplot() +
  # States: fill by status, hatching applied only to the "Reported" class
  geom_sf_pattern(
    data = states_cls,
    aes(fill = status, pattern = status),
    colour                    = "grey25",
    linewidth                 = 0.15,
    pattern_colour            = HATCH_COL,
    pattern_fill              = HATCH_COL,
    pattern_angle             = HATCH_ANGLE,
    pattern_spacing           = HATCH_SPACING,
    pattern_density           = HATCH_DENSITY,
    pattern_size              = HATCH_SIZE,
    pattern_key_scale_factor  = KEY_SCALE
  ) +
  scale_fill_manual(values = pal_endemic, breaks = status_levels, name = NULL) +
  scale_pattern_manual(
    values = c("Core endemic"             = "none",
               "Reported cases/foci"      = "stripe",   # or "circle" for dots
               "No documented cases/foci" = "none"),
    breaks = status_levels, name = NULL
  ) +
  # Country outline on top, heavier stroke
  geom_sf(data = sudan, fill = NA, colour = "black", linewidth = 0.3) +
  # State labels with a white halo so they read over the solid dark fills
  geom_text_repel(
    data = lab_df, aes(x = X, y = Y, label = label),
    size = 2.5, colour = "grey10",
    bg.color = "white", bg.r = 0.14,
    segment.colour = "grey40", segment.size = 0.2,
    min.segment.length = 0.3, box.padding = 0.16, point.padding = 0,
    force = 1.2, max.overlaps = Inf, seed = 1238
  ) +
  add_scalebar(location = "br", width_hint = 0.15,
               pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm")) +
  coord_sf(xlim = map_xlim, ylim = map_ylim, crs = 4326, expand = FALSE) +
  theme_map() +
  theme(
    legend.position      = "bottom",
    legend.justification = "center",
    legend.text          = element_text(size = 7),
    legend.background    = element_blank(),
    legend.margin        = margin(0, 0, 0, 0)
  ) +
  guides(
    fill    = guide_legend(nrow = 1, direction = "horizontal"),
    pattern = guide_legend(nrow = 1, direction = "horizontal")
  )
 
# ------------------------------ Africa inset ---------------------------------

 
p_inset <- ggplot() +
  geom_sf(data = africa, fill = "grey88", colour = "white", linewidth = 0.2) +
  geom_sf(data = sudan, fill = "#B2182B", colour = "black", linewidth = 0.3) +
  coord_sf(xlim = c(-18, 55), ylim = c(-5, 38), crs = 4326) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA,
                                    linewidth = 0.4)
  )
 
# ------------------------------ Combine --------------------------------------
# Inset sits over the empty north-west corner. Raise `bottom` to keep it clear
# of Northern state; lower it to make the inset larger.
 
fig_status <- p_main +
  inset_element(p_inset,
                left = 0.02, bottom = 0.70,
                right = 0.26, top = 0.98)
 
save_fig(file.path(DIR_FIGS, "fig_endemic_status.png"), fig_status,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)
save_fig(file.path(DIR_FIGS, "fig_endemic_status.pdf"), fig_status,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)
 
cat("Saved fig_endemic_status\n")

# ============================================================================
# THREE ALGORITHM MAP
# Rebuilds the three-way suitability comparison from cached surfaces.
#
# Inputs:  outputs/surfaces/{maxent,rf,gbt}_suitability.tif
# Outputs: outputs/figures/suitability_three_models.png
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2)
  library(patchwork); library(ggspatial); library(geodata)
})

suit_mx  <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))
suit_rf  <- rast(file.path(DIR_SURFACES, "rf_suitability.tif"))
suit_gbt <- rast(file.path(DIR_SURFACES, "gbt_suitability.tif"))

sudan  <- st_as_sf(gadm("SDN", level = 0, path = here::here("data", "raw")))
states <- st_as_sf(gadm("SDN", level = 1, path = here::here("data", "raw")))

to_df <- function(r) {
  d <- as.data.frame(mask(r, vect(sudan)), xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "suitability"
  d
}

XLIM <- c(21.5, 39); YLIM <- c(8.5, 22.5)

make_panel <- function(df, label) {
  ggplot() +
    geom_raster(data = df, aes(x, y, fill = suitability)) +
    scale_fill_suitability(guide = guide_colorbar(title.position = "top",
                                                 title.hjust = 0)) +
    layer_admin1(states) +
    layer_country(sudan) +
    coord_sf(xlim = XLIM, ylim = YLIM) +
    labs(title = label) +
    theme_map() +
    theme(legend.position = "none")
}

p_maps <- make_panel(to_df(suit_mx),  "(a) MaxEnt") +
          make_panel(to_df(suit_rf),  "(b) Random Forest") +
          make_panel(to_df(suit_gbt), "(c) Gradient Boosted Trees") +
          add_scalebar(location = "br") +
  plot_layout(guides = "collect") &
  theme(legend.position    = "bottom",
        plot.title         = element_text(size = 9, hjust = 0,
                                          margin = margin(b = 2)),
        legend.key.width   = unit(2, "cm"),
        legend.key.height  = unit(0.25, "cm"),
        legend.title       = element_text(size = 8),
        legend.text        = element_text(size = 7),
        legend.margin      = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        plot.margin        = margin(2, 2, 2, 2, unit = "pt"))

panel_aspect <- diff(YLIM) / (diff(XLIM) * cos(mean(YLIM) * pi / 180))
fig_h <- ((FIG_WIDTH_FULL - 0.5) / 3) * panel_aspect + 1.9

ggsave(file.path(DIR_FIGS, "suitability_three_models.png"), p_maps,
       width = FIG_WIDTH_FULL, height = FIG_WIDTH_HALF, dpi = 300, bg = "white")