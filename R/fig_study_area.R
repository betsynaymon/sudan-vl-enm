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

save_fig(file.path(DIR_FIGS, "fig_study_area.png"), fig_study,
         width = FIG_WIDTH_FULL, height = 14)
save_fig(file.path(DIR_FIGS, "fig_study_area.pdf"), fig_study,
         width = FIG_WIDTH_FULL, height = 14)
cat("Saved fig_study_area\n")