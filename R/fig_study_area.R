# ============================================================================
# fig_study_area.R
# Two maps of Sudan as administered (01):
#   fig_study_area      states and the occurrence records by evidence type
#   fig_endemic_status  states by literature-classified VL endemic status
#                       (STATE_STATUS_FILE, as 16), with a locator inset
# The Halaib Triangle is grey and unlabelled.
#
# Inputs:  ADM1_FILE, EXCLUDED_FILE, DISPLAY_ADM1_FILE, DISPLAY_ADM0_FILE, OCC_FILE, TRAIN_FILE,
#          STATE_STATUS_FILE; Natural Earth countries (rnaturalearth) for the inset
# Outputs: outputs/figures/fig_study_area.png / .pdf,
#          fig_endemic_status.png / .pdf
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(ggplot2)
  library(ggpattern); library(ggrepel); library(patchwork); library(rnaturalearth)
})

# ------------------------------ Load inputs ---------------------------------
display     <- st_as_sf(vect(DISPLAY_ADM0_FILE))
excluded    <- st_as_sf(vect(EXCLUDED_FILE))
states      <- st_as_sf(vect(ADM1_FILE))           # modelled states (study-area map)
states_full <- st_as_sf(vect(DISPLAY_ADM1_FILE))   # GADM states, Triangle in Red Sea (status map)
occ         <- read.csv(OCC_FILE)
train       <- readRDS(TRAIN_FILE)
lit         <- read.csv(STATE_STATUS_FILE)
africa      <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")
lim         <- display_limits(display)

stopifnot(
  "coordinate_id is not unique in OCC_FILE" = !anyDuplicated(occ$coordinate_id),
  "OCC_FILE and the training presences differ (rerun 06)" =
    nrow(occ) == nrow(train$occ_clean) &&
    setequal(occ$coordinate_id, train$occ_clean$coordinate_id),
  "A presence type has no entry in occ_type_labels (plotting_theme.R)" =
    all(occ$presence_type %in% names(occ_type_labels)),
  "State names differ across STATE_STATUS_FILE, ADM1_FILE and DISPLAY_ADM1_FILE" =
    setequal(lit$state, states_full$NAME_1) && !anyDuplicated(lit$state) &&
    setequal(states_full$NAME_1, states$NAME_1),
  "Unknown status in STATE_STATUS_FILE" = all(lit$status %in% names(endemic_labels))
)

# ---------------------------- Study area map --------------------------------
# Records by evidence type. Repeat years share a location, so points
# overplot: the most common type is drawn first and the rarest on top.
# Counts for the caption are printed.

n_type   <- table(factor(occ$presence_type, levels = names(occ_type_labels)))
type_lab <- unname(sprintf("%s (%d)", occ_type_labels, as.integer(n_type)))
loc      <- paste(occ$longitude, occ$latitude)
cat("Records:", nrow(occ), "at", n_distinct(loc), "locations |",
    paste(sprintf("%s %d", names(n_type), n_type), collapse = ", "), "\n")
cat("Locations with more than one evidence type:",
    sum(tapply(occ$presence_type, loc, n_distinct) > 1), "\n")

occ_sf <- occ[order(-as.integer(n_type[occ$presence_type])), ] |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

p_area <- ggplot() +
  geom_sf(data = states, fill = "grey95", colour = NA) +
  layer_excluded(excluded) +
  layer_admin1(states) +
  layer_country(display) +
  geom_sf(data = occ_sf, aes(fill = presence_type), shape = 21, size = 2,
          stroke = 0.3, colour = "black") +
  scale_fill_manual(values = pal_occurrences, breaks = names(occ_type_labels),
                    labels = type_lab, name = NULL) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3))) +
  add_scalebar(location = "tl", width_hint = 0.15) +
  coord_display(lim) +
  theme_map() +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.text = element_text(size = 7), legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0))

# ------------------------- Endemic status map --------------------------------
# Status as in 16 (STATE_STATUS_FILE); map names from state_labels.

cat("States by status:",
    paste(sprintf("%s %d", names(endemic_labels),
                  as.integer(table(factor(lit$status, levels = names(endemic_labels))))),
          collapse = " | "),
    "| without a source:", sum(lit$source %in% c("", NA)), "\n")

st_cls <- states_full |>
  left_join(lit, by = c("NAME_1" = "state")) |>
  mutate(status = factor(status, levels = names(endemic_labels)),
         label  = coalesce(unname(state_labels[NAME_1]), NAME_1))

# point_on_surface sits inside concave polygons; the lon/lat warning is
# irrelevant for label placement
xy     <- suppressWarnings(st_coordinates(st_point_on_surface(st_geometry(st_cls))))
lab_df <- data.frame(label = st_cls$label, X = xy[, "X"], Y = xy[, "Y"])

p_status <- ggplot() +
  geom_sf_pattern(
    data = st_cls, aes(fill = status, pattern = status),
    colour = "grey25", linewidth = 0.15,
    pattern_colour = HATCH_COL, pattern_fill = HATCH_COL, pattern_angle = HATCH_ANGLE,
    pattern_spacing = HATCH_SPACING, pattern_density = HATCH_DENSITY,
    pattern_size = HATCH_SIZE, pattern_key_scale_factor = HATCH_KEY_SCALE) +
  scale_fill_manual(values = pal_endemic, breaks = names(endemic_labels),
                    labels = unname(endemic_labels), name = NULL) +
  scale_pattern_manual(values = endemic_patterns, breaks = names(endemic_labels),
                       labels = unname(endemic_labels), name = NULL) +
  layer_country(display) +
  geom_text_repel(
    data = lab_df, aes(X, Y, label = label), size = 2.5, colour = "grey10",
    bg.color = "white", bg.r = 0.14, segment.colour = "grey40", segment.size = 0.2,
    min.segment.length = 0.3, box.padding = 0.16, point.padding = 0,
    force = 1.2, max.overlaps = Inf, seed = SEED) +
  add_scalebar(location = "br", width_hint = 0.15) +
  coord_display(lim) +
  theme_map() +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.text = element_text(size = 7), legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0)) +
  guides(fill = guide_legend(nrow = 1), pattern = guide_legend(nrow = 1))

p_inset <- ggplot() +
  geom_sf(data = africa, fill = "grey88", colour = "white", linewidth = 0.2) +
  geom_sf(data = display, fill = col_locator, colour = "black", linewidth = 0.3) +
  coord_sf(xlim = INSET_XLIM, ylim = INSET_YLIM, crs = 4326) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4))

fig_status <- p_status +
  inset_element(p_inset, left = INSET_BOX[["left"]], bottom = INSET_BOX[["bottom"]],
                right = INSET_BOX[["right"]], top = INSET_BOX[["top"]])

# --------------------------------- Save -------------------------------------

for (ext in c("png", "pdf")) {
  save_fig(file.path(DIR_FIGS, paste0("fig_study_area.", ext)), p_area,
           width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)
  save_fig(file.path(DIR_FIGS, paste0("fig_endemic_status.", ext)), fig_status,
           width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)
}
cat("fig_study_area.R complete\n")