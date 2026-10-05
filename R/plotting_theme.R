# ============================================================================
# plotting_theme.R
# Centralized styling for dissertation figures
# ============================================================================

library(ggplot2)
library(scales)

# ----------------------------------------------------------------------------
# 1. EXPORT DIMENSIONS
#    All figures saved at the same base width → consistent text sizes.    
# ----------------------------------------------------------------------------

FIG_WIDTH_FULL   <- 16      # cm, full text width (single-column A4)
FIG_WIDTH_HALF   <- 7.5     # cm, for side-by-side panels
FIG_HEIGHT_MAP   <- 14      # cm, portrait map 
FIG_HEIGHT_PLOT  <- 10      # cm, standard non-map figure
FIG_DPI          <- 300     # print quality

# Usage: save_fig("figures/suitability_map.png", p, width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)
save_fig <- function(filename, plot = last_plot(),
                     width = FIG_WIDTH_FULL,
                     height = FIG_HEIGHT_PLOT,
                     dpi = FIG_DPI,
                     bg = "white",
                     ...) {
  # PDFs through Cairo: R's pdf() device mis-clipped ggpattern's hatching on
  # multi-part polygons (fig_endemic_status, Red Sea); Cairo matches the PNGs.
  dev <- if (grepl("\\.pdf$", filename, ignore.case = TRUE)) grDevices::cairo_pdf else NULL
  ggsave(filename = filename, plot = plot, width = width, height = height,
         units = "cm", dpi = dpi, bg = bg, device = dev, ...)
}


# ----------------------------------------------------------------------------
# 2. TYPOGRAPHY
#    Base size calibrated so axis text is legible at FIG_WIDTH_FULL / 300 dpi.
# ----------------------------------------------------------------------------

FONT_FAMILY <- ""           # empty string = ggplot2 default (sans)
BASE_SIZE   <- 10           # axis tick labels
TITLE_SIZE  <- 12           # panel/facet strip titles
CAPTION_SIZE <- 8           # source notes, figure captions


# ----------------------------------------------------------------------------
# 3. COLOUR PALETTES
# ----------------------------------------------------------------------------

# — Suitability surface 
pal_suitability <- c(
    "#2166AC",   # deep blue  (0.00)
    "#67A9CF",   # mid blue   (0.25)
    "#F7F7F7",   # near-white (0.50) 
    "#EF8A62",   # salmon     (0.75)
    "#B2182B"    # deep red   (1.00)
)

# Use in maps:
#   scale_fill_gradientn(colours = pal_suitability, limits = c(0, 1),
#                        name = "Habitat\nsuitability")
scale_fill_suitability <- function(name = "Habitat\nsuitability",
                                   limits = c(0, 1), ...) {
  scale_fill_gradientn(
    colours = pal_suitability,
    limits  = limits,
    name    = name,
    ...
  )
}

# — MESS: diverging palette centered at 0 (interpolation vs. extrapolation) —
pal_mess <- c(
    "#B2182B",   # strong extrapolation (negative)
    "#FDDBC7",   # mild extrapolation
    "#F7F7F7",   # zero
    "#D1E5F0",   # mild interpolation
    "#2166AC"    # strong interpolation (positive)
)

# — MESS binary: just two classes —
pal_mess_binary <- c(
    "Extrapolation"  = "#D6604D",
    "Interpolation"  = "#4393C3"
)

# — Algorithms (17) —
pal_models <- c(
    maxent = "#4575B4",    # steel blue
    rf     = "#A50026",    # dark red
    gbt    = "#E66101"     # orange
)
model_labels <- c(maxent = "MaxEnt", rf = "Random forest", gbt = "Boosted trees")

# — Occurrence records by evidence type (fig_study_area), keyed by presence_type —
pal_occurrences <- c(
    vl_case         = "#332288",   # indigo
    facility_report = "#44AA99",   # teal
    vector_pos      = "#CC6677"    # rose
)
occ_type_labels <- c(vl_case = "Human VL case", facility_report = "Facility report",
                     vector_pos = "Vector-positive site")

# — Endemic status by state (fig_study_area), keyed by STATE_STATUS_FILE's status —
# Both non-core classes are white so grey means only "not modelled"; hatching
# marks "Reported".
endemic_labels   <- c("Core endemic"  = "Core endemic",
                      "Reported"      = "Reported cases or foci",
                      "No documented" = "No documented cases or foci")
pal_endemic      <- c("Core endemic" = "#4393C3", "Reported" = "white", "No documented" = "white")
endemic_patterns <- c("Core endemic" = "none", "Reported" = "stripe", "No documented" = "none")
HATCH_COL       <- "#4393C3"   # stripes for "Reported"
HATCH_ANGLE     <- 45
HATCH_SPACING   <- 0.005       # npc units: smaller = denser
HATCH_DENSITY   <- 0.28        # share of area covered by stripes
HATCH_SIZE      <- 0.01        # stripe line width
HATCH_KEY_SCALE <- 0.7         # pattern scale in the legend key

# — State names on maps where they differ from GADM's NAME_1 —
state_labels <- c("Al Qadarif" = "Gedaref", "Al Jazirah" = "Gezira",
                  "North Kurdufan" = "North Kordofan", "South Kurdufan" = "South Kordofan",
                  "West Kurdufan" = "West Kordofan")

# — Locator inset —
col_locator <- "#B2182B"
INSET_XLIM  <- c(-18, 55)
INSET_YLIM  <- c(-5, 38)
INSET_BOX   <- c(left = 0.02, bottom = 0.70, right = 0.26, top = 0.98)   # share of the figure

# — Accessibility comparison (11) —
pal_access <- c(
    "Presences"                     = "#C4666D",   # muted rose
    "Background"                    = "#969696",   # grey
    "Background, in presence range" = "#4575B4"    # steel blue
)

# — Accessibility-corrected backgrounds (12) —
pal_bias <- c(
    uniform = "#969696",   # grey
    half    = "#91BFDB",   # light blue
    matched = "#4575B4"    # steel blue
)

# — Threshold reference lines —
col_p10    <- "#4575B4"   # blue, liberal threshold
col_maxsss <- "#A50026"   # dark red, conservative threshold
col_observed <- "#A50026"   # observed statistic against a null distribution (13)

# — Covariate labels (07, 24, 15) —
cov_labels <- c(
  slope         = "Slope (degrees)",
  river_dist    = "Distance to river (m)",
  vertisols     = "Vertisols (0/1)",
  lst_night     = "LST night (\u00b0C)",
  rainfall      = "Rainfall (mm/yr)",
  travel_time   = "Travel time (min)",
  elevation     = "Elevation (m)",
  lst_day       = "LST day (\u00b0C)",
  ndvi          = "NDVI",
  treecover     = "Tree cover (%)",
  lst_night_dry = "LST night, dry season (\u00b0C)",
  lst_night_wet = "LST night, wet season (\u00b0C)"
)

# — Covariate variants (15) —
pal_variants <- c(
  primary        = "grey30",
  lst_dry        = "#B15928",
  lst_wet        = "#5E3C99",
  ndvi_for_rain  = "#1B9E77",
  plus_lst_day   = "#D95F02",
  plus_treecover = "#7570B3",
  plus_elevation = "#E7298A",
  minus_river    = "#A6761D"
)

# — Accessibility covariate test (24) —
pal_tt <- c(primary = "#969696", augmented = "#4575B4")

# — Darfur vs the rest (20) —
pal_region <- c("Rest of Sudan" = "#B2182B", "Darfur" = "#2166AC")

# ----------------------------------------------------------------------------
# 4. BASE THEME (non-map figures)
#    Histogram, density, line plots, PDPs, threshold curve, etc.
# ----------------------------------------------------------------------------

theme_dissertation <- function(base_size = BASE_SIZE,
                               base_family = FONT_FAMILY,
                               gridlines = "y") {
  # Start from theme_minimal for a clean canvas
  t <- theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # — Text hierarchy —
      plot.title         = element_text(size = TITLE_SIZE, face = "plain",
                                        hjust = 0, margin = margin(b = 4)),
      plot.subtitle      = element_text(size = base_size, colour = "grey40",
                                        hjust = 0, margin = margin(b = 8)),
      plot.caption       = element_text(size = CAPTION_SIZE, colour = "grey50",
                                        hjust = 1),
      axis.title         = element_text(size = base_size, colour = "grey20"),
      axis.text          = element_text(size = base_size - 1, colour = "grey30"),
      strip.text         = element_text(size = base_size, face = "plain",
                                        hjust = 0),

      # — Panel —
      panel.background   = element_rect(fill = "white", colour = NA),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.border       = element_blank(),

      # — Gridlines: default to horizontal only (y); override per-figure —
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.major.x = if (gridlines == "both") {
                              element_line(colour = "grey90", linewidth = 0.3)
                            } else {
                              element_blank()
                            },
      panel.grid.minor   = element_blank(),

      # — Axes —
      axis.ticks         = element_line(colour = "grey70", linewidth = 0.3),
      axis.ticks.length  = unit(2, "pt"),
      axis.line          = element_blank(),

      # — Legend —
      legend.position    = "right",
      legend.title       = element_text(size = base_size, face = "bold"),
      legend.text        = element_text(size = base_size - 1),
      legend.key.size    = unit(0.9, "lines"),
      legend.background  = element_blank(),

      # — Spacing —
      plot.margin        = margin(8, 8, 8, 8, unit = "pt"),
      panel.spacing      = unit(0.8, "lines")
    )
  t
}


# ----------------------------------------------------------------------------
# 5. MAP THEME
#    Use with ggspatial::annotation_scale() and annotation_north_arrow().
# ----------------------------------------------------------------------------

theme_map <- function(base_size = BASE_SIZE,
                      base_family = FONT_FAMILY) {
  theme_dissertation(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Kill all axes — maps get scale bar + north arrow instead
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Legend: inside the map panel to save space
      legend.position      = c(0.98, 0.25),
      legend.justification = c(1, 0),
      legend.background    = element_rect(fill = alpha("white", 0.85),
                                          colour = NA),
      legend.key.height    = unit(1.2, "lines"),
      legend.key.width     = unit(0.6, "lines"),

      # Tighter margins for maps
      plot.margin = margin(4, 4, 4, 4, unit = "pt")
    )
}


# ----------------------------------------------------------------------------
# 6. MAP FURNITURE HELPERS
#    Wrappers around ggspatial to not repeat style args.
#    Require: library(ggspatial)
# ----------------------------------------------------------------------------
INSET_BOX   <- c(left = 0.01, bottom = 0.86, right = 0.19, top = 0.99)   # share of the map panel
# Scale bar — bottom-left by default
add_scalebar <- function(location = "bl", width_hint = 0.2, ...) {
  ggspatial::annotation_scale(
    location   = location,
    width_hint = width_hint,
    style      = "ticks",
    text_cex   = 0.7,
    line_width = 0.4,
    ...
  )
}

# North arrow — top-right by default, minimal style
add_north_arrow <- function(location = "tr", ...) {
  ggspatial::annotation_north_arrow(
    location = location,
    which_north = "true",
    height = unit(0.8, "cm"),
    width  = unit(0.6, "cm"),
    style  = ggspatial::north_arrow_minimal(),
    ...
  )
}


# ----------------------------------------------------------------------------
# 7. BOUNDARY LAYERS
#    Helpers for the admin-boundary line weights.
#    These return geom_sf layers.
# ----------------------------------------------------------------------------

# Country outline: heavier stroke
layer_country <- function(data, colour = "black", linewidth = 0.3, ...) {
  geom_sf(data = data, fill = NA, colour = colour, linewidth = linewidth, ...)
}

# State/province boundaries: thinner
layer_admin1 <- function(data, colour = "black", linewidth = 0.15, ...) {
  geom_sf(data = data, fill = NA, colour = colour, linewidth = linewidth, ...)
}

# Occurrence points: white fill, black stroke
layer_occurrences <- function(data, mapping = aes(), size = 1.8,
                              stroke = 0.4, shape = 21,
                              fill = "white", colour = "black", ...) {
  geom_sf(data = data, mapping = mapping,
          shape = shape, size = size, stroke = stroke,
          fill = fill, colour = colour, ...)
}

# Area outside the modelled domain (Halaib Triangle): grey, unlabelled
# (figure captions identify it)
layer_excluded <- function(data) {
  ggplot2::geom_sf(data = data, fill = "grey80", colour = NA)
}

# Map extent: the full display outline plus a margin in degrees (07, 12, 17,
# fig_study_area and the remaining map scripts).
MAP_MARGIN_DEG <- 0.5
display_limits <- function(display, margin = MAP_MARGIN_DEG) {
  bb <- sf::st_bbox(display)
  list(x = c(bb[["xmin"]], bb[["xmax"]]) + c(-margin, margin),
       y = c(bb[["ymin"]], bb[["ymax"]]) + c(-margin, margin))
}
coord_display <- function(lim) coord_sf(xlim = lim$x, ylim = lim$y, crs = 4326, expand = FALSE)
# ----------------------------------------------------------------------------
# 8. SET DEFAULTS
#    Automatically applied when this file is sourced.
# ----------------------------------------------------------------------------

theme_set(theme_dissertation())

# Continuous fill default → suitability palette
# (Override per-figure when you need MESS or something else)
update_geom_defaults("bar",  list(fill = "grey60", colour = NA))
update_geom_defaults("col",  list(fill = "grey60", colour = NA))
update_geom_defaults("line", list(linewidth = 0.7))
update_geom_defaults("point", list(size = 1.5))

message("dissertation theme loaded — theme_dissertation() set as default")
