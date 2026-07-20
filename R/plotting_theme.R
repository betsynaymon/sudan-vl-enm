# ============================================================================
# plotting_theme.R
# Centralized styling for dissertation figures
# Source this at the top of every figure script:
#   source("plotting_theme.R")
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

# Wrapper: enforces consistent dimensions and device settings.
# Usage: save_fig("figures/suitability_map.png", p, width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)
save_fig <- function(filename, plot = last_plot(),
                     width = FIG_WIDTH_FULL,
                     height = FIG_HEIGHT_PLOT,
                     dpi = FIG_DPI,
                     bg = "white",
                     ...) {
  ggsave(
    filename = filename,
    plot     = plot,
    width    = width,
    height   = height,
    units    = "cm",
    dpi      = dpi,
    bg       = bg,
    ...
  )
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
    "#F7F7F7",   # near-white (0.50) — keeps the "low = cool" reading
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

# — Model comparison (PDP, robustness) —
# Named so you can use scale_colour_manual(values = pal_models)
pal_models <- c(
    "MaxEnt" = "#4575B4",    # steel blue
    "RF"     = "#A50026",    # dark red
    "GBT"    = "#E66101"     # orange
)

# — Occurrence types (for study area map) —
pal_occurrences <- c(
    "Human VL Case"      = "#332288",   # indigo
    "Facility Report"   = "#44AA99",   # teal
    "Vector"        = "#CC6677"   # rose
)

# — Diagnostic / two-group plots (accessibility bias, etc.) —
pal_two_group <- c(
    "Background"   = "#969696",    # grey
    "Occurrences"  = "#C4666D"     # muted rose 
)

# — Threshold reference lines —
col_p10    <- "#4575B4"   # blue, liberal threshold
col_maxsss <- "#A50026"   # dark red, conservative threshold


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
      plot.title         = element_text(size = TITLE_SIZE, face = "bold",
                                        hjust = 0, margin = margin(b = 4)),
      plot.subtitle      = element_text(size = base_size, colour = "grey40",
                                        hjust = 0, margin = margin(b = 8)),
      plot.caption       = element_text(size = CAPTION_SIZE, colour = "grey50",
                                        hjust = 1),
      axis.title         = element_text(size = base_size, colour = "grey20"),
      axis.text          = element_text(size = base_size - 1, colour = "grey30"),
      strip.text         = element_text(size = base_size, face = "bold",
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
#    Strips lat/lon axes and graticule. Keeps legend and essentials.
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
