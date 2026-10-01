# ============================================================================
# 11_accessibility_bias_diagnostic.R
# How much closer to cities are the presences than the places they are
# compared with? Travel time to the nearest city (Weiss et al. 2018) at the
# presences, by record type, against three references: the background (the
# model's comparison set, which 12 and 13 reweight), background within the
# presences' environmental range (MESS >= 0, which removes the remote
# desert), and the people living within that range. Also reports gaps in the
# travel-time surface, and how suitability co-varies with accessibility (the
# confound 12 and 13 address).
#
# Inputs:  TRAIN_FILE, TT_FILE (downloaded once from TT_URL), ADM0_FILE,
#          ADM1_FILE, DOMAIN_FILE, SUIT_FILE, MESS_FILE, POP_ALIGNED_FILE
# Outputs: outputs/tables/accessibility_comparison.csv
#          outputs/tables/travel_time_gaps_by_state.csv
#          outputs/tables/suitability_travel_time_correlation.csv
#          outputs/figures/accessibility_bias_diagnostic.png
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(ggplot2); library(httr)
})

# ------------------------------ Load inputs ---------------------------------

train  <- readRDS(TRAIN_FILE)
occ    <- train$occ_clean
bg     <- train$bg_clean
suit_r <- rast(SUIT_FILE)
mess_r <- rast(MESS_FILE)
pop    <- rast(POP_ALIGNED_FILE)
cat("Presences:", nrow(occ), "| background:", nrow(bg), "\n")

# ------------------------- Travel-time surface ------------------------------
# Raw input, downloaded once (not a cache of a derived output).

if (!file.exists(TT_FILE)) {
  response <- GET(TT_URL, user_agent("R - MSc dissertation, e.p.naymon@lse.ac.uk"),
                  write_disk(TT_FILE, overwrite = TRUE), progress())
  stopifnot("Travel-time download failed" = status_code(response) == 200)
}

tt_raw <- rast(TT_FILE)
e_tt   <- as.vector(ext(tt_raw))
e_dom  <- as.vector(ext(vect(ADM0_FILE)))
stopifnot(
  "Travel-time CRS differs from the covariates" = same.crs(tt_raw, suit_r),
  "Travel-time raster does not cover the domain" =
    e_tt[["xmin"]] <= e_dom[["xmin"]] && e_tt[["xmax"]] >= e_dom[["xmax"]] &&
    e_tt[["ymin"]] <= e_dom[["ymin"]] && e_tt[["ymax"]] >= e_dom[["ymax"]],
  "MESS or population grid differs from the surface" =
    compareGeom(suit_r, mess_r, stopOnError = FALSE) &&
    compareGeom(suit_r, pop,    stopOnError = FALSE)
)
cat("Travel time:", ncol(tt_raw), "x", nrow(tt_raw), "at", res(tt_raw)[1], "deg | range",
    paste(round(unlist(global(tt_raw, "range", na.rm = TRUE))), collapse = " to "), "min\n")

tt <- travel_time(suit_r)

# --------------------------------- Gaps -------------------------------------
# Domain cells with a prediction but no travel time. Bilinear resampling
# returns NA next to a missing cell; nearest-neighbour shows which gaps are
# in the source grid.

zones   <- state_zones(suit_r)
gap     <- !is.na(suit_r) & is.na(tt)
gap_src <- !is.na(suit_r) & is.na(travel_time(suit_r, method = "near"))
gap_df  <- zonal(c(gap, gap_src, pop * gap), zones, fun = "sum", na.rm = TRUE)
names(gap_df) <- c("state", "cells", "cells_in_source", "population")
gap_df  <- gap_df |> filter(cells > 0) |> arrange(desc(population))

cat("\nDomain cells without travel time:", sum(gap_df$cells),
    "| in the source grid:", sum(gap_df$cells_in_source),
    "| population:", fmt(sum(gap_df$population)), "\n")
print(mutate(gap_df, population = fmt(population)), row.names = FALSE)

# -------------------------------- Points ------------------------------------

xy <- function(d) as.matrix(d[, c("longitude", "latitude")])
occ$tt      <- terra::extract(tt, xy(occ))[, 1]
bg$tt       <- terra::extract(tt, xy(bg))[, 1]
bg$in_range <- terra::extract(mess_r, xy(bg))[, 1] >= 0
stopifnot("A presence has no travel time"         = !anyNA(occ$tt),
          "A background point has no MESS value"  = !anyNA(bg$in_range))
cat("\nBackground points without travel time (excluded below):",
    sum(is.na(bg$tt)), "of", nrow(bg), "\n")

v <- data.frame(tt = values(tt, mat = FALSE), suit = values(suit_r, mat = FALSE),
                pop = values(pop, mat = FALSE), in_range = values(mess_r, mat = FALSE) >= 0)
v <- v[!is.na(v$tt) & !is.na(v$suit), ]
stopifnot("A predicted cell has no MESS value" = !anyNA(v$in_range))

# ------------------------------ Comparison ----------------------------------
# Reading fixed before the run. p_presence_closer is the chance a presence is
# closer to a city than a reference point (0.5 = no difference); 0.56 / 0.64
# / 0.71 are the conventional small / medium / large benchmarks (Cohen's d
# 0.2 / 0.5 / 0.8). Against background within the presence range, a medium
# or larger effect means recording favours accessible places beyond the
# desert's remoteness; a small one means the full-background contrast is
# mostly the desert. A presence median near the population median means
# records sit where people live.

summ <- function(group, x, w = rep(1, length(x)), compare = FALSE) {
  keep <- !is.na(x) & !is.na(w); x <- x[keep]; w <- w[keep]
  o  <- order(x); cw <- cumsum(w[o]) / sum(w)
  q  <- sapply(c(0.25, 0.5, 0.75), function(p) x[o][which(cw >= p)[1]])
  data.frame(group = group, n = round(sum(w)), q25 = q[1], median = q[2], q75 = q[3],
             pct_beyond_remote = 100 * sum(w[x > TT_REMOTE_MIN]) / sum(w),
             p_presence_closer = if (compare) auc_ties(-occ$tt, -x) else NA_real_)
}

types <- sort(unique(occ$presence_type))
comp <- bind_rows(
  summ("Presences", occ$tt),
  bind_rows(lapply(types, function(t)
    summ(paste("Presences:", t), occ$tt[occ$presence_type == t]))),
  summ("Background", bg$tt, compare = TRUE),
  summ("Background, in presence range", bg$tt[bg$in_range], compare = TRUE),
  summ("Population, in presence range", v$tt[v$in_range], w = v$pop[v$in_range])
)

cat("\nTravel time to the nearest city (minutes); pct_beyond_remote = share beyond",
    TT_REMOTE_MIN, "min; n = people for the population row\n")
comp |> mutate(across(c(q25, median, q75), round),
               pct_beyond_remote = round(pct_beyond_remote, 1),
               p_presence_closer = round(p_presence_closer, 3)) |> print(row.names = FALSE)

# --------------------- Suitability and accessibility ------------------------
# Not a test of bias: the suitable belt is also the populated, connected one.
# Describes the confound 12 and 13 address. Every cell, no sampling.

rho <- c(domain            = cor(v$suit, v$tt, method = "spearman"),
         in_presence_range = cor(v$suit[v$in_range], v$tt[v$in_range], method = "spearman"))
cat("\nSpearman correlation, suitability vs travel time:\n"); print(round(rho, 3))

# -------------------------------- Figure ------------------------------------

plot_df <- bind_rows(
  data.frame(group = "Presences", tt = occ$tt),
  data.frame(group = "Background", tt = bg$tt),
  data.frame(group = "Background, in presence range", tt = bg$tt[bg$in_range])
) |> filter(!is.na(tt)) |> mutate(group = factor(group, levels = names(pal_access)))

p_acc <- ggplot(plot_df, aes(tt, fill = group)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(xintercept = TT_REMOTE_MIN, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 10, 60, 300, 1440, 5000)) +
  scale_fill_manual(values = pal_access) +
  labs(x = "Travel time to nearest city (minutes, log scale)", y = "Density", fill = NULL) +
  theme(legend.position = "bottom")

save_fig(file.path(DIR_FIGS, "accessibility_bias_diagnostic.png"), p_acc,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT)

# --------------------------------- Save -------------------------------------

write.csv(comp,   file.path(DIR_TABLES, "accessibility_comparison.csv"), row.names = FALSE)
write.csv(gap_df, file.path(DIR_TABLES, "travel_time_gaps_by_state.csv"), row.names = FALSE)
write.csv(data.frame(scope = names(rho), spearman = rho),
          file.path(DIR_TABLES, "suitability_travel_time_correlation.csv"), row.names = FALSE)
cat("11_accessibility_bias_diagnostic.R complete\n")



