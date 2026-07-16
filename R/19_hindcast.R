# ============================================================================
# 19_hindcast.R
# Projects the MaxEnt model onto 2005 covariates with 2005 WorldPop population
# to compare against Alvar et al. (2006)'s expert-derived Gedaref ARP of
# 0.98M — an entirely independent evidence stream from environmental niche
# modelling.
#
# Inputs:  outputs/models/maxent_final.rds
#          outputs/models/selected_tuning.rds
#          outputs/models/training_data.rds
#          outputs/models/retained_vars.rds
#          outputs/surfaces/maxent_suitability.tif
#          outputs/tables/arp_summary.csv
#          data/raw/ (2005 covariate rasters + population)
# Outputs: outputs/surfaces/maxent_suitability_2005.tif
#          outputs/tables/arp_hindcast_comparison.csv
#          outputs/figures/suitability_2005_vs_ltm.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(maxnet)
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(geodata)
  library(patchwork)
  library(rnaturalearth)
})

set.seed(SEED)

# ------------------------------ Load inputs ---------------------------------

mod   <- readRDS(file.path(DIR_MODELS, "maxent_final.rds"))
sel   <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))
train <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
vars  <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))

pred_ltm <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))

cat("Model:", sel$fc, "rm =", sel$rm, "\n")

# Thresholds from training data
pred_occ <- predict(mod, train$occ_env, type = "cloglog")[, 1]
pred_bg  <- predict(mod, train$bg_env, type = "cloglog")[, 1]

p10 <- unname(quantile(pred_occ, 0.10))
candidates <- sort(unique(c(pred_occ, pred_bg)))
sens <- sapply(candidates, function(t) mean(pred_occ >= t))
spec <- sapply(candidates, function(t) mean(pred_bg < t))
maxsss <- candidates[which.max(sens + spec)]

cat("p10:   ", round(p10, 4), "\n")
cat("maxSSS:", round(maxsss, 4), "\n")

# ---------------------- Build 2005 covariate stack --------------------------

cov_files_2005 <- c(
  slope      = "slope_1km.tif",
  river_dist = "river_distance_1km.tif",
  vertisols  = "vertisols_1km.tif",
  lst_night  = "lst_night_annual_2005_1km.tif",
  rainfall   = "rainfall_2005_1km.tif"
)

covs_2005 <- rast(file.path(DIR_COVARIATES, cov_files_2005[vars]))
names(covs_2005) <- vars

covs_ltm <- rast(file.path(DIR_COVARIATES, COV_FILES[vars]))
names(covs_ltm) <- vars

stopifnot(
  crs(covs_2005) == crs(covs_ltm),
  all(res(covs_2005) == res(covs_ltm)),
  ext(covs_2005) == ext(covs_ltm)
)
cat("Alignment check passed\n")

for (v in c("lst_night", "rainfall")) {
  diff <- global(covs_2005[[v]] - covs_ltm[[v]], c("min", "max", "mean"), na.rm = TRUE)
  cat(v, "\u2014 2005 minus LTM: mean", round(diff$mean, 2),
      "| range [", round(diff$min, 2), ",", round(diff$max, 2), "]\n")
}

# ----------------------- Predict 2005 surface -------------------------------

sudan  <- ne_countries(country = "Sudan", scale = 50, returnclass = "sf")
mask_r <- rast(file.path(DIR_COVARIATES, "ecological_mask_150mm.tif"))

# Full-Sudan prediction (for ARP overlay)
covs_2005_sudan <- mask(covs_2005, vect(sudan))
pred_2005_full <- terra::predict(covs_2005_sudan, mod, type = "cloglog", na.rm = TRUE)

cat("\n2005 surface (full): mean",
    round(global(pred_2005_full, "mean", na.rm = TRUE)[[1]], 4), "\n")
cat("LTM surface:         mean",
    round(global(pred_ltm, "mean", na.rm = TRUE)[[1]], 4), "\n")

writeRaster(pred_2005_full, file.path(DIR_SURFACES, "maxent_suitability_2005.tif"),
            overwrite = TRUE)

# -------------------- 2005 population and ARP -------------------------------

pop_2005_raw <- rast(file.path(DIR_COVARIATES, "population", "sdn_ppp_2005_UNadj.tif"))

pop_2005 <- terra::project(pop_2005_raw, pred_2005_full, method = "sum")
pop_2005 <- mask(pop_2005, pred_2005_full)

total_pop_2005 <- global(pop_2005, "sum", na.rm = TRUE)[[1]]
cat("Total population (2005):", format(round(total_pop_2005), big.mark = ","), "\n")

# National ARP
suit_2005_p10    <- pred_2005_full >= p10
suit_2005_maxsss <- pred_2005_full >= maxsss

arp_p10_2005      <- global(pop_2005 * suit_2005_p10, "sum", na.rm = TRUE)[[1]]
arp_maxsss_2005   <- global(pop_2005 * suit_2005_maxsss, "sum", na.rm = TRUE)[[1]]
arp_weighted_2005 <- global(pop_2005 * pred_2005_full, "sum", na.rm = TRUE)[[1]]

cat("\nAt-Risk Population (2005 surface, 2005 population):\n")
cat("  p10:           ", format(round(arp_p10_2005), big.mark = ","), "\n")
cat("  maxSSS:        ", format(round(arp_maxsss_2005), big.mark = ","), "\n")
cat("  Risk-weighted: ", format(round(arp_weighted_2005), big.mark = ","), "\n")

# Gedaref extraction (direct comparison to Alvar 0.98M)
adm1 <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))
gedaref <- adm1[adm1$NAME_1 == "Al Qadarif", ]

arp_gedaref_weighted <- terra::extract(pop_2005 * pred_2005_full, gedaref,
                                       fun = "sum", na.rm = TRUE, ID = FALSE)[[1]]
arp_gedaref_maxsss   <- terra::extract(pop_2005 * suit_2005_maxsss, gedaref,
                                       fun = "sum", na.rm = TRUE, ID = FALSE)[[1]]
pop_gedaref <- terra::extract(pop_2005, gedaref, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]]

cat("\nGedaref (Al Qadarif):\n")
cat("  Population:    ", format(round(pop_gedaref), big.mark = ","), "\n")
cat("  ARP weighted:  ", format(round(arp_gedaref_weighted), big.mark = ","), "\n")
cat("  ARP maxSSS:    ", format(round(arp_gedaref_maxsss), big.mark = ","), "\n")
cat("  Alvar (2006):   980,000\n")
cat("  Ratio (model/Alvar):", round(arp_gedaref_weighted / 980000, 2), "\n")

# -------------------- Comparison table --------------------------------------

arp_ltm <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))

comparison <- data.frame(
  source = c("Alvar 2006 (expert-derived)",
             "This model (2005 covariates)",
             "This model (LTM)",
             "Pigott 2014 (global BRT)"),
  year_conditions = c("2005", "2005", "2000\u20132024", "~2010"),
  population_year = c("~2005", "2005", "2025", "2010"),
  boundaries = c("Pre-split", "Current Sudan", "Current Sudan", "Current Sudan"),
  arp_estimate = c("2.78M (incl. S. Sudan)",
                   format(round(arp_weighted_2005), big.mark = ","),
                   format(round(arp_ltm$arp[arp_ltm$metric == "risk_weighted"]), big.mark = ","),
                   "16,259,580"),
  method = c("Expert/case-based", "ENM risk-weighted", "ENM risk-weighted",
             "BRT binary (0.19)")
)

cat("\nCross-study ARP comparison:\n")
print(comparison, right = FALSE, row.names = FALSE)

write.csv(comparison, file.path(DIR_TABLES, "arp_hindcast_comparison.csv"),
          row.names = FALSE)

# ========================= FIGURE ===========================================

suit_colours <- c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7",
                  "#EF8A62", "#B2182B")

make_suit_map <- function(r, title) {
  df <- as.data.frame(r, xy = TRUE)
  names(df) <- c("x", "y", "suitability")
  df <- df[!is.na(df$suitability), ]

  ggplot() +
    geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
    geom_raster(data = df, aes(x = x, y = y, fill = suitability)) +
    scale_fill_gradientn(colours = suit_colours, limits = c(0, 1),
                         na.value = "transparent", name = "Suitability") +
    geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
    coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 23), crs = 4326) +
    labs(title = title) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.title = element_blank())
}

p_2005 <- make_suit_map(pred_2005_full, "2005 (hindcast)")
p_ltm  <- make_suit_map(pred_ltm, "Long-term mean (2000\u20132024)")

p_hindcast <- (p_2005 | p_ltm) +
  plot_annotation(title = "MaxEnt suitability: 2005 hindcast vs. long-term mean")

ggsave(file.path(DIR_FIGS, "suitability_2005_vs_ltm.png"), p_hindcast,
       width = 14, height = 7, dpi = 300, bg = "white")
cat("Saved suitability_2005_vs_ltm.png\n")

cat("\n19_hindcast.R complete\n")