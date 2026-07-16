# ============================================================================
# 18_2025_surface.R
# Projects the fitted MaxEnt model onto 2025 single-year covariates to test
# how much interannual climate variation inflates the ARP estimate. The
# comparison against the long-term mean surface diagnoses whether the LST
# night step-function amplifies modest temperature shifts into large ARP
# swings — empirical evidence for the LTM choice.
#
# Inputs:  outputs/models/maxent_final.rds
#          outputs/models/selected_tuning.rds
#          outputs/models/training_data.rds
#          outputs/models/retained_vars.rds
#          outputs/surfaces/maxent_suitability.tif
#          outputs/surfaces/worldpop_2025_aligned.tif
#          outputs/tables/arp_summary.csv
#          outputs/tables/arp_by_state.csv
#          data/raw/ (covariate rasters incl. 2025 annuals)
# Outputs: outputs/surfaces/maxent_suitability_2025.tif
#          outputs/surfaces/maxent_suitability_2025_full.tif
#          outputs/surfaces/binary_p10_2025.tif
#          outputs/surfaces/binary_maxsss_2025.tif
#          outputs/tables/arp_comparison_2025.csv
#          outputs/tables/arp_state_comparison_2025.csv
#          outputs/figures/suitability_2025_vs_ltm.png
#          outputs/figures/lst_night_shift_diagnostic.png
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
  library(ggspatial)
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
cat("Covariates:", paste(vars, collapse = ", "), "\n")

# ---------------------- Build 2025 covariate stack --------------------------

cov_files_2025 <- c(
  slope      = "slope_1km.tif",
  river_dist = "river_distance_1km.tif",
  vertisols  = "vertisols_1km.tif",
  lst_night  = "lst_night_annual_2025_1km.tif",
  rainfall   = "rainfall_2025_1km.tif"
)

covs_2025 <- rast(file.path(DIR_COVARIATES, cov_files_2025[vars]))
names(covs_2025) <- vars

covs_ltm <- rast(file.path(DIR_COVARIATES, COV_FILES[vars]))
names(covs_ltm) <- vars

stopifnot(
  crs(covs_2025) == crs(covs_ltm),
  all(res(covs_2025) == res(covs_ltm)),
  ext(covs_2025) == ext(covs_ltm)
)
cat("Alignment check passed\n")

for (v in c("lst_night", "rainfall")) {
  diff <- global(covs_2025[[v]] - covs_ltm[[v]], c("min", "max", "mean"), na.rm = TRUE)
  cat(v, "\u2014 2025 minus LTM: mean", round(diff$mean, 2),
      "| range [", round(diff$min, 2), ",", round(diff$max, 2), "]\n")
}

# ----------------------- Predict 2025 surface -------------------------------

adm0        <- gadm(country = "SDN", level = 0, path = here::here("data", "raw"))
sudan       <- st_as_sf(adm0)
mask_r <- rast(file.path(DIR_COVARIATES, "ecological_mask_150mm.tif"))

# Masked version (within ecological mask)
covs_2025_masked <- mask(mask(covs_2025, vect(sudan)), mask_r, maskvalues = 0)
pred_2025 <- terra::predict(covs_2025_masked, mod, type = "cloglog", na.rm = TRUE)

cat("\n2025 surface (masked): mean",
    round(global(pred_2025, "mean", na.rm = TRUE)[[1]], 4), "\n")

# Full-Sudan version (no ecological mask)
covs_2025_sudan <- mask(covs_2025, vect(sudan))
pred_2025_full <- terra::predict(covs_2025_sudan, mod, type = "cloglog", na.rm = TRUE)

cat("2025 surface (full):   mean",
    round(global(pred_2025_full, "mean", na.rm = TRUE)[[1]], 4), "\n")
cat("LTM surface:           mean",
    round(global(pred_ltm, "mean", na.rm = TRUE)[[1]], 4), "\n")

diff_full <- pred_2025_full - pred_ltm
cat("\nDifference (2025 minus LTM): mean",
    round(global(diff_full, "mean", na.rm = TRUE)[[1]], 4), "\n")

writeRaster(pred_2025, file.path(DIR_SURFACES, "maxent_suitability_2025.tif"),
            overwrite = TRUE)
writeRaster(pred_2025_full, file.path(DIR_SURFACES, "maxent_suitability_2025_full.tif"),
            overwrite = TRUE)

# ---------------------- Thresholds and binary surfaces ----------------------

pred_occ <- predict(mod, train$occ_env, type = "cloglog")[, 1]
pred_bg  <- predict(mod, train$bg_env, type = "cloglog")[, 1]

p10 <- unname(quantile(pred_occ, 0.10))
candidates <- sort(unique(c(pred_occ, pred_bg)))
sens <- sapply(candidates, function(t) mean(pred_occ >= t))
spec <- sapply(candidates, function(t) mean(pred_bg < t))
maxsss <- candidates[which.max(sens + spec)]

suit_2025_p10    <- pred_2025_full >= p10
suit_2025_maxsss <- pred_2025_full >= maxsss

writeRaster(suit_2025_p10, file.path(DIR_SURFACES, "binary_p10_2025.tif"),
            overwrite = TRUE)
writeRaster(suit_2025_maxsss, file.path(DIR_SURFACES, "binary_maxsss_2025.tif"),
            overwrite = TRUE)

# ----------------------------- ARP 2025 ------------------------------------

pop_aligned <- rast(file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"))
total_pop   <- global(pop_aligned, "sum", na.rm = TRUE)[[1]]

arp_p10_2025      <- global(pop_aligned * suit_2025_p10, "sum", na.rm = TRUE)[[1]]
arp_maxsss_2025   <- global(pop_aligned * suit_2025_maxsss, "sum", na.rm = TRUE)[[1]]
arp_weighted_2025 <- global(pop_aligned * pred_2025_full, "sum", na.rm = TRUE)[[1]]

arp_ltm <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))

arp_compare <- data.frame(
  metric    = c("p10", "maxSSS", "Risk-weighted"),
  threshold = c(round(p10, 4), round(maxsss, 4), NA),
  arp_ltm   = arp_ltm$arp,
  arp_2025  = c(round(arp_p10_2025), round(arp_maxsss_2025), round(arp_weighted_2025))
) |>
  mutate(
    change     = arp_2025 - arp_ltm,
    pct_change = round(100 * change / arp_ltm, 1)
  )

cat("\n--- ARP comparison: LTM vs 2025 ---\n")
arp_compare |>
  mutate(across(c(arp_ltm, arp_2025, change),
                ~ format(., big.mark = ","))) |>
  print(right = FALSE)

write.csv(arp_compare, file.path(DIR_TABLES, "arp_comparison_2025.csv"),
          row.names = FALSE)

# -------------------- State-level comparison --------------------------------

adm1 <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))

risk_weighted_2025_r <- pop_aligned * pred_2025_full

state_2025 <- data.frame(
  state      = adm1$NAME_1,
  total_pop  = terra::extract(pop_aligned, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]],
  arp_w_2025 = terra::extract(risk_weighted_2025_r, adm1, fun = "sum", na.rm = TRUE, ID = FALSE)[[1]]
)

state_ltm <- read.csv(file.path(DIR_TABLES, "arp_by_state.csv"))

state_comp <- state_2025 |>
  left_join(state_ltm |> select(state, arp_w_ltm = arp_weighted), by = "state") |>
  mutate(
    change     = round(arp_w_2025 - arp_w_ltm),
    pct_change = round(100 * change / arp_w_ltm, 1)
  ) |>
  arrange(desc(arp_w_2025))

cat("\nState-level risk-weighted ARP: 2025 vs LTM:\n")
state_comp |>
  mutate(across(c(total_pop, arp_w_2025, arp_w_ltm, change),
                ~ format(round(.), big.mark = ","))) |>
  print(right = FALSE)

write.csv(state_comp, file.path(DIR_TABLES, "arp_state_comparison_2025.csv"),
          row.names = FALSE)

# ====================== FIGURES =============================================

# Side-by-side suitability maps
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

p_ltm  <- make_suit_map(pred_ltm, "Long-term mean (2000\u20132024)")
p_2025 <- make_suit_map(pred_2025_full, "2025 annual")

diff_df <- as.data.frame(diff_full, xy = TRUE)
names(diff_df) <- c("x", "y", "diff")
diff_df <- diff_df[!is.na(diff_df$diff), ]

p_diff <- ggplot() +
  geom_sf(data = sudan, fill = "grey90", colour = "grey30", linewidth = 0.5) +
  geom_raster(data = diff_df, aes(x = x, y = y, fill = diff)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "\u0394 Suitability",
                       limits = c(-0.5, 0.5), oob = scales::squish) +
  geom_sf(data = sudan, fill = NA, colour = "grey30", linewidth = 0.5) +
  coord_sf(xlim = c(21.5, 39), ylim = c(8.5, 23), crs = 4326) +
  labs(title = "Difference (2025 minus LTM)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank())

p_comparison <- (p_ltm | p_2025 | p_diff) +
  plot_annotation(title = "MaxEnt suitability: long-term mean vs. 2025 projection")

ggsave(file.path(DIR_FIGS, "suitability_2025_vs_ltm.png"), p_comparison,
       width = 18, height = 7, dpi = 300)
cat("Saved suitability_2025_vs_ltm.png\n")

# LST night shift diagnostic
ltm_vals   <- values(covs_ltm[["lst_night"]], na.rm = TRUE)
v2025_vals <- values(covs_2025[["lst_night"]], na.rm = TRUE)

hist_df <- rbind(
  data.frame(lst = as.numeric(ltm_vals), surface = "Long-term mean"),
  data.frame(lst = as.numeric(v2025_vals), surface = "2025")
)

p_lst <- ggplot(hist_df, aes(x = lst, fill = surface)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = c(18, 25), linetype = "dashed", colour = "grey40") +
  annotate("rect", xmin = 18, xmax = 25, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "red") +
  annotate("text", x = 21.5, y = Inf, vjust = 2,
           label = "Steep response zone", size = 3.5) +
  scale_fill_manual(values = c("Long-term mean" = "steelblue",
                               "2025" = "firebrick")) +
  labs(x = "LST night (\u00b0C)", y = "Density", fill = NULL,
       title = "LST night distribution shift: LTM vs 2025",
       subtitle = "Shaded band = steep part of suitability response curve") +
  theme_minimal()

ggsave(file.path(DIR_FIGS, "lst_night_shift_diagnostic.png"), p_lst,
       width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved lst_night_shift_diagnostic.png\n")

cat("\n18_2025_surface.R complete\n")