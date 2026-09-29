# ============================================================================
# 20_east_west_diagnostic.R
# Extracts predicted suitability at each thinned presence location, classifies
# east (non-Darfur) vs west (Darfur), and compares distributions. Tests
# whether the western gap persists across all covariate variants and after
# accessibility correction. 
#
# Inputs:  outputs/models/maxent_final.rds
#          outputs/models/training_data.rds
#          outputs/surfaces/maxent_suitability.tif
#          outputs/surfaces/maxent_suitability_{variant}.tif  (9 variants)
#          outputs/surfaces/maxent_suitability_{variant}_fixed.tif
#          outputs/surfaces/maxent_suitability_biased_bg.tif
#          data/processed/occurrences_thinned.csv
# Outputs: outputs/tables/east_west_diagnostic.csv
#          outputs/figures/east_west_suitability_diagnostic.png
#          outputs/figures/fig_east_west.png
#          outputs/figures/fig_east_west.pdf
#          outputs/figures/presence_suitability_by_longitude.png
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "plotting_theme.R"))


suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(geodata)
  library(maxnet)
})

set.seed(SEED)

# ------------------------------ Load inputs ---------------------------------

pred_ltm <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))
train    <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
mod      <- readRDS(file.path(DIR_MODELS, "maxent_final.rds"))

# Thresholds
pred_occ <- predict(mod, train$occ_env, type = "cloglog")[, 1]
pred_bg  <- predict(mod, train$bg_env, type = "cloglog")[, 1]

p10 <- unname(quantile(pred_occ, 0.10))
candidates <- sort(unique(c(pred_occ, pred_bg)))
sens <- sapply(candidates, function(t) mean(pred_occ >= t))
spec <- sapply(candidates, function(t) mean(pred_bg < t))
maxsss <- candidates[which.max(sens + spec)]

cat("p10:", round(p10, 4), "| maxSSS:", round(maxsss, 4), "\n")

# ---------------------- Classify east vs west -------------------------------

occ_all <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))

# Keep the 98 training presences; print the thinned record(s) not in training
occ_keys   <- paste(round(occ_all$longitude, 5), round(occ_all$latitude, 5), occ_all$year)
train_keys <- paste(round(train$occ_clean$longitude, 5),
                    round(train$occ_clean$latitude, 5), train$occ_clean$year)

occ <- occ_all[occ_keys %in% train_keys, ]
cat("Presences used:", nrow(occ), "of", nrow(occ_all), "thinned records\n")
cat("Thinned record(s) not in the training set:\n")
print(occ_all[!occ_keys %in% train_keys, ])

adm1 <- gadm(country = "SDN", level = 1, path = here::here("data", "raw"))

occ_sf <- st_as_sf(occ, coords = c("longitude", "latitude"), crs = 4326)
occ_states <- st_join(occ_sf, st_as_sf(adm1)[, "NAME_1"])

darfur_states <- c("Central Darfur", "East Darfur", "North Darfur",
                   "South Darfur", "West Darfur")

occ_states$region <- ifelse(occ_states$NAME_1 %in% darfur_states, "West", "East")

# Assign border points (NA state) to nearest state
na_idx <- is.na(occ_states$NAME_1)
if (any(na_idx)) {
  nearest <- st_nearest_feature(occ_states[na_idx, ], st_as_sf(adm1))
  occ_states$NAME_1[na_idx] <- adm1$NAME_1[nearest]
  occ_states$region[na_idx] <- ifelse(
    adm1$NAME_1[nearest] %in% darfur_states, "West", "East"
  )
  cat("Assigned", sum(na_idx), "border points to nearest state\n")
}

cat("\nRecords by region:\n")
print(table(occ_states$region))

cat("\nWestern records by state:\n")
occ_states |> filter(region == "West") |>
  st_drop_geometry() |> count(NAME_1) |> print()

# ---------------------- Extract and compare ---------------------------------

occ_coords <- st_coordinates(occ_states)
occ_states$suitability <- terra::extract(pred_ltm, occ_coords)[, 1]

cat("\nSuitability by region:\n")
region_summary <- occ_states |>
  st_drop_geometry() |>
  group_by(region) |>
  summarise(
    n           = n(),
    mean_suit   = round(mean(suitability, na.rm = TRUE), 3),
    median_suit = round(median(suitability, na.rm = TRUE), 3),
    min_suit    = round(min(suitability, na.rm = TRUE), 3),
    max_suit    = round(max(suitability, na.rm = TRUE), 3),
    .groups     = "drop"
  )
print(region_summary)

# ---------------------- Strip plot ---------------------------

plot_df <- occ_states |>
  st_drop_geometry() |>
  select(region, suitability, NAME_1)

east_med <- region_summary$median_suit[region_summary$region == "East"]
west_med <- region_summary$median_suit[region_summary$region == "West"]
n_east   <- region_summary$n[region_summary$region == "East"]
n_west   <- region_summary$n[region_summary$region == "West"]

p_ew <- ggplot(plot_df, aes(x = region, y = suitability, colour = region)) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = maxsss, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = p10, linetype = "dotted", colour = "grey60") +
  annotate("text", x = 2.4, y = p10 + 0.03, label = "p10 threshold",
           size = 3, colour = "grey60") +
  annotate("text", x = 2.4, y = maxsss + 0.03, label = "maxSSS threshold",
           size = 3, colour = "grey40") +
  scale_colour_manual(values = c("East" = "#B2182B", "West" = "#2166AC")) +
  labs(x = NULL, y = "Predicted suitability (LTM surface)",
       title = "Model performance at known VL locations") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(DIR_FIGS, "east_west_suitability_diagnostic.png"), p_ew,
       width = 7, height = 6, dpi = 300, bg = "white")
cat("Saved east_west_suitability_diagnostic.png\n")

# ------------------- Cross-variant sensitivity ---------------------------

variant_names <- c("B_annual", "C_annual", "D_annual", "A_dry", "A_wet",
                   "E_elevation", "F_lstday", "G_treecover", "H_noriver")
surface_names <- c(variant_names, paste0(variant_names, "_fixed"))

variant_surfaces <- c(
  list(A_annual = rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))),
  setNames(lapply(surface_names, function(v) {
    rast(file.path(DIR_SURFACES, paste0("maxent_suitability_", v, ".tif")))
  }), surface_names)
)

west_pts    <- occ_states |> filter(region == "West")
west_coords <- st_coordinates(west_pts)

west_suit <- sapply(variant_surfaces, function(r) {
  terra::extract(r, west_coords)[, 1]
})

west_suit_df <- as.data.frame(west_suit)
west_suit_df$state <- west_pts$NAME_1
west_suit_df <- west_suit_df |> select(state, everything())

cat("\nSuitability at western presences across variants:\n")
print(round(west_suit_df[, -1], 4))

cat("\nMedian suitability per variant:\n")
sapply(west_suit_df[, -1], median) |> round(4) |> print()

cat("\nMax suitability per variant:\n")
sapply(west_suit_df[, -1], max) |> round(4) |> print()

# ------------------- Accessibility-corrected check --------------------------

pred_biased <- rast(file.path(DIR_SURFACES, "maxent_suitability_biased_bg.tif"))
west_suit_bias <- terra::extract(pred_biased, west_coords)[, 1]

cat("\nWestern presences \u2014 accessibility-corrected surface:\n")
print(round(west_suit_bias, 4))
cat("Median:", round(median(west_suit_bias), 4),
    "| Max:", round(max(west_suit_bias), 4), "\n")

# ------------- Breakdown by state, below-p10, and by longitude -------

occ_states$cell <- cellFromXY(pred_ltm, occ_coords)
occ_states$lon  <- occ_coords[, 1]

cat("\nPresence suitability by state:\n")
occ_states |>
  st_drop_geometry() |>
  group_by(NAME_1) |>
  summarise(n = n(), distinct_cells = n_distinct(cell),
            median_suit = round(median(suitability), 3),
            min_suit    = round(min(suitability), 3), .groups = "drop") |>
  arrange(desc(median_suit)) |>
  print(n = Inf)

cat("\nPresences below p10:\n")
occ_states |>
  st_drop_geometry() |>
  filter(suitability < p10) |>
  select(NAME_1, region, year, source, lon, suitability) |>
  arrange(suitability) |>
  print()

p_lon <- ggplot(st_drop_geometry(occ_states),
                aes(x = lon, y = suitability, colour = region)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = p10, linetype = "dotted", colour = "grey50") +
  scale_colour_manual(values = c("East" = "#B2182B", "West" = "#2166AC")) +
  labs(x = "Longitude (\u00b0E)", y = "Predicted suitability (LTM surface)",
       colour = NULL) +
  theme_minimal()

ggsave(file.path(DIR_FIGS, "presence_suitability_by_longitude.png"), p_lon,
       width = 7, height = 5, dpi = 300, bg = "white")
cat("Saved presence_suitability_by_longitude.png\n")

# -------------- Breakdown by source type -----------------
occ_states$source_type <- case_when(
  grepl("radio_dabanga|IFRC|reliefweb|sudan_tribune", occ_states$source,
        ignore.case = TRUE)          ~ "Humanitarian/news",
  grepl("Pigott", occ_states$source) ~ "Pigott database",
  TRUE                               ~ "Peer-reviewed"
)

occ_states |>
  st_drop_geometry() |>
  group_by(source_type) |>
  summarise(n = n(), below_p10 = sum(suitability < p10),
            median_suit = round(median(suitability), 3), .groups = "drop") |>
  print()

# --------- Below-p10 share by source type: humanitarian/news vs others -------
# Exploratory (grouping chosen after inspecting the below-p10 records).

src_df <- occ_states |>
  st_drop_geometry() |>
  mutate(news      = source_type == "Humanitarian/news",
         below_p10 = suitability < p10)

tab_all <- table(news = src_df$news, below_p10 = src_df$below_p10)
print(tab_all)

n_news      <- sum(src_df$news)
n_other     <- sum(!src_df$news)
k_news      <- sum(src_df$below_p10 & src_df$news)
k_other     <- sum(src_df$below_p10 & !src_df$news)
share_news  <- k_news / n_news
share_other <- k_other / n_other

ft_all <- fisher.test(tab_all)

cat("\nBelow p10 — humanitarian/news:", k_news, "of", n_news,
    "(", round(100 * share_news, 1), "% ) | all others:", k_other, "of", n_other,
    "(", round(100 * share_other, 1), "% )\n")
cat("Relative risk:", round(share_news / share_other, 1),
    "| Fisher's exact p =", signif(ft_all$p.value, 2), "\n")

# Same comparison against peer-reviewed records only
pr_df <- src_df |> filter(source_type %in% c("Humanitarian/news", "Peer-reviewed"))
ft_pr <- fisher.test(table(pr_df$news, pr_df$below_p10))
cat("Against peer-reviewed only: Fisher's exact p =", signif(ft_pr$p.value, 2), "\n")

# --------------------------------- Save -------------------------------------

east_west_out <- occ_states |>
  cbind(st_coordinates(occ_states)) |>
  st_drop_geometry() |>
  rename(longitude = X, latitude = Y) |>
  select(longitude, latitude, year, source,
         state = NAME_1, region, suitability)

write.csv(east_west_out, file.path(DIR_TABLES, "east_west_diagnostic.csv"),
          row.names = FALSE)

cat("\n--- East-west diagnostic summary ---\n")
cat("Eastern presences (n=", n_east, "): median suitability", east_med, "\n")
cat("Western presences (n=", n_west, "): median suitability", west_med, "\n")


# ================== DISSERTATION FIGURE =======================================
# Objects needed: plot_df, region_summary, east_med, west_med,
#                 n_east, n_west, p10, maxsss
# =============================================================================

p_ew_diss <- ggplot(plot_df, aes(x = region, y = suitability, colour = region)) +
  geom_point(position = position_jitter(width = 0.15, seed = 42),
             size = 2.5, alpha = 0.7) +
  geom_hline(yintercept = p10, linetype = "dotted", colour = "grey60") +
  geom_hline(yintercept = maxsss, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 2.35, y = p10 + 0.03, label = "p10 threshold",
           size = 3, colour = "grey50", hjust = 1) +
  annotate("text", x = 2.35, y = maxsss + 0.03, label = "maxSSS threshold",
           size = 3, colour = "grey40", hjust = 1) +
  scale_colour_manual(values = c("East" = "#B2182B", "West" = "#2166AC")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL,
       y = "Predicted suitability" #,
       #title = "Predicted suitability at documented VL locations"
       ) +
  theme_dissertation(gridlines = "both") +
  theme(legend.position = "none",
        panel.grid.major = element_line(colour = "grey92"),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 12))

save_fig(file.path(DIR_FIGS, "fig_east_west.png"), p_ew_diss,
         width = FIG_WIDTH_HALF * 1.5, height = FIG_HEIGHT_PLOT)
save_fig(file.path(DIR_FIGS, "fig_east_west.pdf"), p_ew_diss,
         width = FIG_WIDTH_HALF * 1.5, height = FIG_HEIGHT_PLOT)
cat("Saved fig_east_west\n")

cat("\n20_east_west_diagnostic.R complete\n")
