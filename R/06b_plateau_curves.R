# ============================================================================
# 06b_plateau_curves.R
# Compares response curves of the plateau candidates from 06, so the choice
# among statistically equivalent configurations can be made on ecological
# plausibility (params.R). Candidates are refitted on all data with the same
# settings as 06.
#
# Inputs:  TRAIN_FILE, TUNING_FILE, retained_vars.rds, enmeval_results.csv
# Outputs: outputs/figures/plateau_response_curves.png
#          outputs/figures/plateau_rain_lst_surface.png
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(terra)
})

train  <- readRDS(TRAIN_FILE)
tuning <- readRDS(TUNING_FILE)
retained_vars <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))

occ_env <- train$occ_env[, retained_vars]
bg_env  <- train$bg_env[, retained_vars]
pa      <- c(rep(1, nrow(occ_env)), rep(0, nrow(bg_env)))
env_all <- rbind(occ_env, bg_env)

# ------------------------------ Candidates ----------------------------------
# Ends of the plateau for the selected feature classes, plus LQH at the lower
# rm as a reference without product features.

stopifnot("TUNING_FILE has no plateau: rerun 06" = !is.null(tuning$plateau))
ends <- tuning$plateau |> filter(fc == tuning$fc) |>
  slice(c(which.min(rm), which.max(rm))) |> select(fc, rm)
candidates <- bind_rows(ends, data.frame(fc = "LQH", rm = min(ends$rm))) |> distinct()

res <- read.csv(file.path(DIR_TABLES, "enmeval_results.csv"))
cat("Candidates, cross-validated:\n")
res |> semi_join(candidates, by = c("fc", "rm")) |>
  select(fc, rm, cbi.val.avg, cbi.val.sd, cbi_wet.avg, cbi_wet.sd, or.10p.avg) |>
  mutate(across(where(is.numeric), ~ round(., 3))) |> print()

mods <- lapply(seq_len(nrow(candidates)), function(i)
  fit_maxnet(pa, env_all, candidates$fc[i], candidates$rm[i]))
names(mods) <- paste0(candidates$fc, " rm ", candidates$rm)

for (nm in names(mods)) {
  cat("\n", nm, ":", length(mods[[nm]]$betas), "features\n")
  print(signif(sort(mods[[nm]]$betas), 3))
}

# --------------------------- Reference values -------------------------------
# Others held at the presence median. With the desert in the background, the
# background median is a desert cell and would flatten every curve.

ref <- sapply(occ_env, median)
cat("\nReference (presence median):\n"); print(signif(ref, 3))

ref_frame <- function(n) as.data.frame(matrix(rep(ref, each = n), ncol = length(ref),
                                              dimnames = list(NULL, names(ref))))

# --------------------------- Response curves --------------------------------

curve_df <- bind_rows(lapply(names(mods), function(nm) {
  bind_rows(lapply(retained_vars, function(v) {
    x  <- if (v == "vertisols") c(0, 1) else
      seq(min(env_all[[v]]), max(env_all[[v]]), length.out = 200)
    nd <- ref_frame(length(x)); nd[[v]] <- x
    data.frame(model = nm, variable = v, value = x,
               suit = as.numeric(predict(mods[[nm]], nd, type = "cloglog", clamp = TRUE)))
  }))
}))

cat("\nVertisols (suitability off vs on clay):\n")
curve_df |> filter(variable == "vertisols") |>
  select(model, value, suit) |> mutate(suit = round(suit, 3)) |>
  pivot_wider(names_from = value, names_prefix = "vertisols_", values_from = suit) |> print()

occ_long   <- occ_env |> pivot_longer(everything(), names_to = "variable")
pres_range <- occ_long |> group_by(variable) |> summarise(lo = min(value), hi = max(value))

p_curves <- ggplot(filter(curve_df, variable != "vertisols"),
                   aes(value, suit, colour = model)) +
  geom_rect(data = filter(pres_range, variable != "vertisols"),
            aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey92") +
  geom_line(linewidth = 0.8) +
  geom_rug(data = filter(occ_long, variable != "vertisols"),
           aes(x = value), inherit.aes = FALSE, sides = "b", alpha = 0.4) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(x = NULL, y = "Suitability (cloglog)", colour = NULL,
       title = "Plateau candidates: response curves",
       subtitle = "Others at the presence median; grey band and rug = presences") +
  theme_minimal() + theme(legend.position = "top")

ggsave(file.path(DIR_FIGS, "plateau_response_curves.png"), p_curves,
       width = 10, height = 6, dpi = 300)

# --------------------- Rainfall x night temperature -------------------------

grid <- expand.grid(
  rainfall  = seq(min(env_all$rainfall),  max(env_all$rainfall),  length.out = 100),
  lst_night = seq(min(env_all$lst_night), max(env_all$lst_night), length.out = 100))

surf_df <- bind_rows(lapply(names(mods), function(nm) {
  nd <- ref_frame(nrow(grid)); nd$rainfall <- grid$rainfall; nd$lst_night <- grid$lst_night
  cbind(grid, model = nm,
        suit = as.numeric(predict(mods[[nm]], nd, type = "cloglog", clamp = TRUE)))
}))

p_surf <- ggplot(surf_df, aes(rainfall, lst_night, fill = suit)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Suitability", limits = c(0, 1)) +
  geom_point(data = bg_env, aes(rainfall, lst_night), inherit.aes = FALSE,
             colour = "grey80", size = 0.2, alpha = 0.15) +
  geom_point(data = occ_env, aes(rainfall, lst_night), inherit.aes = FALSE,
             shape = 21, fill = "white", size = 1.3) +
  facet_wrap(~ model) +
  labs(x = "Rainfall (mm/yr)", y = "LST night (\u00b0C)",
       title = "Rainfall x night temperature, plateau candidates",
       subtitle = "Others at the presence median; white = presences, grey = background") +
  theme_minimal()

ggsave(file.path(DIR_FIGS, "plateau_rain_lst_surface.png"), p_surf,
       width = 12, height = 4.5, dpi = 300)

# ------------------ Candidates on the prediction surface --------------------
# Do the candidates differ where it matters: across Sudan's long-term surface?

dom  <- rast(DOMAIN_FILE)
covs <- rast(file.path(DIR_COVARIATES, COV_FILES[retained_vars]))
names(covs) <- retained_vars
covs <- mask(covs, dom, maskvalues = 0)
cat("\nLong-term range across the domain:\n")
print(global(covs[[c("rainfall", "lst_night")]], c("min", "max"), na.rm = TRUE))

surf <- do.call(c, lapply(mods, function(m)
  terra::predict(covs, m, type = "cloglog", na.rm = TRUE)))
names(surf) <- names(mods)
writeRaster(surf, file.path(DIR_SURFACES, "plateau_candidate_surfaces.tif"), overwrite = TRUE)

v <- values(surf)
cat("\nMean suitability across the domain:\n"); print(round(colMeans(v, na.rm = TRUE), 3))
cat("\nCorrelation between candidate surfaces:\n"); print(round(cor(v, use = "complete.obs"), 3))

d48 <- v[, "LQHP rm 4"] - v[, "LQHP rm 8"]
cat("\nLQHP rm 4 minus rm 8 (1st / 50th / 99th percentiles):\n")
print(round(quantile(d48, c(0.01, 0.5, 0.99), na.rm = TRUE), 3))

big <- which(abs(d48) >= quantile(abs(d48), 0.99, na.rm = TRUE))
cat("\nCovariates where the difference is largest (top 1%):\n")
print(summary(values(covs)[big, c("rainfall", "lst_night")]))

cat("\n06b_plateau_curves.R complete\n")