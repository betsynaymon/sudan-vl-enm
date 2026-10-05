# ============================================================================
# 14_uncertainty_surface.R
# How much do the surface and the estimate depend on which block of presences
# the model sees? Takes 06's fold models (the selected configuration fitted
# without each spatial CV fold), predicts each over the domain, and reports
# the spread across the four surfaces and estimates, nationally and by state.
# This is sensitivity to leaving out one geographic block (and a quarter of
# the presences), not a confidence interval: there are four models, and an
# SD on a 0-1 scale is near zero wherever suitability is near 0 or 1.
#
# Inputs:  TRAIN_FILE, TUNING_FILE, MODEL_FILE, SUIT_FILE, POP_ALIGNED_FILE,
#          SENS_MASK_FILE, DOMAIN_FILE, ADM1_FILE, DISPLAY_ADM0_FILE,
#          EXCLUDED_FILE, retained_vars.rds, COV_FILES,
#          tuning_by_fold.csv (06), arp_summary.csv and
#          arp_plateau_candidates.csv (08)
# Outputs: FOLD_SURFACES_FILE, FOLD_SD_FILE
#          outputs/tables/arp_fold_excluded.csv, fold_excluded_by_state.csv
#          outputs/figures/fig_fold_sd.png / .pdf
# Runtime: about 3-6 min (four fits, four domain predictions).
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(ggplot2); library(maxnet)
})

# ------------------------------ Load inputs ---------------------------------

train   <- readRDS(TRAIN_FILE)
tuning  <- readRDS(TUNING_FILE)
mod     <- readRDS(MODEL_FILE)
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
suit_r  <- rast(SUIT_FILE)
pop     <- rast(POP_ALIGNED_FILE)
belt_r  <- rast(SENS_MASK_FILE)
fold_06 <- read.csv(file.path(DIR_TABLES, "tuning_by_fold.csv"))
arp_08  <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))
cand_08 <- read.csv(file.path(DIR_TABLES, "arp_plateau_candidates.csv"))

occ     <- train$occ_clean
bg      <- train$bg_clean
occ_env <- train$occ_env[, vars]
bg_env  <- train$bg_env[, vars]

# Each row carries its fold and >= 150 mm flag from 06.
pa   <- c(rep(1, nrow(occ_env)), rep(0, nrow(bg_env)))
env  <- rbind(occ_env, bg_env)
fold <- c(occ$fold, bg$fold)
wet  <- c(occ$wet, bg$wet)
stopifnot(
  "Presence rows and values differ"            = nrow(occ) == nrow(occ_env),
  "Background rows and values differ"          = nrow(bg)  == nrow(bg_env),
  "A training value is missing"                = !anyNA(env),
  "A row has no valid fold"                    = all(fold %in% seq_len(K_FOLDS)),
  "A row has no >= 150 mm flag"                = is.logical(wet) && !anyNA(wet),
  "River distance is not a retained covariate" = "river_dist" %in% vars,
  "Population grid differs from the surface"   = compareGeom(pop,    suit_r, stopOnError = FALSE),
  "Belt mask grid differs from the surface"    = compareGeom(belt_r, suit_r, stopOnError = FALSE)
)
zones     <- state_zones(suit_r)
occ$state <- as.character(terra::extract(zones, as.matrix(occ[, c("longitude", "latitude")]))[, 1])
stopifnot("A presence has no state" = !anyNA(occ$state))
cat("Presences:", sum(pa == 1), "| background:", sum(pa == 0),
    "| test presences per fold:", table(fold[pa == 1]), "\n")
cat("\nPresences by state and fold (each fold-excluded model drops one column):\n")
print(table(occ$state, fold = occ$fold))

# --------------- Fold models reproduce 06; full model reproduces 08 ----------

cat("\nFold models,", tuning$fc, "rm", tuning$rm, "...\n")
cv        <- cv_maxnet(pa, env, fold, tuning$fc, tuning$rm, within = wet, keep_models = TRUE)
fold_mods <- attr(cv, "models")
thr_full  <- thresholds_from(maxnet_prob(mod, occ_env), maxnet_prob(mod, bg_env))
arp_full  <- arp_estimates(suit_r, pop, thr_full)
rw_full   <- arp_full[["risk_weighted"]]

m_cols  <- c("auc", "cbi", "or_10p", "auc_wet", "cbi_wet", "or_10p_wet")
f06     <- fold_06 |> filter(fc == tuning$fc, rm == tuning$rm) |> arrange(fold)
thr_08  <- setNames(arp_08$threshold, arp_08$metric)
arp_08v <- setNames(arp_08$arp,       arp_08$metric)
stopifnot(
  "Fold-model CV metrics differ from tuning_by_fold.csv (06)" =
    nrow(f06) == K_FOLDS && length(fold_mods) == K_FOLDS &&
    isTRUE(all.equal(as.matrix(cv[, m_cols]), as.matrix(f06[, m_cols]),
                     check.attributes = FALSE)),
  "Thresholds differ from arp_summary.csv (08)" = all(abs(thr_full - thr_08[names(thr_full)]) < 1e-6),
  "Estimates differ from arp_summary.csv (08)"  = all(round(arp_full) == arp_08v[names(arp_full)])
)
cat("Fold models reproduce 06 (CV by fold); the full model reproduces 08 (thresholds, estimates)\n")
cat("Features per fold model:", sapply(fold_mods, function(m) length(m$betas)),
    "| full model:", length(mod$betas), "\n")

# ------------------------ Fold-excluded surfaces ------------------------------

covs   <- domain_covs(vars)
fold_r <- rast(lapply(seq_len(K_FOLDS), function(k) {
  cat("  predicting without fold", k, "\n")
  terra::predict(covs, fold_mods[[k]], type = "cloglog", na.rm = TRUE)
}))
names(fold_r) <- paste0("without_fold_", seq_len(K_FOLDS))
writeRaster(fold_r, FOLD_SURFACES_FILE, overwrite = TRUE)
fold_r <- rast(FOLD_SURFACES_FILE)   # saved (FLT4S) values
stopifnot("A fold surface covers different cells from SUIT_FILE" =
  all(sapply(seq_len(K_FOLDS), function(k)
    global(is.na(fold_r[[k]]) != is.na(suit_r), "sum")[[1]] == 0)))

# ------------------------------ Estimates -------------------------------------
# Each model's thresholds from its own training rows (year-matched), as 08
# does for the full model.

plateau <- max(abs(cand_08$change_vs_selected[startsWith(cand_08$candidate,
                                                         paste0(tuning$fc, " "))]))
by_st <- lapply(seq_len(K_FOLDS), function(k) rw_by_scale(fold_r[[k]], suit_r, pop, zones))

est <- bind_rows(lapply(seq_len(K_FOLDS), function(k) {
  tr  <- fold != k
  m   <- fold_mods[[k]]
  thr <- thresholds_from(maxnet_prob(m, env[tr & pa == 1, ]), maxnet_prob(m, env[tr & pa == 0, ]))
  a   <- arp_estimates(fold_r[[k]], pop, thr)
  data.frame(model = names(fold_r)[k], train_presences = sum(tr & pa == 1),
             features = length(m$betas), p10 = thr[["p10"]], maxsss = thr[["maxsss"]],
             arp_weighted = a[["risk_weighted"]], arp_p10 = a[["p10"]],
             arp_maxsss = a[["maxsss"]], common = sum(by_st[[k]]$common))
})) |>
  mutate(own_pct    = 100 * (arp_weighted / rw_full - 1),
         common_pct = 100 * (common / rw_full - 1),
         reading    = case_when(abs(common_pct) > plateau ~ "geography",
                                abs(own_pct)    > plateau ~ "scale",
                                TRUE                      ~ "within"))
stopifnot(
  "rw_by_scale and arp_estimates disagree" =
    all(abs(sapply(by_st, function(d) sum(d$own)) / est$arp_weighted - 1) < 1e-9),
  "rw_by_scale's reference differs from the full estimate" =
    all(abs(sapply(by_st, function(d) sum(d$ref)) / rw_full - 1) < 1e-9)
)

cat("\nPopulation at risk without each fold | full model: risk-weighted", fmt(rw_full),
    "| maxSSS", fmt(arp_full[["maxsss"]]), "| p10", fmt(arp_full[["p10"]]),
    "| plateau spread +/-", round(plateau, 1), "%\n")
est |>
  mutate(across(c(arp_weighted, common, arp_maxsss, arp_p10), fmt),
         across(c(own_pct, common_pct), ~ round(., 1)),
         across(c(p10, maxsss), ~ round(., 3))) |>
  select(model, train_presences, features, arp_weighted, own_pct, common, common_pct,
         reading, arp_maxsss, maxsss, arp_p10, p10) |>
  print(row.names = FALSE)
cat("Range, risk-weighted: own scale", fmt(min(est$arp_weighted)), "to", fmt(max(est$arp_weighted)),
    "| full model's scale", fmt(min(est$common)), "to", fmt(max(est$common)), "\n")

# ------------------------------- By state -------------------------------------

stopifnot("States differ between fold models" =
  all(sapply(by_st, function(d) identical(d$state, by_st[[1]]$state))))
common_pct <- sapply(by_st, function(d) 100 * (d$common / d$ref - 1))
own_pct    <- sapply(by_st, function(d) 100 * (d$own    / d$ref - 1))
colnames(common_pct) <- paste0("common_wo", seq_len(K_FOLDS))
colnames(own_pct)    <- paste0("own_wo",    seq_len(K_FOLDS))

sd_r  <- app(fold_r, "sd");    names(sd_r) <- "sd"
mm    <- app(fold_r, "range"); rng_r <- mm[[2]] - mm[[1]]
writeRaster(sd_r, FOLD_SD_FILE, overwrite = TRUE)

pop_m  <- mask(pop, sd_r)
s_suit <- zonal(suit_r, zones, fun = "mean", na.rm = TRUE)
s_sd   <- zonal(sd_r,   zones, fun = "mean", na.rm = TRUE)
s_pw   <- zonal(c(sd_r * pop_m, pop_m), zones, fun = "sum", na.rm = TRUE)
names(s_suit) <- c("state", "suit_area")
names(s_sd)   <- c("state", "sd_area")
names(s_pw)   <- c("state", "sd_pop_sum", "pop")

st_tab <- data.frame(state = by_st[[1]]$state, full = by_st[[1]]$ref, common_pct, own_pct) |>
  left_join(s_suit, by = "state") |>
  left_join(s_sd,   by = "state") |>
  left_join(transmute(s_pw, state, sd_pop = sd_pop_sum / pop), by = "state") |>
  arrange(desc(full))
stopifnot("A state lacks a summary" = !anyNA(st_tab))

cat("\nBy state: change without each fold on the full model's scale (%), and SD across",
    "the four surfaces (area mean; population-weighted) beside mean suitability:\n")
st_tab |>
  mutate(full = fmt(full), across(starts_with("common_wo"), ~ round(., 1)),
         across(c(suit_area, sd_area, sd_pop), ~ round(., 3))) |>
  select(state, full, starts_with("common_wo"), suit_area, sd_area, sd_pop) |>
  print(row.names = FALSE)

v_sd    <- values(sd_r, mat = FALSE)
v_riv   <- values(covs[["river_dist"]], mat = FALSE)
ok      <- !is.na(v_sd)
belt    <- values(belt_r, mat = FALSE) %in% 1
rho_riv <- c(domain = cor(v_sd[ok], v_riv[ok], method = "spearman"),
             belt   = cor(v_sd[ok & belt], v_riv[ok & belt], method = "spearman"))
cat(sprintf("\nSD across the four surfaces: domain mean %.3f, max %.3f | max - min: mean %.3f, max %.3f\n",
            mean(v_sd[ok]), max(v_sd[ok]),
            global(rng_r, "mean", na.rm = TRUE)[[1]], global(rng_r, "max", na.rm = TRUE)[[1]]))
cat(sprintf("Spearman rho, SD vs river distance: domain %.3f | within >= 150 mm %.3f -> drainage reading %s\n",
            rho_riv[["domain"]], rho_riv[["belt"]],
            if (rho_riv[["domain"]] <= -0.1) "kept" else "dropped"))

# -------------------------------- Figure --------------------------------------

display  <- st_as_sf(vect(DISPLAY_ADM0_FILE))
excluded <- st_as_sf(vect(EXCLUDED_FILE))
states   <- st_as_sf(vect(ADM1_FILE))
lim      <- display_limits(display)
occ_sf   <- st_as_sf(occ, coords = c("longitude", "latitude"), crs = 4326)
sd_df    <- as.data.frame(sd_r, xy = TRUE, na.rm = TRUE)

p_sd <- ggplot() +
  geom_sf(data = states, fill = "grey95", colour = NA) +
  geom_raster(data = sd_df, aes(x, y, fill = sd)) +
  layer_excluded(excluded) +
  scale_fill_spread() +
  layer_admin1(states) +
  layer_country(display) +
  layer_occurrences(occ_sf, size = 1.2, stroke = 0.3) +
  add_scalebar() +
  coord_display(lim) +
  theme_map() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
for (ext in c("png", "pdf"))
  save_fig(file.path(DIR_FIGS, paste0("fig_fold_sd.", ext)), p_sd,
           width = FIG_WIDTH_FULL, height = FIG_HEIGHT_MAP)

# --------------------------------- Save -------------------------------------

write.csv(est,    file.path(DIR_TABLES, "arp_fold_excluded.csv"),     row.names = FALSE)
write.csv(st_tab, file.path(DIR_TABLES, "fold_excluded_by_state.csv"), row.names = FALSE)
cat("\n14_uncertainty_surface.R complete\n")