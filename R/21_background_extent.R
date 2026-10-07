# ============================================================================
# 21_background_extent.R
# Supplement. Does drawing the background only from the >= 150 mm region
# (the dissertation's ecological mask) change the surface and the estimate?
# The primary draws background from the whole domain. Here the model is
# refitted on the primary's background points inside the region. Presences
# (all of them, including those below 150 mm), their folds and year-matched
# values, the prediction stack and the population are unchanged. The
# restricted background is a subset of the primary's, so within-belt CV
# scores both models on the same test rows.
#
# Inputs:  TRAIN_FILE, TUNING_FILE, MODEL_FILE, SUIT_FILE, POP_ALIGNED_FILE,
#          SENS_MASK_FILE, DOMAIN_FILE, ADM1_FILE, DISPLAY_ADM0_FILE,
#          EXCLUDED_FILE, retained_vars.rds, enmeval_results.csv (06),
#          arp_summary.csv and arp_plateau_candidates.csv (08), COV_FILES
# Outputs: BELT_SURFACES_FILE (layers fixed, tuned)
#          outputs/tables/belt_background_grid.csv
#          outputs/tables/belt_background_summary.csv
#          outputs/tables/belt_background_by_state.csv
#          outputs/tables/belt_background_curves.csv
#          outputs/figures/response_curves_belt_background.png
#          outputs/figures/fig_belt_background_difference.png / .pdf
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(maxnet)
  library(ggplot2); library(patchwork)
})
t_start <- Sys.time()

# ------------------------------ Load inputs ---------------------------------

train   <- readRDS(TRAIN_FILE)
tuning  <- readRDS(TUNING_FILE)
mod     <- readRDS(MODEL_FILE)
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
enm_06  <- read.csv(file.path(DIR_TABLES, "enmeval_results.csv"))
arp_08  <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))
cand_08 <- read.csv(file.path(DIR_TABLES, "arp_plateau_candidates.csv"))
suit_r  <- rast(SUIT_FILE)
pop     <- rast(POP_ALIGNED_FILE)
dom_r   <- rast(DOMAIN_FILE)
occ     <- train$occ_clean
bg      <- train$bg_clean
occ_env <- train$occ_env[, vars]
bg_env  <- train$bg_env[, vars]
region  <- paste0(">= ", SENS_MASK_MM, " mm")
cat("Model:", tuning$fc, "rm =", tuning$rm, "| presences:", nrow(occ),
    "| background:", nrow(bg), "\n")

# --------------------------- >= 150 mm region -------------------------------
# SENS_MASK_FILE covers the GEE export rectangle, not Sudan, so the region is
# its 1 cells inside the domain. It is the same CHIRPS computation as the
# rainfall covariate (python/01_study_area.ipynb, 02_covariates.ipynb), so the
# two must agree on every domain cell with a rainfall value.

dom_v  <- values(dom_r, mat = FALSE) %in% 1
sens_v <- values(rast(SENS_MASK_FILE), mat = FALSE)
stopifnot("SENS_MASK_FILE is not coded 0 / 1 over the domain" = all(sens_v[dom_v] %in% c(0, 1)))
belt_v <- dom_v & sens_v %in% 1

rain_v <- values(rast(file.path(DIR_COVARIATES, COV_FILES[["rainfall"]])), mat = FALSE)
has_r  <- dom_v & !is.na(rain_v)
n_off  <- sum((rain_v[has_r] >= SENS_MASK_MM) != belt_v[has_r])
if (n_off > 0) stop(n_off, " domain cells where SENS_MASK_FILE disagrees with rainfall >= SENS_MASK_MM")
cat("Domain cells in the", region, "region:", fmt(sum(belt_v)), "of", fmt(sum(dom_v)),
    "| disagreements with the rainfall covariate:", n_off, "\n")

# ------------------------- Restricted background ----------------------------
# The primary's background points inside the region, which 06 flags as `wet`.
# Nothing is redrawn: every point keeps its 06 values and fold, and this
# script uses no random numbers.

cell_of <- function(d) cellFromXY(dom_r, as.matrix(d[, c("longitude", "latitude")]))
stopifnot("TRAIN_FILE's region flags differ from SENS_MASK_FILE: rerun 06" =
  identical(bg$wet, belt_v[cell_of(bg)]) && identical(occ$wet, belt_v[cell_of(occ)]))

keep <- list(primary = rep(TRUE, nrow(bg)), belt = bg$wet)
bg_per_fold <- table(factor(bg$fold[keep$belt], levels = seq_len(K_FOLDS)))
stopifnot("A fold has no background in the region" = all(bg_per_fold > 0))
cat("Background in the region:", sum(keep$belt), "of", nrow(bg), "| by fold:", bg_per_fold, "\n")
cat("Presences below", SENS_MASK_MM, "mm (kept; no background there under the restriction):",
    sum(!occ$wet), "of", nrow(occ), "\n")
print(occ[!occ$wet, c("coordinate_id", "presence_type", "source", "year")], row.names = FALSE)

pa_of   <- function(k) c(rep(1, nrow(occ_env)), rep(0, sum(k)))
env_of  <- function(k) rbind(occ_env, bg_env[k, ])
fold_of <- function(k) c(occ$fold, bg$fold[k])
wet_of  <- function(k) c(occ$wet, bg$wet[k])

# ---------------------- Grid on the restricted background -------------------
# tune_grid with within-belt scores. Selection by 06's automatic rule, without
# its plausibility override: highest mean CV CBI among configurations scored
# on every fold; ties to fewer feature classes, then higher rm. The rule is
# checked on 06's own table first. CBI here is against the restricted
# background, so it is not comparable with 06's. If the best sits at the top
# of ENM_RM with CBI still rising by more than one SE, it is flagged and the
# tuned row is read as a direction only (the grid is not extended).

select_best <- function(tune)
  tune |> filter(n_cbi == K_FOLDS) |> arrange(desc(cbi), nchar(fc), desc(rm)) |> slice(1)

b06 <- select_best(transmute(enm_06, fc, rm, cbi = cbi.val.avg,
                             n_cbi = ifelse(n_folds == K_FOLDS, n_cbi, 0L)))
stopifnot("select_best() does not reproduce 06's highest-CBI configuration" =
  b06$fc == tuning$plateau$fc[1] && b06$rm == tuning$plateau$rm[1])

cat("\nGrid on the restricted background:", length(ENM_FC) * length(ENM_RM),
    "configurations x", K_FOLDS, "folds\n")
t0     <- Sys.time()
grid_b <- tune_grid(pa_of(keep$belt), env_of(keep$belt), fold_of(keep$belt),
                    within = wet_of(keep$belt), progress = TRUE)
cat(sprintf("Grid: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

best <- select_best(grid_b)
sel  <- filter(grid_b, fc == tuning$fc, rm == tuning$rm)
stopifnot("The primary settings were not scored on every fold of the restricted background" =
  nrow(sel) == 1 && sel$n_cbi == K_FOLDS)
se_best <- best$cbi_sd / sqrt(K_FOLDS)
lower   <- filter(grid_b, fc == best$fc, rm < best$rm, n_cbi == K_FOLDS)
at_edge <- best$rm == max(ENM_RM) && nrow(lower) > 0 && best$cbi - max(lower$cbi) > se_best

cat("Scored on every fold:", sum(grid_b$n_cbi == K_FOLDS), "of", nrow(grid_b), "\n")
cat("Best:", best$fc, "rm", best$rm, "| CBI", round(best$cbi, 3), "+/-", round(best$cbi_sd, 3),
    if (at_edge) "| AT THE GRID EDGE (CBI still rising)" else "", "\n")
cat("Primary settings: CBI ", round(sel$cbi, 3), " | rank ",
    sum(grid_b$cbi > sel$cbi, na.rm = TRUE) + 1, " of ", sum(!is.na(grid_b$cbi)),
    " | within one SE of the best: ", sel$cbi >= best$cbi - se_best, "\n", sep = "")
cat("Within one SE of the best:\n")
grid_b |> filter(n_cbi == K_FOLDS, cbi >= best$cbi - se_best) |> arrange(desc(cbi)) |>
  select(fc, rm, cbi, cbi_sd, auc, cbi_wet) |>
  mutate(across(where(is.double), ~ round(., 3))) |> print(row.names = FALSE)

# -------------------------------- Refits ------------------------------------
# Each model is fitted on all its rows, cross-validated with within-belt
# scores and predicted over the domain with the long-term means. Thresholds
# are taken against the full background for every model, so binary estimates
# share one reference (as in 12); p10 depends on the presences only.

covs <- domain_covs(vars)
fit_one <- function(k, fc, rm) {
  r <- refit_maxnet(occ_env, bg_env[k, ], occ$fold, bg$fold[k], fc, rm, within = wet_of(k))
  r$thr  <- thresholds_from(maxnet_prob(r$mod, occ_env), maxnet_prob(r$mod, bg_env))
  r$surf <- terra::predict(covs, r$mod, type = "cloglog", na.rm = TRUE)
  r$arp  <- arp_estimates(r$surf, pop, r$thr)
  r$cfg  <- list(fc = fc, rm = rm, keep = k)
  cat(sprintf("  %s rm %s, %s background points: %d features\n",
              fc, rm, fmt(sum(k)), length(r$mod$betas)))
  r
}

cat("\nRefits:\n")
same_cfg <- best$fc == tuning$fc && best$rm == tuning$rm
fits <- list(primary = fit_one(keep$primary, tuning$fc, tuning$rm),
             fixed   = fit_one(keep$belt,    tuning$fc, tuning$rm))
fits$tuned <- if (same_cfg) fits$fixed else fit_one(keep$belt, best$fc, best$rm)
if (same_cfg) cat("  tuned = fixed: the primary settings are also the best here\n")

# ------------------------- Primary reproduces 06-08 -------------------------

p       <- fits$primary
wet_06  <- enm_06$cbi_wet.avg[enm_06$fc == tuning$fc & enm_06$rm == tuning$rm]
arp_08v <- setNames(arp_08$arp, arp_08$metric)
stopifnot(
  "Coefficients differ from MODEL_FILE (06)" = isTRUE(all.equal(p$mod$betas, mod$betas)),
  "CV CBI differs from 06"                   = abs(mean(p$cv$cbi) - tuning$cbi) < 1e-8,
  "Within-belt CV CBI differs from 06"       = abs(mean(p$cv$cbi_wet, na.rm = TRUE) - wet_06) < 1e-8,
  "Surface differs from SUIT_FILE (07)" =
    global(abs(p$surf - suit_r), "max", na.rm = TRUE)[[1]] < 1e-6,
  "Estimates differ from 08" =
    all(round(arp_estimates(suit_r, pop, p$thr)) == arp_08v[names(p$arp)]),
  "Fixed refit's CV differs from its grid row" =
    abs(mean(fits$fixed$cv$cbi, na.rm = TRUE) - sel$cbi) < 1e-8
)
cat("Primary reproduces 06 (coefficients, CV and within-belt CBI), 07 (surface) and 08 (estimates)\n")

# ------------------------------ Comparison ----------------------------------
# Readings fixed before the run:
#  1. National, fixed row: a risk-weighted change within the plateau spread
#     from 08 means the background extent matters less than the plateau
#     choice; beyond it, the estimate depends on the background extent, and
#     the range is reported beside the other sensitivity ranges.
#  2. Within-belt discrimination, fixed row: both models are scored on the
#     same within-belt test rows. Within-belt CV CBI above the primary's by
#     more than one SE (the primary's within-belt SE) means the desert in the
#     background dilutes discrimination within the belt; below it by more
#     than one SE, restricting the background loses discrimination there;
#     otherwise no detectable difference.
#  3. Descriptive: the tuned row, CV against each model's own background,
#     rank agreement, the split at the region boundary, the share of each
#     estimate where the model clamps, and the binary estimates.

plateau <- max(abs(cand_08$change_vs_selected[startsWith(cand_08$candidate, paste0(tuning$fc, " "))]))
se_wet  <- sd(p$cv$cbi_wet, na.rm = TRUE) / sqrt(K_FOLDS)

pop_v  <- values(pop, mat = FALSE)
v_prim <- values(p$surf, mat = FALSE)
cov_v  <- values(covs)
rw_sum <- function(v, cond) sum(pop_v[cond] * v[cond], na.rm = TRUE)

# Cells where a model clamps: any covariate beyond its training range
# (maxnet's varmin / varmax, presences plus background).
clamped <- function(m) {
  out <- rep(FALSE, nrow(cov_v))
  for (x in vars) out <- out | cov_v[, x] < m$varmin[[x]] | cov_v[, x] > m$varmax[[x]]
  out %in% TRUE
}

summary_df <- bind_rows(lapply(names(fits), function(nm) {
  r <- fits[[nm]]; v <- values(r$surf, mat = FALSE); ok <- !is.na(v) & !is.na(v_prim)
  data.frame(
    model = nm, background = if (all(r$cfg$keep)) "domain" else region,
    fc = r$cfg$fc, rm = r$cfg$rm, n_background = sum(r$cfg$keep),
    n_coef = length(r$mod$betas),
    cv_cbi_own_bg = mean(r$cv$cbi, na.rm = TRUE), cv_auc_own_bg = mean(r$cv$auc),
    cv_cbi_belt = mean(r$cv$cbi_wet, na.rm = TRUE), cv_auc_belt = mean(r$cv$auc_wet),
    rho_domain = cor(v[ok], v_prim[ok], method = "spearman"),
    rho_belt   = cor(v[ok & belt_v], v_prim[ok & belt_v], method = "spearman"),
    arp_weighted    = r$arp[["risk_weighted"]],
    rw_in_region    = rw_sum(v, belt_v),
    rw_below_region = rw_sum(v, dom_v & !belt_v),
    rw_clamped      = rw_sum(v, clamped(r$mod)),
    p10 = r$thr[["p10"]], arp_p10 = r$arp[["p10"]],
    maxsss = r$thr[["maxsss"]], arp_maxsss = r$arp[["maxsss"]])
})) |>
  mutate(rw_change_pct    = 100 * (arp_weighted / arp_weighted[model == "primary"] - 1),
         beyond_plateau   = abs(rw_change_pct) > plateau,
         in_region_share  = 100 * rw_in_region / arp_weighted,
         in_region_pct    = 100 * (rw_in_region / rw_in_region[model == "primary"] - 1),
         below_region_pct = 100 * (rw_below_region / rw_below_region[model == "primary"] - 1),
         clamped_share    = 100 * rw_clamped / arp_weighted,
         cbi_belt_diff    = cv_cbi_belt - cv_cbi_belt[model == "primary"],
         belt_reading     = case_when(model == "primary"     ~ NA_character_,
                                      cbi_belt_diff >  se_wet ~ "higher by > 1 SE",
                                      cbi_belt_diff < -se_wet ~ "lower by > 1 SE",
                                      TRUE                    ~ "within 1 SE"))
stopifnot("Region split does not sum to the national estimate" =
  all(abs((summary_df$rw_in_region + summary_df$rw_below_region) / summary_df$arp_weighted - 1) < POP_TOL))

cat("\nFit. Own background: not comparable across rows. Within the region: same test rows.\n")
summary_df |>
  select(model, background, fc, rm, n_background, n_coef, cv_cbi_own_bg, cv_auc_own_bg,
         cv_cbi_belt, cv_auc_belt, cbi_belt_diff, belt_reading) |>
  mutate(across(cv_cbi_own_bg:cbi_belt_diff, ~ round(., 3))) |>
  print(row.names = FALSE)
cat("Primary within-belt SE:", round(se_wet, 3), "\n")

cat("\nSurface and population at risk | plateau spread: +/-", round(plateau, 1), "%\n")
summary_df |>
  select(model, rho_domain, rho_belt, arp_weighted, rw_change_pct, beyond_plateau,
         in_region_share, in_region_pct, below_region_pct, clamped_share) |>
  mutate(across(c(rho_domain, rho_belt), ~ round(., 3)), arp_weighted = fmt(arp_weighted),
         across(c(rw_change_pct, in_region_share, in_region_pct, below_region_pct, clamped_share),
                ~ round(., 1))) |>
  print(row.names = FALSE)
cat("Risk-weighted in / below the region:",
    paste(summary_df$model, fmt(summary_df$rw_in_region), "/", fmt(summary_df$rw_below_region),
          collapse = " | "), "\n")

cat("\nBinary (thresholds against the full background):\n")
summary_df |> select(model, p10, arp_p10, maxsss, arp_maxsss) |>
  mutate(across(c(p10, maxsss), ~ round(., 3)), across(c(arp_p10, arp_maxsss), fmt)) |>
  print(row.names = FALSE)

# -------------------------------- States ------------------------------------

zones    <- state_zones(suit_r)
surf_all <- do.call(c, unname(lapply(fits, `[[`, "surf")))
names(surf_all) <- names(fits)
rw_all   <- pop * surf_all
names(rw_all) <- names(fits)

state_df <- zonal(rw_all, zones, fun = "sum", na.rm = TRUE)
names(state_df)[1] <- "state"
stopifnot("States do not sum to the national estimates" =
  all(abs(colSums(state_df[names(fits)]) / summary_df$arp_weighted - 1) < POP_TOL))
state_df <- state_df |>
  mutate(across(c(fixed, tuned), ~ round(100 * (. / primary - 1), 1), .names = "{.col}_pct")) |>
  arrange(desc(primary))

cat("\nRisk-weighted estimate by state (% change from the primary):\n")
state_df |> select(state, primary, fixed, fixed_pct, tuned_pct) |>
  mutate(across(c(primary, fixed), fmt)) |> print(right = FALSE, row.names = FALSE)

# ---------------------------- Response curves -------------------------------
# Others at the presence median (as 07), over the primary's data range so the
# models share one axis. Beyond a restricted model's own range its curve is
# flat (clamped).

ref     <- sapply(occ_env, median)
env_all <- rbind(occ_env, bg_env)
shown   <- if (same_cfg) c("primary", "fixed") else names(fits)
curves  <- bind_rows(lapply(shown, function(nm)
  data.frame(model = nm, response_curves(fits[[nm]]$mod, ref, env_all, vars))))
feat <- bind_rows(lapply(shown, function(nm)
  data.frame(model = nm, curve_table(filter(curves, model == nm)))))

cat("\nResponse-curve positions (others at the presence median):\n")
feat |> filter(variable %in% c("lst_night", "rainfall")) |>
  select(variable, model, peak_x, rise_50, rise_90, fall_90, fall_50, min_suit) |>
  mutate(across(where(is.double), ~ signif(., 3))) |> arrange(variable) |>
  print(row.names = FALSE)

labels_rt <- belt_labels
labels_rt[["tuned"]] <- paste0(belt_labels[["tuned"]], " (", best$fc, " rm ", best$rm, ")")
curves <- curves |> mutate(model = factor(model, levels = names(pal_belt)),
                           label = factor(cov_labels[variable], levels = cov_labels))
p_curves <- ggplot(curves, aes(value, suit, colour = model)) +
  geom_line(data = filter(curves, variable != "vertisols"), linewidth = 0.7) +
  geom_point(data = filter(curves, variable == "vertisols"), size = 2) +
  facet_wrap(~ label, scales = "free_x", nrow = 2) +
  scale_colour_manual(values = pal_belt, labels = labels_rt, name = NULL) +
  labs(x = NULL, y = "Suitability (cloglog)") +
  theme(legend.position = "bottom", legend.direction = "vertical")
save_fig(file.path(DIR_FIGS, "response_curves_belt_background.png"), p_curves,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT)

# ------------------------------ Difference map ------------------------------
# Each restricted surface minus the primary, on one symmetric scale, with the
# region outlined.

display  <- st_as_sf(vect(DISPLAY_ADM0_FILE))
excluded <- st_as_sf(vect(EXCLUDED_FILE))
states   <- st_as_sf(vect(ADM1_FILE))
lim      <- display_limits(display)
belt_r   <- rast(dom_r); values(belt_r) <- ifelse(belt_v, 1, NA)
outline  <- st_as_sf(as.polygons(belt_r))

alt    <- setdiff(shown, "primary")
diff_r <- surf_all[[alt]] - surf_all[["primary"]]
names(diff_r) <- alt
d_lim  <- max(global(abs(diff_r), "max", na.rm = TRUE)[[1]])

diff_panel <- function(nm, title) {
  df <- as.data.frame(diff_r[[nm]], xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "diff")
  ggplot() +
    geom_raster(data = df, aes(x, y, fill = diff)) +
    scale_fill_diff(lim = d_lim) +
    layer_excluded(excluded) +
    layer_admin1(data = states) +
    layer_country(data = display) +
    layer_region(outline) +
    labs(title = title) +
    coord_display(lim) +
    theme_map() +
    theme(plot.title = element_text(size = 9, hjust = 0))
}
panels <- lapply(seq_along(alt), function(i)
  diff_panel(alt[i], paste0("(", letters[i], ") ", labels_rt[[alt[i]]], ", minus primary")))
fig_diff <- wrap_plots(panels, nrow = 1) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
for (ext in c("png", "pdf"))
  save_fig(file.path(DIR_FIGS, paste0("fig_belt_background_difference.", ext)), fig_diff,
           width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT)

# --------------------------------- Save -------------------------------------

writeRaster(surf_all[[c("fixed", "tuned")]], BELT_SURFACES_FILE, overwrite = TRUE)
write.csv(grid_b,     file.path(DIR_TABLES, "belt_background_grid.csv"),     row.names = FALSE)
write.csv(summary_df, file.path(DIR_TABLES, "belt_background_summary.csv"),  row.names = FALSE)
write.csv(state_df,   file.path(DIR_TABLES, "belt_background_by_state.csv"), row.names = FALSE)
write.csv(feat,       file.path(DIR_TABLES, "belt_background_curves.csv"),   row.names = FALSE)
cat(sprintf("21_background_extent.R complete (%.1f min)\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))