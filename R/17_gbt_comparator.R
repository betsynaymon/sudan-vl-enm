# ============================================================================
# 17_gbt_comparator.R
# Does the result depend on the algorithm? Adds gradient boosted trees (GBT)
# beside MaxEnt and the random forest (09), fitted to the same presences,
# background, covariates and spatial folds, and compares all three on
# cross-validated fit, surface ranks, population at risk, response-curve
# positions and held-out prediction at the presences. MaxEnt and the forest
# run through the same code first and must reproduce 06-09. Draws the
# three-algorithm figures.
#
# Readings, fixed before the run (as 09 where they overlap):
#   Selection  highest mean CV CBI among configurations scored on every fold;
#              ties to fewer, then shallower, trees. Stops if the selection
#              is at GBT_NTREES_MAX: the grid didn't bracket it.
#   Fit        an algorithm fits as well as MaxEnt if its CV CBI is within
#              one SE of MaxEnt's (fold SD / sqrt(K_FOLDS)).
#   Estimate   beyond the plateau if its risk-weighted change from MaxEnt
#              exceeds 08's plateau spread. The algorithm range is the
#              lowest to highest of the three.
#   Transfer   in each state whose presences MaxEnt fails to predict held
#              out (held-out median below MaxEnt's p10), an algorithm whose
#              map median is at or above its own p10 but whose held-out
#              median is below it maps the state because its records are in
#              training, not because the environment transfers. (09's check
#              on the forest was post hoc; this applies it in advance.)
#   Curves     descriptive.
#
# Inputs:  TRAIN_FILE, TUNING_FILE, MODEL_FILE, SUIT_FILE, POP_ALIGNED_FILE,
#          SENS_MASK_FILE, DOMAIN_FILE, ADM0_FILE, ADM1_FILE,
#          DISPLAY_ADM0_FILE, EXCLUDED_FILE, retained_vars.rds, COV_FILES,
#          tuning_by_fold.csv and heldout_presence_predictions.csv (06),
#          response_curve_features.csv (07), arp_summary.csv and
#          arp_plateau_candidates.csv (08), RF_MODEL_FILE, RF_SUIT_FILE,
#          COMPARATOR_FILE, comparator_by_state.csv,
#          comparator_curve_features.csv (09)
# Outputs: GBT_MODEL_FILE, GBT_SUIT_FILE, THREE_MODEL_FILE
#          outputs/tables/gbt_tuning.csv, three_model_by_state.csv,
#          three_model_curve_features.csv, three_model_heldout.csv,
#          three_model_heldout_by_state.csv
#          outputs/figures/suitability_three_models.png / .pdf,
#          response_curves_three_models.png / .pdf
# Runtime: about 30-45 min, mostly the tuning loop (one line per setting).
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(ggplot2); library(patchwork)
  library(maxnet); library(ranger); library(gbm)
})

# ------------------------------ Load inputs ---------------------------------

train   <- readRDS(TRAIN_FILE)
tuning  <- readRDS(TUNING_FILE)
mod     <- readRDS(MODEL_FILE)
rf      <- readRDS(RF_MODEL_FILE)
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
suit_r  <- rast(SUIT_FILE)
rf_surf <- rast(RF_SUIT_FILE)
pop     <- rast(POP_ALIGNED_FILE)
belt_r  <- rast(SENS_MASK_FILE)
fold_06 <- read.csv(file.path(DIR_TABLES, "tuning_by_fold.csv"))
ho_06   <- read.csv(file.path(DIR_TABLES, "heldout_presence_predictions.csv"))
feat_07 <- read.csv(file.path(DIR_TABLES, "response_curve_features.csv"))
arp_08  <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))
cand_08 <- read.csv(file.path(DIR_TABLES, "arp_plateau_candidates.csv"))
cmp_09  <- read.csv(COMPARATOR_FILE)
st_09   <- read.csv(file.path(DIR_TABLES, "comparator_by_state.csv"))
feat_09 <- read.csv(file.path(DIR_TABLES, "comparator_curve_features.csv"))

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
  "Presence rows and values differ"          = nrow(occ) == nrow(occ_env),
  "Background rows and values differ"        = nrow(bg)  == nrow(bg_env),
  "A training value is missing"              = !anyNA(env),
  "A row has no valid fold"                  = all(fold %in% seq_len(K_FOLDS)),
  "A row has no >= 150 mm flag"              = is.logical(wet) && !anyNA(wet),
  "Population grid differs from the surface" = compareGeom(pop,     suit_r, stopOnError = FALSE),
  "Belt mask grid differs from the surface"  = compareGeom(belt_r,  suit_r, stopOnError = FALSE),
  "RF surface grid differs from the surface" = compareGeom(rf_surf, suit_r, stopOnError = FALSE),
  "RF out-of-bag predictions don't match the training rows" = nrow(rf$predictions) == length(pa),
  "COMPARATOR_FILE (09) should hold MaxEnt and RF only" =
    identical(sort(cmp_09$model), c("maxent", "rf"))
)
cat("Presences:", sum(pa == 1), "| background:", sum(pa == 0),
    "| test presences per fold:", table(fold[pa == 1]), "\n")
cat("gbm", as.character(packageVersion("gbm")),
    "| ranger", as.character(packageVersion("ranger")), "\n")

# ------------------------ MaxEnt reproduces 06-08 ---------------------------
# Through the code the other algorithms use below. Estimates from the saved
# surface (FLT4S); 09 checked the surface itself.

cat("\nMaxEnt", tuning$fc, "rm", tuning$rm, ": fit, CV, held-out...\n")
ref       <- sapply(occ_env, median)   # presence median, as 07 and 09
mx_mod    <- fit_maxnet(pa, env, tuning$fc, tuning$rm)
mx_cv     <- cv_maxnet(pa, env, fold, tuning$fc, tuning$rm, within = wet)
mx_ho     <- heldout_preds(pa, env, fold,
                           function(p, d) fit_maxnet(p, d, tuning$fc, tuning$rm), maxnet_prob)
mx_thr    <- thresholds_from(maxnet_prob(mx_mod, occ_env), maxnet_prob(mx_mod, bg_env))
mx_arp    <- arp_estimates(suit_r, pop, mx_thr)
mx_curves <- response_curves(mx_mod, ref, env, vars)

m_cols  <- c("auc", "cbi", "or_10p", "auc_wet", "cbi_wet", "or_10p_wet")
f06     <- fold_06 |> filter(fc == tuning$fc, rm == tuning$rm) |> arrange(fold)
ho_06   <- ho_06[match(occ$coordinate_id, ho_06$coordinate_id), ]
thr_08  <- setNames(arp_08$threshold, arp_08$metric)
arp_08v <- setNames(arp_08$arp,       arp_08$metric)
stopifnot(
  "Coefficients differ from MODEL_FILE (06)" = isTRUE(all.equal(mx_mod$betas, mod$betas)),
  "CV metrics by fold differ from tuning_by_fold.csv (06)" =
    nrow(f06) == K_FOLDS &&
    isTRUE(all.equal(as.matrix(mx_cv[, m_cols]), as.matrix(f06[, m_cols]),
                     check.attributes = FALSE)),
  "CV CBI differs from 06" = abs(mean(mx_cv$cbi) - tuning$cbi) < 1e-8,
  "Held-out predictions differ from heldout_presence_predictions.csv (06)" =
    !anyNA(ho_06$pred) && isTRUE(all.equal(mx_ho$pred, ho_06$pred)) &&
    all(mx_ho$below_threshold == ho_06$below_threshold),
  "Curve features differ from response_curve_features.csv (07)" =
    isTRUE(all.equal(as.data.frame(curve_table(mx_curves)), feat_07, check.attributes = FALSE)),
  "Thresholds differ from arp_summary.csv (08)" = all(abs(mx_thr - thr_08[names(mx_thr)]) < 1e-6),
  "Estimates differ from arp_summary.csv (08)"  = all(round(mx_arp) == arp_08v[names(mx_arp)])
)
cat("MaxEnt (", length(mx_mod$betas), " features) reproduces 06 (coefficients, CV by fold,",
    " held-out predictions), 07 (curve features) and 08 (thresholds, estimates)\n", sep = "")

# --------------------- Random forest reproduces 09 --------------------------
# The saved forest, its out-of-bag thresholds and its saved surface; CV and
# held-out predictions refitted (fit_rf is seeded). The row is checked
# against 09's table once all three are built.

cat("\nRandom forest, mtry", rf$mtry, ": CV, held-out...\n")
rf_cv     <- cv_fit(pa, env, fold, function(p, d) fit_rf(p, d, rf$mtry), rf_prob, within = wet)
rf_ho     <- heldout_preds(pa, env, fold, function(p, d) fit_rf(p, d, rf$mtry), rf_prob)
rf_thr    <- thresholds_from(rf$predictions[pa == 1, "1"], rf$predictions[pa == 0, "1"])
rf_arp    <- arp_estimates(rf_surf, pop, rf_thr)
rf_curves <- response_curves(rf, ref, env, vars, pred = rf_prob)
stopifnot("RF surface covers different cells from SUIT_FILE" =
  global(is.na(rf_surf) != is.na(suit_r), "sum")[[1]] == 0)

# ----------------------- Boosted trees: tune by spatial CV -------------------
# Same rows, folds and metrics as MaxEnt and the forest. For each learning
# rate and depth, every fold model is grown to GBT_NTREES_MAX trees and
# scored at every GBT_TREE_STEP; fit_gbt is seeded, so that equals growing
# each number of trees separately (checked against cv_fit at the selection).
# OR10 uses each fold model's own training presences, which boosting partly
# reproduces, so it is printed but not compared.

steps <- seq(GBT_TREE_STEP, GBT_NTREES_MAX, by = GBT_TREE_STEP)
grid  <- expand.grid(depth = GBT_DEPTH, lr = GBT_LR)
cat("\nBoosted trees:", nrow(grid), "settings x", K_FOLDS, "folds, grown to",
    GBT_NTREES_MAX, "trees, scored every", GBT_TREE_STEP, "\n")

score_steps <- function(m, k) {   # one fold model, scored at every step
  tr   <- fold != k
  y    <- pa[!tr]
  w    <- wet[!tr]
  p_te <- matrix(predict(m, env[!tr, ], n.trees = steps, type = "response"),
                 ncol = length(steps))
  p_tr <- matrix(predict(m, env[tr & pa == 1, ], n.trees = steps, type = "response"),
                 ncol = length(steps))
  bind_rows(lapply(seq_along(steps), function(s) {
    belt_s <- score_preds(p_te[y == 1 & w, s], p_te[y == 0 & w, s], p_tr[, s])
    names(belt_s) <- paste0(names(belt_s), "_wet")
    cbind(n_trees = steps[s], fold = k, n_test_pres = sum(y == 1),
          score_preds(p_te[y == 1, s], p_te[y == 0, s], p_tr[, s]), belt_s)
  }))
}

gbt_folds <- bind_rows(lapply(seq_len(nrow(grid)), function(g) {
  t0  <- Sys.time()
  out <- bind_rows(lapply(seq_len(K_FOLDS), function(k) {
    m <- fit_gbt(pa[fold != k], env[fold != k, ], grid$lr[g], grid$depth[g], GBT_NTREES_MAX)
    score_steps(m, k)
  }))
  cbi_t <- tapply(out$cbi, out$n_trees, mean)
  cat(sprintf("  lr %.3f depth %d | best %5s trees | CBI %.3f | %.1f min\n",
              grid$lr[g], grid$depth[g], names(which.max(cbi_t)), max(cbi_t, na.rm = TRUE),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  cbind(lr = grid$lr[g], depth = grid$depth[g], out)
}))

gbt_tune <- gbt_folds |>
  group_by(lr, depth, n_trees) |>
  summarise(n_cbi = sum(!is.na(cbi)), cbi_sd = sd(cbi, na.rm = TRUE),
            cbi = mean(cbi, na.rm = TRUE), cbi_wet = mean(cbi_wet, na.rm = TRUE),
            auc = mean(auc), or_10p = mean(or_10p), .groups = "drop") |>
  as.data.frame()
scored <- filter(gbt_tune, n_cbi == K_FOLDS)
stopifnot("No configuration was scored on every fold" = nrow(scored) > 0)
best <- scored |> arrange(desc(cbi), n_trees, depth) |> slice(1)
cat(sprintf("\nSelected: lr %s, depth %d, %d trees | CV CBI %.3f (SD %.3f) | within-belt %.3f | AUC %.3f\n",
            best$lr, best$depth, best$n_trees, best$cbi, best$cbi_sd, best$cbi_wet, best$auc))
cat(sprintf("CV CBI over the %d configurations: median %.3f, max %.3f; %d within 0.01 of the max\n",
            nrow(scored), median(scored$cbi), max(scored$cbi),
            sum(scored$cbi >= max(scored$cbi) - 0.01)))
stopifnot("Selection is at GBT_NTREES_MAX: the grid didn't bracket it; raise GBT_NTREES_MAX" =
  best$n_trees < GBT_NTREES_MAX)

gbt_cfg <- function(p, d) fit_gbt(p, d, best$lr, best$depth, best$n_trees)
gbt_cv  <- cv_fit(pa, env, fold, gbt_cfg, gbt_prob, within = wet)
sel     <- gbt_folds |>
  filter(lr == best$lr, depth == best$depth, n_trees == best$n_trees) |> arrange(fold)
stopifnot("The tuning loop and cv_fit disagree at the selection" =
  isTRUE(all.equal(as.matrix(gbt_cv[, m_cols]), as.matrix(sel[, m_cols]),
                   check.attributes = FALSE)))

# ----------------------- Boosted trees: final model --------------------------

gbt  <- gbt_cfg(pa, env)
infl <- summary(gbt, plotit = FALSE)
cat("\nRelative influence (%):",
    paste(sprintf("%s %.1f", as.character(infl$var), infl$rel.inf), collapse = " | "), "\n")
saveRDS(gbt, GBT_MODEL_FILE)

# Thresholds from cross-fitted predictions: every training point predicted by
# a model fitted without it (GBT_THR_FOLDS random folds within each class),
# the boosting counterpart of the forest's out-of-bag predictions. At its own
# training points boosting predicts close to what it was fitted to, so
# in-sample thresholds are printed as the diagnostic.
set.seed(SEED)
cf_fold <- integer(length(pa))
for (cl in c(1, 0))
  cf_fold[pa == cl] <- sample(rep_len(seq_len(GBT_THR_FOLDS), sum(pa == cl)))
cf <- rep(NA_real_, length(pa))
for (k in seq_len(GBT_THR_FOLDS))
  cf[cf_fold == k] <- gbt_prob(gbt_cfg(pa[cf_fold != k], env[cf_fold != k, ]),
                               env[cf_fold == k, ])
stopifnot("A training point has no cross-fitted prediction" = !anyNA(cf))
gbt_thr    <- thresholds_from(cf[pa == 1], cf[pa == 0])
gbt_thr_in <- thresholds_from(gbt_prob(gbt, occ_env), gbt_prob(gbt, bg_env))
cat(sprintf("GBT thresholds | cross-fitted: p10 %.3f, maxSSS %.3f | in-sample: p10 %.3f, maxSSS %.3f\n",
            gbt_thr[["p10"]], gbt_thr[["maxsss"]], gbt_thr_in[["p10"]], gbt_thr_in[["maxsss"]]))

cat("Predicting the boosted trees over the domain...\n")
gbt_surf <- terra::predict(domain_covs(vars), gbt, fun = gbt_prob, na.rm = TRUE)
names(gbt_surf) <- "suitability"
writeRaster(gbt_surf, GBT_SUIT_FILE, overwrite = TRUE)
gbt_surf <- rast(GBT_SUIT_FILE)   # saved (FLT4S) values
stopifnot("GBT surface covers different cells from SUIT_FILE" =
  global(is.na(gbt_surf) != is.na(suit_r), "sum")[[1]] == 0)
cat("GBT surface range:", paste(round(minmax(gbt_surf)[, 1], 3), collapse = " to "), "\n")
gbt_arp <- arp_estimates(gbt_surf, pop, gbt_thr)
arp_in  <- arp_estimates(gbt_surf, pop, gbt_thr_in)
cat("Diagnostic | GBT binary estimates at in-sample thresholds: p10", fmt(arp_in[["p10"]]),
    "| maxSSS", fmt(arp_in[["maxsss"]]), "\n")

# ------------------------- Surfaces and states -------------------------------
# Spearman rank agreement over every domain cell and within >= 150 mm, as 09.

surfs <- c(suit_r, rf_surf, gbt_surf)
names(surfs) <- c("maxent", "rf", "gbt")
v    <- sapply(names(surfs), function(m) values(surfs[[m]], mat = FALSE))
ok   <- !is.na(v[, "maxent"])
belt <- values(belt_r, mat = FALSE) %in% 1
rho_of <- function(a, b) c(domain = cor(v[ok, a], v[ok, b], method = "spearman"),
                           belt   = cor(v[ok & belt, a], v[ok & belt, b], method = "spearman"))
rho <- list(rf = rho_of("maxent", "rf"), gbt = rho_of("maxent", "gbt"),
            gbt_rf = rho_of("rf", "gbt"))

zones     <- state_zones(suit_r)
rw        <- c(pop * suit_r, pop * rf_surf, pop * gbt_surf)
names(rw) <- c("maxent", "rf", "gbt")
state_df  <- zonal(rw, zones, fun = "sum", na.rm = TRUE)
names(state_df)[1] <- "state"
st_chk <- st_09[match(state_df$state, st_09$state), ]
stopifnot(
  "States do not sum to the national estimates" =
    all(abs(colSums(state_df[c("maxent", "rf", "gbt")]) /
            c(mx_arp[["risk_weighted"]], rf_arp[["risk_weighted"]],
              gbt_arp[["risk_weighted"]]) - 1) < POP_TOL),
  "MaxEnt or RF state totals differ from comparator_by_state.csv (09)" =
    nrow(st_09) == nrow(state_df) && !anyNA(st_chk$state) &&
    isTRUE(all.equal(state_df[c("maxent", "rf")], st_chk[c("maxent", "rf")],
                     check.attributes = FALSE))
)
state_df <- state_df |>
  mutate(rf_pct = 100 * (rf / maxent - 1), gbt_pct = 100 * (gbt / maxent - 1),
         share_maxent = 100 * maxent / sum(maxent), share_rf = 100 * rf / sum(rf),
         share_gbt = 100 * gbt / sum(gbt)) |>
  arrange(desc(maxent))
rho_state <- c(rf  = cor(state_df$maxent, state_df$rf,  method = "spearman"),
               gbt = cor(state_df$maxent, state_df$gbt, method = "spearman"))

# ---------------------------- Response curves -------------------------------
# All three with the others at the presence median, as 07 and 09.
# Descriptive. The rainfall decline is the fragile finding (10, 12, 15, 24):
# suitability at the wettest presence as a percentage of each curve's peak.

gbt_curves <- response_curves(gbt, ref, env, vars, pred = gbt_prob)
curves_l   <- list(maxent = mx_curves, rf = rf_curves, gbt = gbt_curves)
feat <- bind_rows(lapply(names(curves_l), function(m)
  data.frame(model = m, curve_table(curves_l[[m]])))) |>
  arrange(variable, model)
stopifnot("MaxEnt or RF curve features differ from comparator_curve_features.csv (09)" =
  isTRUE(all.equal(filter(feat, model != "gbt"), feat_09, check.attributes = FALSE)))
rain_max <- max(occ_env$rainfall)
at_wet   <- sapply(curves_l, pct_of_peak, var = "rainfall", x = rain_max)
vert     <- sapply(curves_l, function(cv) cv$suit[cv$variable == "vertisols"])

# ------------------------------ Comparison ----------------------------------

se_mx   <- sd(mx_cv$cbi) / sqrt(K_FOLDS)
plateau <- max(abs(cand_08$change_vs_selected[startsWith(cand_08$candidate,
                                                         paste0(tuning$fc, " "))]))
three <- bind_rows(
  comparator_row("maxent", paste(tuning$fc, "rm", tuning$rm), mx_cv, mx_thr, mx_arp),
  comparator_row("rf",  paste("mtry", rf$mtry), rf_cv, rf_thr, rf_arp, rho$rf),
  comparator_row("gbt", sprintf("lr %s depth %d trees %d", best$lr, best$depth, best$n_trees),
                 gbt_cv, gbt_thr, gbt_arp, rho$gbt)
) |>
  vs_primary(se_mx, plateau) |>
  mutate(rho_states = c(NA_real_, rho_state[["rf"]], rho_state[["gbt"]]),
         rain_at_wettest_pct = unname(at_wet[model]))
stopifnot("MaxEnt or RF row differs from COMPARATOR_FILE (09)" =
  isTRUE(all.equal(three[match(cmp_09$model, three$model), names(cmp_09)], cmp_09,
                   check.attributes = FALSE)))
cat("\nMaxEnt and RF reproduce 09 (rows, state totals, curve features)\n")

cat("\nCV CBI by fold (all of Sudan):\n")
data.frame(fold = mx_cv$fold, test_presences = mx_cv$n_test_pres,
           maxent = mx_cv$cbi, rf = rf_cv$cbi, gbt = gbt_cv$cbi) |>
  mutate(across(maxent:gbt, ~ round(., 3))) |> print(row.names = FALSE)

cat("\nFit | one SE below MaxEnt's CBI:", round(three$cv_cbi[1] - se_mx, 3), "\n")
three |>
  select(model, config, cv_cbi, cv_cbi_sd, cv_cbi_wet, cv_auc, cv_or_10p, cbi_within_se) |>
  mutate(across(where(is.double), ~ round(., 3))) |> print(row.names = FALSE)

cat("\nPopulation at risk | plateau spread: +/-", round(plateau, 1), "%\n")
three |>
  select(model, arp_weighted, rw_change_pct, beyond_plateau, arp_maxsss, maxsss, arp_p10, p10) |>
  mutate(across(c(arp_weighted, arp_maxsss, arp_p10), fmt),
         rw_change_pct = round(rw_change_pct, 1),
         across(c(maxsss, p10), ~ round(., 3))) |>
  print(row.names = FALSE)
cat("Algorithm range, risk-weighted:", fmt(min(three$arp_weighted)), "to",
    fmt(max(three$arp_weighted)), "\n")

cat(sprintf(paste0("\nSurface rank agreement (Spearman), domain | within >= 150 mm:\n",
                   "  MaxEnt-RF  %.3f | %.3f\n  MaxEnt-GBT %.3f | %.3f\n  RF-GBT     %.3f | %.3f\n"),
            rho$rf[["domain"]], rho$rf[["belt"]], rho$gbt[["domain"]], rho$gbt[["belt"]],
            rho$gbt_rf[["domain"]], rho$gbt_rf[["belt"]]))
cat("State rank agreement with MaxEnt: RF", round(rho_state[["rf"]], 3),
    "| GBT", round(rho_state[["gbt"]], 3), "\n")

cat("\nRisk-weighted by state (change from MaxEnt, %; share of national, %):\n")
state_df |>
  mutate(across(c(maxent, rf, gbt), fmt), across(rf_pct:share_gbt, ~ round(., 1))) |>
  print(right = FALSE, row.names = FALSE)

cat("\nResponse-curve positions (others at the presence median):\n")
feat |> mutate(across(where(is.numeric), ~ signif(., 3))) |> print(row.names = FALSE)
cat(sprintf("Rainfall at the wettest presence (%.0f mm), %% of peak: MaxEnt %.0f | RF %.0f | GBT %.0f\n",
            rain_max, at_wet[["maxent"]], at_wet[["rf"]], at_wet[["gbt"]]))
cat("Vertisols, suitability off / on:",
    paste(sprintf("%s %.3f / %.3f", colnames(vert), vert[1, ], vert[2, ]), collapse = " | "), "\n")

# ------------------- Held-out prediction at the presences --------------------
# Map suitability (each algorithm's surface) and held-out prediction (the
# fold model fitted without the record's fold) at every presence, by state.
# The transfer reading is in the header.

gbt_ho <- heldout_preds(pa, env, fold, gbt_cfg, gbt_prob)
xy   <- as.matrix(occ[, c("longitude", "latitude")])
map  <- terra::extract(surfs, xy)
pres <- data.frame(coordinate_id = occ$coordinate_id, fold = occ$fold,
                   state = as.character(terra::extract(zones, xy)[, 1]),
                   map_maxent = map$maxent, map_rf = map$rf, map_gbt = map$gbt,
                   ho_maxent = mx_ho$pred, ho_rf = rf_ho$pred, ho_gbt = gbt_ho$pred)
stopifnot("A presence lacks a state or a value" = !anyNA(pres))

p10 <- setNames(three$p10, three$model)
by_state <- pres |>
  group_by(state) |>
  summarise(n = n(), across(c(starts_with("map_"), starts_with("ho_")), median),
            .groups = "drop") |>
  mutate(maxent_fails    = ho_maxent < p10[["maxent"]],
         rf_fitted_only  = maxent_fails & map_rf  >= p10[["rf"]]  & ho_rf  < p10[["rf"]],
         gbt_fitted_only = maxent_fails & map_gbt >= p10[["gbt"]] & ho_gbt < p10[["gbt"]]) |>
  arrange(desc(n)) |>
  as.data.frame()

cat("\nMedian suitability at the presences by state, map and held out | p10:",
    paste(sprintf("%s %.3f", names(p10), p10), collapse = ", "), "\n")
by_state |> mutate(across(where(is.double), ~ round(., 3))) |> print(row.names = FALSE)

fail    <- by_state$state[by_state$maxent_fails]
in_fail <- state_df$state %in% fail
d_fail  <- sapply(c(rf = "rf", gbt = "gbt"), function(m) sum((state_df[[m]] - state_df$maxent)[in_fail]))
d_all   <- sapply(c(rf = "rf", gbt = "gbt"), function(m) sum(state_df[[m]] - state_df$maxent))
cat("States MaxEnt fails to predict held out:", paste(fail, collapse = ", "), "\n")
cat("Change from MaxEnt (risk-weighted) in those states | nationally: RF",
    fmt(d_fail[["rf"]]), "|", fmt(d_all[["rf"]]), "; GBT",
    fmt(d_fail[["gbt"]]), "|", fmt(d_all[["gbt"]]), "\n")

# -------------------------------- Figures -----------------------------------

display  <- st_as_sf(vect(DISPLAY_ADM0_FILE))
excluded <- st_as_sf(vect(EXCLUDED_FILE))
sudan    <- st_as_sf(vect(ADM0_FILE))
states   <- st_as_sf(vect(ADM1_FILE))
lim      <- display_limits(display)

map_panel <- function(r, title) {
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "suitability"
  ggplot() +
    geom_sf(data = sudan, fill = "grey95", colour = NA) +
    geom_raster(data = d, aes(x, y, fill = suitability)) +
    layer_excluded(excluded) +
    scale_fill_suitability(guide = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    layer_admin1(states) +
    layer_country(display) +
    coord_display(lim) +
    labs(title = title) +
    theme_map() +
    theme(plot.title = element_text(size = 9, hjust = 0))
}
panels <- lapply(seq_along(model_labels), function(i) {
  m <- names(model_labels)[i]
  map_panel(surfs[[m]], paste0("(", letters[i], ") ", model_labels[[m]]))
})
p_maps <- (panels[[1]] + panels[[2]] + panels[[3]] + add_scalebar(location="tl")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.25, "cm"), legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))
panel_h <- (FIG_WIDTH_FULL / 3) * diff(lim$y) /
  (diff(lim$x) * cos(mean(lim$y) * pi / 180))
for (ext in c("png", "pdf"))
  save_fig(file.path(DIR_FIGS, paste0("suitability_three_models.", ext)), p_maps,
           width = FIG_WIDTH_FULL, height = panel_h + 2.5)

lab <- cov_labels[vars]
curves <- bind_rows(lapply(names(curves_l), function(m) data.frame(model = m, curves_l[[m]]))) |>
  mutate(model = factor(model, levels = names(model_labels)),
         label = factor(lab[variable], levels = lab))
pres_range <- data.frame(variable = vars, lo = sapply(occ_env, min), hi = sapply(occ_env, max)) |>
  filter(variable != "vertisols") |>
  mutate(label = factor(lab[variable], levels = lab))
p_curves <- ggplot(curves, aes(value, suit, colour = model)) +
  geom_rect(data = pres_range, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey92") +
  geom_line(data = filter(curves, variable != "vertisols")) +
  geom_point(data = filter(curves, variable == "vertisols"),
             position = position_dodge(width = 0.3), size = 2) +
  facet_wrap(~ label, scales = "free_x", nrow = 2) +
  scale_x_continuous(breaks = function(lim)
     if (diff(lim) < 2) c(0, 1) else scales::extended_breaks()(lim)) +
  scale_colour_manual(values = pal_models, labels = model_labels, name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Predicted suitability") +
  theme_dissertation(gridlines = "both") +
  theme(legend.position = "bottom")
for (ext in c("png", "pdf"))
  save_fig(file.path(DIR_FIGS, paste0("response_curves_three_models.", ext)), p_curves,
           width = FIG_WIDTH_FULL, height = 12)
cat("Saved suitability_three_models and response_curves_three_models\n")

# --------------------------------- Save -------------------------------------

write.csv(three,    THREE_MODEL_FILE, row.names = FALSE)
write.csv(gbt_tune, file.path(DIR_TABLES, "gbt_tuning.csv"), row.names = FALSE)
write.csv(mutate(state_df, across(rf_pct:share_gbt, ~ round(., 1))),
          file.path(DIR_TABLES, "three_model_by_state.csv"), row.names = FALSE)
write.csv(feat,     file.path(DIR_TABLES, "three_model_curve_features.csv"), row.names = FALSE)
write.csv(pres,     file.path(DIR_TABLES, "three_model_heldout.csv"), row.names = FALSE)
write.csv(by_state, file.path(DIR_TABLES, "three_model_heldout_by_state.csv"), row.names = FALSE)
cat("\n17_gbt_comparator.R complete\n")