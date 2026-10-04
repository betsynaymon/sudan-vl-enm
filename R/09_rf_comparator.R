# ============================================================================
# 09_rf_comparator.R
# Does the surface or the estimate depend on the algorithm? Fits a random
# forest to the same presences, background, covariates and spatial folds as
# MaxEnt and compares the two on cross-validated fit, surface ranks,
# response-curve positions and population at risk. The forest is
# down-sampled: every tree draws as many background points as presences
# (fit_rf, helpers.R; Valavi et al. 2021). MaxEnt runs through the same code
# first and must reproduce 06-08.
#
# Inputs:  TRAIN_FILE, TUNING_FILE, MODEL_FILE, SUIT_FILE, POP_ALIGNED_FILE,
#          SENS_MASK_FILE, DOMAIN_FILE, ADM1_FILE, retained_vars.rds,
#          tuning_by_fold.csv (06), response_curve_features.csv (07),
#          arp_summary.csv and arp_plateau_candidates.csv (08), COV_FILES
# Outputs: RF_MODEL_FILE, RF_SUIT_FILE, COMPARATOR_FILE
#          outputs/tables/rf_tuning.csv, comparator_by_state.csv,
#          comparator_curve_features.csv
# Figures are drawn in 17, with the gradient boosted trees.
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(maxnet); library(ranger)
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
feat_07 <- read.csv(file.path(DIR_TABLES, "response_curve_features.csv"))
arp_08  <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))
cand_08 <- read.csv(file.path(DIR_TABLES, "arp_plateau_candidates.csv"))

occ     <- train$occ_clean
bg      <- train$bg_clean
occ_env <- train$occ_env[, vars]
bg_env  <- train$bg_env[, vars]

# Values and points were saved together by 06; each row carries its own fold
# and its >= 150 mm flag from there.
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
  "Population grid differs from the surface" = compareGeom(pop, suit_r, stopOnError = FALSE),
  "Belt mask grid differs from the surface"  = compareGeom(belt_r, suit_r, stopOnError = FALSE)
)
cat("Presences:", sum(pa == 1), "| background:", sum(pa == 0),
    "| test presences per fold:", table(fold[pa == 1]), "\n")
cat("ranger", as.character(packageVersion("ranger")), "\n")

# ------------------------ MaxEnt reproduces 06-08 ---------------------------
# Through the same code as the forest below (cv_fit, score_preds,
# response_curves, curve_features). If this stops, the comparison code
# differs from 06-08 and nothing below is comparable.

cat("\nMaxEnt", tuning$fc, "rm", tuning$rm, ": fit, CV and surface...\n")
covs      <- domain_covs(vars)
ref       <- sapply(occ_env, median)
mx_mod    <- fit_maxnet(pa, env, tuning$fc, tuning$rm)
mx_cv     <- cv_maxnet(pa, env, fold, tuning$fc, tuning$rm, within = wet)
mx_thr    <- thresholds_from(maxnet_prob(mx_mod, occ_env), maxnet_prob(mx_mod, bg_env))
mx_surf   <- terra::predict(covs, mx_mod, type = "cloglog", na.rm = TRUE)
mx_curves <- response_curves(mx_mod, ref, env, vars)
mx_arp    <- arp_estimates(suit_r, pop, mx_thr)

m_cols  <- c("auc", "cbi", "or_10p", "auc_wet", "cbi_wet", "or_10p_wet")
f06     <- fold_06 |> filter(fc == tuning$fc, rm == tuning$rm) |> arrange(fold)
thr_08  <- setNames(arp_08$threshold, arp_08$metric)
arp_08v <- setNames(arp_08$arp,       arp_08$metric)
stopifnot(
  "Coefficients differ from MODEL_FILE (06)" = isTRUE(all.equal(mx_mod$betas, mod$betas)),
  "CV metrics by fold differ from tuning_by_fold.csv (06)" =
    nrow(f06) == K_FOLDS &&
    isTRUE(all.equal(as.matrix(mx_cv[, m_cols]), as.matrix(f06[, m_cols]),
                     check.attributes = FALSE)),
  "CV CBI differs from 06" = abs(mean(mx_cv$cbi) - tuning$cbi) < 1e-8,
  "Surface differs from SUIT_FILE (07)" =
    global(abs(mx_surf - suit_r), "max", na.rm = TRUE)[[1]] < 1e-6 &&
    global(is.na(mx_surf) != is.na(suit_r), "sum")[[1]] == 0,
  "Curve features differ from response_curve_features.csv (07)" =
    isTRUE(all.equal(as.data.frame(curve_table(mx_curves)), feat_07, check.attributes = FALSE)),
  "Thresholds differ from arp_summary.csv (08)" = all(abs(mx_thr - thr_08[names(mx_thr)]) < 1e-6),
  "Estimates differ from arp_summary.csv (08)"  = all(round(mx_arp) == arp_08v[names(mx_arp)])
)
cat("MaxEnt (", length(mx_mod$betas), " features) reproduces 06 (coefficients, CV by fold),",
    " 07 (surface, curve features) and 08 (thresholds, estimates)\n", sep = "")

# ------------------------- Random forest: tune mtry -------------------------
# Same rows, folds and metrics as MaxEnt. mtry from 1 to the number of
# covariates; selected: highest mean CV CBI among settings scored on every
# fold, ties to the smaller mtry. OR10 uses each fold model's predictions at
# its own training presences, which a forest partly reproduces, so the
# forest's OR10 is printed but not compared.

cat("\nRandom forest CV,", RF_NTREES, "trees per forest:\n")
rf_cvs <- lapply(seq_along(vars), function(m) {
  cv <- cv_fit(pa, env, fold, function(p, d) fit_rf(p, d, m), rf_prob, within = wet)
  cat(sprintf("  mtry %d | CBI %.3f (SD %.3f) | within-belt CBI %.3f | AUC %.3f\n",
              m, mean(cv$cbi), sd(cv$cbi), mean(cv$cbi_wet), mean(cv$auc)))
  cv
})
rf_tune <- bind_rows(lapply(seq_along(rf_cvs), function(m) {
  cv <- rf_cvs[[m]]
  data.frame(mtry = m, n_cbi = sum(!is.na(cv$cbi)),
             cbi = mean(cv$cbi, na.rm = TRUE), cbi_sd = sd(cv$cbi, na.rm = TRUE),
             cbi_wet = mean(cv$cbi_wet, na.rm = TRUE), auc = mean(cv$auc),
             or_10p = mean(cv$or_10p))
}))
scored <- filter(rf_tune, n_cbi == K_FOLDS)
stopifnot("No mtry was scored on every fold" = nrow(scored) > 0)
best_mtry <- scored$mtry[which.max(scored$cbi)]
rf_cv     <- rf_cvs[[best_mtry]]
cat("Selected mtry:", best_mtry, "\n")

# ------------------------ Random forest: final model ------------------------
# Every tree must draw exactly as many presences and as many background
# points as there are presences: the sampling the comparison assumes.

rf    <- fit_rf(pa, env, best_mtry, keep.inbag = TRUE)
inbag <- do.call(cbind, rf$inbag.counts)
n1    <- sum(pa == 1)
stopifnot("A tree's sample is not balanced" =
  all(colSums(inbag[pa == 1, ]) == n1) && all(colSums(inbag[pa == 0, ]) == n1))
cat("Every tree drew", n1, "presences and", n1, "background points;",
    "distinct presences per tree, median", median(colSums(inbag[pa == 1, ] > 0)), "\n")

# Thresholds from out-of-bag predictions: at its own training points a forest
# predicts close to what it was grown on. In-sample thresholds are printed as
# the diagnostic.
oob <- rf$predictions[, "1"]
stopifnot("A training point was never out of bag" = !anyNA(oob))
rf_thr    <- thresholds_from(oob[pa == 1], oob[pa == 0])
rf_thr_in <- thresholds_from(rf_prob(rf, occ_env), rf_prob(rf, bg_env))
cat(sprintf("RF thresholds | out-of-bag: p10 %.3f, maxSSS %.3f | in-sample: p10 %.3f, maxSSS %.3f\n",
            rf_thr[["p10"]], rf_thr[["maxsss"]], rf_thr_in[["p10"]], rf_thr_in[["maxsss"]]))

rf$inbag.counts <- NULL
saveRDS(rf, RF_MODEL_FILE)

cat("Predicting the forest over the domain...\n")
rf_surf <- terra::predict(covs, rf, fun = rf_prob, na.rm = TRUE)
names(rf_surf) <- "suitability"
writeRaster(rf_surf, RF_SUIT_FILE, overwrite = TRUE)
rf_surf <- rast(RF_SUIT_FILE)   # saved (FLT4S) values, as 17 reads them
stopifnot("RF surface covers different cells from SUIT_FILE" =
  global(is.na(rf_surf) != is.na(suit_r), "sum")[[1]] == 0)
cat("RF surface range:", paste(round(minmax(rf_surf)[, 1], 3), collapse = " to "), "\n")

# ------------------------------ Comparison ----------------------------------
se_mx   <- sd(mx_cv$cbi) / sqrt(K_FOLDS)
plateau <- max(abs(cand_08$change_vs_selected[startsWith(cand_08$candidate, paste0(tuning$fc, " "))]))
rf_arp  <- arp_estimates(rf_surf, pop, rf_thr)

v_mx <- values(suit_r,  mat = FALSE)
v_rf <- values(rf_surf, mat = FALSE)
ok   <- !is.na(v_mx)
belt <- values(belt_r, mat = FALSE) %in% 1
rho  <- c(domain = cor(v_mx[ok], v_rf[ok], method = "spearman"),
          belt   = cor(v_mx[ok & belt], v_rf[ok & belt], method = "spearman"))

summary_df <- bind_rows(
  comparator_row("maxent", paste(tuning$fc, "rm", tuning$rm), mx_cv, mx_thr, mx_arp),
  comparator_row("rf", paste("mtry", best_mtry), rf_cv, rf_thr, rf_arp, rho)
) |> vs_primary(se_mx, plateau)

cat("\nCV CBI by fold (all of Sudan; within >= 150 mm):\n")
data.frame(fold = mx_cv$fold, test_presences = mx_cv$n_test_pres,
           maxent = mx_cv$cbi, rf = rf_cv$cbi,
           maxent_belt = mx_cv$cbi_wet, rf_belt = rf_cv$cbi_wet) |>
  mutate(across(maxent:rf_belt, ~ round(., 3))) |> print(row.names = FALSE)

cat("\nFit | one SE below MaxEnt's CBI:", round(summary_df$cv_cbi[1] - se_mx, 3), "\n")
summary_df |>
  select(model, config, cv_cbi, cv_cbi_sd, cv_cbi_wet, cv_auc, cv_or_10p, cbi_within_se) |>
  mutate(across(where(is.double), ~ round(., 3))) |> print(row.names = FALSE)

cat("\nPopulation at risk | plateau spread: +/-", round(plateau, 1), "%\n")
summary_df |>
  select(model, arp_weighted, rw_change_pct, beyond_plateau, arp_maxsss, maxsss, arp_p10, p10) |>
  mutate(across(c(arp_weighted, arp_maxsss, arp_p10), fmt),
         rw_change_pct = round(rw_change_pct, 1),
         across(c(maxsss, p10), ~ round(., 3))) |>
  print(row.names = FALSE)

cat(sprintf("\nSurface rank agreement (Spearman): domain %.3f | within >= 150 mm %.3f\n",
            rho[["domain"]], rho[["belt"]]))

# -------------------------------- States ------------------------------------

zones     <- state_zones(suit_r)
rw        <- c(pop * suit_r, pop * rf_surf)
names(rw) <- c("maxent", "rf")
state_df  <- zonal(rw, zones, fun = "sum", na.rm = TRUE)
names(state_df)[1] <- "state"
stopifnot("States do not sum to the national estimates" =
  all(abs(colSums(state_df[c("maxent", "rf")]) / summary_df$arp_weighted - 1) < POP_TOL))
state_df <- state_df |>
  mutate(rf_pct       = round(100 * (rf / maxent - 1), 1),
         share_maxent = round(100 * maxent / sum(maxent), 1),
         share_rf     = round(100 * rf / sum(rf), 1)) |>
  arrange(desc(maxent))
rho_state <- cor(state_df$maxent, state_df$rf, method = "spearman")

cat("\nRisk-weighted by state (RF change, %; share of national, %) | state rank agreement",
    round(rho_state, 3), "\n")
state_df |> mutate(across(c(maxent, rf), fmt)) |> print(right = FALSE, row.names = FALSE)

# ---------------------------- Response curves -------------------------------
# Both models with the others at the presence median, as 07. Descriptive. The
# rainfall decline is the fragile finding (10, 12, 15, 24): for each model,
# suitability at the wettest presence as a percentage of its peak.

rf_curves <- response_curves(rf, ref, env, vars, pred = rf_prob)
feat <- bind_rows(data.frame(model = "maxent", curve_table(mx_curves)),
                  data.frame(model = "rf",     curve_table(rf_curves))) |>
  arrange(variable, model)
cat("\nResponse-curve positions (others at the presence median):\n")
feat |> mutate(across(where(is.numeric), ~ signif(., 3))) |> print(row.names = FALSE)

rain_max <- max(occ_env$rainfall)
at_wet <- sapply(list(maxent = mx_curves, rf = rf_curves), pct_of_peak,
                 var = "rainfall", x = rain_max)
cat(sprintf("Rainfall at the wettest presence (%.0f mm), %% of peak: MaxEnt %.0f | RF %.0f\n",
            rain_max, at_wet[["maxent"]], at_wet[["rf"]]))
vert <- sapply(list(maxent = mx_curves, rf = rf_curves),
               function(cv) cv$suit[cv$variable == "vertisols"])
cat("Vertisols, suitability off / on: MaxEnt", round(vert[, "maxent"], 3),
    "| RF", round(vert[, "rf"], 3), "\n")

# --------------------------------- Save -------------------------------------

summary_df$rho_states          <- c(NA_real_, rho_state)
summary_df$rain_at_wettest_pct <- unname(at_wet[summary_df$model])
write.csv(summary_df, COMPARATOR_FILE, row.names = FALSE)
write.csv(rf_tune,  file.path(DIR_TABLES, "rf_tuning.csv"), row.names = FALSE)
write.csv(state_df, file.path(DIR_TABLES, "comparator_by_state.csv"), row.names = FALSE)
write.csv(feat,     file.path(DIR_TABLES, "comparator_curve_features.csv"), row.names = FALSE)
cat("\n09_rf_comparator.R complete\n")