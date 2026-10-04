# ============================================================================
# 10_precision_sensitivity.R
# Does the population at risk rest on the weakest records? Refits the selected
# model without each set of questionable presences (background, folds and
# configuration unchanged) and compares fit, surface and population at risk
# with the primary model, nationally and by state.
#   precision:      map-digitised records with positional error > 4 km
#   facility:       records locating a treatment facility, not a site
#   khartoum:       Khartoum case records
#   central_cities: Khartoum and Wad Madani case records
#   id_121:         ID 121 alone (in the precision set; the only dry-riverine
#                   presence), to separate its effect from positional error
# The primary model runs through the same code and must reproduce 06-08.
#
# Inputs:  TRAIN_FILE, TUNING_FILE, MODEL_FILE, SUIT_FILE, POP_ALIGNED_FILE,
#          DOMAIN_FILE, ADM1_FILE, SENS_MASK_FILE, retained_vars.rds,
#          outputs/tables/arp_summary.csv, arp_plateau_candidates.csv (08)
# Outputs: DQ_SURFACES_FILE
#          outputs/tables/data_quality_refits.csv
#          outputs/tables/data_quality_refits_by_state.csv
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(maxnet)
})

# ------------------------------ Load inputs ---------------------------------

train  <- readRDS(TRAIN_FILE)
tuning <- readRDS(TUNING_FILE)
mod    <- readRDS(MODEL_FILE)
suit_r <- rast(SUIT_FILE)
pop    <- rast(POP_ALIGNED_FILE)
vars   <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
arp_08 <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))
cand   <- read.csv(file.path(DIR_TABLES, "arp_plateau_candidates.csv"))

occ <- train$occ_clean
bg  <- train$bg_clean
stopifnot(
  "Presence rows and values differ"          = nrow(occ) == nrow(train$occ_env),
  "Background rows and values differ"        = nrow(bg)  == nrow(train$bg_env),
  "Population grid differs from the surface" = compareGeom(pop, suit_r, stopOnError = FALSE),
  "Belt mask grid differs from the surface" = compareGeom(rast(SENS_MASK_FILE), suit_r, stopOnError = FALSE)
)
cat("Model:", tuning$fc, "rm =", tuning$rm, "| presences:", nrow(occ),
    "| background:", nrow(bg), "\n")

zones <- state_zones(suit_r)
occ$state_zone <- as.character(
  terra::extract(zones, as.matrix(occ[, c("longitude", "latitude")]))[, 1])

# ------------------------------ Refit sets ----------------------------------
# Sets of coordinate_id, selected by ID. "primary" drops
# nothing and must reproduce 06-08.

drop_sets <- list(
  primary        = integer(0),
  precision      = occ$coordinate_id[occ$source %in% IMPRECISE_SOURCES],
  facility       = occ$coordinate_id[occ$presence_type == FACILITY_TYPE],
  khartoum       = KHARTOUM_CASE_IDS,
  central_cities = c(KHARTOUM_CASE_IDS, CENTRAL_CITY_CASE_IDS),
  id_121         = KHARTOUM_VECTOR_ID
)

stopifnot(
  "A named ID is not a training presence" =
    all(c(KHARTOUM_CASE_IDS, CENTRAL_CITY_CASE_IDS, KHARTOUM_VECTOR_ID) %in% occ$coordinate_id),
  "An imprecise source matches no training presence" = all(IMPRECISE_SOURCES %in% occ$source),
  "No facility records"                = length(drop_sets$facility) > 0,
  "ID 121 is not in the precision set" = KHARTOUM_VECTOR_ID %in% drop_sets$precision
)

for (s in names(drop_sets)[-1]) {
  cat("\n", s, ": ", length(drop_sets[[s]]), " record(s)\n", sep = "")
  occ |> filter(coordinate_id %in% drop_sets[[s]]) |>
    select(coordinate_id, source, presence_type, georef_method, year, state_zone, fold) |>
    print(row.names = FALSE)
}
cat("\nIn both the precision and facility sets:",
    length(intersect(drop_sets$precision, drop_sets$facility)), "\n")

pres_left <- sapply(drop_sets, function(ids)
  table(factor(occ$fold[!occ$coordinate_id %in% ids], levels = seq_len(K_FOLDS))))
cat("\nTest presences per fold after each drop (floor ", MIN_TEST_PRES, "):\n", sep = "")
print(pres_left)

# --------------------------- Refit and score --------------------------------
# Background, folds (each row's own) and configuration unchanged; predicted
# as in 07.

covs_dom <- domain_covs(vars)

refits <- setNames(lapply(names(drop_sets), function(s) {
  keep <- !occ$coordinate_id %in% drop_sets[[s]]
  r <- refit_maxnet(train$occ_env[keep, vars], train$bg_env[, vars],
                    occ$fold[keep], bg$fold, tuning$fc, tuning$rm)
  r$surf <- terra::predict(covs_dom, r$mod, type = "cloglog", na.rm = TRUE)
  r$arp  <- arp_estimates(r$surf, pop, r$thr)
  r$n    <- sum(keep)
  cat("  ", s, ": ", r$n, " presences\n", sep = "")
  r
}), names(drop_sets))

# ------------------------ Primary reproduces 06-08 --------------------------
# If this stops, the refit code differs from 06-08 and nothing below is
# comparable.

p       <- refits$primary
thr_08  <- setNames(arp_08$threshold, arp_08$metric)
arp_08v <- setNames(arp_08$arp,       arp_08$metric)
stopifnot(
  "Coefficients differ from MODEL_FILE (06)" = isTRUE(all.equal(p$mod$betas, mod$betas)),
  "CV CBI differs from 06"                   = abs(mean(p$cv$cbi) - tuning$cbi) < 1e-8,
  "Surface differs from SUIT_FILE (07)" =
    global(abs(p$surf - suit_r), "max", na.rm = TRUE)[[1]] < 1e-6 &&
    global(is.na(p$surf) != is.na(suit_r), "sum")[[1]] == 0,
  "Thresholds differ from arp_summary.csv (08)" = all(abs(p$thr - thr_08[names(p$thr)]) < 1e-6),
  "Estimates differ from arp_summary.csv (08)" =
    all(round(arp_estimates(suit_r, pop, p$thr)) == arp_08v[names(p$arp)])
)
cat("\nPrimary reproduces 06 (coefficients, CV CBI), 07 (surface), 08 (thresholds, estimates)\n")

# ------------------------------ Comparison ----------------------------------
# Rules fixed before the refits were run: a fold's CBI counts only with at
# least MIN_TEST_PRES test presences; CBI changes are read against the
# primary's SE; national changes against the plateau spread from 08 (other
# candidates with the selected feature classes), the model-choice
# uncertainty already reported.

se_cbi  <- sd(p$cv$cbi) / sqrt(K_FOLDS)
plateau <- max(abs(cand$change_vs_selected[startsWith(cand$candidate, paste0(tuning$fc, " "))]))
v_prim  <- values(p$surf, mat = FALSE)
belt    <- values(rast(SENS_MASK_FILE), mat = FALSE) %in% 1

summary_df <- bind_rows(lapply(names(refits), function(s) {
  r  <- refits[[s]]
  u  <- r$cv$n_test_pres >= MIN_TEST_PRES
  v  <- values(r$surf, mat = FALSE)
  ok <- !is.na(v) & !is.na(v_prim)
  data.frame(
    refit = s, n_dropped = nrow(occ) - r$n, n_presences = r$n,
    n_coef = length(r$mod$betas), folds_scored = sum(u),
    cbi = mean(r$cv$cbi[u], na.rm = TRUE), cbi_sd = sd(r$cv$cbi[u], na.rm = TRUE),
    auc = mean(r$cv$auc[u]),
    rho_domain = cor(v[ok], v_prim[ok], method = "spearman"),
    rho_belt   = cor(v[ok & belt], v_prim[ok & belt], method = "spearman"),
    p10 = r$thr[["p10"]], maxsss = r$thr[["maxsss"]],
    arp_weighted = r$arp[["risk_weighted"]],
    arp_maxsss = r$arp[["maxsss"]], arp_p10 = r$arp[["p10"]])
})) |>
  mutate(rw_change_pct  = 100 * (arp_weighted / arp_weighted[refit == "primary"] - 1),
         cbi_change     = cbi - cbi[refit == "primary"],
         beyond_plateau = abs(rw_change_pct) > plateau)

cat("\nFit and surface | primary CBI SE:", round(se_cbi, 3), "\n")
summary_df |>
  select(refit, n_dropped, folds_scored, cbi, cbi_change, cbi_sd, auc, n_coef,
         rho_domain, rho_belt) |>
  mutate(across(where(is.double), ~ round(., 3))) |> print(row.names = FALSE)

cat("\nCBI by fold (NA below the floor):\n")
print(sapply(refits, function(r)
  round(ifelse(r$cv$n_test_pres >= MIN_TEST_PRES, r$cv$cbi, NA), 3)))

cat("\nPopulation at risk, own thresholds | plateau spread: +/-", round(plateau, 1), "%\n")
summary_df |>
  select(refit, arp_weighted, rw_change_pct, beyond_plateau, arp_maxsss, maxsss, arp_p10, p10) |>
  mutate(across(starts_with("arp"), fmt), rw_change_pct = round(rw_change_pct, 1),
         across(c(maxsss, p10), ~ round(., 3))) |>
  print(row.names = FALSE)

# -------------------------------- States ------------------------------------

surf_all <- do.call(c, unname(lapply(refits, `[[`, "surf")))
names(surf_all) <- names(refits)
rw_all <- pop * surf_all
names(rw_all) <- names(refits)

state_df <- zonal(rw_all, zones, fun = "sum", na.rm = TRUE)
names(state_df)[1] <- "state"
stopifnot("States do not sum to the national estimates" =
  all(abs(colSums(state_df[names(refits)]) / summary_df$arp_weighted - 1) < POP_TOL))

state_df <- state_df |>
  mutate(n_presences = as.integer(table(factor(occ$state_zone, levels = state))),
         across(all_of(names(refits)[-1]), ~ round(100 * (. / primary - 1), 1),
                .names = "{.col}_pct")) |>
  arrange(desc(primary))

cat("\nRisk-weighted estimate by state (primary; % change under each refit):\n")
state_df |> select(state, n_presences, primary, ends_with("_pct")) |>
  mutate(primary = fmt(primary)) |> print(right = FALSE, row.names = FALSE)

# --------------------------------- Save -------------------------------------

writeRaster(surf_all, DQ_SURFACES_FILE,
            overwrite = TRUE)
write.csv(summary_df, file.path(DIR_TABLES, "data_quality_refits.csv"), row.names = FALSE)
write.csv(state_df,   file.path(DIR_TABLES, "data_quality_refits_by_state.csv"), row.names = FALSE)
cat("10_precision_sensitivity.R complete\n")

