# ============================================================================
# 19_hindcast.R
# The model's figure for Gedaref beside the only earlier population-at-risk
# figure for the state: Alvar et al. (2006), ALVAR_STATE_ARP in HINDCAST_YEAR,
# from unpublished Federal Ministry of Health data with "at risk" undefined.
# Both surfaces carry HINDCAST_YEAR population: the long-term surface (07),
# which is the model as reported, and the HINDCAST_YEAR surface (the
# dissertation's comparison), which adds one year's weather (18). National
# figures are context only.
#
# Population: WorldPop has no constrained surface before 2015, so
# POP_HINDCAST_FILE is the unconstrained, UN-adjusted product, which spreads
# people over all land; POP_FILE (08) is constrained. Within a population year
# the two surfaces compare; across population years, growth, redistribution
# and the change of product are mixed.
#
# Inputs:  MODEL_FILE, TUNING_FILE, TRAIN_FILE, SUIT_FILE, POP_ALIGNED_FILE,
#          POP_HINDCAST_FILE (downloaded once from POP_HINDCAST_URL),
#          DOMAIN_FILE, ADM1_FILE, retained_vars.rds, COV_FILES, COV_ANNUAL,
#          arp_summary.csv (08), single_year_estimates.csv (18: run 18 first)
# Outputs: outputs/tables/hindcast_national.csv, hindcast_gedaref.csv,
#          arp_literature_context.csv
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(maxnet)
})

# ------------------------------ Load inputs ---------------------------------

mod     <- readRDS(MODEL_FILE)
tuning  <- readRDS(TUNING_FILE)
train   <- readRDS(TRAIN_FILE)
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
suit_r  <- rast(SUIT_FILE)
pop_now <- rast(POP_ALIGNED_FILE)
arp_08  <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))
yr_18   <- read.csv(file.path(DIR_TABLES, "single_year_estimates.csv"))

hc_lab   <- as.character(HINDCAST_YEAR)
hc_files <- year_files(HINDCAST_YEAR)[vars]
stopifnot(
  "HINDCAST_YEAR falls in Terra's drift period" = HINDCAST_YEAR < TERRA_DRIFT_FROM,
  "Missing annual rasters for HINDCAST_YEAR"    = all(file.exists(file.path(DIR_COVARIATES, hc_files))),
  "HINDCAST_YEAR is not in 18's table (run 18)" = sum(yr_18$year %in% HINDCAST_YEAR) == 1,
  "Population grid differs from the surface"    = compareGeom(pop_now, suit_r, stopOnError = FALSE)
)
cat("Model:", tuning$fc, "rm", tuning$rm, "| hindcast year:", HINDCAST_YEAR, "| annual files:",
    paste(hc_files[intersect(vars, names(COV_ANNUAL))], collapse = ", "), "\n")

# ----------- The prediction code reproduces 07 (surface) and 08 -------------

cat("\nLong-term means: predicting...\n")
ltm_r   <- terra::predict(domain_covs(vars), mod, type = "cloglog", na.rm = TRUE)
thr     <- thresholds_from(maxnet_prob(mod, train$occ_env[, vars]),
                           maxnet_prob(mod, train$bg_env[, vars]))
thr_08  <- setNames(arp_08$threshold, arp_08$metric)
arp_08v <- setNames(arp_08$arp,       arp_08$metric)
arp_lt  <- arp_estimates(suit_r, pop_now, thr)
stopifnot(
  "Surface differs from SUIT_FILE (07)" =
    global(abs(ltm_r - suit_r), "max", na.rm = TRUE)[[1]] < 1e-6 &&
    global(is.na(ltm_r) != is.na(suit_r), "sum")[[1]] == 0,
  "Thresholds differ from arp_summary.csv (08)" = all(abs(thr - thr_08[names(thr)]) < 1e-6),
  "Estimates differ from arp_summary.csv (08)"  = all(round(arp_lt) == arp_08v[names(arp_lt)])
)
cat("Prediction code reproduces 07 (surface) and 08 (thresholds, estimates)\n")
cat("p10:", round(thr[["p10"]], 3), "| maxSSS:", round(thr[["maxsss"]], 3), "\n")

# ---------------- The hindcast-year surface reproduces 18 --------------------
# 18 predicts the same model onto each year with POP_FILE population, on the
# cells the long-term surface covers. With that population, 19's surface must
# give 18's row exactly.

cat("\n", hc_lab, " covariates: predicting...\n", sep = "")
hc_r <- mask(terra::predict(domain_covs(vars, year_files(HINDCAST_YEAR)), mod,
                            type = "cloglog", na.rm = TRUE), suit_r)
row_18 <- yr_18[yr_18$year %in% HINDCAST_YEAR, ]
lt_18  <- yr_18[yr_18$group == "long-term mean", ]
r18    <- unlist(row_18[c("risk_weighted", "p10", "maxsss")])
arp_hc_now <- arp_estimates(hc_r, pop_now, thr)
stopifnot(
  "Hindcast surface lacks cells the long-term surface covers" =
    global(is.na(hc_r) != is.na(suit_r), "sum")[[1]] == 0,
  "Hindcast surface with POP_FILE population differs from 18's row" =
    all(abs(arp_hc_now[names(r18)] / r18 - 1) < 1e-9)
)
cat(hc_lab, "surface reproduces 18 (cells; estimates with", POP_YEAR, "population)\n")
cat(sprintf("%s against the long-term means (18, domain): night LST %+.2f C, rainfall %+.0f mm\n",
            hc_lab, row_18$lst_night - lt_18$lst_night, row_18$rainfall - lt_18$rainfall))

# ---------------------------- Hindcast population ----------------------------
# Outside the domain is the Halaib Triangle, by design. Inside the domain
# without a prediction are cells lacking a covariate (11: Red Sea coast).

cat("\nWorldPop", hc_lab, ": aligning...\n")
download_once(POP_HINDCAST_FILE, POP_HINDCAST_URL)
al_hc  <- align_pop(rast(POP_HINDCAST_FILE), suit_r)
pop_hc <- al_hc$pop
tot_hc <- global(pop_hc, "sum", na.rm = TRUE)[[1]]
in_dom <- global(mask(pop_hc, rast(DOMAIN_FILE), maskvalues = 0), "sum", na.rm = TRUE)[[1]]
in_srf <- global(mask(pop_hc, suit_r), "sum", na.rm = TRUE)[[1]]
cat(sprintf(paste0("WorldPop %s (unconstrained): %s | outside the domain %s | in the domain ",
                   "without a prediction %s | with a prediction %s\n"),
            hc_lab, fmt(al_hc$raw_total), fmt(tot_hc - in_dom), fmt(in_dom - in_srf), fmt(in_srf)))

# ------------------------- National (context only) --------------------------

nat_row <- function(surface, pop_year, surf_r, pop_r) {
  a   <- arp_estimates(surf_r, pop_r, thr)
  cov <- global(mask(pop_r, surf_r), "sum", na.rm = TRUE)[[1]]
  data.frame(surface = surface, population = pop_year, pop_covered = cov,
             risk_weighted = a[["risk_weighted"]], maxsss = a[["maxsss"]], p10 = a[["p10"]],
             rw_pct_of_pop = 100 * a[["risk_weighted"]] / cov)
}
national <- bind_rows(
  nat_row("long-term", POP_YEAR,      suit_r, pop_now),
  nat_row("long-term", HINDCAST_YEAR, suit_r, pop_hc),
  nat_row(hc_lab,      POP_YEAR,      hc_r,   pop_now),
  nat_row(hc_lab,      HINDCAST_YEAR, hc_r,   pop_hc)) |>
  group_by(population) |>
  mutate(rw_pct_vs_long_term = 100 * (risk_weighted / risk_weighted[surface == "long-term"] - 1)) |>
  ungroup() |> as.data.frame()

cat("\nNational (context only; across population years, growth, redistribution and",
    "the change of WorldPop product are mixed):\n")
national |>
  mutate(across(c(pop_covered, risk_weighted, maxsss, p10), fmt),
         across(c(rw_pct_of_pop, rw_pct_vs_long_term), ~ round(., 1))) |>
  print(row.names = FALSE)

# ---------------------------------- Gedaref ----------------------------------

zones <- state_zones(suit_r)
stopifnot("ALVAR_STATE is not a state in ADM1_FILE" = ALVAR_STATE %in% levels(zones)[[1]][[2]])
occ_state <- as.character(terra::extract(
  zones, as.matrix(train$occ_clean[, c("longitude", "latitude")]))[, 1])
cat("\nTraining presences in", ALVAR_STATE, ":", sum(occ_state == ALVAR_STATE, na.rm = TRUE),
    "of", nrow(train$occ_clean), "\n")

lay <- c(pop_hc, pop_hc * suit_r, pop_hc * (suit_r >= thr[["maxsss"]]),
         pop_hc * hc_r, pop_hc * (hc_r >= thr[["maxsss"]]), pop_now)
names(lay) <- c("pop_hc", "rw_ltm", "maxsss_ltm", "rw_hc", "maxsss_hc", "pop_now")
st <- zonal(lay, zones, fun = "sum", na.rm = TRUE)
names(st)[1] <- "state"
nat_hc <- national[national$population == HINDCAST_YEAR, ]
stopifnot("State estimates do not sum to the national estimates" =
  abs(sum(st$rw_ltm) / nat_hc$risk_weighted[nat_hc$surface == "long-term"] - 1) < POP_TOL &&
  abs(sum(st$rw_hc)  / nat_hc$risk_weighted[nat_hc$surface == hc_lab]      - 1) < POP_TOL)
g <- st[st$state == ALVAR_STATE, ]

gedaref <- data.frame(
  source   = c(rep("This model, long-term surface", 2),
               rep(paste0("This model, ", hc_lab, " surface"), 2), "Alvar et al. 2006"),
  measure  = c(rep(c("risk-weighted index", "people in cells >= maxSSS"), 2),
               "people at risk (undefined)"),
  estimate = c(g$rw_ltm, g$maxsss_ltm, g$rw_hc, g$maxsss_hc, ALVAR_STATE_ARP)) |>
  mutate(pct_of_state_pop = 100 * estimate / g$pop_hc,
         ratio_to_alvar   = estimate / ALVAR_STATE_ARP)

cat(sprintf("%s population (WorldPop): %s %s | %d %s\n", ALVAR_STATE,
            hc_lab, fmt(g$pop_hc), POP_YEAR, fmt(g$pop_now)))
cat("With", hc_lab, "population (descriptive; reading in the header):\n")
gedaref |>
  mutate(estimate = fmt(estimate), pct_of_state_pop = round(pct_of_state_pop, 1),
         ratio_to_alvar = round(ratio_to_alvar, 2)) |>
  print(row.names = FALSE)

# --------------------------------- Save -------------------------------------

write.csv(national, file.path(DIR_TABLES, "hindcast_national.csv"),      row.names = FALSE)
write.csv(gedaref,  file.path(DIR_TABLES, "hindcast_gedaref.csv"),       row.names = FALSE)
write.csv(lit,      file.path(DIR_TABLES, "arp_literature_context.csv"), row.names = FALSE)
cat("\n19_hindcast.R complete\n")