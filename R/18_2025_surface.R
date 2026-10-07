# ============================================================================
# 18_2025_surface.R
# How much does the estimate move with a single year's conditions? The model
# is fitted on year-matched annual covariates and predicted on the 2000-2024
# long-term means (07). This predicts the same model onto each single year's
# rasters, every occurrence year plus PROJ_YEAR, with 2025 population
# throughout, so only the covariates change. (Not for the common-scale
# step: same model, so a change in level is the model's response.)
#
# MODIS Terra, the source of night LST (MOD11A2), has drifted to an earlier
# overpass since 2020 (TERRA_DRIFT_FROM, params.R). An earlier night pass
# reads warmer, so night LST from those years is not comparable with
# earlier years.
#
#
# Inputs:  MODEL_FILE, TUNING_FILE, TRAIN_FILE, SUIT_FILE, POP_ALIGNED_FILE,
#          DOMAIN_FILE, ADM1_FILE, OCC_RAW_FILE (years), retained_vars.rds,
#          COV_FILES, COV_ANNUAL, arp_summary.csv (08)
# Outputs: outputs/tables/single_year_estimates.csv, single_year_by_state.csv
#          outputs/figures/fig_single_year.png / .pdf
# Runtime: about 10-15 min (one domain prediction per year).
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(ggplot2); library(patchwork); library(maxnet)
})

# ------------------------------ Load inputs ---------------------------------

mod    <- readRDS(MODEL_FILE)
tuning <- readRDS(TUNING_FILE)
train  <- readRDS(TRAIN_FILE)
vars   <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
suit_r <- rast(SUIT_FILE)
pop    <- rast(POP_ALIGNED_FILE)
arp_08 <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))

dyn   <- intersect(vars, names(COV_ANNUAL))
years <- sort(unique(c(read.csv(OCC_RAW_FILE)$year, PROJ_YEAR)))

need    <- unlist(lapply(years, function(y) year_files(y)[dyn]))
missing <- need[!file.exists(file.path(DIR_COVARIATES, need))]
stopifnot(
  "Population grid differs from the surface" = compareGeom(pop, suit_r, stopOnError = FALSE),
  "Night LST is not a retained covariate"    = "lst_night" %in% vars,
  "A retained dynamic covariate was not exported for PROJ_YEAR" = all(dyn %in% names(PROJ_ANNUAL))
)
if (length(missing) > 0) stop("Missing annual rasters: ", paste(missing, collapse = ", "))
cat("Model:", tuning$fc, "rm", tuning$rm, "| dynamic covariates:", paste(dyn, collapse = ", "),
    "\nYears:", years, "| drift-affected from", TERRA_DRIFT_FROM, "\n")

# ----------- The prediction code reproduces 07 (surface) and 08 -------------

cat("\nLong-term means: predicting...\n")
ltm_r <- terra::predict(domain_covs(vars), mod, type = "cloglog", na.rm = TRUE)
thr   <- thresholds_from(maxnet_prob(mod, train$occ_env[, vars]),
                         maxnet_prob(mod, train$bg_env[, vars]))
arp_lt  <- arp_estimates(suit_r, pop, thr)
thr_08  <- setNames(arp_08$threshold, arp_08$metric)
arp_08v <- setNames(arp_08$arp,       arp_08$metric)
stopifnot(
  "Surface differs from SUIT_FILE (07)" =
    global(abs(ltm_r - suit_r), "max", na.rm = TRUE)[[1]] < 1e-6 &&
    global(is.na(ltm_r) != is.na(suit_r), "sum")[[1]] == 0,
  "Thresholds differ from arp_summary.csv (08)" = all(abs(thr - thr_08[names(thr)]) < 1e-6),
  "Estimates differ from arp_summary.csv (08)"  = all(round(arp_lt) == arp_08v[names(arp_lt)])
)
cat("Prediction code reproduces 07 (surface) and 08 (thresholds, estimates)\n")

# ------------------------- Cells every year covers ---------------------------
# A single year's composite can lack a cell the long-term mean has (no
# clear-sky night all year). Comparisons use cells every surface covers.

common <- !is.na(suit_r)
for (y in years) common <- common & noNA(domain_covs(vars, year_files(y)))
pop_all  <- global(mask(pop, suit_r), "sum", na.rm = TRUE)[[1]]
pop_comm <- global(mask(pop, common, maskvalues = 0), "sum", na.rm = TRUE)[[1]]
n_lost   <- global(!is.na(suit_r) & !common, "sum", na.rm = TRUE)[[1]]
cat("Cells without a value in at least one year:", fmt(n_lost),
    "| population there:", fmt(pop_all - pop_comm), "\n")
stopifnot("Cells missing in some year hold more than POP_TOL of the population" =
  (pop_all - pop_comm) / pop_all <= POP_TOL)

# -------------------------------- Each year ---------------------------------

zones   <- state_zones(suit_r)
lst_max <- max(train$occ_env$lst_night)   # warmest presence night (year-matched), as 08

summ_year <- function(covs_c, pred_c, year) {
  a   <- arp_estimates(pred_c, pop, thr)
  hot <- global(pop * pred_c * (covs_c[["lst_night"]] > lst_max), "sum", na.rm = TRUE)[[1]]
  s   <- zonal(pop * pred_c, zones, fun = "sum", na.rm = TRUE)
  names(s) <- c("state", if (is.na(year)) "ltm" else as.character(year))
  list(row = data.frame(
         year = year,
         lst_night = global(covs_c[["lst_night"]], "mean", na.rm = TRUE)[[1]],
         rainfall  = if ("rainfall" %in% vars)
                       global(covs_c[["rainfall"]], "mean", na.rm = TRUE)[[1]] else NA_real_,
         risk_weighted = a[["risk_weighted"]], maxsss = a[["maxsss"]], p10 = a[["p10"]],
         hot_pct = 100 * hot / a[["risk_weighted"]]),
       state = s)
}

res <- list(summ_year(mask(domain_covs(vars), common, maskvalues = 0),
                      mask(suit_r, common, maskvalues = 0), NA_integer_))
for (y in years) {
  t0     <- Sys.time()
  covs_y <- mask(domain_covs(vars, year_files(y)), common, maskvalues = 0)
  pred_y <- terra::predict(covs_y, mod, type = "cloglog", na.rm = TRUE)
  res[[length(res) + 1]] <- summ_year(covs_y, pred_y, y)
  cat(sprintf("  %d | %.1f min\n", y, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

est <- bind_rows(lapply(res, `[[`, "row")) |>
  mutate(group = case_when(is.na(year)              ~ "long-term mean",
                           year < TERRA_DRIFT_FROM  ~ "before drift",
                           TRUE                     ~ "drift"))
lt  <- filter(est, group == "long-term mean")
pre <- filter(est, group == "before drift")
est <- est |>
  mutate(lst_anom      = lst_night - lt$lst_night,
         rw_pct        = 100 * (risk_weighted / lt$risk_weighted - 1),
         above_pre_lst = group == "drift" & lst_night > max(pre$lst_night))

cat("\nBy year (cells every year covers; 2025 population throughout):\n")
est |>
  mutate(year = ifelse(is.na(year), "long-term", year),
         across(c(risk_weighted, maxsss, p10), fmt),
         across(c(lst_night, lst_anom, rw_pct, hot_pct), ~ round(., 1)),
         rainfall = round(rainfall)) |>
  select(year, group, lst_night, lst_anom, rainfall, risk_weighted, rw_pct,
         maxsss, p10, hot_pct, above_pre_lst) |>
  print(row.names = FALSE)

cat(sprintf(paste0("\nBefore drift (%d years): risk-weighted %s to %s (%+.1f%% to %+.1f%% of the ",
                   "long-term estimate %s); median %s; mean of the annual estimates %s\n"),
            nrow(pre), fmt(min(pre$risk_weighted)), fmt(max(pre$risk_weighted)),
            100 * (min(pre$risk_weighted) / lt$risk_weighted - 1),
            100 * (max(pre$risk_weighted) / lt$risk_weighted - 1),
            fmt(lt$risk_weighted), fmt(median(pre$risk_weighted)), fmt(mean(pre$risk_weighted))))
conf <- est$year[est$above_pre_lst]
cat("Drift years with mean night LST above every earlier year (confounded, not read as climate):",
    if (length(conf)) paste(conf, collapse = ", ") else "none", "\n")
cat("Rank of", PROJ_YEAR, "among", length(years), "years (1 = highest risk-weighted):",
    rank(-est$risk_weighted[!is.na(est$year)])[est$year[!is.na(est$year)] == PROJ_YEAR], "\n")

# ----------------------- Presences from drift years --------------------------

occ      <- train$occ_clean
ltm_lst  <- terra::extract(domain_covs("lst_night"),
                           as.matrix(occ[, c("longitude", "latitude")]))[, 1]
gap      <- train$occ_env$lst_night - ltm_lst
in_drift <- occ$year >= TERRA_DRIFT_FROM
med_gap  <- c(before = median(gap[!in_drift]),
              drift  = if (any(in_drift)) median(gap[in_drift]) else NA_real_)
cat(sprintf(paste0("\nPresences from drift years: %d (%s) | median year-matched minus ",
                   "long-term night LST: drift years %.2f C, earlier %.2f C\n"),
            sum(in_drift),
            paste(names(table(occ$year[in_drift])), table(occ$year[in_drift]),
                  sep = ": ", collapse = ", "),
            med_gap[["drift"]], med_gap[["before"]]))
cat("Raise for 06:", any(in_drift) && med_gap[["drift"]] - med_gap[["before"]] > 0.5, "\n")

# -------------------------------- By state ----------------------------------

st <- Reduce(function(a, b) merge(a, b, by = "state"), lapply(res, `[[`, "state"))
pre_cols <- as.character(pre$year)
dr_cols  <- as.character(est$year[est$group == "drift"])
st_tab <- data.frame(state = st$state, ltm = st$ltm,
                     pre_min_pct = 100 * (apply(st[pre_cols], 1, min) / st$ltm - 1),
                     pre_max_pct = 100 * (apply(st[pre_cols], 1, max) / st$ltm - 1),
                     setNames(100 * (st[dr_cols] / st$ltm - 1), paste0("pct_", dr_cols)),
                     check.names = FALSE) |>
  arrange(desc(ltm))
cat("\nBy state: risk-weighted, long-term surface; before-drift range and drift years (% change):\n")
st_tab |> mutate(ltm = fmt(ltm), across(-c(state, ltm), ~ round(., 1))) |>
  print(row.names = FALSE)

# --------------------------------- Figure -----------------------------------

yr   <- filter(est, group != "long-term mean")
yr_p <- function(var, ref, ylab, legend = FALSE) {
  ggplot(yr, aes(year, .data[[var]], colour = group)) +
    geom_hline(yintercept = ref, linetype = "dashed", colour = "grey50") +
    geom_point(size = 2) +
    scale_colour_manual(values = pal_drift, labels = drift_labels, name = NULL,
                        guide = if (legend) "legend" else "none") +
    labs(x = NULL, y = ylab)
}
p_rw <- yr_p("risk_weighted", lt$risk_weighted, "Population at risk,\nrisk-weighted", legend = TRUE) +
  scale_y_continuous(labels = function(x) paste0(x / 1e6, "M")) +
  annotate("text", x = min(yr$year), y = lt$risk_weighted, label = "Long-term mean surface",
           hjust = 0, vjust = -0.6, size = 2.8, colour = "grey40") +
  theme(legend.position = "top")
p_lst  <- yr_p("lst_night", lt$lst_night, "Mean night LST,\nSudan (\u00b0C)")
p_rain <- yr_p("rainfall",  lt$rainfall,  "Mean rainfall,\nSudan (mm)")
for (ext in c("png", "pdf"))
  save_fig(file.path(DIR_FIGS, paste0("fig_single_year.", ext)),
           (p_rw / p_lst / p_rain) + plot_annotation(tag_levels = "a"),
           width = FIG_WIDTH_FULL, height = 16)

# --------------------------------- Save -------------------------------------

write.csv(est,    file.path(DIR_TABLES, "single_year_estimates.csv"), row.names = FALSE)
write.csv(st_tab, file.path(DIR_TABLES, "single_year_by_state.csv"),  row.names = FALSE)
cat("\n18_2025_surface.R complete\n")