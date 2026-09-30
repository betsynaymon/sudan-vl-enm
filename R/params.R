# params.R - source all modelling decisions
# Sourced at top of each script

library("here")

# ----- Reproducibility -----
SEED <- 1234L        # Used at every stochastic step

# ----- Paths -----
DIR_RAW        <- here("data", "raw")        # source data; pipeline writes only downloads here
DIR_COVARIATES <- DIR_RAW                    # covariate rasters live in raw/
DIR_POP        <- file.path(DIR_RAW, "population")
DIR_PROCESSED  <- here("data", "processed")
DIR_OUTPUTS    <- here("outputs")
DIR_FIGS       <- file.path(DIR_OUTPUTS, "figures")
DIR_TABLES     <- file.path(DIR_OUTPUTS, "tables")
DIR_MODELS     <- file.path(DIR_OUTPUTS, "models")
DIR_SURFACES   <- file.path(DIR_OUTPUTS, "surfaces")
DIR_SENS       <- file.path(DIR_OUTPUTS, "sensitivity")   # script 21
for (d in c(DIR_POP, DIR_PROCESSED, DIR_OUTPUTS, DIR_FIGS, DIR_TABLES,
            DIR_MODELS, DIR_SURFACES, DIR_SENS))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ----- Files -----
OCC_RAW_FILE <- file.path(DIR_RAW, "compiled_vl_presences.csv")
OCC_FILE     <- file.path(DIR_PROCESSED, "occurrences_thinned.csv")   # written by 02
BG_FILE      <- file.path(DIR_PROCESSED, "background_points.csv")     # written by 04
TT_FILE      <- file.path(DIR_RAW, "weiss_travel_time.tif")           # downloaded by 11
POP_FILE     <- file.path(DIR_POP, "sdn_pop_2025_100m_constrained.tif") # downloaded by 08


# ----- Covariate filename lookup -----
# All candidates screened in 03 (static layers and 2000-2024 long-term means).
# 03 writes the retained subset to retained_vars.rds; downstream scripts index
# COV_FILES by it and never use the full set.
COV_FILES <- c(
  elevation  = "elevation_1km.tif",
  slope      = "slope_1km.tif",
  river_dist = "river_distance_1km.tif",
  vertisols  = "vertisols_1km.tif",
  lst_day    = "lst_day_annual_mean_2000_2024_1km.tif",
  lst_night  = "lst_night_annual_mean_2000_2024_1km.tif",
  ndvi       = "ndvi_annual_mean_2000_2024_1km.tif",
  rainfall   = "rainfall_mean_2000_2024_1km.tif",
  treecover  = "treecover_mean_2000_2024_1km.tif"
)

# Dynamic covariates: one raster per year with occurrence records, matched to
# each record's year (exported in python/02_covariates.ipynb). COV_FILES holds
# their long-term means (2000-2024).
COV_ANNUAL <- c(
  lst_night = "lst_night_annual_{year}_1km.tif",
  rainfall  = "rainfall_{year}_1km.tif"
)

# ----- Study domain -----
# Calibration extent = prediction extent: Sudan (GADM level 0) rasterised to
# the covariate grid. Coded 1 inside, 0 outside.
DOMAIN_FILE <- file.path(DIR_PROCESSED, "domain_sudan.tif")

# Supplement only (21): background restricted to >= 150 mm mean annual
# rainfall. Threshold set in python/01_study_area.ipynb.
SENS_MASK_FILE <- file.path(DIR_RAW, "ecological_mask_150mm.tif")

# ----- Study area -----
# Sudan as administered: GADM minus the Egyptian-administered Halaib Triangle,
# taken as land north of the 22nd parallel (the 1899 boundary) east of the
# western vertex of GADM's Halayeb locality there. Matches WorldPop's
# population surface. Built once in 01; other scripts read these files
# instead of calling gadm().
EXCLUDE_ADM2     <- "Halayeb"
EXCLUDE_NORTH_OF <- 22
ADM0_FILE <- file.path(DIR_PROCESSED, "study_adm0.gpkg")
ADM1_FILE <- file.path(DIR_PROCESSED, "study_adm1.gpkg")


# ----- Occurrence filtering -----
# Referral record: a major national referral centre draws patients from the
# whole country, so its location carries no information on where infection
# occurred (02).
#   78: Elnoor et al. 2024 facility report, Tropical Disease Teaching Hospital, Omdurman
REFERRAL_IDS <- c(78)

# Khartoum case records (Pigott et al. 2014, 2001-2003). Retained in the
# primary model: Pigott included only autochthonous cases. Dropped in the
# sensitivity refit (10).
KHARTOUM_CASE_IDS <- c(4, 5, 8, 15)

# Precision check (10): map-digitised sources with positional error > 4 km
# (python/03_digitized_points_check.ipynb)
IMPRECISE_SOURCES <- c(
  "vandebogaart_et_al_2013_fig1", "vandebogaart_et_al_2013_fig2",
  "vandebogaart_et_al_2013_fig3", "vandebogaart_et_al_2013_fig4",
  "vandebogaart_et_al_2013_fig5",
  "hassan_et_al_2020_fig1", "hassan_et_al_2020_fig2",
  "hassan_et_al_2020_fig3", "hassan_et_al_2020_fig4"
)

# Data-quality check (10): records locating a facility rather than a site
FACILITY_TYPE <- "facility_report"

# ----- Spatial thinning -----
THIN_KM <- 5     # see 02_spatial_thinning.R
THIN_REPS    <- 100                      # spThin replicates
THIN_TEST_KM <- c(1, 2, 5, 10, 20, 50)   # retention-curve diagnostic

# ----- Background ----- 
N_BACKGROUND <- 10000     # uniform-random over study domain

# ----- Collinearity Screening -----
# See 03_collinearity.R
COLLIN_SAMPLE_N <- 100000     # random study domain cells for screening
COR_THRESHOLD <- 0.7     # |r| pairwise cut off 
VIF_THRESHOLD <- 10
# RETAINED_VARS is written to disk by 03_collinearity.R and read downstream

# ----- Spatial block CV -----
# Blocks near the covariate autocorrelation range leave most presences in one
# or two folds (Gedaref clustering); 05 prints balance at BLOCK_TEST_KM.
BLOCK_SIZE_M  <- 100000
K_FOLDS       <- 4
BLOCK_TEST_KM <- c(100, 200, 350, 700)
SAC_SAMPLE_N  <- 5000   # cells sampled for cv_spatial_autocor
FOLD_ITER     <- 1000   # cv_spatial keeps the most balanced of these assignments
MIN_TEST_PRES <- 10     # floor on test presences per fold for a usable fold CBI
FOLDS_FILE      <- file.path(DIR_MODELS, "spatial_cv_folds.rds")  # blockCV object (blocks)
FOLD_TABLE_FILE <- file.path(DIR_MODELS, "cv_fold_table.rds")     # pa, id, fold
MODEL_FILE  <- file.path(DIR_MODELS, "maxent_final.rds")
TRAIN_FILE  <- file.path(DIR_MODELS, "training_data.rds")
TUNING_FILE <- file.path(DIR_MODELS, "selected_tuning.rds")


# ----- ENMeval -----
ENM_FC    <- c("L", "LQ", "LQH", "LQHP")   # H-only removed: predict.maxnet fails on binary covariates
ENM_RM <- c(seq(0.5, 4, by = 0.5), 5, 6, 8)   # extends well past the conventional range
BOYCE_RES <- 100                            # ecospat.boyce resolution (06, 15)

# Selection: highest mean validation CBI among configurations with
# ecologically plausible response curves (checked in 07).
# NULL = highest-CBI configuration. To override, set both and record why.
# 06 writes the fitted choice to selected_tuning.rds; downstream reads that.
SELECTED_FC <- NULL
SELECTED_RM <- NULL

# Underfitting check (06): configurations are also scored within the >= 150 mm
# region. If a configuration beats the selected one on within-belt CBI by more
# than one SE, selection moves to within-belt CBI, because the population
# estimate depends on discrimination within the inhabited belt. Both metrics
# are reported.

# ----- Thresholds and importance -----
OMISSION_Q <- 0.10   # p10: 10th percentile of training-presence predictions
N_PERM     <- 50     # permutation-importance repeats (07)

# ----- Accessibility -----
TT_REMOTE_MIN <- 300   # remoteness threshold, minutes to nearest city (22)

# ----- Comparators (supplement) -----
RF_NTREES     <- 1000                    # 09; mtry grid derived from retained vars
GBT_LR        <- c(0.001, 0.005, 0.01)   # 17
GBT_DEPTH     <- c(1, 3, 5)
GBT_NTREES    <- 5000
GBT_TREE_STEP <- 500
GBT_BAG_FRAC  <- 0.75
GBT_MIN_OBS   <- 10


# ----- Null Model Test -----
N_NULL <- 99        # randomizations 
