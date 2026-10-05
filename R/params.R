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
TT_FILE      <- file.path(DIR_RAW, "weiss_travel_time.tif") # downloaded by 11
TT_URL <- paste0(
  "https://data.malariaatlas.org/geoserver/Accessibility/ows?",
  "service=WCS&version=2.0.1&request=GetCoverage&format=image/geotiff&",
  "CoverageId=Accessibility__201501_Global_Travel_Time_to_Cities&",
  "subset=Long(21.5,39)&subset=Lat(8.5,22.5)")
POP_FILE     <- file.path(DIR_POP, "sdn_pop_2025_100m_constrained.tif") # downloaded by 08
POP_URL <- paste0("https://data.worldpop.org/GIS/Population/Global_2015_2030/",
                  "R2025A/2025/SDN/v1/100m/constrained/sdn_pop_2025_CN_100m_R2025A_v1.tif")
POP_ALIGNED_FILE <- file.path(DIR_SURFACES, "worldpop_2025_aligned.tif")
CANDIDATES_FILE  <- file.path(DIR_SURFACES, "plateau_candidate_surfaces.tif")
POP_TOL <- 0.001   # max relative change allowed when aligning or summing population
BIAS_SURFACES_FILE <- file.path(DIR_SURFACES, "bias_correction_surfaces.tif")  # 12; read by 20
VARIANT_SURFACES_FILE <- file.path(DIR_SURFACES, "variant_surfaces.tif")  # 15; read by 20
STATE_STATUS_FILE <- file.path(DIR_RAW, "state_endemic_status.csv")  # literature endemic status by state (16)
DQ_SURFACES_FILE     <- file.path(DIR_SURFACES, "data_quality_refit_surfaces.tif")     # 10; read by 20
TT_COV_SURFACES_FILE <- file.path(DIR_SURFACES, "accessibility_covariate_surfaces.tif") # 24; read by 20

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
# their long-term means (2000-2024). Elevation, slope, river distance and
# vertisols are static.
COV_ANNUAL <- c(
  lst_day   = "lst_day_annual_{year}_1km.tif",
  lst_night = "lst_night_annual_{year}_1km.tif",
  ndvi      = "ndvi_annual_{year}_1km.tif",
  rainfall  = "rainfall_{year}_1km.tif",
  treecover = "treecover_{year}_1km.tif"
)

# Seasonal LST-night composites, tested in 15 only. Kept out of COV_FILES so
# 01's checks and 03's screening are unchanged.
SEASONAL_FILES <- c(
  lst_night_dry = "lst_night_dry_mean_2000_2024_1km.tif",
  lst_night_wet = "lst_night_wet_mean_2000_2024_1km.tif"
)
SEASONAL_ANNUAL <- c(
  lst_night_dry = "lst_night_dry_{year}_1km.tif",
  lst_night_wet = "lst_night_wet_{year}_1km.tif"
)

# ----- Study domain -----
# Calibration extent = prediction extent: Sudan (GADM level 0) rasterised to
# the covariate grid. Coded 1 inside, 0 outside.
DOMAIN_FILE <- file.path(DIR_PROCESSED, "domain_sudan.tif")

# Supplement only (21): background restricted to >= 150 mm mean annual
# rainfall. Threshold set in python/01_study_area.ipynb.
SENS_MASK_FILE <- file.path(DIR_RAW, "ecological_mask_150mm.tif")

# ----- Study area -----
# Modelled domain: Sudan (GADM) excluding the Halaib Triangle, for which
# WorldPop assigns no population and covariate coverage is incomplete. The
# Triangle is taken as land north of the 22nd parallel east of the western
# vertex of GADM's Halayeb locality there. Maps show the full GADM boundary,
# with the Triangle marked as not modelled.
EXCLUDE_ADM2      <- "Halayeb"
EXCLUDE_NORTH_OF  <- 22
ADM0_FILE         <- file.path(DIR_PROCESSED, "study_adm0.gpkg")    # modelled domain
ADM1_FILE         <- file.path(DIR_PROCESSED, "study_adm1.gpkg")
DISPLAY_ADM0_FILE <- file.path(DIR_PROCESSED, "display_adm0.gpkg")  # full GADM outline (maps)
DISPLAY_ADM1_FILE <- file.path(DIR_PROCESSED, "display_adm1.gpkg")  # full GADM states (maps that aren't model output)
EXCLUDED_FILE     <- file.path(DIR_PROCESSED, "excluded_area.gpkg") # Halaib Triangle (maps)


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

# Wad Madani case records (Pigott et al. 2014; one location, 2001 and 2002):
# the same concern as Khartoum. Dropped with them in a widened refit (10).
CENTRAL_CITY_CASE_IDS <- c(6, 13)

# Hassan et al. 2020 fig 4: vector-positive trap site near Khartoum (2013),
# the only dry-riverine presence. In the precision set; also refitted alone
# (10) to separate its effect from positional error.
KHARTOUM_VECTOR_ID <- 121

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

SUIT_FILE      <- file.path(DIR_SURFACES, "maxent_suitability.tif")
MESS_FILE      <- file.path(DIR_SURFACES, "mess_vs_presences.tif")
MESS_VARS_FILE <- file.path(DIR_SURFACES, "mess_vs_presences_by_variable.tif")

# ----- ENMeval -----
ENM_FC    <- c("L", "LQ", "LQH", "LQHP")   # H-only removed: predict.maxnet fails on binary covariates
ENM_RM <- c(seq(0.5, 4, by = 0.5), 5, 6, 8)   # extends well past the conventional range
BOYCE_RES <- 100                            # ecospat.boyce resolution (06, 15)

# Plateau choice (06b): LQHP rm 4. It meets the plausibility criteria within
# the presence range, has the higher within-belt CBI, and is the only
# candidate whose rainfall response peaks and declines within the range of
# the presence data, consistent with the upper rainfall limit reported for
# P. orientalis (Gebre-Michael et al. 2004). Its steeper decline beyond the
# wettest presences is shaped by clamping and is treated as extrapolation.
SELECTED_FC <- "LQHP"
SELECTED_RM <- 4

# 06b also refits these feature classes at the plateau's lowest rm, as a
# reference without product features.
PLATEAU_REF_FC <- "LQH"


# Underfitting check (06): configurations are also scored within the >= 150 mm
# region. If a configuration beats the selected one on within-belt CBI by more
# than one SE, selection moves to within-belt CBI, because the population
# estimate depends on discrimination within the inhabited belt. Both metrics
# are reported.

# Choice within the plateau (07): LQHP rm 4 and rm 8 bracket it and differ in
# within-belt CBI. If both have ecologically plausible response curves, prefer
# the higher within-belt CBI, because the population estimate depends on
# discrimination within the inhabited belt; otherwise prefer the plausible one.
# The other is reported as a sensitivity estimate in 08. (Within-belt CBI was
# introduced as a secondary criterion during analysis, in response to the
# underfitting concern, and is reported as such.)

# ----- Thresholds and importance -----
OMISSION_Q <- 0.10   # p10: 10th percentile of training-presence predictions
N_PERM     <- 50     # permutation-importance repeats (07)

# ----- Accessibility -----
TT_REMOTE_MIN <- 300   # remoteness threshold, minutes to nearest city (22)

# ----- Comparators (09, 17) -----
# Random forest (09): a down-sampled probability forest (fit_rf, helpers.R).
# Every tree draws, with replacement, as many presences and as many background
# points as there are presences (Valavi et al. 2021, Ecography: down-sampling
# improved random forests on presence-background data; class weighting did
# not). mtry is tuned over 1 to the number of retained covariates.
RF_NTREES       <- 1000
RF_MODEL_FILE   <- file.path(DIR_MODELS,   "rf_final.rds")            # read by 17
RF_SUIT_FILE    <- file.path(DIR_SURFACES, "rf_suitability.tif")      # read by 17
COMPARATOR_FILE <- file.path(DIR_TABLES,   "comparator_summary.csv")  # MaxEnt and RF (09); read by 17

# Gradient boosted trees (17; fit_gbt, helpers.R): presences weight 1 and the
# background down-weighted to the same total weight, as for boosted
# regression trees in Valavi et al. 2022.
# Tuned by spatial CV over learning rate x depth x number of trees; the grid
# includes their settings (learning rate 0.001, depth 1 or 5, bag fraction
# 0.75) and the dissertation's, with the number of trees extended beyond the
# dissertation's 5,000. 17 stops if the selection sits at the maximum.
# Thresholds come from cross-fitted predictions, the boosting counterpart of
# the forest's out-of-bag predictions.
GBT_LR           <- c(0.001, 0.005, 0.01)  # shrinkage
GBT_DEPTH        <- c(1, 3, 5)             # gbm interaction.depth: splits per tree
GBT_TREE_STEP    <- 500
GBT_NTREES_MAX   <- 10000
GBT_BAG          <- 0.75                   # share of rows drawn for each tree
GBT_MIN_NODE     <- 10                     # minimum rows per leaf
GBT_THR_FOLDS    <- 10                     # random folds for cross-fitted thresholds
GBT_MODEL_FILE   <- file.path(DIR_MODELS,   "gbt_final.rds")
GBT_SUIT_FILE    <- file.path(DIR_SURFACES, "gbt_suitability.tif")
THREE_MODEL_FILE <- file.path(DIR_TABLES,   "three_model_comparison.csv")  # one row per algorithm (17)

# ----- Fold-excluded refits (14) -----
FOLD_SURFACES_FILE <- file.path(DIR_SURFACES, "fold_excluded_surfaces.tif")  # layers without_fold_k
FOLD_SD_FILE       <- file.path(DIR_SURFACES, "fold_excluded_sd.tif")        # SD across those layers

# ----- Null model test (13) -----
N_NULL      <- 99      # iterations per null; smallest attainable p = 1 / (N_NULL + 1)
N_NULL_POOL <- 50000   # candidate null presences per pool
NULL_ALPHA  <- 0.05    # the accessibility null is beaten if p <= NULL_ALPHA

# ----- Accessibility covariate test (24) -----
N_TT_GRID <- 20   # presence travel-time quantiles averaged over to integrate accessibility out

# ----- Single-year projections (18) -----
# The model is fitted on year-matched annual values and predicted on the
# 2000-2024 long-term means; 18 also predicts it onto each single year's
# rasters: every occurrence year (python/02_covariates.ipynb) and PROJ_YEAR
# (python/04_2025_rasters.ipynb, which exported only these layers).
PROJ_YEAR   <- 2025
PROJ_ANNUAL <- COV_ANNUAL[c("lst_night", "rainfall")]
# MODIS Terra (MOD11A2, night LST) drifts to an earlier overpass: beyond
# 10:30 +/- 1 min from April 2021, about 9:00 by December 2025 (NASA Terra
# mission pages). An earlier night pass reads warmer, so night LST from these
# years is not comparable with earlier years.
TERRA_DRIFT_FROM <- 2021