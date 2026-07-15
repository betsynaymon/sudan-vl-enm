# params.R - source all modelling decisions
# Sourced at top of each script

library("here")

# ----- Reproducibility -----
SEED <- 1234L        # Used at every stochastic step

# ----- Paths -----
# Covariate rasters: exported and grid verified in 01_covariate_setup.R
DIR_COVARIATES <- here("data", "raw")
DIR_OUTPUTS  <- here("outputs")
DIR_FIGS     <- here("outputs", "figures")
DIR_TABLES   <- here("outputs", "tables")
DIR_MODELS   <- here("outputs", "models")
DIR_SURFACES <- here("outputs", "surfaces")
for (d in c(DIR_OUTPUTS, DIR_FIGS, DIR_TABLES, DIR_MODELS, DIR_SURFACES))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ----- Covariate filename lookup -----
# Maps short variable names (used in model formulas) to filenames on disk
COV_FILES <- c(
  slope      = "slope_1km.tif",
  river_dist = "river_distance_1km.tif",
  vertisols  = "vertisols_1km.tif",
  lst_night  = "lst_night_annual_mean_2000_2024_1km.tif",
  rainfall   = "rainfall_mean_2000_2024_1km.tif"
)

# ----- Study area/mask -----
MASK_RAINFALL_MM <- 150     # >= 150mm CHIRPS; conservative outer bound

# ----- Spatial thinning -----
THIN_KM <- 5     # see 02_spatial_thinning.R

# ----- Background ----- 
N_BACKGROUND <- 10000     # uniform-random over masked study area

# ----- Collinearity Screening -----
# See 03_collinearity.R
COLLIN_SAMPLE_N <- 100000     # random study-area cells for screening
COR_THRESHOLD <- 0.7     # |r| pairwise cut off 
VIF_THRESHOLD <- 10
# RETAINED_VARS is written to disk by 03_collinearity.R and read downstream

# ----- Spatial block CV -----
# BLOCK_SIZE_M is estimated in CV step
# Median autocorrelation range = 706km, but blocks at that scale produce imbalanced folds
# due to Gedaref clustering. Reduced to 100km with k=4 for viable per-fold balance.
# See 05_spatial_cv_folds.R
BLOCK_SIZE_M <- 100000
K_FOLDS <- 4

# ----- ENMeval -----
ENM_FC <- c("L", "LQ", "LQH", "LQHP")   # H-only removed: predict.maxnet fails on binary covariates
ENM_RM <- seq(0.5, 4, by = 0.5)           # regularization multipliers
# Model selection: best CBI across the feature class × regularization grid.
# Omission rate (0.20–0.35) does not discriminate under spatial CV with
# clustered presences, so it is not used as a filter. See 06_maxent_tuning.R.
# SELECTED_FC / SELECTED_RM saved to selected_tuning.rds after grid search.

# ----- Evaluation Metrics -----
# Primary: Continuous Boyce Index 
# Reported: AUC (for comparability with prior studies)
# Checked: or.10p (non-discriminating under spatial CV; documented)

# ----- Null Model Test -----
N_NULL <- 99        # randomizations 
N_PRES <- 98        # match observed presence count
N_DRAW <- N_PRES * 2