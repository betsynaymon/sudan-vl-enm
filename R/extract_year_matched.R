# ============================================================================
# extract_year_matched.R
# Year-matched covariate extraction, shared across scripts that build new
# point sets (12, 13). Static covariates come from single rasters; dynamic
# covariates from the annual raster for each point's year. Rasters are masked
# to the ecological mask, as in the primary model (06).
# ============================================================================

static_files <- c(
  slope      = "slope_1km.tif",
  river_dist = "river_distance_1km.tif",
  vertisols  = "vertisols_1km.tif"
)

dynamic_patterns <- c(
  lst_night = "lst_night_annual_{year}_1km.tif",
  rainfall  = "rainfall_{year}_1km.tif"
)

extract_year_matched <- function(pts, vars, mask_r) {
  static_vars  <- intersect(vars, names(static_files))
  dynamic_vars <- intersect(vars, names(dynamic_patterns))

  static_r <- rast(file.path(DIR_COVARIATES, static_files[static_vars]))
  names(static_r) <- static_vars
  static_r <- mask(static_r, mask_r, maskvalues = 0)

  static_vals <- terra::extract(static_r,
    as.matrix(pts[, c("longitude", "latitude")]))

  dynamic_vals <- as.data.frame(matrix(NA, nrow = nrow(pts),
                                       ncol = length(dynamic_vars),
                                       dimnames = list(NULL, dynamic_vars)))

  for (yr in sort(unique(pts$year))) {
    idx <- which(pts$year == yr)
    yr_files <- sapply(dynamic_vars, function(v) {
      file.path(DIR_COVARIATES, gsub("\\{year\\}", yr, dynamic_patterns[v]))
    })
    yr_r <- rast(yr_files)
    names(yr_r) <- dynamic_vars
    yr_r <- mask(yr_r, mask_r, maskvalues = 0)

    dynamic_vals[idx, ] <- terra::extract(yr_r,
      as.matrix(pts[idx, c("longitude", "latitude")]))
  }

  cbind(static_vals, dynamic_vals)
}