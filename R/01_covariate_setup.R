# ============================================================================
# 01_covariate_setup.R
# Downloads covariate rasters from Google Drive and verifies spatial alignment.
#
# Inputs:  Google Drive folder "sudan_enm_covariates" (GEE exports)
# Outputs: Rasters in data/raw/ (downloaded, not modified)
#
# Note: This is a one-time data acquisition step. Once all rasters are
# downloaded, subsequent pipeline runs can skip this script.
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(tidyverse)
})

# ---- Download rasters (skipped if already present) -------------------------

if (!dir.exists(DIR_COVARIATES)) {
  dir.create(DIR_COVARIATES, recursive = TRUE)
}

n_local <- length(list.files(DIR_COVARIATES, pattern = "\\.tif$"))

if (n_local >= 246) {
  cat("All", n_local, "rasters already present in data/raw/ — skipping Drive download\n")
} else {
  cat(n_local, "rasters found locally — connecting to Drive for remaining files\n")

  library(googledrive)
  drive_auth()

  covariate_files <- drive_ls("sudan_enm_covariates")
  cat("Found", nrow(covariate_files), "files in Drive folder\n")

  stopifnot(
    "No files found in Drive folder — check folder name or auth" =
      nrow(covariate_files) > 0
  )

  walk2(covariate_files$id, covariate_files$name, function(file_id, file_name) {
    dest <- file.path(DIR_COVARIATES, file_name)
    if (file.exists(dest)) {
      return(invisible(NULL))
    }
    cat("  Downloading:", file_name, "\n")
    drive_download(
      file   = as_id(file_id),
      path   = dest,
      overwrite = FALSE
    )
  })

  n_local <- length(list.files(DIR_COVARIATES, pattern = "\\.tif$"))
  cat("Files now in data/raw:", n_local, "\n")
}

# ---- Verify spatial alignment ----------------------------------------------
# Checks GEE-exported covariates only

tif_files <- list.files(DIR_COVARIATES, pattern = "\\.tif$", full.names = TRUE)
tif_files <- tif_files[!grepl("weiss_travel_time", tif_files)]

raster_info <- map_dfr(tif_files, function(f) {
  r <- rast(f)
  tibble(
    file = basename(f),
    crs  = crs(r, describe = TRUE)$code,
    xres = res(r)[1],
    yres = res(r)[2],
    nrow = nrow(r),
    ncol = ncol(r),
    xmin = ext(r)$xmin,
    xmax = ext(r)$xmax,
    ymin = ext(r)$ymin,
    ymax = ext(r)$ymax
  )
})

cat("Checked alignment of", nrow(raster_info), "GEE-exported rasters\n")

# Flag any misalignment
n_crs  <- n_distinct(raster_info$crs)
n_res  <- n_distinct(paste(raster_info$xres, raster_info$yres))
n_ext  <- n_distinct(paste(raster_info$xmin, raster_info$xmax,
                           raster_info$ymin, raster_info$ymax))

if (n_crs == 1 & n_res == 1 & n_ext == 1) {
  cat("All rasters aligned: single CRS, resolution, and extent\n")
} else {
  cat("WARNING: Alignment issues detected\n")
  cat("  Unique CRS:", n_crs, "\n")
  cat("  Unique resolutions:", n_res, "\n")
  cat("  Unique extents:", n_ext, "\n")
  # Print only the distinct configurations, not all 249 rows
  raster_info |>
    distinct(crs, xres, yres, nrow, ncol, xmin, xmax, ymin, ymax) |>
    print()
}

cat("01_covariate_setup.R complete\n")