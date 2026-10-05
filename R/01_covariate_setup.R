# ============================================================================
# 01_covariate_setup.R
# Downloads covariate rasters from Google Drive, verifies they share one grid,
# and builds the study-domain raster.
#
# Inputs:  Google Drive folder "sudan_enm_covariates" (GEE exports); GADM;
#          compiled_vl_presences.csv (occurrence years)
# Outputs: Rasters in data/raw/ (downloaded, not modified); DOMAIN_FILE
#
# The download is one-time. The domain step rebuilds every run.
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(tidyverse)
})

# ----- Download rasters (skipped if already present) ------
# Annual rasters exist only for years with occurrence records
# (python/02_covariates.ipynb), so the required years come from the data.
occ_years <- sort(unique(read.csv(OCC_RAW_FILE)$year))
cat("Occurrence years:", occ_years, "\n")

annual_files <- expand_grid(pattern = unname(COV_ANNUAL), year = occ_years) |>
  mutate(file = str_replace(pattern, fixed("{year}"), as.character(year))) |>
  pull(file)

# Single-year rasters for 18's projection year (python/04_2025_rasters.ipynb)
proj_files <- str_replace(unname(PROJ_ANNUAL), fixed("{year}"), as.character(PROJ_YEAR))

required <- c(unname(COV_FILES), annual_files, proj_files)
missing  <- required[!file.exists(file.path(DIR_COVARIATES, required))]

if (length(missing) == 0) {
  cat("All", length(required), "required rasters present - skipping Drive download\n")
} else {
  cat(length(missing), "of", length(required),
      "required rasters missing - connecting to Drive\n")


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

    missing <- required[!file.exists(file.path(DIR_COVARIATES, required))]
  if (length(missing) > 0) {
    print(missing)
    stop("Required rasters still missing after Drive download")
  }
}

# ---- Verify spatial alignment ----------------------------------------------
# Checks GEE-exported covariates only

tif_files <- list.files(DIR_COVARIATES, pattern = "\\.tif$", full.names = TRUE)
tif_files <- tif_files[basename(tif_files) != basename(TT_FILE)]

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
  raster_info |>
    distinct(crs, xres, yres, nrow, ncol, xmin, xmax, ymin, ymax) |>
    print()
  stop("Rasters are not on a single grid; the domain and all extraction assume one")
}

# ---- Study area and domain -------------------------------------------------
# Sudan as administered (params.R). Outline, states and domain raster are all
# built from this one boundary; other scripts read the files.

template  <- rast(file.path(DIR_COVARIATES, COV_FILES[1]))
adm0_gadm <- geodata::gadm(country = "SDN", level = 0, path = DIR_RAW)
adm1_gadm <- geodata::gadm(country = "SDN", level = 1, path = DIR_RAW)
adm2_gadm <- geodata::gadm(country = "SDN", level = 2, path = DIR_RAW)

unit <- adm2_gadm[adm2_gadm$NAME_2 == EXCLUDE_ADM2, ]
stopifnot("EXCLUDE_ADM2 not found in GADM level 2" = nrow(unit) == 1)
tri  <- crop(unit, ext(xmin(unit), xmax(unit), EXCLUDE_NORTH_OF, ymax(unit)))

# Everything north of the parallel, east of the Triangle's western vertex
excluded <- as.polygons(ext(xmin(tri), xmax(adm0_gadm) + 1,
                            EXCLUDE_NORTH_OF, ymax(adm0_gadm) + 1),
                        crs = crs(adm0_gadm))
adm0 <- erase(adm0_gadm, excluded)
adm1 <- erase(adm1_gadm, excluded)

cat("Cut at", EXCLUDE_NORTH_OF, "N east of", round(xmin(tri), 3), "E; removed",
    round(sum(expanse(adm0_gadm, unit = "km")) - sum(expanse(adm0, unit = "km"))),
    "km2\n")

stopifnot("Study area extends beyond the covariate grid" =
  xmin(adm0) >= xmin(template) && xmax(adm0) <= xmax(template) &&
  ymin(adm0) >= ymin(template) && ymax(adm0) <= ymax(template))

writeVector(adm0, ADM0_FILE, overwrite = TRUE)
writeVector(adm1, ADM1_FILE, overwrite = TRUE)

# For maps: the full GADM outline, and the excluded area drawn as not modelled
writeVector(adm0_gadm, DISPLAY_ADM0_FILE, overwrite = TRUE)
writeVector(adm1_gadm, DISPLAY_ADM1_FILE, overwrite = TRUE)
writeVector(terra::intersect(adm0_gadm, excluded), EXCLUDED_FILE, overwrite = TRUE)

# Domain raster: 1 inside, 0 outside. touches = TRUE keeps boundary cells,
# as terra::mask() does with a polygon.
domain <- rasterize(adm0, template, touches = TRUE, background = 0)
names(domain) <- "domain"
stopifnot(compareGeom(domain, template))
writeRaster(domain, DOMAIN_FILE, datatype = "INT1U", overwrite = TRUE)

# Coverage: the GEE export region has curved (geodesic) top and bottom edges,
# so the grid's bounding box overstates coverage. Every domain cell must have
# every static covariate.
static <- rast(file.path(DIR_COVARIATES,
                         COV_FILES[c("elevation", "slope", "river_dist", "vertisols")]))
n_gap <- global(any(is.na(static)) & (domain == 1), "sum")[[1]]
stopifnot("Domain cells outside the covariate export region" = n_gap == 0)


# Diagnostics
print(freq(domain))

dom_km2  <- global(cellSize(domain, unit = "km") * domain, "sum")[[1]]
area_km2 <- sum(expanse(adm0, unit = "km"))
cat("Domain area:", round(dom_km2), "km2 | Study area:", round(area_km2),
    "km2 | ratio:", round(dom_km2 / area_km2, 3), "\n")

covs   <- rast(file.path(DIR_COVARIATES, COV_FILES))
n_na   <- global(any(is.na(covs)) & (domain == 1), "sum")[[1]]
cat("Domain cells with an NA long-term covariate:", n_na, "\n")

cat("01_covariate_setup.R complete\n")