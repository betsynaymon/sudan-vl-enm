# Data

This directory contains all input data for the analysis pipeline. Most files are too large for version control and must be acquired separately. This document describes every dataset, its source, and the expected file path. File names are set in `R/params.R`. 

## Included in Repository

| File | Description |
|------|-------------|
| `raw/compiled_vl_presences.csv` | Compiled VL occurrence records from published literature and digitized study maps (see `digitize_maps/`) |
| `raw/state_endemic_status.csv` | VL endemic status of each state from the literature, used by `R/16_qualitative_state_validation.R` and `R/fig_study_area.R` |
| `processed/occurrences_thinned.csv` | Spatially thinned occurrence records (5 km thinning distance), produced by `R/02_spatial_thinning.R` |

## Requires Download

All environmental covariates were extracted via Google Earth Engine (GEE) and exported to Google Drive. The Python notebooks in `python/` document the exact GEE collections, compositing logic, scale factors, and export parameters for every raster. All layers are exported at 1 km resolution in EPSG:4326 over a boundary rectangle (21.5°E–39°E, 8.5°N–22.5°N). 

All other surfaces are downloaded by the R script that first needs it. 

### Covariates used in the final model

Not all layers below are model predictors. The final MaxEnt model uses five —
**slope**, **distance to rivers**, **vertisols**, **annual nighttime LST**, and
**annual rainfall** — selected on ecological grounds during collinearity screening (`R/03_collinearity.R`).

The rest are extracted but not used in the primary model:

- **NDVI, tree cover, elevation, daytime LST** — retained for collinearity
  screening and as alternative cluster representatives in the variable-sensitivity
  analysis (`R/15_variable_sensitivity.R`).
- **Seasonal (wet/dry) composites** — used only by that sensitivity analysis; the
  wet- and dry-season nighttime LST layers are the seasonal temperature
  representatives it swaps in. 
- **Travel time** — used for the sampling-bias diagnostics and accessibility-weighted
  null (scripts 11–13, 22), and as a sixth covariate in `R/15_variable_sensitivity.R`.
- **≥150 mm rainfall raster** — defines the region for within-belt diagnostics and the background-extent check (script 21).
- **Population and administrative boundaries** — used downstream of the model, not as predictors.

### GEE Export Configuration

Dynamic covariates are extracted as annual composites for each occurrence year (2000–2016 annually, then 2018, 2020, 2022, 2024), plus a long-term mean (2000–2024) for the prediction surface. Seasonal composites use wet season (June–October) and dry season (November–May) windows where applicable.

Extracted in `python/02_covariates.ipynb` unless noted otherwise.

### Ecological Mask

Binary raster: 1 where CHIRPS mean annual rainfall (2000–2024) is at least 150 mm, 0 elsewhere, coded over the export rectangle rather than Sudan. Derived in `python/01_study_area.ipynb`.
 
This was the dissertation's ecological mask. It no longer constrains the study area or the background. Scripts intersect it with the study domain to report results within the ≥150 mm region (for example, within-belt CBI in script 06), and script 21 uses it to restrict the background in a supplementary check, after confirming that it matches the rainfall covariate at `SENS_MASK_MM` (150 mm) cell for cell.

```
raw/ecological_mask_150mm.tif
```

### Static Covariates

**Elevation and Slope (SRTM)**

Source: NASA SRTM 30m (`USGS/SRTMGL1_003`). Elevation resampled to 1 km via mean reducer. Slope derived in degrees using `ee.Terrain.slope()` before resampling.

```
raw/elevation_1km.tif
raw/slope_1km.tif
```

**Distance to Rivers (HydroSHEDS)**

Source: HydroSHEDS 15 arc-second flow accumulation (`WWF/HydroSHEDS/15ACC`). River pixels defined as cells with >500 upstream cells, which includes seasonal flows (khors) relevant to *P. orientalis* habitat. Distance computed via `fastDistanceTransform`, converted to metres. Gaps in the Sahara where HydroSHEDS has no data are filled with 100,000 m (functionally "no river nearby"). Resampled from ~500 m to 1 km via mean reducer.

```
raw/river_distance_1km.tif
```

**Vertisols (HWSD v2.0)**

Source: FAO Harmonized World Soil Database v2.0 (`projects/sat-io/open-datasets/FAO/HWSD_V2_SMU`). Binary layer: WRB2_CODE = 33 (Vertisol). Resampled to 1 km via mode reducer to preserve the binary classification.

```
raw/vertisols_1km.tif
```

### Dynamic Covariates

**Land Surface Temperature (MODIS MOD11A2)**

Source: MODIS Terra LST 8-day composite, 1 km native resolution. Scale factor 0.02, converted from Kelvin to Celsius. Day and night bands extracted separately, each with annual, wet season, and dry season composites.

```
raw/lst_{day|night}_{annual|dry|wet}_{year}_1km.tif
raw/lst_{day|night}_{annual|dry|wet}_mean_2000_2024_1km.tif
```

2025 single-year nighttime LST for the temporal sensitivity diagnostic (script 18), extracted in `python/04_2025_rasters.ipynb`:

```
raw/lst_night_annual_2025_1km.tif
```

**Normalized Difference Vegetation Index (MODIS MOD13Q1)**

Source: MODIS Terra NDVI 16-day composite, 250 m native resolution. Scale factor 0.0001. Annual, wet season, and dry season composites. Resampled from 250 m to 1 km via mean reducer.

```
raw/ndvi_{annual|dry|wet}_{year}_1km.tif
raw/ndvi_{annual|dry|wet}_mean_2000_2024_1km.tif
```

**Rainfall (CHIRPS v2.0)**

Source: CHIRPS Daily v2.0. Annual composites are the sum of daily rainfall within each calendar year (mm/year). No seasonal split — annual totals are the ecologically relevant metric.

```
raw/rainfall_{year}_1km.tif
raw/rainfall_mean_2000_2024_1km.tif
```

2025 single-year rainfall for the temporal sensitivity diagnostic (script 18), extracted in `python/04_2025_rasters.ipynb`:

```
raw/rainfall_2025_1km.tif
```

**Tree Cover (Hansen Global Forest Change v1.13)**

Source: Hansen Global Forest Change 2025 release (`UMD/hansen/global_forest_change_2025_v1_13`), 30 m native resolution. Tree cover for year *t* = year-2000 baseline with pixels that experienced cumulative forest loss through year *t* set to 0. Values are percentage canopy cover (0–100). Resampled from 30 m to 1 km via mean reducer.

```
raw/treecover_{year}_1km.tif
raw/treecover_mean_2000_2024_1km.tif
```

### Accessibility Surface

**Travel Time to Nearest City (Weiss et al. 2018)**

Used for the accessibility bias diagnostic (script 11), the accessibility-weighted background (script 12), the accessibility-matched null model (script 13) and state-level remoteness (script 22), and as a sixth covariate only in the accessibility-covariate test (script 24) Source: Malaria Atlas Project (https://malariaatlas.org/research-project/accessibility-to-cities/). Downloaded programmatically in `R/11_accessibility_bias_diagnostic.R` if not already present.

```
raw/weiss_travel_time.tif
```

### Population

**WorldPop 2025 Constrained (100 m)**

UN-adjusted constrained population estimates for Sudan. Used for the primary ARP estimation (script 08) and all downstream comparisons. Source: WorldPop (https://www.worldpop.org/). Downloaded programmatically in `R/08_pop_estimate.R` if not already present.

```
raw/population/sdn_pop_2025_100m_constrained.tif
```

**WorldPop 2005 Unconstrained, UN-adjusted (~100 m)**

Historical population estimates for the 2005 hindcast comparison against Alvar et al. (2006). WorldPop has no constrained product before 2015, so this is the unconstrained surface, which spreads people over all land, unlike the 2025 file. Source: WorldPop (https://www.worldpop.org/). Must be downloaded manually from the WorldPop archive.

```
raw/population/sdn_ppp_2005_UNadj.tif
```

### Administrative Boundaries

**GADM v4.1**

Downloaded programmatically in `R/01_covariate_setup.R` using `geodata::gadm()`.

```
raw/gadm/
├── gadm41_SDN_0_pk.rds    # National boundary
└── gadm41_SDN_1_pk.rds    # State boundaries (level 1)
└── gadm41_SDN_2_pk.rds    # Localities (level 2), used to locate the Halaib Triangle
```

## Generated by Pipeline

The following files are produced by the R scripts and do not need to be downloaded.

| File | Produced by | Description |
|------|-------------|-------------|
| `processed/domain_sudan.tif` | `R/01_covariate_setup.R` | Study domain on the covariate grid: 1 in the study area, 0 outside. The calibration and prediction extent |
| `processed/study_adm0.gpkg`, `processed/study_adm1.gpkg` | `R/01_covariate_setup.R` | Study area outline and states |
| `processed/display_adm0.gpkg`, `processed/display_adm1.gpkg` | `R/01_covariate_setup.R` | Full GADM outline and states, for maps |
| `processed/excluded_area.gpkg` | `R/01_covariate_setup.R` | The Halaib Triangle |
| `processed/background_points.csv` | `R/04_background_sampling.R` | 10,000 uniform background points from domain cells with complete covariates, each with a year drawn from the occurrence years; keyed by `bg_id` |

## GEE Export Summary

From `python/02_covariates.ipynb`:

| Category | Rasters |
|----------|---------|
| Static covariates | 4 (elevation, slope, river distance, vertisols) |
| NDVI (annual/wet/dry × 21 years) | 63 |
| LST day (annual/wet/dry × 21 years) | 63 |
| LST night (annual/wet/dry × 21 years) | 63 |
| Rainfall (21 years) | 21 |
| Tree cover (21 years) | 21 |
| Long-term means | 11 |
| **Subtotal** | **246** |

From `python/01_study_area.ipynb`: 1 (ecological mask)

From `python/04_2025_rasters.ipynb`: 2 (LST night 2025, rainfall 2025)

**Total: 249 rasters** exported to `Google Drive > sudan_enm_covariates`