# Data

This directory contains all input data for the analysis pipeline. Most files are too large for version control and must be acquired separately. This document describes every dataset, its source, and the expected file path.

## Included in Repository

| File | Description |
|------|-------------|
| `raw/compiled_vl_presences.csv` | Compiled VL occurrence records from published literature and digitized study maps (see `digitize_maps/`) |
| `processed/occurrences_thinned.csv` | Spatially thinned occurrence records (5 km thinning distance), produced by `R/02_spatial_thinning.R` |

## Requires Download

All environmental covariates were extracted via Google Earth Engine (GEE) and exported to Google Drive. The Python notebooks in `python/` document the exact GEE collections, compositing logic, scale factors, and export parameters for every raster. All layers are exported at 1 km resolution in EPSG:4326, clipped to the study extent (21.5°E–39°E, 8.5°N–22.5°N).

### GEE Export Configuration

Dynamic covariates are extracted as annual composites for each occurrence year (2000–2016 annually, then 2018, 2020, 2022, 2024), plus a long-term mean (2000–2024) for the prediction surface. Seasonal composites use wet season (June–October) and dry season (November–May) windows where applicable.

Extracted in `python/02_covariates.ipynb` unless noted otherwise.

### Ecological Mask

Binary raster identifying areas receiving ≥150 mm mean annual rainfall, used to constrain the study area to ecologically plausible VL transmission zones. Derived from CHIRPS long-term mean annual rainfall in `python/01_study_area.ipynb`.

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

**Normalised Difference Vegetation Index (MODIS MOD13Q1)**

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

Used for sampling bias diagnostics (scripts 11–12) and the accessibility-weighted null model (script 13). Not a model covariate. Source: Malaria Atlas Project (https://malariaatlas.org/research-project/accessibility-to-cities/). Downloaded programmatically in `R/11_accessibility_bias_diagnostic.R` if not already present.

```
raw/weiss_travel_time.tif
```

### Population

**WorldPop 2025 Constrained (100 m)**

UN-adjusted constrained population estimates for Sudan. Used for the primary ARP estimation (script 08) and all downstream comparisons. Source: WorldPop (https://www.worldpop.org/). Downloaded programmatically in `R/08_pop_estimate.R` if not already present.

```
raw/population/sdn_pop_2025_100m_constrained.tif
```

**WorldPop 2005 UN-adjusted (~100 m)**

Historical population estimates for the 2005 hindcast comparison against Alvar et al. (2006). Source: WorldPop (https://www.worldpop.org/). Must be downloaded manually from the WorldPop archive.

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
```

## Generated by Pipeline

The following files are produced by the R scripts and do not need to be downloaded.

| File | Produced by | Description |
|------|-------------|-------------|
| `processed/ecological_mask.tif` | `R/01_covariate_setup.R` | Ecological mask aligned to the covariate grid |
| `processed/background_points.csv` | `R/04_background_sampling.R` | 10,000 uniform background points within the ecological mask |

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