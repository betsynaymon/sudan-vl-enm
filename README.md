# Estimating Populations at Risk of Visceral Leishmaniasis in Sudan

Ecological niche model estimating populations at risk of visceral leishmaniasis (VL) in Sudan at 1 km resolution, using publicly available environmental covariates under conditions of data scarcity. MaxEnt primary model with random forest and gradient boosted tree comparators and a 13-test robustness suite.

This repository contains the analysis pipeline for an MSc dissertation submitted to the London School of Economics and Political Science, Department of Methodology (Applied Social Data Science), August 2025.

## Repository Structure

```
sudan-vl-enm/
├── R/                        # Analysis pipeline (run in order)
│   ├── params.R              # Shared parameters and configuration
│   ├── 01_covariate_setup.R
│   ├── 02_spatial_thinning.R
│   ├── ...
│   └── 19_hindcast.R
├── python/                   # Google Earth Engine data acquisition notebooks
│   ├── 01_study_area.ipynb
│   ├── 02_covariates.ipynb
│   ├── 03_digitized_points_check.ipynb
│   └── 04_2025_rasters.ipynb
├── digitize_maps/            # Standalone: occurrence point extraction from published maps
│   ├── digitize_map.py
│   ├── map_config.json
│   ├── maps/                 # Source map images
│   └── outputs/              # Digitized coordinate CSVs
├── data/
│   ├── raw/                  # Source datasets (see Data section below)
│   └── processed/            # Pipeline intermediates
├── outputs/
│   ├── figures/              # Maps and plots
│   ├── tables/               # Summary CSVs
│   ├── models/               # Fitted model objects (.rds)
│   └── surfaces/             # Predicted suitability and binary rasters (.tif)
└── README.md
```

## Pipeline Overview

The analysis proceeds in three stages. The Python notebooks and map digitization are upstream data acquisition steps; the R scripts are the analytical pipeline.

### Data Acquisition

The **Python notebooks** (`python/`) were run in Google Colab to extract environmental covariates from Google Earth Engine. They are included as documentation of data provenance rather than steps to rerun locally.

The **map digitization** workflow (`digitize_maps/`) extracts georeferenced VL occurrence points from published study maps. It is a standalone preprocessing step whose outputs feed into the compiled occurrence dataset.

### R Analysis Pipeline

All R scripts are in `R/` and should be run in numerical order. Shared configuration (file paths, CRS, study area extent, ecological mask threshold, model settings) is stored in `params.R`.

**Data preparation (01–05):** Covariate alignment and stacking, spatial thinning of occurrence records, collinearity screening, background point sampling, and spatial cross-validation fold assignment.

**Model fitting and estimation (06–08):** MaxEnt hyperparameter tuning via ENMeval with spatial block cross-validation, model visualization (suitability surface, response curves, variable importance), and population-at-risk estimation using WorldPop 2025 constrained population.

**Robustness suite (09–19):** Thirteen tests spanning signal validation, data robustness, model robustness, algorithm robustness, external validation, confound assessment, and scope diagnostics:

| Script | Test | Category |
|--------|------|----------|
| 09 | Random forest comparator | Algorithm robustness |
| 10 | Precision sensitivity | Data robustness |
| 11 | Accessibility bias diagnostic | Confound assessment |
| 12 | Sampling bias correction | Data robustness |
| 13 | Null model test (uniform + accessibility-weighted) | Signal validation |
| 14 | Uncertainty surface (fold jackknife) | Model robustness |
| 15 | Variable sensitivity | Model robustness |
| 16 | State-level literature validation | External validation |
| 17 | Gradient boosted tree comparator | Algorithm robustness |
| 18 | 2025 single-year projection | Model robustness |
| 19 | 2005 hindcast | External validation |

## Data

### Included in this repository

| File | Description |
|------|-------------|
| `data/raw/compiled_vl_presences.csv` | Compiled VL occurrence records from published literature and digitized maps |
| `data/processed/occurrences_thinned.csv` | Spatially thinned occurrence records (5 km thinning distance) |

### Requires download

All environmental covariates, administrative boundaries, and population data must be acquired separately. See [`data/README.md`](data/README.md) for sources, download instructions, and expected file paths.

| Dataset | Source | Resolution |
|---------|--------|------------|
| Land surface temperature (day/night, annual/seasonal) | MODIS MOD11A1 via Google Earth Engine | 1 km |
| NDVI (annual/seasonal) | MODIS MOD13A2 via Google Earth Engine | 1 km |
| Rainfall | CHIRPS v2.0 via Google Earth Engine | 1 km (resampled from ~5 km) |
| Slope | SRTM via `geodata::elevation_30s()` | 1 km (resampled from ~90 m) |
| River distance | HydroSHEDS | 1 km |
| Tree cover | MODIS MOD44B via Google Earth Engine | 1 km (resampled from 250 m) |
| Vertisols | HWSD v2.0 (WRB2_CODE = 33) | 1 km |
| Travel time to nearest city | Weiss et al. (2018) | 1 km |
| Population | WorldPop 2025 constrained | 100 m |
| Administrative boundaries | GADM v4.1 (levels 0 and 1) | — |

## Requirements

### R

Core packages: `maxnet`, `ENMeval`, `blockCV`, `ranger`, `gbm`, `terra`, `sf`, `dplyr`, `ggplot2`, `patchwork`, `ggspatial`, `geodata`

### Python

Python notebooks require a Google Earth Engine account and were run in Google Colab. Dependencies: `ee`, `geemap`.

### Map digitization

`digitize_maps/digitize_map.py` requires Python with `opencv-python`, `numpy`, and `pandas`.

## License

Code is licensed under the MIT License. See [LICENSE](LICENSE).

The license covers the code in this repository only. Environmental covariate data are subject to their respective terms of use. The compiled occurrence dataset (`compiled_vl_presences.csv`) is derived from published literature and is provided for reproducibility.

## Citation

Naymon, B. (2025). Estimating populations at risk of visceral leishmaniasis in Sudan: An ecological niche modelling approach under data scarcity. MSc dissertation, London School of Economics and Political Science.
