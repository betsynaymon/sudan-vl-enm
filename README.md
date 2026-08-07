# Estimating Populations Living in Environmental Risk of Visceral Leishmaniasis in Sudan

Ecological niche model estimating populations living in areas of environmental suitability for visceral leishmaniasis (VL) in Sudan at 1 km resolution, using publicly available environmental covariates under conditions of data scarcity. MaxEnt primary model with random forest and gradient boosted tree comparators and a set of robustness and sensitivity tests.

This repository contains the analysis pipeline for an MSc dissertation, *Geographies of Neglect: Bridging Evidence Gaps to Map Visceral Leishmaniasis in Sudan*, submitted to the London School of Economics and Political Science, Department of Methodology (Applied Social Data Science), August 2026.

## Repository Structure

```
sudan-vl-enm/
├── R/                        # Analysis pipeline (run in order)
│   ├── params.R              # Shared parameters and configuration
│   ├── plotting_theme.R      # Shared ggplot theme (sourced by figure scripts)
│   ├── 01_covariate_setup.R
│   ├── 02_spatial_thinning.R
│   ├── 03_collinearity.R
│   ├── 04_background_sampling.R
│   ├── 05_spatial_cv_folds.R
│   ├── 06_maxent_tuning.R
│   ├── 07_model_visualization.R
│   ├── 08_pop_estimate.R
│   ├── 09_rf_comparator.R
│   ├── 10_precision_sensitivity.R
│   ├── 11_accessibility_bias_diagnostic.R
│   ├── 12_sampling_bias_robustness.R
│   ├── 13_null_model_test.R
│   ├── 14_uncertainty_surface.R
│   ├── 15_variable_sensitivity.R
│   ├── 16_qualitative_state_validation.R
│   ├── 17_gbt_comparator.R
│   ├── 18_2025_surface.R
│   ├── 19_hindcast.R
│   ├── 20_east_west_diagnostic.R
│   ├── 21_ecological_mask_check.R
│   ├── 22_state_accessibility.R
│   ├── 23_unverified_low_prediction.R
│   └── fig_study_area.R       # Study-area map (Figure 4); run independently
├── python/                   # Google Earth Engine data acquisition notebooks
│   ├── 01_study_area.ipynb
│   ├── 02_covariates.ipynb
│   ├── 03_digitized_points_check.ipynb
│   └── 04_2025_rasters.ipynb
├── digitize_maps/            # Occurrence extraction from published maps (see note below)
│   ├── digitize_map.py
│   └── map_config.json
├── data/
│   ├── raw/                  # Compiled occurrence records (see Data section)
│   └── processed/            # Pipeline intermediates
├── outputs/
│   ├── figures/              # Maps and plots
│   ├── tables/               # Summary CSVs
│   ├── sensitivity/          # Sensitivity-analysis outputs (e.g. no-mask refit)
│   ├── models/               # Fitted model objects (.rds)
│   └── surfaces/             # Predicted suitability/binary rasters (.tif; generated, not committed)
└── README.md
```

## Pipeline Overview

The analysis proceeds in three stages. The Python notebooks and map digitization are upstream data acquisition steps; the R scripts are the analytical pipeline.

### Data Acquisition

The Python notebooks (`python/`) were run in Google Colab to extract environmental covariates from Google Earth Engine. They are included as documentation of data provenance rather than steps to rerun locally. `04_2025_rasters.ipynb` extracts single-year 2025 LST night and rainfall rasters used by script 18 for the temporal sensitivity diagnostic.

The map digitization workflow (`digitize_maps/`) extracts georeferenced VL occurrence points from published study maps. It is a standalone preprocessing step whose outputs feed into the compiled occurrence dataset. The source map images are extracts of copyrighted journal figures and are not included here; the digitization code and configuration are included, but the workflow is not rerunnable without the original figures.

### R Analysis Pipeline

All numbered R scripts are in `R/` and should be run in numerical order. Shared configuration (file paths, CRS, study area extent, ecological mask threshold, model settings) is stored in `params.R`; `plotting_theme.R` holds the shared figure theme sourced by the plotting scripts. `fig_study_area.R` produces the study-area maps for the dissertation (Figures 3 & 4) and is run independently of the numbered sequence.

**Data preparation (01–05):** Covariate alignment and stacking, spatial thinning of occurrence records, collinearity screening, background point sampling, and spatial cross-validation fold assignment.

**Model fitting and estimation (06–08):** MaxEnt hyperparameter tuning via spatial block cross-validation with year-matched covariate extraction, model visualization (suitability surface, response curves, variable importance), and population-at-risk estimation using WorldPop 2025 constrained population.

**Robustness suite (09–23):** Spans signal validation, data robustness, model robustness, algorithm robustness, external validation, confound assessment, and scope diagnostics:

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
| 20 | East–west diagnostic | Scope diagnostics |
| 21 | Ecological mask sensitivity | Data robustness |
| 22 | State-level accessibility | Confound assessment |
| 23 | Unverified low-prediction extent | Scope diagnostics |

## Data

### Included in this repository

| File | Description |
|------|-------------|
| `data/raw/compiled_vl_presences.csv` | Compiled VL occurrence records from published literature and digitized maps |
| `data/processed/occurrences_thinned.csv` | Spatially thinned occurrence records (5 km thinning distance) |
| `outputs/models/spatial_cv_folds.rds` | Spatial block CV fold assignments (committed for exact reproducibility of CV metrics) |

Fitted figures, summary tables (`outputs/tables/`), and sensitivity outputs (`outputs/sensitivity/`) are also committed. Large environmental rasters and predicted surfaces (`outputs/surfaces/`) are generated by the pipeline and excluded.

### Requires download

All environmental covariates, administrative boundaries, and population data must be acquired separately. See `data/README.md` for sources, download instructions, and expected file paths.

| Dataset | Source | Resolution |
|---------|--------|------------|
| Land surface temperature (day/night, annual/seasonal) | MODIS MOD11A2 via Google Earth Engine | 1 km |
| NDVI (annual/seasonal) | MODIS MOD13Q1 via Google Earth Engine | 1 km (resampled from 250 m) |
| Rainfall | CHIRPS v2.0 via Google Earth Engine | 1 km (resampled from ~5 km) |
| Elevation and slope | SRTM via Google Earth Engine | 1 km (resampled from ~30 m) |
| River distance | HydroSHEDS via Google Earth Engine | 1 km |
| Tree cover | Hansen Global Forest Change v1.13 via Google Earth Engine | 1 km (resampled from 30 m) |
| Vertisols | HWSD v2.0 (WRB2_CODE = 33) via Google Earth Engine | 1 km |
| Travel time to nearest city | Weiss et al. (2018) via Malaria Atlas Project | 1 km |
| Population (2025) | WorldPop 2025 constrained | 100 m |
| Population (2005) | WorldPop 2005 UN-adjusted | ~100 m |
| Administrative boundaries | GADM v4.1 (levels 0 and 1) via geodata R package | — |

### Not included: clinical surveillance data

The dissertation also draws on clinical surveillance records from Oriole Global Health—monthly case counts and operational status for Sudan's VL treatment facilities (2023–2025)—used descriptively to characterize where patients present for treatment and to document infrastructure collapse, not as a model input. These data were accessed under a data-sharing agreement within a secure data environment with no external access. Neither the data nor the scripts used to analyze them are included in this repository, and the facility-level results reported in the dissertation (the state-level comparison table and the facility-functionality figure) are not reproducible from this repo. The pipeline here covers the public-data ENM component only.

## Requirements

### R

Modelling and evaluation: `maxnet`, `blockCV`, `ecospat`, `ranger`, `gbm`

Collinearity: `usdm`

Spatial data and thinning: `terra`, `sf`, `geodata`, `rnaturalearth`, `spThin`

Data handling and I/O: `tidyverse` (provides `dplyr`, `ggplot2`, `tidyr`, `tibble`), `here`, `googledrive`, `httr`

Figures: `patchwork`, `ggspatial`, `ggrepel`, `ggpattern`, `scales`

### Python

Python notebooks require a Google Earth Engine account and were run in Google Colab. Dependencies: `ee`, `geemap`.

Map digitization: `digitize_maps/digitize_map.py` requires Python with `opencv-python`, `numpy`, and `pandas`.

## License

Code is licensed under the MIT License. See `LICENSE`.

The license covers the code in this repository only. Environmental covariate data are subject to their respective terms of use. The compiled occurrence dataset (`compiled_vl_presences.csv`) is derived from published literature and is provided for reproducibility.

## Citation

Naymon, B. (2026). *Geographies of Neglect: Bridging Evidence Gaps to Map Visceral Leishmaniasis in Sudan.* MSc dissertation, London School of Economics and Political Science.