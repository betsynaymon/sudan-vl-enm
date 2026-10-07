# Estimating Populations Living in Environmental Risk of Visceral Leishmaniasis in Sudan

Ecological niche model estimating populations living in areas of environmental suitability for visceral leishmaniasis (VL) in Sudan at 1 km resolution, using publicly available environmental covariates under conditions of data scarcity. MaxEnt primary model with random forest and gradient boosted tree comparators and a set of robustness and sensitivity tests.

This repository contains the analysis pipeline for an MSc dissertation, *Geographies of Neglect: Bridging Evidence Gaps to Map Visceral Leishmaniasis in Sudan*, submitted to the London School of Economics and Political Science, Department of Methodology (Applied Social Data Science), August 2026. The pipeline has since been refactored for a journal manuscript, so its results differ from the dissertation (see [Changes](#changes)).

## Changes
Pipeline was audited and refactored post-submission. Primary change is removing the 150mm rainfall ecological mask that restricted background points. The mask included a judgement about where VL could not occur, but evidence from Oriole Global Health and one vector-positive site outside of the mask suggest making that judgement outright led the estimate to extrapolate beyond the model calibration area. 

| | Dissertation (August 2026) | Current pipeline |
| --- | --- | --- |
| Background| Sampled within the $\geq$ 150mm rainfall mask | Sampled uniformly across all of Sudan |
| Presence records | 98 | 103 (5 Khartoum points that originally fell outside of mask included) |
| MaxEnt configuration | LQH, regularization multiplier 1.5 | LQHP, regularization multiplier 4 |

## Repository Structure

```
sudan-vl-enm/
├── R/                                    # Analysis pipeline (see Run order)
│   ├── params.R                          # Paths, constants and modelling decisions
│   ├── helpers.R                         # Shared functions
│   ├── plotting_theme.R                  # Shared figure styling
│   ├── 01_covariate_setup.R
│   ├── 02_spatial_thinning.R
│   ├── 03_collinearity.R
│   ├── 04_background_sampling.R
│   ├── 05_spatial_cv_folds.R
│   ├── 06_maxent_tuning.R
│   ├── 06b_plateau_curves.R
│   ├── 07_model_visualization.R
│   ├── 08_pop_estimate.R
│   ├── 09_rf_comparator.R
│   ├── 10_data_quality_refits.R
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
│   ├── 21_background_extent.R
│   ├── 22_state_accessibility.R
│   ├── 23_unverified_low_prediction.R    # Retired; not run
│   ├── 24_accessibility_covariate_test.R
│   └── fig_study_area.R                  # Study-area and endemic-status maps
├── python/                               # Earth Engine exports and a digitization check
│   ├── 01_study_area.ipynb
│   ├── 02_covariates.ipynb
│   ├── 03_digitized_points_check.ipynb
│   └── 04_2025_rasters.ipynb
├── digitize_maps/                        # Occurrence extraction from published maps
│   ├── digitize_map.py
│   ├── map_config.json
│   └── README.md
├── data/
│   ├── raw/                              # Occurrence records, state endemic status, downloaded inputs
│   └── processed/                        # Thinned records, background points, study-area files
├── outputs/
│   ├── figures/                          # Maps and plots
│   ├── tables/                           # Summary CSVs
│   ├── models/                           # Fitted models, training data and CV folds (.rds)
│   └── surfaces/                         # Predicted surfaces (.tif; generated, not committed)
└── README.md
```

## Pipeline Overview

The analysis proceeds in three stages. The Python notebooks and map digitization are upstream data acquisition steps; the R scripts are the analytical pipeline.

### Data Acquisition

The Python notebooks (`python/`) use Google Colab to extract environmental covariates from Google Earth Engine. They are included as documentation of data provenance rather than steps to rerun locally. `04_2025_rasters.ipynb` extracts single-year 2025 LST night and rainfall rasters used by script 18 for the temporal sensitivity diagnostic.

The map digitization workflow (`digitize_maps/`) extracts georeferenced VL occurrence points from published study maps. It is a standalone preprocessing step whose outputs feed into the compiled occurrence dataset. The source map images are extracts of copyrighted journal figures and are not included here; the digitization code and configuration are included, but the workflow is not rerunnable without the original figures.

### R Analysis Pipeline

All numbered R scripts are in `R/` and should be run in numerical order. Shared configuration (file paths, constants, model settings and decisions) is stored in `params.R`; `helpers.R` holds shared functions. `plotting_theme.R` holds the shared figure theme sourced by the plotting scripts. `fig_study_area.R` produces the study-area maps for the dissertation (Figures 3 & 4) and is run independently of the numbered sequence.

**Data preparation (01–05):** Covariate alignment and stacking, spatial thinning of occurrence records, collinearity screening, background point sampling, and spatial cross-validation fold assignment.

**Model fitting and estimation (06–08):** MaxEnt hyperparameter tuning via spatial block cross-validation with year-matched covariate extraction, model visualization (suitability surface, response curves, variable importance), and population-at-risk estimation using WorldPop 2025 constrained population.

**Robustness suite (09–24):** Spans signal validation, data robustness, model robustness, algorithm robustness, external validation, confound assessment, and scope diagnostics:

Core: 

| Script | Test | Category |
|---|---|---|
| 10 | Data-quality refits: without imprecise digitized records, facility reports, Khartoum case records (alone and with Wad Madani's), or ID 121 | Data robustness |
| 11 | Accessibility bias diagnostic: travel time at the presences against the background and the population | Confound assessment |
| 12 | Sampling bias correction: background weighted by accessibility at two strengths | Confound assessment |
| 13 | Null model test: uniform and accessibility-matched null presences, cross-validated as in 06 | Signal validation |
| 15 | Variable sensitivity: seasonal LST swaps, NDVI for rainfall, three added covariates, river distance dropped | Model robustness |
| 16 | State-level agreement with literature endemic status | External comparison |
| 20 | East–west diagnostic: Darfur presences against the rest, on the map, held out, and across alternative surfaces | Scope diagnostics |
| 22 | State-level accessibility: share of each state's population more than five hours from a city | Confound assessment |
| 24 | Travel time as a covariate, integrated out for prediction | Confound assessment |


Supplementary:
 
| Script | Test | Category |
|---|---|---|
| 09 | Random forest comparator (down-sampled) | Algorithm robustness |
| 14 | Uncertainty surface: refits without each spatial fold | Model robustness |
| 17 | Gradient boosted tree comparator, with the three-algorithm comparison | Algorithm robustness |
| 18 | Single-year projections: each occurrence year and 2025 | Temporal sensitivity |
| 19 | 2005 hindcast: Gedaref against Alvar et al. (2006) | External comparison |
| 21 | Background extent: background restricted to the ≥150 mm region | Data robustness |
 
Script 23 (unverified low-prediction extent) is retired.

## Run order
 
Paths are set in `params.R` with `here`, relative to the repository root. Run the numbered scripts in order, with two exceptions: 06b runs between 06 and 07, and 24 runs before 20, which reads its surfaces.
 
```
01 → 02 → 03 → 04 → 05 → 06 → 06b → 07 → 08
09 → 10 → 11 → 12 → 13 → 14 → 15 → 16 → 17 → 18 → 19
24 → 20 → 21 → 22
```
  
For the core analyses alone, run 01–08 (with 06b), then 10, 11, 12, 13, 15, 16, 24, 20 and 22. The supplementary scripts need only core outputs, except that 17 needs 09 and 19 needs 18.
 
The MaxEnt configuration is fixed in `params.R` (`SELECTED_FC`, `SELECTED_RM`). To select it again, for example after changing the data or the tuning grid, set both to `NULL` so that 06 applies its automatic rule, run 06 and 06b, set the chosen configuration, and run 06 and 06b again before 07.


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

Modelling and evaluation: `maxnet`, `blockCV`, `ecospat`, `ranger`, `gbm`, `dismo`

Spatial data and thinning: `terra`, `sf`, `raster`, `geodata`, `rnaturalearth`, `spThin`, `rnaturalearthdata`

Data handling and I/O: `tidyverse` (provides `dplyr`, `ggplot2`, `tidyr`, `tibble`), `here`, `googledrive`, `httr`

Figures: `patchwork`, `ggspatial`, `ggrepel`, `ggpattern`, `scales`. PDFs written with `grDevices::cairo_pdf`. 

### Python

Python notebooks require a Google Earth Engine account and were run in Google Colab. Dependencies: `ee`, `geemap`.

Map digitization: `digitize_maps/digitize_map.py` requires Python with `opencv-python`, `numpy`, and `pandas`.

## License

Code is licensed under the MIT License. See `LICENSE`.

The license covers the code in this repository only. Environmental covariate data are subject to their respective terms of use. The compiled occurrence dataset (`compiled_vl_presences.csv`) is derived from published literature and is provided for reproducibility.

## Citation

Naymon, B. (2026). *Geographies of Neglect: Bridging Evidence Gaps to Map Visceral Leishmaniasis in Sudan.* MSc dissertation, London School of Economics and Political Science.