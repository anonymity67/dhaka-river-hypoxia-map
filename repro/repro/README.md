# Hypoxia risk in Dhaka's rivers: reproducible analysis

Code for a Bayesian spatial (INLA-SPDE) model of the probability that dissolved oxygen (DO) falls below 2 mg/L in Dhaka's peripheral rivers, by season. Hypoxia is a habitat condition; no fish or fisheries data are used. A public map built from the results is at the repository `dhaka-river-hypoxia-map`.

## Data (not included; download each file)

| File to place in `data/` | Source |
|---|---|
| `Water quality.csv` | Sinha, P., Debnath, N., & Bhattacharjee, P. (2026). *Water Quality and Hydroclimatic Dataset for Rivers and Streams of Dhaka, Bangladesh* (Version 1). Mendeley Data. https://doi.org/10.17632/xftpv7jvx3.1 (CC BY 4.0). |
| `waterways.geojson` | OpenStreetMap waterways for Bangladesh (HOT export, "hotosm_bgd_waterways_osm_geojson", via the Humanitarian Data Exchange). (c) OpenStreetMap contributors, ODbL. |
| `pollution_int.tif` (script 07 only) | Kaggle, Bangladesh-River-Pollution-GIS (MIT licence). |
| `Pollution_Intensity_Coordinates_Updated.csv` (Python scripts only) | Kaggle, Bangladesh-River-Pollution-GIS (MIT licence), the point-grid file. Place in `data/`. |

Paths are set in `config.R`. Scripts must be run from the repository root.

## Software

- R 4.x with `INLA` (R-INLA 26.08.07 was used), `fmesher`, `sf`, `Matrix`, `terra` (script 07), `jsonlite`.
  Install INLA with `install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)`.
- Python 3 with `pandas`, `numpy`, `scipy`, `statsmodels`, `rasterio` (Python scripts only).
- Google Earth Engine account (Sentinel-2 step only).

## Run order

| Script | What it does | Time |
|---|---|---|
| `R/01_main_model.R` | Descriptives, model A (no map) and model B (with map), season odds ratios, spatial range, WAIC; files for Earth Engine. `RUN_LOSO <- TRUE` in `config.R` adds leave-one-station-out AUC. | minutes (hours with LOSO) |
| `R/02_river_risk_map.R` | Risk every 500 m along OpenStreetMap rivers, river tables, calibration check, and the data for the map website (`results/site/`). | minutes |
| `R/03_robustness.R` | Mesh size, station effect, Laplace strategy, and model C (wet-season field). | 10-20 min |
| `R/04_without_2020.R` | Refit without 2020: do the conclusions hold? | minutes |
| `R/05_year_2020_checks.R` | What the data say about the 2020 dip (sampling, rain, other variables). | minutes |
| `R/06_sentinel2.R` | Sentinel-2 indices added to the model (exploratory). Needs `results/s2_indices_per_sample.csv` (below). | 10-20 min, longer with LOSO |
| `R/07_pollution_raster_check.R` | Kaggle raster vs station risk. Needs script 02 and `data/pollution_int.tif`. | minutes |
| `python/01_descriptives_hotspots_grid.py` | Beta-binomial river and season tables, logistic regression, Moran's I, Getis-Ord Gi* hotspots, Kaggle point-grid check. | seconds |
| `python/02_raster_check.py` | Raster vs station check (Python version). | seconds |
| `extras/osm_sources_test.R` | Not in the paper: OpenStreetMap industrial sources vs hypoxia. Null result, kept for transparency. | minutes |

`run_all.R` runs scripts 01-05 and then 06-07 if their inputs exist.

### Sentinel-2 step (script 06)

1. Run script 01 (it writes `results/samples_for_gee.csv` and `results/stations_for_gee.csv`).
2. In Earth Engine, upload `samples_for_gee.csv` as a table asset (X = `lon`, Y = `lat`), paste `gee/gee_samples.js` into the Code Editor, set the asset path at the top, run it and start the export task. Save the result as `results/s2_indices_per_sample.csv`.
3. Run `R/06_sentinel2.R`.

## Values to expect (paper)

Use these to check your run. Small differences can come from INLA or package versions.

| Quantity | Value |
|---|---|
| Samples, stations, rivers | 1,495; 60; 12 |
| Share of samples DO < 5 mg/L, DO < 2 mg/L | 0.82; 0.57 |
| Mesh nodes | 1,136 |
| WAIC, model A and model B | 1450.9; 1060.4 |
| Season odds ratios vs winter (pre-monsoon, monsoon, post-monsoon) | 0.47, 0.006, 0.029 |
| Spatial range (median, km) and spatial sd | about 26-28; about 2.7 |
| Leave-one-station-out AUC, A and B | 0.813; 0.908 |
| Without 2020: monsoon and post-monsoon odds ratios | 0.006; 0.033 |
| Without 2020: field correlation and shared top-15 hotspots | 0.97; 12 of 15 |
| Calibration (observed and fitted hypoxia share) | 0.567; 0.567 |

## Differences from the exploratory scripts

This repository consolidates scripts that were first run interactively. Two small changes: river-line points are assigned to their own line (earlier: nearest line), and each point takes the river of its nearest station. Neither changes the model fit.

## Limits

Hypoxia is a habitat condition (DO < 2 mg/L), not a measure of fish or other life. Sampling was about once per season, 2020 differs from other years (see script 05), COD appears filled in and is not used, and the pollution raster has no documented method.

## Licences and citation

Code: MIT (`LICENSE`). Derived data: CC BY 4.0. Cite the Mendeley dataset above; add the paper reference when available.
