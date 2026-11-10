# =====================================================================
# Configuration: 
# =====================================================================
# Run every script from the repository root (the folder that holds this file).
# Defaults are relative to the repository root. Put the input files in data/ (see README),
# or override any path with the environment variables shown.

OUT_DIR   <- Sys.getenv("HYPOXIA_OUT_DIR", "results")

# Input files (download them into data/; see README)
DATA_FILE <- Sys.getenv("HYPOXIA_DATA_FILE", "data/Water quality.csv")   # Sinha et al. (2026), Mendeley Data
WATERWAYS <- Sys.getenv("HYPOXIA_WATERWAYS", "data/waterways.geojson")   # OpenStreetMap waterways for Bangladesh (HOT export)
RASTER    <- Sys.getenv("HYPOXIA_RASTER", "data/pollution_int.tif")   # Kaggle Bangladesh-River-Pollution-GIS raster
S2_FILE   <- file.path(OUT_DIR, "s2_indices_per_sample.csv")   # exported from Google Earth Engine (gee/)

# Model settings (these are the values used in the paper)
HYPOXIA_MGL    <- 2      # DO below this (mg/L) counts as hypoxic
RANGE_PRIOR_KM <- 15     # PC prior: P(spatial range < 15 km) = 0.5
SIGMA_PRIOR    <- c(3, 0.01)   # PC prior: P(spatial sd > 3) = 0.01
YEAR_PRIOR     <- c(3, 0.01)   # PC prior: P(year sd > 3) = 0.01

# Slow optional steps (leave-one-station-out cross-validation)
RUN_LOSO <- FALSE
