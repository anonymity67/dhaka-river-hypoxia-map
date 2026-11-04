# =====================================================================
# Configuration: the only file you need to edit
# =====================================================================
# Run every script from the repository root (the folder that holds this file).
# Input and output locations are set below (absolute paths for this computer).

OUT_DIR   <- Sys.getenv("HYPOXIA_OUT_DIR", "C:/Antar/Munim Project/River Pollution Research/results")

# Input files (absolute paths on the author's computer; change them if the files move)
DATA_FILE <- Sys.getenv("HYPOXIA_DATA_FILE", "C:/Antar/Munim Project/Fresh Water/ClimatePlus-Collaborator-Kit/ClimatePlus-Collaborator-Kit/research/water/Water quality.csv")   # Sinha et al. (2026), Mendeley Data
WATERWAYS <- "C:/Antar/Munim Project/Fresh Water/hotosm_bgd_waterways_osm_geojson/waterways.geojson"   # OpenStreetMap waterways for Bangladesh (HOT export)
RASTER    <- "C:/Antar/Munim Project/River Pollution Research/archive/Pollution_intensity_Map/Pollution_intensity_Map/pollution_int.tif"   # Kaggle Bangladesh-River-Pollution-GIS raster
S2_FILE   <- file.path(OUT_DIR, "s2_indices_per_sample.csv")   # exported from Google Earth Engine (gee/)

# Model settings (these are the values used in the paper)
HYPOXIA_MGL    <- 2      # DO below this (mg/L) counts as hypoxic
RANGE_PRIOR_KM <- 15     # PC prior: P(spatial range < 15 km) = 0.5
SIGMA_PRIOR    <- c(3, 0.01)   # PC prior: P(spatial sd > 3) = 0.01
YEAR_PRIOR     <- c(3, 0.01)   # PC prior: P(year sd > 3) = 0.01

# Slow optional steps (leave-one-station-out cross-validation)
RUN_LOSO <- FALSE
