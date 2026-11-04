# =====================================================================
# run_all.R - run the whole analysis in order (from the repository root)
# =====================================================================
# Each script can also be run alone: source("R/01_main_model.R"), etc.
# Slow optional steps (leave-one-station-out AUC) are switched on in config.R (RUN_LOSO <- TRUE).
# Script 06 needs results/s2_indices_per_sample.csv from Google Earth Engine; script 07 needs data/pollution_int.tif.
steps <- c("R/01_main_model.R", "R/02_river_risk_map.R", "R/03_robustness.R", "R/04_without_2020.R", "R/05_year_2020_checks.R", "R/06_mesh_figure.R")
optional <- c("R/06_sentinel2.R", "R/07_pollution_raster_check.R")
source("config.R"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (f in steps) { cat("\n==========", f, "==========\n"); source(f, local = FALSE) }
source("config.R")
for (f in optional) {
  need <- if (grepl("sentinel2", f)) S2_FILE else RASTER
  if (file.exists(need)) { cat("\n==========", f, "==========\n"); source(f, local = FALSE) }
  else cat("\nSkipping", f, "(missing", need, ")\n")
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\nAll done. Outputs are in", OUT_DIR, "\n")
