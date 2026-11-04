# =====================================================================
# EXTRA (not in the paper): do pollution-source locations from OpenStreetMap explain the hotspots?
# Result: no association (null). Kept here for transparency.
# =====================================================================
# Self-contained. Downloads industrial areas, wastewater plants and outfalls
# for the study area from OpenStreetMap (Overpass), measures how close each of
# the 60 stations is to them, and tests whether closer stations are more often
# hypoxic.
#
# Run script 01 first. Needs: stations_for_gee.csv in OUT_DIR (it has station, river, lon, lat, p_hyp).
# Needs internet. The download can take a few minutes.
#
# DECIDED BEFORE LOOKING AT THE DATA (so we do not go fishing):
#   Main test  : industrial area within 2 km of the station vs hypoxia share
#   Extras     : distance to nearest industrial feature, wastewater plant, outfall
#   All use Spearman correlation, a within-river permutation test, and the
#   correlation after removing distance from central Dhaka.
# Only 60 stations, so treat a result as suggestive, not proof.
#
# HOW TO RUN: section by section, Ctrl+Enter. Send me what sections 1-3 print.
# =====================================================================

## ---- 0. SETTINGS -------------------------------------------------------
source("config.R"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
BBOX      <- c(90.00, 23.45, 90.75, 24.25)   # xmin, ymin, xmax, ymax around the 60 stations
BUFFER_KM <- 2
N_PERM    <- 5000
CENTRE    <- c(90.4074, 23.7104)             # central Dhaka (lon, lat), used as a control

for (p in c("osmdata", "sf")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(osmdata)
library(sf)
sf_use_s2(FALSE)

## ---- 1. DOWNLOAD THE SOURCE LOCATIONS ------------------------------------
# Each query returns points and/or polygons. We keep the geometry only.
get_osm <- function(key, value) {
  q <- try(opq(bbox = BBOX, timeout = 300) |> add_osm_feature(key = key, value = value) |>
             osmdata_sf(), silent = TRUE)
  if (inherits(q, "try-error")) { cat("  query failed:", key, "=", value, "\n"); return(NULL) }
  g <- list()
  for (nm in c("osm_points", "osm_polygons", "osm_multipolygons")) {
    x <- q[[nm]]
    if (!is.null(x) && nrow(x) > 0) g[[nm]] <- st_geometry(st_transform(x, 32646))
  }
  if (length(g) == 0) return(NULL)
  do.call(c, g)
}
classes <- list(
  industrial = list(c("landuse", "industrial"), c("building", "industrial"), c("man_made", "works")),
  wwtp       = list(c("man_made", "wastewater_plant")),
  outfall    = list(c("man_made", "outfall"), c("waterway", "outfall"))
)
src <- lapply(names(classes), function(cl) {
  cat("Downloading", cl, "...\n")
  parts <- lapply(classes[[cl]], function(kv) get_osm(kv[1], kv[2]))
  parts <- parts[!sapply(parts, is.null)]
  if (length(parts) == 0) return(NULL)
  do.call(c, parts)
})
names(src) <- names(classes)
cat("\nFeatures found (points + polygons):\n")
print(sapply(src, function(x) if (is.null(x)) 0L else length(x)))
saveRDS(src, file.path(OUT_DIR, "osm_sources.rds"))
# If a class has fewer than about 5 features, its test is not meaningful and is skipped below.

## ---- 2. DISTANCES AND AREAS AT EACH STATION --------------------------------
st <- read.csv(file.path(OUT_DIR, "stations_for_gee.csv"), stringsAsFactors = FALSE)
stp <- st_transform(st_as_sf(st, coords = c("lon", "lat"), crs = 4326), 32646)
ctr <- st_transform(st_sfc(st_point(CENTRE), crs = 4326), 32646)
st$d_centre_km <- as.numeric(st_distance(stp, ctr)) / 1000

nearest_km <- function(g) {
  if (is.null(g) || length(g) < 5) return(rep(NA_real_, nrow(stp)))
  apply(st_distance(stp, g), 1, min) / 1000
}
st$d_industrial_km <- nearest_km(src$industrial)
st$d_wwtp_km       <- nearest_km(src$wwtp)
st$d_outfall_km    <- nearest_km(src$outfall)

# Industrial area (km2) inside the buffer around each station
st$ind_area_km2 <- NA_real_
if (!is.null(src$industrial) && length(src$industrial) >= 5) {
  poly <- src$industrial[st_geometry_type(src$industrial) %in% c("POLYGON", "MULTIPOLYGON")]
  poly <- st_union(st_make_valid(poly))
  buf  <- st_buffer(stp, BUFFER_KM * 1000)
  st$ind_area_km2 <- sapply(seq_len(nrow(buf)), function(i) {
    x <- suppressWarnings(st_intersection(st_geometry(buf)[i], poly))
    if (length(x) == 0) 0 else as.numeric(sum(st_area(x))) / 1e6
  })
}
write.csv(st, file.path(OUT_DIR, "osm_source_distances.csv"), row.names = FALSE)
cat("\nSummary of station measures:\n"); print(summary(st[, c("d_industrial_km", "d_wwtp_km", "d_outfall_km", "ind_area_km2")]))

## ---- 3. TESTS ----------------------------------------------------------------
spear <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
partial_spear <- function(x, y, z) {                      # remove the effect of z by ranks
  ok <- complete.cases(x, y, z)
  rx <- resid(lm(rank(x[ok]) ~ rank(z[ok]))); ry <- resid(lm(rank(y[ok]) ~ rank(z[ok])))
  cor(rx, ry)
}
perm_p <- function(x, y, river) {                          # shuffle y within each river
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; river <- river[ok]
  obs <- spear(x, y); set.seed(1)
  sim <- replicate(N_PERM, {
    y2 <- y; for (r in unique(river)) { i <- which(river == r); if (length(i) > 1) y2[i] <- sample(y[i]) }
    spear(x, y2) })
  mean(abs(sim) >= abs(obs))
}
tests <- data.frame(
  measure = c("ind_area_km2 (MAIN)", "d_industrial_km", "d_wwtp_km", "d_outfall_km"),
  col     = c("ind_area_km2", "d_industrial_km", "d_wwtp_km", "d_outfall_km"),
  expected_sign = c("+", "-", "-", "-"), stringsAsFactors = FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(tests)), function(i) {
  x <- st[[tests$col[i]]]
  if (all(is.na(x))) return(data.frame(measure = tests$measure[i], n = 0, rho = NA, partial_rho_centre = NA,
                                       p_within_river = NA, note = "too few features, skipped"))
  data.frame(measure = tests$measure[i], n = sum(!is.na(x)),
             rho = round(spear(x, st$p_hyp), 2),
             partial_rho_centre = round(partial_spear(x, st$p_hyp, st$d_centre_km), 2),
             p_within_river = round(perm_p(x, st$p_hyp, st$river), 3), note = "")
}))
cat("\nHypoxia share vs source measures (expected sign: + for area, - for distances):\n"); print(res)
cat("\nFor reference, distance from central Dhaka vs hypoxia share: rho =",
    round(spear(st$d_centre_km, st$p_hyp), 2), "\n")

# Hotspot view: the 15 stations with the highest hypoxia share vs the rest
st$hot <- rank(-st$p_hyp, ties.method = "first") <= 15
for (cl in c("d_industrial_km", "d_wwtp_km", "d_outfall_km")) if (!all(is.na(st[[cl]])))
  cat(cl, ": median hotspot", round(median(st[[cl]][st$hot], na.rm = TRUE), 2), "km, others",
      round(median(st[[cl]][!st$hot], na.rm = TRUE), 2), "km\n")
write.csv(res, file.path(OUT_DIR, "osm_source_tests.csv"), row.names = FALSE)
cat("\nDone. Send me everything printed above.\n")
