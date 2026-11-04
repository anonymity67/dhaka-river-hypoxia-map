# =====================================================================
# 07_pollution_raster_check.R - does the Kaggle pollution raster track hypoxia risk?
# =====================================================================
# Needs data/pollution_int.tif and results/river_objects.rds (from script 02).
# PART A  station level (60 stations): raster value vs observed hypoxia share, median DO and mean WQI
# PART B  river network: raster value vs the model's modelled risk, by season
# The raster has no date, unit or documented method: this checks agreement, it does not show the raster measures pollution.
source("R/00_setup.R")
CENTRE  <- c(lon = 90.40, lat = 23.78)    # central Dhaka, to test whether the raster just follows distance from the centre
WINDOWS <- c(1, 5, 11)                    # raster averaging window in pixels; 11 is about 2 km each side
MAIN_W  <- 11
N_PERM  <- 2000
stopifnot(file.exists(RASTER), file.exists(file.path(OUT_DIR, "river_objects.rds")))
obj <- readRDS(file.path(OUT_DIR, "river_objects.rds")); rp <- obj$rp; risk_B <- obj$risk
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
library(terra)

## ---- 1. READ THE RASTER ------------------------------------------------
r <- rast(RASTER)[[1]]
if (!is.lonlat(r)) r <- project(r, "EPSG:4326")
cat("Raster size:", nrow(r), "x", ncol(r), "| pixel about", round(res(r)[1] * 111.32, 2), "km\n")
cat("Value range:", paste(round(unlist(global(r, "range", na.rm = TRUE)), 2), collapse = " to "),
    "| share of pixels with data:", round(mean(!is.na(values(r))), 3), "\n")

# One layer per averaging window (a moving-average of the raster)
lay <- lapply(WINDOWS, function(w) if (w == 1) r else focal(r, w = w, fun = "mean", na.rm = TRUE))
names(lay) <- paste0("w", WINDOWS)
get_vals <- function(lon, lat) {
  v <- vect(cbind(lon, lat), crs = "EPSG:4326")
  out <- sapply(lay, function(l) extract(l, v)[, 2]); as.data.frame(out)
}
dist_centre_km <- function(lon, lat)
  sqrt(((lon - CENTRE["lon"]) * cos(CENTRE["lat"] * pi / 180))^2 + (lat - CENTRE["lat"])^2) * 111.32
rho  <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
prho <- function(x, y, z) {               # Spearman correlation of x and y after removing z
  ok <- complete.cases(x, y, z)
  cor(resid(lm(rank(x[ok]) ~ rank(z[ok]))), resid(lm(rank(y[ok]) ~ rank(z[ok]))))
}


## ---- 2. PART A: STATION LEVEL ----------------------------------------------
first <- !duplicated(d$station)
stn <- data.frame(station = d$station[first], river = d$river[first],
                  lon = d$lon[first], lat = d$lat[first],
                  p_hyp  = tapply(d$y,   d$station, mean)[d$station[first]],
                  do_med = tapply(d$do,  d$station, median, na.rm = TRUE)[d$station[first]],
                  wqi_mean = tapply(d$wqi, d$station, mean, na.rm = TRUE)[d$station[first]],
                  row.names = NULL)
stn <- cbind(stn, get_vals(stn$lon, stn$lat))
stn$dist_c <- dist_centre_km(stn$lon, stn$lat)
cat("\nStations with a raster value (window", MAIN_W, "):", sum(!is.na(stn[[paste0("w", MAIN_W)]])), "of", nrow(stn), "\n")

resA <- do.call(rbind, lapply(WINDOWS, function(w) {
  x <- stn[[paste0("w", w)]]
  data.frame(window_px = w,
             rho_hypoxia = round(rho(x, stn$p_hyp), 2),
             rho_DO_median = round(rho(x, stn$do_med), 2),
             rho_WQI_mean = round(rho(x, stn$wqi_mean), 2),
             partial_hypoxia_given_dist = round(prho(x, stn$p_hyp, stn$dist_c), 2))
}))
cat("\nPART A: Spearman rho between raster value and station measurements\n")
cat("(hypoxia share: positive = good; median DO: negative = good; WQI runs backwards: positive = good)\n")
print(resA)

xm <- stn[[paste0("w", MAIN_W)]]
# leave one river out
lor <- sapply(unique(stn$river), function(rv) rho(xm[stn$river != rv], stn$p_hyp[stn$river != rv]))
cat("\nLeave-one-river-out rho with hypoxia share: min", round(min(lor), 2), "| max", round(max(lor), 2), "\n")
# within-river permutation: shuffle raster values among stations of the same river
set.seed(1)
obs <- rho(xm, stn$p_hyp)
perm <- replicate(N_PERM, rho(ave(xm, stn$river, FUN = function(v) v[sample.int(length(v))]), stn$p_hyp))
cat("Within-river permutation: observed rho", round(obs, 2), "| p =", round(mean(abs(perm) >= abs(obs)), 3),
    "(share of shuffles at least as strong)\n")


## ---- 3. PART B: ALONG THE RIVER NETWORK ---------------------------------------
ll   <- st_coordinates(st_transform(rp, 4326))
rv   <- cbind(data.frame(lon = ll[, 1], lat = ll[, 2]), get_vals(ll[, 1], ll[, 2]))
rv$dist_c <- dist_centre_km(rv$lon, rv$lat)
n_r  <- nrow(rv)
cat("\nPART B: river points with a raster value (window", MAIN_W, "):",
    sum(!is.na(rv[[paste0("w", MAIN_W)]])), "of", n_r, "\n")
xr <- rv[[paste0("w", MAIN_W)]]
resB <- do.call(rbind, lapply(seasons, function(sn) {
  p <- risk_B$p_mean[risk_B$season == sn]
  data.frame(season = sn, rho_with_model_risk = round(rho(xr, p), 2),
             partial_given_dist_centre = round(prho(xr, p, rv$dist_c), 2))
}))
print(resB)
cat("Spatial correlation between neighbouring points inflates significance, so use these as descriptions, not as tests.\n")


## ---- 4. SAVE ---------------------------------------------------------------------
write.csv(stn,  file.path(OUT_DIR, "raster_check_stations.csv"), row.names = FALSE)
write.csv(resA, file.path(OUT_DIR, "raster_check_stations_summary.csv"), row.names = FALSE)
write.csv(resB, file.path(OUT_DIR, "raster_check_river_summary.csv"), row.names = FALSE)

png(file.path(OUT_DIR, "raster_vs_stations.png"), width = 1000, height = 800, res = 130)
cols <- as.integer(factor(stn$river))
plot(xm, stn$p_hyp, pch = 19, col = cols, xlab = paste0("Raster value (", MAIN_W, "-pixel mean)"),
     ylab = "Observed hypoxia share (DO below 2 mg/L)",
     main = paste0("Raster vs station hypoxia, rho = ", round(obs, 2)))
legend("bottomright", legend = levels(factor(stn$river)), col = seq_along(levels(factor(stn$river))),
       pch = 19, cex = 0.6, bg = "white")
dev.off()
cat("Saved to", OUT_DIR, "\n")
