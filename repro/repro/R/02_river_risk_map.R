# =====================================================================
# 02_river_risk_map.R - predict risk along the rivers and build the website data
# =====================================================================
# Re-fits model B with prediction points every 500 m along OpenStreetMap
# river/canal/stream lines within 5 km of a station, for each season (typical
# year: the year effect is left NA). Writes tables, GeoJSON for the map website
# and an .rds used by script 07.
source("R/00_setup.R")
BBOX_BUFFER_KM <- 6                          # read waterways inside the station box + this margin
MAX_DIST_KM    <- 5                          # keep river points within this distance of a station
SPACING_M      <- 500                        # one prediction point every 500 m of river
KEEP_TYPES     <- c("river", "canal", "stream")
BALU_SYSTEM    <- c("Balu", "Begunbari Canal", "Pagar Canal", "Bashundhara Canal")  # grouped in tables

## ---- 1. River lines near the stations --------------------------------------------
stopifnot(file.exists(WATERWAYS))
bb  <- st_bbox(st_transform(st_buffer(st_union(pts), BBOX_BUFFER_KM * 1000), 4326))
ww  <- st_read(WATERWAYS, wkt_filter = st_as_text(st_as_sfc(bb)), quiet = TRUE)   # reads only the box (large file)
ww  <- ww[st_geometry_type(ww) %in% c("LINESTRING", "MULTILINESTRING"), ]
ww  <- ww[ww$waterway %in% KEEP_TYPES, ]
ww  <- st_transform(ww, 32646)
ww  <- suppressWarnings(st_cast(st_cast(ww, "MULTILINESTRING"), "LINESTRING"))
pick <- function(col) if (col %in% names(ww)) as.character(ww[[col]]) else rep(NA_character_, nrow(ww))
lab  <- pick("name_en")
for (alt in c("name_latin", "name")) { a <- pick(alt); lab <- ifelse(is.na(lab) | lab == "", a, lab) }
ww$osm_name <- ifelse(is.na(lab) | lab == "", NA_character_, lab)
ww <- ww[lengths(st_is_within_distance(ww, pts, dist = MAX_DIST_KM * 1000)) > 0, ]
ww$line_id <- seq_len(nrow(ww))
cat("River lines kept:", nrow(ww), "\n")

## ---- 2. Prediction points every 500 m ------------------------------------------------
smp <- st_line_sample(ww, density = 1 / SPACING_M)
rp  <- st_sf(line_id = ww$line_id, osm_name = ww$osm_name, waterway = ww$waterway, geometry = smp)
rp  <- suppressWarnings(st_cast(rp[!st_is_empty(rp), ], "POINT"))
nn  <- st_nearest_feature(rp, pts)                                   # nearest sampling point
rp$river   <- d$river[nn]                                            # river name = river of the nearest station
rp$dist_km <- as.numeric(st_distance(rp, pts[nn, ], by_element = TRUE)) / 1000
rp <- rp[rp$dist_km <= MAX_DIST_KM, ]
cat("Prediction points along rivers:", nrow(rp), "(x 4 seasons)\n")

## ---- 3. Fit model B with the river points as prediction sites ---------------------------------
xy_r <- st_coordinates(rp) / 1000; ll_r <- st_coordinates(st_transform(rp, 4326)); n_r <- nrow(xy_r)
cov_r <- do.call(rbind, lapply(seasons, function(sn) data.frame(
  Intercept = 1, s_pre = as.integer(sn == "Pre-monsoon"), s_mon = as.integer(sn == "Monsoon"),
  s_post = as.integer(sn == "Post-monsoon"), year_id = NA_integer_)))[rep(seq_along(seasons), each = n_r), ]
A_r  <- fmesher::fm_basis(mesh, loc = xy_r)
stk_r <- inla.stack(tag = "pred", data = list(y = rep(NA_real_, n_r * 4), Ntrials = rep(1, n_r * 4)),
                    A = list(do.call(rbind, rep(list(A_r), 4)), 1),
                    effects = list(s = seq_len(spde$n.spde), cov_r))
stk_all <- inla.stack(stk_est(d$y, A_obs, cov_obs), stk_r)
fit_R <- run_inla(form_B, stk_all, pred = TRUE)

## ---- 4. Probabilities (logit-scale mean converted with the usual normal correction) ----------------
ip <- inla.stack.index(stk_all, "pred")$data
lp <- fit_R$summary.linear.predictor[ip, ]
risk <- data.frame(season = rep(seasons, each = n_r), line_id = rep(rp$line_id, 4), river = rep(rp$river, 4),
                   lon = rep(ll_r[, 1], 4), lat = rep(ll_r[, 2], 4), dist_station_km = rep(rp$dist_km, 4),
                   p_mean = plogis(lp$mean / sqrt(1 + pi * lp$sd^2 / 8)),
                   p_lo = plogis(lp$`0.025quant`), p_hi = plogis(lp$`0.975quant`),
                   p_gt_50 = 1 - pnorm(0, mean = lp$mean, sd = lp$sd))
risk$class <- cut(risk$p_mean, c(-Inf, 0.25, 0.5, 0.75, Inf), labels = c("Low", "Moderate", "High", "Very high"))
risk$river_group <- ifelse(risk$river %in% BALU_SYSTEM, "Balu system", risk$river)
write.csv(risk, file.path(OUT_DIR, "risk_river_points.csv"), row.names = FALSE)

# Calibration check: fitted and observed hypoxia share should agree overall (paper: 0.567 and 0.567)
ie <- inla.stack.index(stk_all, "est")$data
cat("\nCalibration: observed", round(mean(d$y), 3), "| fitted", round(mean(fit_R$summary.fitted.values$mean[ie]), 3), "\n")

## ---- 5. River tables ------------------------------------------------------------------------------
tab <- aggregate(p_mean ~ river_group + season, data = risk, FUN = mean)
wide <- tapply(tab$p_mean, list(tab$river_group, tab$season), mean)[, seasons]
wide <- round(wide[order(-wide[, "Winter"]), ], 2)
st <- aggregate(cbind(lon, lat) ~ station + river, data = d, FUN = mean)
st$river_group <- ifelse(st$river %in% BALU_SYSTEM, "Balu system", st$river)
out <- data.frame(group = rownames(wide), stations = as.integer(table(st$river_group)[rownames(wide)]),
                  points = as.integer(table(risk$river_group[risk$season == "Winter"])[rownames(wide)]), wide,
                  row.names = NULL, check.names = FALSE)
cat("\nMean modelled probability of hypoxia by river group:\n"); print(out)
write.csv(out, file.path(OUT_DIR, "risk_by_river_group.csv"), row.names = FALSE)

## ---- 6. Data for the map website (site/data/) ---------------------------------------------------------
dir.create(file.path(OUT_DIR, "site"), showWarnings = FALSE)
mode_name <- function(v) names(sort(table(v), decreasing = TRUE))[1]
ln <- ww[ww$line_id %in% risk$line_id, c("line_id", "osm_name", "waterway")]
ln$river  <- as.character(tapply(rp$river, rp$line_id, mode_name)[as.character(ln$line_id)])
ln$g      <- ifelse(ln$river %in% BALU_SYSTEM, "Balu system", ln$river)
ln$k      <- as.integer(table(factor(rp$line_id, levels = ln$line_id)))
key <- c(Winter = "w", `Pre-monsoon` = "p", Monsoon = "m", `Post-monsoon` = "o")
for (sn in seasons) ln[[key[sn]]] <- round(as.numeric(tapply(risk$p_mean[risk$season == sn], risk$line_id[risk$season == sn],
                                                              mean)[as.character(ln$line_id)]), 2)
names(ln)[names(ln) == "osm_name"] <- "n"; names(ln)[names(ln) == "waterway"] <- "t"
ln <- ln[, c("n", "g", "w", "p", "m", "o", "t", "k")]
f_rivers <- file.path(OUT_DIR, "site", "rivers.geojson")
st_write(st_transform(ln, 4326), f_rivers, delete_dsn = TRUE, quiet = TRUE, layer_options = "COORDINATE_PRECISION=5")

if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
hs <- aggregate(y ~ station, data = d, FUN = mean)
stn <- data.frame(s = st$station, g = st$river_group, x = round(st$lon, 4), y = round(st$lat, 4),
                  h = round(hs$y[match(st$station, hs$station)], 2))
jsonlite::write_json(stn, file.path(OUT_DIR, "site", "stations.json"), digits = NA, auto_unbox = TRUE)

saveRDS(list(ww = ww, rp = rp, risk = risk), file.path(OUT_DIR, "river_objects.rds"))
cat("\nSaved to", OUT_DIR, "(site files in", file.path(OUT_DIR, "site"), ")\n")
