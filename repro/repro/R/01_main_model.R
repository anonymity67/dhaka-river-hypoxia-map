# =====================================================================
# 01_main_model.R - descriptives, the main INLA-SPDE model, validation
# =====================================================================
# Model: hypoxia (DO < 2 mg/L) ~ Bernoulli, logit link
#        logit(p) = intercept + season (Winter = reference) + year (iid) + Matern spatial field (SPDE)
# Outputs (results/): descriptive_*.csv, odds_ratios_season.csv, model_hyperparameters.csv,
#                     model_waic.csv, [loso_auc.csv], samples_for_gee.csv, stations_for_gee.csv
source("R/00_setup.R")

## ---- 1. Descriptives ------------------------------------------------------------
cat("\nShare of samples below 5 mg/L:", round(mean(d$do < 5), 3), "| below 2 mg/L:", round(mean(d$y), 3), "\n")
by_season <- do.call(rbind, lapply(seasons, function(s) { g <- d[d$season == s, ]
  data.frame(season = s, samples = nrow(g), lt5 = round(mean(g$do < 5), 2), lt2 = round(mean(g$y), 2)) }))
print(by_season)
by_river <- do.call(rbind, lapply(sort(unique(d$river)), function(r) { g <- d[d$river == r, ]
  data.frame(river = r, samples = nrow(g), stations = length(unique(g$station)),
             lt5 = round(mean(g$do < 5), 2), lt2 = round(mean(g$y), 2), median_do = round(median(g$do), 2)) }))
print(by_river[order(-by_river$lt2), ], row.names = FALSE)
write.csv(by_season, file.path(OUT_DIR, "descriptive_season.csv"), row.names = FALSE)
write.csv(by_river,  file.path(OUT_DIR, "descriptive_river.csv"),  row.names = FALSE)

## ---- 2. Files for Google Earth Engine (script 06) -------------------------------------
txt <- substr(trimws(d$date), 1, 10); dt <- as.Date(txt, format = "%Y-%m-%d")
for (f in c("%d/%m/%Y", "%m/%d/%Y", "%d-%m-%Y", "%Y/%m/%d")) {
  bad <- is.na(dt); if (!any(bad)) break; dt[bad] <- as.Date(txt[bad], format = f) }
stopifnot(!anyNA(dt))
write.csv(data.frame(sample_id = d$sample_id, station = d$station, river = d$river, lon = d$lon, lat = d$lat,
                     date = format(dt, "%Y-%m-%d"), season = d$season, do = d$do, hyp = d$y),
          file.path(OUT_DIR, "samples_for_gee.csv"), row.names = FALSE)
st <- aggregate(cbind(lon, lat, p_hyp = y) ~ station + river, data = d, FUN = mean)
write.csv(st, file.path(OUT_DIR, "stations_for_gee.csv"), row.names = FALSE)

## ---- 3. Fit model A (no map) and model B (with map) ---------------------------------------
stk <- stk_est(d$y, A_obs, cov_obs)
fit_A <- run_inla(form_A, stk)
fit_B <- run_inla(form_B, stk)
waic <- data.frame(model = c("A: season + year", "B: season + year + spatial field"),
                   WAIC = round(c(fit_A$waic$waic, fit_B$waic$waic), 1))
cat("\nWAIC (lower is better):\n"); print(waic)          # paper: 1450.9 and 1060.4
write.csv(waic, file.path(OUT_DIR, "model_waic.csv"), row.names = FALSE)

orr <- odds_ratios(fit_B)
cat("\nSeason odds ratios against winter (model B):\n"); print(orr)   # paper: 0.47, 0.006, 0.029
write.csv(orr, file.path(OUT_DIR, "odds_ratios_season.csv"), row.names = FALSE)

res_f   <- inla.spde2.result(fit_B, "s", spde, do.transf = TRUE)
q3      <- c(0.025, 0.5, 0.975)
range_q <- inla.qmarginal(q3, res_f$marginals.range.nominal[[1]])
sd_q    <- sqrt(inla.qmarginal(q3, res_f$marginals.variance.nominal[[1]]))
year_sd <- inla.qmarginal(q3, inla.tmarginal(function(x) 1 / sqrt(x), fit_B$marginals.hyperpar[["Precision for year_id"]]))
hp <- data.frame(quantity = c("spatial range (km)", "spatial sd (logit)", "year sd (logit)"),
                 lower95 = round(c(range_q[1], sd_q[1], year_sd[1]), 2),
                 median  = round(c(range_q[2], sd_q[2], year_sd[2]), 2),
                 upper95 = round(c(range_q[3], sd_q[3], year_sd[3]), 2))
cat("\nSpatial range, spatial sd, year sd:\n"); print(hp)
write.csv(hp, file.path(OUT_DIR, "model_hyperparameters.csv"), row.names = FALSE)
writeLines(capture.output(summary(fit_B)), file.path(OUT_DIR, "model_B_summary.txt"))

## ---- 4. Leave-one-station-out AUC (slow; set RUN_LOSO <- TRUE in config.R) ------------------
if (RUN_LOSO) {
  cv <- data.frame(model = c("A: no map", "B: with map"),
                   AUC = round(c(loso_auc(form_A), loso_auc(form_B)), 3))
  cat("\nLeave-one-station-out AUC:\n"); print(cv)        # paper: 0.813 and 0.908
  write.csv(cv, file.path(OUT_DIR, "loso_auc.csv"), row.names = FALSE)
} else cat("\n(Leave-one-station-out AUC skipped: set RUN_LOSO <- TRUE in config.R.)\n")
