# =====================================================================
# 05_year_2020_checks.R - what does the data say about the 2020 dip?
# =====================================================================
# 1. Was the sampling different in 2020 (stations, months)?
# 2. Does the rainfall recorded with each sample explain it?  (INLA, with and without rain)
# 3. Did other variables (pH, phosphate, turbidity, ...) shift at the same stations?
# 4. Is the 2020 DO distribution odd?
# Result in the paper: the dip is not explained by station mix or recorded rain; several variables shift together,
# so it is reported as a limitation.
source("R/00_setup.R")
WET <- c("Monsoon", "Post-monsoon")
wet <- d[d$season %in% WET, ]

## ---- 1. Sampling by year and season -----------------------------------------------
cat("\nSamples by year and season:\n"); print(table(d$year, d$season))
tab <- do.call(rbind, lapply(split(wet, list(wet$year, wet$season), drop = TRUE), function(g)
  data.frame(year = g$year[1], season = g$season[1], samples = nrow(g), stations = length(unique(g$station)),
             months = paste(sort(unique(g$month)), collapse = ","), hypoxic_share = round(mean(g$y), 2),
             mean_do = round(mean(g$do, na.rm = TRUE), 2), mean_rain = round(mean(g$rain, na.rm = TRUE), 1))))
rownames(tab) <- NULL; print(tab[order(tab$season, tab$year), ])
write.csv(tab, file.path(OUT_DIR, "year2020_wet_season_by_year.csv"), row.names = FALSE)

# Same stations: wet-season hypoxic share in 2020 and in other years
paired <- function(v, data) {
  m <- aggregate(data[[v]], list(station = data$station, is2020 = data$year == 2020), mean, na.rm = TRUE)
  w <- reshape(m, idvar = "station", timevar = "is2020", direction = "wide"); w <- w[complete.cases(w), ]
  data.frame(variable = v, stations = nrow(w), mean_2020 = round(mean(w$x.TRUE), 2), mean_other = round(mean(w$x.FALSE), 2),
             stations_higher_in_2020 = sum(w$x.TRUE > w$x.FALSE), stations_lower_in_2020 = sum(w$x.TRUE < w$x.FALSE),
             p_paired = signif(suppressWarnings(wilcox.test(w$x.TRUE, w$x.FALSE, paired = TRUE, exact = FALSE)$p.value), 2))
}
cat("\nHypoxic share at the same stations, wet seasons:\n"); print(paired("y", wet), row.names = FALSE)

## ---- 2. Rain anomaly and the year effect (INLA) ----------------------------------------------
d$rain_anom <- d$rain - ave(d$rain, d$season, FUN = function(x) mean(x, na.rm = TRUE))
ok <- !is.na(d$rain_anom)
u <- d[ok, ]; u$z_rain <- as.numeric(scale(u$rain_anom))
A_u <- A_obs[ok, , drop = FALSE]; cov_u <- cbind(cov_obs[ok, , drop = FALSE], z_rain = u$z_rain)
form_rain <- as.formula(paste("y ~ 0 + Intercept + s_pre + s_mon + s_post +", YEAR_TERM, "+ z_rain + f(s, model = spde)"))
f0 <- run_inla(form_B,    stk_est(u$y, A_u, cov_u))
f1 <- run_inla(form_rain, stk_est(u$y, A_u, cov_u))
cat("\nWAIC without rain =", round(f0$waic$waic, 1), "| with rain anomaly =", round(f1$waic$waic, 1), "\n")
or <- exp(as.numeric(unlist(f1$summary.fixed["z_rain", c("mean", "0.025quant", "0.975quant")])))
cat("Odds ratio of hypoxia per +1 SD rain anomaly:", round(or[1], 2), "(95%", round(or[2], 2), "-", round(or[3], 2), ")\n")
yr <- function(f) data.frame(year = sort(unique(u$year)), mean = round(f$summary.random$year_id$mean, 2),
                             lo = round(f$summary.random$year_id$`0.025quant`, 2), hi = round(f$summary.random$year_id$`0.975quant`, 2))
cat("\nYear effects without rain:\n"); print(yr(f0)); cat("\nYear effects with rain anomaly:\n"); print(yr(f1))

## ---- 3. Other variables, 2020 vs other years, same stations ------------------------------------------
VARS <- c("temp", "ph", "do", "cod", "no3", "turb", "po4")
cmp <- do.call(rbind, lapply(VARS, paired, data = wet))
cat("\nWet seasons, station means, 2020 vs other years:\n"); print(cmp, row.names = FALSE)
write.csv(cmp, file.path(OUT_DIR, "year2020_other_variables_wet.csv"), row.names = FALSE)
cat("Note: COD is probably filled in (mean exactly 70 in 2020), so do not interpret it.\n")
cat("A dry-season comparison is confounded: 2020 has only winter samples, other years also pre-monsoon.\n")

## ---- 4. DO distribution by year ---------------------------------------------------------------------------
print(do.call(rbind, lapply(split(d, d$year), function(g) { x <- g$do
  data.frame(year = g$year[1], n = length(x), min = min(x), q1 = quantile(x, .25), median = median(x),
             q3 = quantile(x, .75), max = max(x)) })), row.names = FALSE)
