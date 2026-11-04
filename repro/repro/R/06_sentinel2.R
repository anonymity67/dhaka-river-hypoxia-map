# =====================================================================
# 06_sentinel2.R - do Sentinel-2 water indices improve the model? (exploratory)
# =====================================================================
# Needs results/s2_indices_per_sample.csv, exported from Google Earth Engine with the scripts in gee/
# (see README). Fits four models on the samples that have usable satellite values, so they are
# compared like for like:
#   A0 season + year                 AX = A0 + 4 satellite indices (MNDWI, NDCI, NDTI, red)
#   B0 season + year + spatial map   BX = B0 + 4 satellite indices
# and scores each by WAIC and (optionally) by leave-one-station-out AUC with a station bootstrap.
source("R/00_setup.R")
MIN_WATER_PX <- 5      # keep a sample only if at least this many clear water pixels
SR_ONLY      <- FALSE  # TRUE = keep only surface-reflectance images (drops 2017-2018 top-of-atmosphere ones); outputs get suffix "_sr"
RUN_CV       <- RUN_LOSO   # slow: about 10-20 minutes per run
N_BOOT       <- 500
stopifnot(file.exists(S2_FILE))

## ---- 1. ATTACH THE SATELLITE VALUES TO THE SAMPLES -----------------------
s2 <- read.csv(S2_FILE, stringsAsFactors = FALSE)
d$sample_id <- seq_len(nrow(d))                       # same order as the export
stopifnot(nrow(d) == nrow(s2), all(d$sample_id %in% s2$sample_id))
m <- match(d$sample_id, s2$sample_id)
stopifnot(all(abs(d$do - s2$do[m]) < 1e-6, na.rm = TRUE))   # rows line up: DO must match

for (v in c("mndwi", "ndci", "ndti", "red", "n_water_px", "source"))
  d[[v]] <- s2[[v]][m]
ok <- !is.na(d$n_water_px) & d$n_water_px >= MIN_WATER_PX & !is.na(d$mndwi) &
      !is.na(d$ndci) & !is.na(d$ndti) & !is.na(d$red)
cat("Image source before the SR filter:\n"); print(table(d$source[ok], useNA = "ifany"))
if (SR_ONLY) {
  is_sr <- grepl("SR", toupper(d$source)) & !grepl("TOA", toupper(d$source))
  ok <- ok & !is.na(is_sr) & is_sr
  cat("SR_ONLY = TRUE: kept", sum(ok), "samples with surface-reflectance images\n")
}
SFX <- if (SR_ONLY) "_sr" else ""
cat("Samples with usable satellite values:", sum(ok), "of", nrow(d), "\n")
cat("Stations:", length(unique(d$station[ok])), "| share hypoxic:", round(mean(d$y[ok]), 3), "\n")
print(table(d$season[ok]))
print(table(d$source[ok]))

u    <- d[ok, ]
xy_u <- xy[ok, , drop = FALSE]
z    <- function(x) as.numeric(scale(x))              # 1 unit = 1 standard deviation
u$z_mndwi <- z(u$mndwi); u$z_ndci <- z(u$ndci); u$z_ndti <- z(u$ndti); u$z_red <- z(u$red)


## ---- 2. MODELS ---------------------------------------------------------
spde_u <- spde   # same spatial model as the main analysis
A_u   <- fmesher::fm_basis(mesh, loc = xy_u)
cov_u <- data.frame(Intercept = 1, s_pre = u$s_pre, s_mon = u$s_mon, s_post = u$s_post,
                    year_id = u$year_id, z_mndwi = u$z_mndwi, z_ndci = u$z_ndci,
                    z_ndti = u$z_ndti, z_red = u$z_red)

YEAR <- "f(year_id, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(3, 0.01))))"
SAT  <- "z_mndwi + z_ndci + z_ndti + z_red"
SP   <- "f(s, model = spde_u)"
mk   <- function(sat, map) as.formula(paste(
  "y ~ 0 + Intercept + s_pre + s_mon + s_post +", YEAR,
  if (sat) paste("+", SAT) else "", if (map) paste("+", SP) else ""))

fit_u <- function(form, y_vec, ctrl = list(), waic = TRUE, pred = FALSE) {
  stk <- inla.stack(tag = "est", data = list(y = y_vec, Ntrials = rep(1, length(y_vec))),
                    A = list(A_u, 1), effects = list(s = seq_len(spde_u$n.spde), cov_u))
  dd  <- inla.stack.data(stk)
  fit <- inla(form, family = "binomial", Ntrials = dd$Ntrials, data = dd,
              control.predictor = list(A = inla.stack.A(stk), compute = pred),
              control.compute = list(waic = waic), control.inla = ctrl)
  list(fit = fit, stk = stk)
}

forms <- list(A0 = mk(FALSE, FALSE), AX = mk(TRUE, FALSE), B0 = mk(FALSE, TRUE), BX = mk(TRUE, TRUE))
fits  <- list()
for (nm in names(forms)) {
  cat("Fitting", nm, "...\n")
  fits[[nm]] <- fit_u(forms[[nm]], u$y)$fit
}


## ---- 3. RESULTS: WAIC AND THE SATELLITE EFFECTS ----------------------------
tab <- data.frame(model = names(fits), WAIC = sapply(fits, function(f) round(f$waic$waic, 1)),
                  row.names = NULL)
cat("\nWAIC on the same", nrow(u), "samples (lower is better):\n"); print(tab)

# Odds ratio of hypoxia for +1 standard deviation of each index
or_sat <- do.call(rbind, lapply(c("AX", "BX"), function(nm) {
  fx <- fits[[nm]]$summary.fixed
  fx <- fx[grepl("^z_", rownames(fx)), c("mean", "0.025quant", "0.975quant")]
  data.frame(model = nm, index = sub("^z_", "", rownames(fx)),
             odds_ratio = round(exp(fx$mean), 2), lower95 = round(exp(fx$`0.025quant`), 2),
             upper95 = round(exp(fx$`0.975quant`), 2), row.names = NULL)
}))
cat("\nOdds ratio of hypoxia per +1 SD of each index (a CI excluding 1 = clear effect):\n")
print(or_sat)

# Does the spatial map shrink once the indices are in?
sp_sum <- do.call(rbind, lapply(c("B0", "BX"), function(nm) {
  r <- inla.spde2.result(fits[[nm]], "s", spde_u, do.transf = TRUE)
  data.frame(model = nm,
             range_km = round(inla.qmarginal(0.5, r$marginals.range.nominal[[1]]), 1),
             sd_space = round(sqrt(inla.qmarginal(0.5, r$marginals.variance.nominal[[1]])), 2))
}))
cat("\nSpatial map (median range in km, sd on the logit scale):\n"); print(sp_sum)

write.csv(tab,    file.path(OUT_DIR, paste0("s2_waic", SFX, ".csv")),         row.names = FALSE)
write.csv(or_sat, file.path(OUT_DIR, paste0("s2_odds_ratios", SFX, ".csv")),  row.names = FALSE)
write.csv(sp_sum, file.path(OUT_DIR, paste0("s2_spatial_sd", SFX, ".csv")),   row.names = FALSE)


## ---- 4. LEAVE-ONE-STATION-OUT AUC (slow) --------------------------------------
if (RUN_CV) {
  auc <- function(y, p) {
    r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  cv_stations <- names(which(table(u$station) >= 5))
  loso <- function(nm) {
    held_p <- rep(NA_real_, nrow(u))
    for (k in seq_along(cv_stations)) {
      held <- u$station == cv_stations[k]
      y_cv <- u$y; y_cv[held] <- NA
      r    <- fit_u(forms[[nm]], y_cv, ctrl = list(int.strategy = "eb"), waic = FALSE, pred = TRUE)
      pr   <- r$fit$summary.fitted.values$mean[inla.stack.index(r$stk, "est")$data]
      held_p[held] <- pr[held]
      if (k %% 10 == 0) cat(nm, ": station", k, "of", length(cv_stations), "\n")
    }
    held_p
  }
  held <- list()
  for (nm in names(forms)) { cat("Leave-one-station-out:", nm, "\n"); held[[nm]] <- loso(nm) }
  saveRDS(held, file.path(OUT_DIR, paste0("s2_loso_predictions", SFX, ".rds")))

  ok2 <- !is.na(held$B0)
  cv  <- data.frame(model = names(held),
                    AUC = sapply(held, function(p) round(auc(u$y[ok2], p[ok2]), 3)), row.names = NULL)
  cat("\nLeave-one-station-out AUC on", sum(ok2), "samples:\n"); print(cv)

  # Gain from adding the indices, with a station bootstrap interval
  boot_gain <- function(base, plus) {
    st <- unique(u$station[ok2]); idx <- split(which(ok2), u$station[ok2])
    g <- replicate(N_BOOT, {
      ix <- unlist(idx[sample(st, length(st), replace = TRUE)], use.names = FALSE)
      if (length(unique(u$y[ix])) < 2) return(NA_real_)
      auc(u$y[ix], held[[plus]][ix]) - auc(u$y[ix], held[[base]][ix])
    })
    quantile(g, c(0.025, 0.975), na.rm = TRUE)
  }
  set.seed(1)
  cat("\nAUC gain from the indices (95% interval, station bootstrap):\n")
  cat("  without map  (AX - A0):", round(auc(u$y[ok2], held$AX[ok2]) - auc(u$y[ok2], held$A0[ok2]), 3),
      " [", round(boot_gain("A0", "AX"), 3), "]\n")
  cat("  with map     (BX - B0):", round(auc(u$y[ok2], held$BX[ok2]) - auc(u$y[ok2], held$B0[ok2]), 3),
      " [", round(boot_gain("B0", "BX"), 3), "]\n")

  # Within each season: do the indices help rank samples inside a season?
  by_season <- do.call(rbind, lapply(unique(u$season), function(sn) {
    w <- ok2 & u$season == sn
    if (length(unique(u$y[w])) < 2) return(NULL)
    data.frame(season = sn, n = sum(w),
               B0 = round(auc(u$y[w], held$B0[w]), 3), BX = round(auc(u$y[w], held$BX[w]), 3))
  }))
  cat("\nAUC within each season:\n"); print(by_season)

  write.csv(cv,        file.path(OUT_DIR, paste0("s2_loso_auc", SFX, ".csv")),      row.names = FALSE)
  write.csv(by_season, file.path(OUT_DIR, paste0("s2_loso_by_season", SFX, ".csv")), row.names = FALSE)
}

cat("\nDone. Send me everything printed above.\n")
