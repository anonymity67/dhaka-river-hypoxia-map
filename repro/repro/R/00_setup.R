# =====================================================================
# 00_setup.R - shared code: read and clean the data, build the mesh and the
# spatial model, define helper functions. Every other script starts with
#   source("R/00_setup.R")
# =====================================================================
if (!file.exists("config.R")) stop("Run from the repository root (the folder that contains config.R).")
source("config.R")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(INLA); library(Matrix); library(sf) })
if (!requireNamespace("fmesher", quietly = TRUE)) stop("Install the fmesher package.")
options(timeout = 3600)
set.seed(1)

## ---- Data ----------------------------------------------------------------
stopifnot(file.exists(DATA_FILE))
raw <- read.csv(DATA_FILE, fileEncoding = "UTF-8-BOM", check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(ncol(raw) == 16)
names(raw) <- c("date", "river", "location", "lat", "lon", "temp", "ph", "do",
                "cod", "no3", "turb", "po4", "rain", "season", "wqi", "label")
d <- raw
d$river <- trimws(d$river)
d$river[d$river == "Tongi khal"] <- "Tongi Khal"              # two spellings of one river
key <- tolower(gsub("[- ]", "", trimws(d$season)))
d$season <- unname(c(winter = "Winter", premonsoon = "Pre-monsoon",
                     monsoon = "Monsoon", postmonsoon = "Post-monsoon")[key])
stopifnot(!anyNA(d$season))
for (v in c("temp", "ph", "do", "cod", "no3", "turb", "po4", "rain", "wqi"))
  d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
d$sample_id <- seq_len(nrow(d))                                   # row order of the file
d$year    <- as.integer(substr(d$date, 1, 4))
d$month   <- suppressWarnings(as.integer(substr(d$date, 6, 7)))
d$year_id <- d$year - min(d$year) + 1
d$y       <- as.integer(d$do < HYPOXIA_MGL)                       # 1 = hypoxic sample
d$station <- paste(d$river, d$location, sep = " | ")
d$s_pre   <- as.integer(d$season == "Pre-monsoon")
d$s_mon   <- as.integer(d$season == "Monsoon")
d$s_post  <- as.integer(d$season == "Post-monsoon")               # Winter is the reference
seasons   <- c("Winter", "Pre-monsoon", "Monsoon", "Post-monsoon")

## ---- Coordinates (UTM 46N, kilometres), mesh, spatial model -------------
pts <- st_transform(st_as_sf(d, coords = c("lon", "lat"), crs = 4326), 32646)
xy  <- st_coordinates(pts) / 1000
mesh <- fmesher::fm_mesh_2d_inla(loc = xy, max.edge = c(3, 12), cutoff = 1, offset = c(5, 20))
cat("Samples:", nrow(d), "| stations:", length(unique(d$station)), "| mesh nodes:", mesh$n, "\n")
if (nrow(d) != 1495 || mesh$n != 1136)
  warning("Expected 1,495 samples and 1,136 mesh nodes (as in the paper); check the data file and fmesher version.")

spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2,
                            prior.range = c(RANGE_PRIOR_KM, 0.5), prior.sigma = SIGMA_PRIOR)
A_obs   <- fmesher::fm_basis(mesh, loc = xy)
cov_obs <- data.frame(Intercept = 1, s_pre = d$s_pre, s_mon = d$s_mon, s_post = d$s_post, year_id = d$year_id)

## ---- Model formulas ---------------------------------------------------------
YEAR_TERM <- sprintf("f(year_id, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(%g, %g))))",
                     YEAR_PRIOR[1], YEAR_PRIOR[2])
form_A <- as.formula(paste("y ~ 0 + Intercept + s_pre + s_mon + s_post +", YEAR_TERM))               # season + year
form_B <- as.formula(paste("y ~ 0 + Intercept + s_pre + s_mon + s_post +", YEAR_TERM, "+ f(s, model = spde)"))  # + spatial field

## ---- Helper functions ---------------------------------------------------------
stk_est <- function(y, A, cov, tag = "est")                       # estimation stack for any subset of rows
  inla.stack(tag = tag, data = list(y = y, Ntrials = rep(1, length(y))), A = list(A, 1),
             effects = list(s = seq_len(spde$n.spde), cov))

run_inla <- function(form, stk, waic = TRUE, pred = FALSE, ctrl = list()) {
  dd <- inla.stack.data(stk)
  inla(form, family = "binomial", Ntrials = dd$Ntrials, data = dd,
       control.predictor = list(A = inla.stack.A(stk), compute = pred),
       control.compute = list(waic = waic), control.inla = ctrl)
}

auc <- function(y, p) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
                        (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

# Leave-one-station-out predictions: hide each station (with >= min_n samples), predict it from the rest.
# 'keep' selects which rows of d are used at all. Returns the predicted probability for every kept row.
loso_pred <- function(form, keep = rep(TRUE, nrow(d)), min_n = 5) {
  idx <- which(keep); u <- d[idx, ]; A <- A_obs[idx, , drop = FALSE]; cov <- cov_obs[idx, , drop = FALSE]
  held_p <- rep(NA_real_, nrow(u)); sts <- names(which(table(u$station) >= min_n))
  for (k in seq_along(sts)) {
    held <- u$station == sts[k]; y_cv <- u$y; y_cv[held] <- NA
    stk <- stk_est(y_cv, A, cov)
    fit <- run_inla(form, stk, waic = FALSE, pred = TRUE, ctrl = list(int.strategy = "eb"))
    pr  <- fit$summary.fitted.values$mean[inla.stack.index(stk, "est")$data]
    held_p[held] <- pr[held]
    if (k %% 10 == 0) cat("  station", k, "of", length(sts), "\n")
  }
  list(idx = idx, p = held_p)
}
loso_auc <- function(form, keep = rep(TRUE, nrow(d))) {
  r <- loso_pred(form, keep); ok <- !is.na(r$p); auc(d$y[r$idx][ok], r$p[ok]) }

odds_ratios <- function(fit, nm = c("s_pre", "s_mon", "s_post")) {
  q <- as.matrix(fit$summary.fixed[nm, c("mean", "0.025quant", "0.975quant")])
  data.frame(effect = nm, OR = round(exp(q[, 1]), 3), lower95 = round(exp(q[, 2]), 3),
             upper95 = round(exp(q[, 3]), 3), row.names = NULL)
}
