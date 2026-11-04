# =====================================================================
# 03_robustness.R - how stable is model B?  Plus model C (wet-season field)
# =====================================================================
# Refits model B with: a finer and a coarser mesh, a station-level random
# effect, the full Laplace strategy; and model C = B + one extra spatial field
# that acts only in the wet seasons (Monsoon, Post-monsoon).
# Output: results/robustness_summary.csv.   Each fit takes about 1-5 minutes.
source("R/00_setup.R")
RUN_MESH <- TRUE; RUN_STATION <- TRUE; RUN_LAPLACE <- TRUE; RUN_WET <- TRUE
WET <- c("Monsoon", "Post-monsoon")
d$station_id <- as.integer(factor(d$station))
d$wet <- as.integer(d$season %in% WET)

YEAR_PC <- sprintf("c(%g, %g)", YEAR_PRIOR[1], YEAR_PRIOR[2])
PC_STATION <- paste0(" + f(station_id, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = ", YEAR_PC, ")))")
WET_FIELD  <- " + f(s2, model = spde_tmp)"
make_form <- function(extra = "") as.formula(paste0(
  "y ~ 0 + Intercept + s_pre + s_mon + s_post + ",
  "f(year_id, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = ", YEAR_PC, "))) + ",
  "f(s, model = spde_tmp)", extra))

build_stk <- function(yvec, A_i, sp, wet = FALSE, tag = "est") {
  idx <- seq_len(sp$n.spde); cov_i <- cbind(cov_obs, station_id = d$station_id)
  dat <- list(y = yvec, Ntrials = rep(1, length(yvec)))
  if (wet) inla.stack(tag = tag, data = dat, A = list(A_i, Diagonal(x = d$wet) %*% A_i, 1),
                      effects = list(s = idx, s2 = idx, cov_i))
  else inla.stack(tag = tag, data = dat, A = list(A_i, 1), effects = list(s = idx, cov_i))
}
fit_est <- function(mesh_i, extra = "", wet = FALSE, ctrl = list()) {
  sp <- inla.spde2.pcmatern(mesh_i, alpha = 2, prior.range = c(RANGE_PRIOR_KM, 0.5), prior.sigma = SIGMA_PRIOR)
  assign("spde_tmp", sp, envir = .GlobalEnv)            # the formula looks the model up by this name
  stk <- build_stk(d$y, fmesher::fm_basis(mesh_i, loc = xy), sp, wet); dd <- inla.stack.data(stk)
  fit <- inla(make_form(extra), family = "binomial", Ntrials = dd$Ntrials, data = dd,
              control.predictor = list(A = inla.stack.A(stk), compute = FALSE),
              control.compute = list(waic = TRUE), control.inla = ctrl)
  list(fit = fit, spde = sp, nodes = mesh_i$n)
}
med_sd <- function(m) inla.qmarginal(0.5, inla.tmarginal(function(x) 1 / sqrt(x), m))
grab <- function(res, label) {
  fit <- res$fit; fx <- exp(fit$summary.fixed$mean); names(fx) <- rownames(fit$summary.fixed)
  hp <- fit$marginals.hyperpar; r1 <- inla.spde2.result(fit, "s", res$spde, do.transf = TRUE)
  out <- data.frame(model = label, nodes = res$nodes, WAIC = round(fit$waic$waic, 1),
    OR_pre = round(fx["s_pre"], 3), OR_mon = round(fx["s_mon"], 3), OR_post = round(fx["s_post"], 3),
    range_km = round(inla.qmarginal(0.5, r1$marginals.range.nominal[[1]]), 1),
    sd_space = round(sqrt(inla.qmarginal(0.5, r1$marginals.variance.nominal[[1]])), 2),
    sd_year = round(med_sd(hp[["Precision for year_id"]]), 2), sd_station = NA_real_, sd_wet = NA_real_, row.names = NULL)
  if ("Precision for station_id" %in% names(hp)) out$sd_station <- round(med_sd(hp[["Precision for station_id"]]), 2)
  if ("s2" %in% names(fit$summary.random)) {
    r2 <- inla.spde2.result(fit, "s2", res$spde, do.transf = TRUE)
    out$sd_wet <- round(sqrt(inla.qmarginal(0.5, r2$marginals.variance.nominal[[1]])), 2) }
  out
}

rows <- list(base = grab(fit_est(mesh), "B: baseline (max.edge 3/12)"))
if (RUN_MESH) {
  m_fine   <- fmesher::fm_mesh_2d_inla(loc = xy, max.edge = c(2, 8),  cutoff = 0.7, offset = c(5, 20))
  m_coarse <- fmesher::fm_mesh_2d_inla(loc = xy, max.edge = c(5, 20), cutoff = 1.5, offset = c(5, 20))
  rows$fine   <- grab(fit_est(m_fine),   "B: finer mesh (2/8)")
  rows$coarse <- grab(fit_est(m_coarse), "B: coarser mesh (5/20)")
}
if (RUN_STATION) rows$station <- grab(fit_est(mesh, extra = PC_STATION), "B + station effect")
if (RUN_LAPLACE) rows$laplace <- grab(fit_est(mesh, ctrl = list(strategy = "laplace", int.strategy = "ccd")), "B: full Laplace strategy")
if (RUN_WET)     rows$wet     <- grab(fit_est(mesh, extra = WET_FIELD, wet = TRUE), "C: B + wet-season field")

tab <- do.call(rbind, rows); rownames(tab) <- NULL
print(tab)
write.csv(tab, file.path(OUT_DIR, "robustness_summary.csv"), row.names = FALSE)
cat("\nHow to read: odds ratios and range should stay close to the baseline row; WAIC differences of a few points are small.\n")
