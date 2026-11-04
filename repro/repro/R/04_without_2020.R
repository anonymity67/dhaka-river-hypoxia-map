# =====================================================================
# 04_without_2020.R - does the main result hold without the 2020 campaign?
# =====================================================================
# Fits model B on all years and on 2017-2019 + 2021 (same mesh, same priors) and compares season
# odds ratios, spatial range/sd, and the spatial field at the 60 stations.
# Output: results/without_2020_season_or.csv, without_2020_stations.csv
source("R/00_setup.R")
keep <- d$year != 2020
fit_all <- run_inla(form_B, stk_est(d$y, A_obs, cov_obs))
fit_no  <- run_inla(form_B, stk_est(d$y[keep], A_obs[keep, ], cov_obs[keep, ]))
cat("Samples: all =", nrow(d), "| without 2020 =", sum(keep), "\n")

orr <- rbind(cbind(sample = "all years", odds_ratios(fit_all)), cbind(sample = "without 2020", odds_ratios(fit_no)))
print(orr, row.names = FALSE)
write.csv(orr, file.path(OUT_DIR, "without_2020_season_or.csv"), row.names = FALSE)
hp <- function(f) round(f$summary.hyperpar[grep("Range|Stdev", rownames(f$summary.hyperpar)), c("mean", "0.025quant", "0.975quant")], 2)
cat("\nSpatial range (km) and sd, all years:\n");    print(hp(fit_all))
cat("\nSpatial range (km) and sd, without 2020:\n"); print(hp(fit_no))

idx  <- match(unique(d$station), d$station)
A_st <- fmesher::fm_basis(mesh, loc = xy[idx, , drop = FALSE])
f_all <- as.numeric(A_st %*% fit_all$summary.random$s$mean)
f_no  <- as.numeric(A_st %*% fit_no$summary.random$s$mean)
top <- function(x, k = 15) order(-x)[seq_len(k)]
cat("\nSpatial field at the stations, all years vs without 2020:\n")
cat("  Spearman correlation =", round(cor(f_all, f_no, method = "spearman"), 3), "\n")
cat("  Top-15 hotspot stations in common =", length(intersect(top(f_all), top(f_no))), "of 15\n")
write.csv(data.frame(station = unique(d$station), field_all = round(f_all, 3), field_no2020 = round(f_no, 3)),
          file.path(OUT_DIR, "without_2020_stations.csv"), row.names = FALSE)

if (RUN_LOSO) cat("\nLeave-one-station-out AUC: all years =", round(loso_auc(form_B), 3),
                  "| without 2020 =", round(loso_auc(form_B, keep), 3), "\n")
