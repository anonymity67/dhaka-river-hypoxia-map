# =====================================================================
# 06_mesh_figure.R - evidence for the SPDE mesh: figure, summary table, saved object
# =====================================================================
# Rebuilds the mesh with the exact settings used in the paper (defined in 00_setup.R),
# draws it with the monitoring stations, and records its size and edge lengths.
# Output (in OUT_DIR): fig_mesh.png, fig_mesh.pdf, mesh_summary.csv, mesh.rds, mesh_session_info.txt
source("R/00_setup.R")

tv  <- mesh$graph$tv                 # triangle vertex indices
loc <- mesh$loc[, 1:2]               # node coordinates (UTM 46N, km)
edge_len <- function(a, b) sqrt(rowSums((loc[a, ] - loc[b, ])^2))
el <- c(edge_len(tv[, 1], tv[, 2]), edge_len(tv[, 2], tv[, 3]), edge_len(tv[, 3], tv[, 1]))

st_idx <- match(unique(d$station), d$station)         # one row per station
hull   <- chull(xy[st_idx, ])                          # study-area outline

draw_mesh <- function() {
  plot(NA, xlim = range(loc[, 1]), ylim = range(loc[, 2]), asp = 1, xlab = "Easting (km, UTM 46N)",
       ylab = "Northing (km, UTM 46N)", main = "")
  segments(loc[tv[, 1], 1], loc[tv[, 1], 2], loc[tv[, 2], 1], loc[tv[, 2], 2], col = "grey60", lwd = 0.3)
  segments(loc[tv[, 2], 1], loc[tv[, 2], 2], loc[tv[, 3], 1], loc[tv[, 3], 2], col = "grey60", lwd = 0.3)
  segments(loc[tv[, 3], 1], loc[tv[, 3], 2], loc[tv[, 1], 1], loc[tv[, 1], 2], col = "grey60", lwd = 0.3)
  polygon(xy[st_idx, ][hull, ], border = "steelblue", lty = 2, lwd = 1.2)
  points(xy[st_idx, ], pch = 21, bg = "firebrick", col = "white", cex = 0.8)
  legend("topright", bty = "n", cex = 0.8, pch = c(21, NA), lty = c(NA, 2), col = c("white", "steelblue"),
         pt.bg = c("firebrick", NA), legend = c("Monitoring stations (60)", "Station outline"))
}
png(file.path(OUT_DIR, "fig_mesh.png"), width = 2400, height = 2400, res = 300); draw_mesh(); dev.off()
pdf(file.path(OUT_DIR, "fig_mesh.pdf"), width = 8, height = 8); draw_mesh(); dev.off()

summ <- data.frame(
  quantity = c("Nodes", "Triangles", "max.edge inner / outer (km)", "cutoff (km)", "offset inner / outer (km)",
               "Edge length min (km)", "Edge length median (km)", "Edge length max (km)", "Samples", "Stations",
               "fmesher version", "INLA version"),
  value = c(mesh$n, nrow(tv), "3 / 12", "1", "5 / 20",
            round(min(el), 2), round(median(el), 2), round(max(el), 2), nrow(d), length(st_idx),
            as.character(packageVersion("fmesher")), as.character(packageVersion("INLA"))))
print(summ, row.names = FALSE)
write.csv(summ, file.path(OUT_DIR, "mesh_summary.csv"), row.names = FALSE)
saveRDS(mesh, file.path(OUT_DIR, "mesh.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "mesh_session_info.txt"))
cat("\nMesh nodes:", mesh$n, "(paper: 1,136). Figure and tables written to", OUT_DIR, "\n")
