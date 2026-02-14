rm(list = ls())

# Build the six TSP figure panels used in the manuscript (fig:TSPrejection_probabilities)
# Reads simulation tables from output/simulations/TSP and writes PNGs to output/figures/TSP.

library(fda)

source(file.path("R", "genData.R"))
source(file.path("R", "generalMethods.R"))

plotPowerCurve <- function(base_path, Delta, aVector, q) {
  ntables <- length(aVector)
  n50m50list <- readTables(file.path(base_path, "n50m50"), ntables)
  n50m100list <- readTables(file.path(base_path, "n50m100"), ntables)
  n100m100list <- readTables(file.path(base_path, "n100m100"), ntables)
  n100m200list <- readTables(file.path(base_path, "n100m200"), ntables)

  plot(ylim = c(0, 1), x = aVector, y = empRejProb(n50m50list, Delta = Delta, q = q),
       xlab = expression(a), ylab = "Empirical rejection probability",
       type = "l", lty = 1, font.lab = 2, lwd = 0.5, cex.lab = 1,
       bty = "n", axes = FALSE)

  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)

  lines(x = aVector, y = empRejProb(n50m100list, Delta = Delta, q = q), lty = 2, lwd = 0.5)
  lines(x = aVector, y = empRejProb(n100m100list, Delta = Delta, q = q), lty = 3, lwd = 0.5)
  lines(x = aVector, y = empRejProb(n100m200list, Delta = Delta, q = q), lty = 4, lwd = 0.5)

  lines(x = c(0, aVector[ntables]), y = c(0.05, 0.05), lwd = 0.5)
  abline(v = sqrt(Delta * 30), lty = 2, lwd = 0.5)

  legend(0, 1, c("n = m = 50", "n = 50, m = 100", "n = m = 100", "n = 100, m = 200"),
         lty = c(1, 2, 3, 4), lwd = 0.5, cex = 1, box.lwd = 0.5)
}

plotPowerCurve_hong <- function(base_path, Delta, aVector, q) {
  ntables <- length(aVector)
  n50m50list <- readTables_hong(file.path(base_path, "n50m50"), ntables)
  n50m100list <- readTables_hong(file.path(base_path, "n50m100"), ntables)
  n100m100list <- readTables_hong(file.path(base_path, "n100m100"), ntables)
  n100m200list <- readTables_hong(file.path(base_path, "n100m200"), ntables)

  plot(ylim = c(0, 1), x = aVector, y = empRejProb(n50m50list, Delta = Delta, q = q),
       xlab = expression(a), ylab = "Empirical rejection probability",
       type = "l", lty = 1, font.lab = 2, lwd = 0.5, cex.lab = 1,
       bty = "n", axes = FALSE)

  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)

  lines(x = aVector, y = empRejProb(n50m100list, Delta = Delta, q = q), lty = 2, lwd = 0.5)
  lines(x = aVector, y = empRejProb(n100m100list, Delta = Delta, q = q), lty = 3, lwd = 0.5)
  lines(x = aVector, y = empRejProb(n100m200list, Delta = Delta, q = q), lty = 4, lwd = 0.5)

  lines(x = c(0, aVector[ntables]), y = c(0.05, 0.05), lwd = 0.5)
  abline(v = sqrt(Delta * 30), lty = 2, lwd = 0.5)

  legend(0, 1, c("n = m = 50", "n = 50, m = 100", "n = m = 100", "n = 100, m = 200"),
         lty = c(1, 2, 3, 4), lwd = 0.5, cex = 1, box.lwd = 0.5)
}

# -----------------------------
# Critical values
# -----------------------------
simInfo_cv <- list(normalizerType = "V_n")

set.seed(2)
BMsample <- BM(N = 300, n = 1000)
q_shao <- quant.fct(BMsample, prob = 0.95, simInfo = simInfo_cv)
q_hong <- quant.fct_hong(BMsample, prob = 0.95, simInfo = simInfo_cv)

# -----------------------------
# Figure settings (from manuscript)
# -----------------------------
Delta <- 0.2^2 / 30

base_sim <- file.path("output", "simulations", "TSP")
out_dir <- file.path("output", "figures", "TSP")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# fIID
aVector <- 0:50/100
png(file.path(out_dir, "TSPfIID_hong.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(file.path(base_sim, "fIID"), Delta = Delta, aVector = aVector, q = q_hong)
dev.off()

png(file.path(out_dir, "TSPfIID.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(file.path(base_sim, "fIID"), Delta = Delta, aVector = aVector, q = q_shao)
dev.off()

# fMA1 (independent)
aVector <- 0:50/100
png(file.path(out_dir, "TSPfMA1_independent_hong.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(file.path(base_sim, "fMA1_independent"), Delta = Delta, aVector = aVector, q = q_hong)
dev.off()

png(file.path(out_dir, "TSPfMA1_independent.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(file.path(base_sim, "fMA1_independent"), Delta = Delta, aVector = aVector, q = q_shao)
dev.off()

# Brownian bridge
aVector <- 0:(14-1)/10
png(file.path(out_dir, "TSPBB_hong.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(file.path(base_sim, "BB"), Delta = Delta, aVector = aVector, q = q_hong)
dev.off()

png(file.path(out_dir, "TSPBB.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(file.path(base_sim, "BB"), Delta = Delta, aVector = aVector, q = q_shao)
dev.off()

message("TSP figures written to: ", normalizePath(out_dir))
