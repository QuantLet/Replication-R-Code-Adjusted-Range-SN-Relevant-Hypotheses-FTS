rm(list = ls())

# Build the six OSP figure panels used in the manuscript (fig:OSPrejection_probabilities)
# Reads simulation tables from output/simulations/OSP and writes PNGs to output/figures/OSP.

library(fda)

source(file.path("R", "genData.R"))
source(file.path("R", "generalMethods.R"))

# -----------------------------
# Plot helpers
# -----------------------------
plotPowerCurve <- function(base_path, Delta, deltaVector, q) {
  ntables <- length(deltaVector)
  n25list <- readTables(file.path(base_path, "n25"), ntables)
  n50list <- readTables(file.path(base_path, "n50"), ntables)
  n100list <- readTables(file.path(base_path, "n100"), ntables)

  plot(ylim = c(0, 1), x = deltaVector, y = empRejProb(n25list, Delta = Delta, q = q),
       xlab = expression(delta), ylab = "Empirical rejection probability",
       type = "l", lty = 1, font.lab = 2, lwd = 0.5, cex.lab = 1,
       bty = "n", axes = FALSE)

  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)

  lines(x = deltaVector, y = empRejProb(n50list, Delta = Delta, q = q), lty = 2, lwd = 0.5)
  lines(x = deltaVector, y = empRejProb(n100list, Delta = Delta, q = q), lty = 3, lwd = 0.5)

  lines(x = c(0, deltaVector[ntables]), y = c(0.05, 0.05), lwd = 0.5)
  abline(v = Delta, lty = 2, lwd = 0.5)

  legend(0.01, 1, c("n = 25", "n = 50", "n = 100"),
         lty = c(1, 2, 3), lwd = 0.5, cex = 1, box.lwd = 0.5)
}

plotPowerCurve_hong <- function(base_path, Delta, deltaVector, q) {
  ntables <- length(deltaVector)
  n25list <- readTables_hong(file.path(base_path, "n25"), ntables)
  n50list <- readTables_hong(file.path(base_path, "n50"), ntables)
  n100list <- readTables_hong(file.path(base_path, "n100"), ntables)

  plot(ylim = c(0, 1), x = deltaVector, y = empRejProb(n25list, Delta = Delta, q = q),
       xlab = expression(delta), ylab = "Empirical rejection probability",
       type = "l", lty = 1, font.lab = 2, lwd = 0.5, cex.lab = 1,
       bty = "n", axes = FALSE)

  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)

  lines(x = deltaVector, y = empRejProb(n50list, Delta = Delta, q = q), lty = 2, lwd = 0.5)
  lines(x = deltaVector, y = empRejProb(n100list, Delta = Delta, q = q), lty = 3, lwd = 0.5)

  lines(x = c(0, deltaVector[ntables]), y = c(0.05, 0.05), lwd = 0.5)
  abline(v = Delta, lty = 2, lwd = 0.5)

  legend(0.01, 1, c("n = 25", "n = 50", "n = 100"),
         lty = c(1, 2, 3), lwd = 0.5, cex = 1, box.lwd = 0.5)
}

# -----------------------------
# Critical values (Brownian-motion Monte‑Carlo)
# -----------------------------
set.seed(2)
BMsample <- BM(N = 300, n = 1000)

q_shao <- quant.fct(BMsample, prob = 0.95)
q_hong <- quant.fct_hong(BMsample, prob = 0.95)

# -----------------------------
# Figure settings (from manuscript)
# -----------------------------
Delta <- 0.02
deltaVector <- 10:40/1000

base_sim <- file.path("output", "simulations", "OSP")
out_dir <- file.path("output", "figures", "OSP")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Produce the six panels used in the manuscript
# -----------------------------

# fIID
png(file.path(out_dir, "OSPfIID_hong.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(file.path(base_sim, "fIID"), Delta = Delta, deltaVector = deltaVector, q = q_hong)
dev.off()

png(file.path(out_dir, "OSPfIID.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(file.path(base_sim, "fIID"), Delta = Delta, deltaVector = deltaVector, q = q_shao)
dev.off()

# fMA1
png(file.path(out_dir, "OSPfMA1_hong.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(file.path(base_sim, "fMA1"), Delta = Delta, deltaVector = deltaVector, q = q_hong)
dev.off()

png(file.path(out_dir, "OSPfMA1.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(file.path(base_sim, "fMA1"), Delta = Delta, deltaVector = deltaVector, q = q_shao)
dev.off()

# Brownian bridge
png(file.path(out_dir, "OSPBB_hong.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(file.path(base_sim, "BB"), Delta = Delta, deltaVector = deltaVector, q = q_hong)
dev.off()

png(file.path(out_dir, "OSPBB.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(file.path(base_sim, "BB"), Delta = Delta, deltaVector = deltaVector, q = q_shao)
dev.off()

message("OSP figures written to: ", normalizePath(out_dir))
