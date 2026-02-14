rm(list = ls())

# Build the CPP power-curve panels (fig:CPP_rejection_probabilities)
# Reads simulation tables from output/simulations/CPP/fIID and writes PNGs to output/figures/CPP.

library(fda)

source(file.path("R", "genData.R"))
source(file.path("R", "generalMethods.R"))

plotPowerCurve <- function(base_path, Delta, aVector, q) {
  ntables <- length(aVector)
  n200list <- readTables(file.path(base_path, "n200"), ntables)
  n500list <- readTables(file.path(base_path, "n500"), ntables)

  plot(ylim = c(0, 1), x = aVector, y = empRejProb(n200list, Delta = Delta, q = q),
       xlab = expression(a), ylab = "Empirical rejection probability",
       type = "l", lty = 1, font.lab = 2, lwd = 0.5, cex.lab = 1,
       bty = "n", axes = FALSE)

  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)

  lines(x = aVector, y = empRejProb(n500list, Delta = Delta, q = q), lty = 2, lwd = 0.5)

  lines(x = c(0, aVector[ntables]), y = c(0.05, 0.05), lwd = 0.5)
  abline(v = sqrt(Delta * 30), lty = 2, lwd = 0.5)

  legend(0, 1, c("N = 200", "N = 500"), lty = c(1, 2), lwd = 0.5, cex = 1,
         box.lwd = 0.5)
}

plotPowerCurve_hong <- function(base_path, Delta, aVector, q) {
  ntables <- length(aVector)
  n200list <- readTables_hong(file.path(base_path, "n200"), ntables)
  n500list <- readTables_hong(file.path(base_path, "n500"), ntables)

  plot(ylim = c(0, 1), x = aVector, y = empRejProb(n200list, Delta = Delta, q = q),
       xlab = expression(a), ylab = "Empirical rejection probability",
       type = "l", lty = 1, font.lab = 2, lwd = 0.5, cex.lab = 1,
       bty = "n", axes = FALSE)

  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)

  lines(x = aVector, y = empRejProb(n500list, Delta = Delta, q = q), lty = 2, lwd = 0.5)

  lines(x = c(0, aVector[ntables]), y = c(0.05, 0.05), lwd = 0.5)
  abline(v = sqrt(Delta * 30), lty = 2, lwd = 0.5)

  legend(0, 1, c("N = 200", "N = 500"), lty = c(1, 2), lwd = 0.5, cex = 1,
         box.lwd = 0.5)
}

# -----------------------------
# Critical values
# -----------------------------
set.seed(2)
BMsample <- BM(N = 300, n = 1000)
q_shao <- quant.fct(BMsample, prob = 0.95)

set.seed(3)
BMsample2 <- BM(N = 300, n = 1000)
q_hong <- quant.fct_hong(BMsample2, prob = 0.95)

# -----------------------------
# Build the six panels
# -----------------------------
aVector <- 0:25/50
base_sim <- file.path("output", "simulations", "CPP", "fIID")
out_dir <- file.path("output", "figures", "CPP")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Δ = 0.1^2 / 30
Delta <- 0.1^2 / 30
png(file.path(out_dir, "a01power_hong.png"), width = 5.25, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(base_sim, Delta = Delta, aVector = aVector, q = q_hong)
dev.off()

png(file.path(out_dir, "a01power.png"), width = 5.25, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(base_sim, Delta = Delta, aVector = aVector, q = q_shao)
dev.off()

# Δ = 0.2^2 / 30
Delta <- 0.2^2 / 30
png(file.path(out_dir, "a02power_hong.png"), width = 5.25, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(base_sim, Delta = Delta, aVector = aVector, q = q_hong)
dev.off()

png(file.path(out_dir, "a02power.png"), width = 5.25, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(base_sim, Delta = Delta, aVector = aVector, q = q_shao)
dev.off()

# Δ = 0.3^2 / 30
Delta <- 0.3^2 / 30
png(file.path(out_dir, "a03power_hong.png"), width = 5.25, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(base_sim, Delta = Delta, aVector = aVector, q = q_hong)
dev.off()

png(file.path(out_dir, "a03power.png"), width = 5.25, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(base_sim, Delta = Delta, aVector = aVector, q = q_shao)
dev.off()

message("CPP power-curve panels written to: ", normalizePath(out_dir))
