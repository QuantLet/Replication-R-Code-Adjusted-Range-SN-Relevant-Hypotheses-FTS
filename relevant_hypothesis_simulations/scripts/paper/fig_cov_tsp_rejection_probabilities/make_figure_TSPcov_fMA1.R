rm(list = ls())

# Build the TSP covariance-operator panels (fig:cov_tsp_rejection_probabilities)
# Reads simulation tables from output/simulations/TSPcov/fMA1 and writes PNGs to output/figures/Covariance.

library(fda)

source(file.path("R", "genData.R"))
source(file.path("R", "generalMethods.R"))

plotPowerCurve <- function(base_path, Delta, factorVector, q) {
  ntables <- length(factorVector)
  n100m100list <- readTables(file.path(base_path, "n100m100"), ntables)
  n200m200list <- readTables(file.path(base_path, "n200m200"), ntables)
  n500m500list <- readTables(file.path(base_path, "n500m500"), ntables)

  plot(ylim = c(0, 1), x = factorVector, y = empRejProb(n100m100list, Delta = Delta, q = q),
       xlab = expression(a), ylab = "Empirical rejection probability",
       type = "l", lty = 1, font.lab = 2, lwd = 0.5, cex.lab = 1,
       bty = "n", axes = FALSE)

  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)

  lines(x = factorVector, y = empRejProb(n200m200list, Delta = Delta, q = q), lty = 2, lwd = 0.5)
  lines(x = factorVector, y = empRejProb(n500m500list, Delta = Delta, q = q), lty = 3, lwd = 0.5)

  lines(x = c(factorVector[1], factorVector[ntables]), y = c(0.05, 0.05), lwd = 0.5)

  legend(1, 1, c("n = m = 100", "n = m = 200", "n = m = 500"),
         lty = c(1, 2, 3), lwd = 0.5, cex = 1, box.lwd = 0.5)
}

plotPowerCurve_hong <- function(base_path, Delta, factorVector, q) {
  ntables <- length(factorVector)
  n100m100list <- readTables_hong(file.path(base_path, "n100m100"), ntables)
  n200m200list <- readTables_hong(file.path(base_path, "n200m200"), ntables)
  n500m500list <- readTables_hong(file.path(base_path, "n500m500"), ntables)

  plot(ylim = c(0, 1), x = factorVector, y = empRejProb(n100m100list, Delta = Delta, q = q),
       xlab = expression(a), ylab = "Empirical rejection probability",
       type = "l", lty = 1, font.lab = 2, lwd = 0.5, cex.lab = 1,
       bty = "n", axes = FALSE)

  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)

  lines(x = factorVector, y = empRejProb(n200m200list, Delta = Delta, q = q), lty = 2, lwd = 0.5)
  lines(x = factorVector, y = empRejProb(n500m500list, Delta = Delta, q = q), lty = 3, lwd = 0.5)

  lines(x = c(factorVector[1], factorVector[ntables]), y = c(0.05, 0.05), lwd = 0.5)

  legend(1, 1, c("n = m = 100", "n = m = 200", "n = m = 500"),
         lty = c(1, 2, 3), lwd = 0.5, cex = 1, box.lwd = 0.5)
}

# Critical values
set.seed(2)
BMsample <- BM(N = 300, n = 1000)
q_shao <- quant.fct(BMsample, prob = 0.95)

set.seed(3)
BMsample2 <- BM(N = 300, n = 1000)
q_hong <- quant.fct_hong(BMsample2, prob = 0.95)

factorVector <- seq(from = 1, to = 2.5, by = 0.025)

# Manuscript settings
a0 <- 1.5
D <- 21

# varType A
SigmaA <- 1 / seq(1, D)
DeltaA <- (1 - a0^2)^2 * sum(SigmaA^4) * (1 + 0.7^2)^2

# varType B
SigmaB <- 1.2^(-seq(1, D))
DeltaB <- (1 - a0^2)^2 * sum(SigmaB^4) * (1 + 0.7^2)^2

base_sim <- file.path("output", "simulations", "TSPcov", "fMA1")
out_dir <- file.path("output", "figures", "Covariance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# varType A panels
png(file.path(out_dir, "TSPcov_fMA1_varTypeA_hong.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(file.path(base_sim, "varType_A"), Delta = DeltaA, factorVector = factorVector, q = q_hong)
abline(v = a0, lty = 2, lwd = 0.5)
dev.off()

png(file.path(out_dir, "TSPcov_fMA1_varTypeA.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(file.path(base_sim, "varType_A"), Delta = DeltaA, factorVector = factorVector, q = q_shao)
abline(v = a0, lty = 2, lwd = 0.5)
dev.off()

# varType B panels
png(file.path(out_dir, "TSPcov_fMA1_varTypeB_hong.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(file.path(base_sim, "varType_B"), Delta = DeltaB, factorVector = factorVector, q = q_hong)
abline(v = a0, lty = 2, lwd = 0.5)
dev.off()

png(file.path(out_dir, "TSPcov_fMA1_varTypeB.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(file.path(base_sim, "varType_B"), Delta = DeltaB, factorVector = factorVector, q = q_shao)
abline(v = a0, lty = 2, lwd = 0.5)
dev.off()

message("TSP covariance panels written to: ", normalizePath(out_dir))
