rm(list = ls())

# Build the CPP covariance-operator panels (fig:cov_cpp_rejection_probabilities)
# Reads simulation tables from output/simulations/CPPcov/fMA1 and writes PNGs to output/figures/Covariance.

library(fda)

source(file.path("R", "genData.R"))
source(file.path("R", "generalMethods.R"))

plotPowerCurve <- function(base_path, Delta, factorVector, q) {
  ntables <- length(factorVector)
  N200list <- readTables(file.path(base_path, "N200"), ntables)
  N400list <- readTables(file.path(base_path, "N400"), ntables)
  N1000list <- readTables(file.path(base_path, "N1000"), ntables)

  plot(ylim = c(0, 1), x = factorVector, y = empRejProb(N200list, Delta = Delta, q = q),
       xlab = expression(a), ylab = "Empirical rejection probability",
       type = "l", lty = 1, font.lab = 2, lwd = 0.5, cex.lab = 1,
       bty = "n", axes = FALSE)

  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)

  lines(x = factorVector, y = empRejProb(N400list, Delta = Delta, q = q), lty = 2, lwd = 0.5)
  lines(x = factorVector, y = empRejProb(N1000list, Delta = Delta, q = q), lty = 3, lwd = 0.5)

  lines(x = c(factorVector[1], factorVector[ntables]), y = c(0.05, 0.05), lwd = 0.5)

  legend(1, 1, c("N = 200", "N = 400", "N = 1000"),
         lty = c(1, 2, 3), lwd = 0.5, cex = 1, box.lwd = 0.5)
}

plotPowerCurve_hong <- function(base_path, Delta, factorVector, q) {
  ntables <- length(factorVector)
  N200list <- readTables_hong(file.path(base_path, "N200"), ntables)
  N400list <- readTables_hong(file.path(base_path, "N400"), ntables)
  N1000list <- readTables_hong(file.path(base_path, "N1000"), ntables)

  plot(ylim = c(0, 1), x = factorVector, y = empRejProb(N200list, Delta = Delta, q = q),
       xlab = expression(a), ylab = "Empirical rejection probability",
       type = "l", lty = 1, font.lab = 2, lwd = 0.5, cex.lab = 1,
       bty = "n", axes = FALSE)

  box(lwd = 0.5)
  axis(side = 1, lwd = 0.5)
  axis(side = 2, lwd = 0.5)

  lines(x = factorVector, y = empRejProb(N400list, Delta = Delta, q = q), lty = 2, lwd = 0.5)
  lines(x = factorVector, y = empRejProb(N1000list, Delta = Delta, q = q), lty = 3, lwd = 0.5)

  lines(x = c(factorVector[1], factorVector[ntables]), y = c(0.05, 0.05), lwd = 0.5)

  legend(1, 1, c("N = 200", "N = 400", "N = 1000"),
         lty = c(1, 2, 3), lwd = 0.5, cex = 1, box.lwd = 0.5)
}

# Critical values
set.seed(2)
BMsample <- BM(N = 300, n = 1000)
q_shao <- quant.fct(BMsample, prob = 0.95)

set.seed(3)
BMsample2 <- BM(N = 300, n = 1000)
q_hong <- quant.fct_hong(BMsample2, prob = 0.95)

factorVector <- seq(from = 1, to = 2.5, by = 0.1)

# Manuscript settings
a0 <- 1.5
D <- 21

SigmaA <- 1 / seq(1, D)
DeltaA <- (1 - a0^2)^2 * sum(SigmaA^4) * (1 + 0.7^2)^2

SigmaB <- 1.2^(-seq(1, D))
DeltaB <- (1 - a0^2)^2 * sum(SigmaB^4) * (1 + 0.7^2)^2

base_sim <- file.path("output", "simulations", "CPPcov", "fMA1")
out_dir <- file.path("output", "figures", "Covariance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# varType A panels
png(file.path(out_dir, "CPPcov_fMA1_varTypeA_hong.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(file.path(base_sim, "varType_A"), Delta = DeltaA, factorVector = factorVector, q = q_hong)
abline(v = a0, lty = 2, lwd = 0.5)
dev.off()

png(file.path(out_dir, "CPPcov_fMA1_varTypeA.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(file.path(base_sim, "varType_A"), Delta = DeltaA, factorVector = factorVector, q = q_shao)
abline(v = a0, lty = 2, lwd = 0.5)
dev.off()

# varType B panels
png(file.path(out_dir, "CPPcov_fMA1_varTypeB_hong.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve_hong(file.path(base_sim, "varType_B"), Delta = DeltaB, factorVector = factorVector, q = q_hong)
abline(v = a0, lty = 2, lwd = 0.5)
dev.off()

png(file.path(out_dir, "CPPcov_fMA1_varTypeB.png"), width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)
plotPowerCurve(file.path(base_sim, "varType_B"), Delta = DeltaB, factorVector = factorVector, q = q_shao)
abline(v = a0, lty = 2, lwd = 0.5)
dev.off()

message("CPP covariance panels written to: ", normalizePath(out_dir))
