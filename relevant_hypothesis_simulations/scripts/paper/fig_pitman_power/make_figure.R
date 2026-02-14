rm(list = ls())

# Pitman local power curves (m = 20)
# This script is self-contained but assumes your working directory is the repository root.

library(fda)

source(file.path("R", "genData.R"))  # provides BM(...)
# We do not use quant.fct()/W() here; we compute V and H directly for speed/clarity.

# -----------------------------
# User-tunable parameters
# -----------------------------
alpha <- 0.05
m <- 20              # number of grid points used in the normalizer (matches the manuscript figure label "m=20")
N_BM <- 300          # Brownian motion discretization (grid size inside BM())
n_mc <- 20000        # Monte‑Carlo size (increase for smoother curves; e.g. 50000)
kappa_grid <- seq(0, 6, by = 0.1)  # kappa = c/tau in the manuscript

set.seed(1)

# -----------------------------
# Simulate Brownian motion sample
# -----------------------------
BMsample <- BM(N = N_BM, n = n_mc)

# Evaluate B(t) at the m grid points and at t=1
lambda <- 1:m / m
B_lambda <- eval.fd(lambda, BMsample)          # (m x n_mc) matrix
B_1 <- as.numeric(eval.fd(1, BMsample))        # length n_mc vector

# Build the "bridge-like" values: lambda * (B(lambda) - lambda * B(1))
values_mat <- matrix(0, nrow = m, ncol = n_mc)
for (k in 1:m) {
  values_mat[k, ] <- lambda[k] * (B_lambda[k, ] - lambda[k] * B_1)
}

# Quadratic normalizer (baseline)
V_vec <- sqrt(colSums(values_mat^2) / (m - 1))

# Adjusted-range normalizer (proposed)
H_vec <- apply(values_mat, 2, function(x) max(x) - min(x))

# Limit statistics under H0
S_quadratic <- B_1 / V_vec
S_range <- B_1 / H_vec

# Critical values (Monte‑Carlo)
q_quadratic <- as.numeric(quantile(S_quadratic, probs = 1 - alpha, na.rm = TRUE))
q_range <- as.numeric(quantile(S_range, probs = 1 - alpha, na.rm = TRUE))

# Local power curves: P( (B(1)+kappa)/normalizer > q )
power_quadratic <- sapply(kappa_grid, function(kappa) mean((B_1 + kappa) / V_vec > q_quadratic, na.rm = TRUE))
power_range <- sapply(kappa_grid, function(kappa) mean((B_1 + kappa) / H_vec > q_range, na.rm = TRUE))

# -----------------------------
# Plot
# -----------------------------
out_dir <- file.path("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png(file.path(out_dir, "pitman_power_curves_m20.png"),
    width = 5.65, height = 5.25, units = "cm", res = 1200, pointsize = 4)

plot(kappa_grid, power_range, type = "l", ylim = c(0, 1),
     xlab = expression(kappa), ylab = "Local power",
     lty = 1, lwd = 0.5, bty = "n", axes = FALSE)

box(lwd = 0.5)
axis(side = 1, lwd = 0.5)
axis(side = 2, lwd = 0.5)

lines(kappa_grid, power_quadratic, lty = 2, lwd = 0.5)

legend("bottomright",
       legend = c("Adjusted-range (proposed)", "Quadratic (Dette et al., 2020)"),
       lty = c(1, 2), lwd = 0.5, cex = 1, box.lwd = 0.5)

dev.off()
