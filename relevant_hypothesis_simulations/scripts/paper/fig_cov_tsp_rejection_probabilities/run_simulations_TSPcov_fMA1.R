rm(list = ls())

# Two-sample covariance operator simulations for Figure fig:cov_tsp_rejection_probabilities
# Assumes working directory is the repository root.

library(fda)
library(parallel)

source(file.path("R", "genData.R"))
source(file.path("R", "TSPcovMethods.R"))

ensure_dir <- function(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# Grid of covariance scaling factors
factorVector <- seq(from = 1, to = 2.5, by = 0.025)

simInfo <- list(
  nsim = 1000,
  cores = max(1, detectCores() - 1)
)

base_out <- file.path("output", "simulations", "TSPcov", "fMA1")

run_one <- function(varType_tag) {
  dataInfo <- list(
    n = 100, m = 100,
    type = "fMA1", varType = varType_tag, basisType = "fourier",
    kappa = 0.7,
    factor = NULL,
    Psi = diag(21)
  )

  # n = m = 100
  out_dir <- file.path(base_out, paste0("varType_", varType_tag), "n100m100"); ensure_dir(out_dir)
  set.seed(101); writeStatistics(out_dir, factorVector, dataInfo, simInfo); writeStatistics_hong(out_dir, factorVector, dataInfo, simInfo)

  # n = m = 200
  dataInfo$n <- 200; dataInfo$m <- 200
  out_dir <- file.path(base_out, paste0("varType_", varType_tag), "n200m200"); ensure_dir(out_dir)
  set.seed(100); writeStatistics(out_dir, factorVector, dataInfo, simInfo); writeStatistics_hong(out_dir, factorVector, dataInfo, simInfo)

  # n = m = 500
  dataInfo$n <- 500; dataInfo$m <- 500
  out_dir <- file.path(base_out, paste0("varType_", varType_tag), "n500m500"); ensure_dir(out_dir)
  set.seed(99); writeStatistics(out_dir, factorVector, dataInfo, simInfo); writeStatistics_hong(out_dir, factorVector, dataInfo, simInfo)
}

run_one("A")
run_one("B")

message("TSPcov simulations finished. Raw tables written to: ", normalizePath(base_out))
