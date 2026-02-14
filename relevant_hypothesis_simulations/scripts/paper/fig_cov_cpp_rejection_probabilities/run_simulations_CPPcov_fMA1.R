rm(list = ls())

# Change-point covariance operator simulations for Figure fig:cov_cpp_rejection_probabilities
# Assumes working directory is the repository root.

library(fda)
library(parallel)

source(file.path("R", "genData.R"))
source(file.path("R", "CPPcovMethods.R"))

ensure_dir <- function(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)

factorVector <- seq(from = 1, to = 2.5, by = 0.1)

simInfo <- list(
  nsim = 1000,
  cores = max(1, detectCores() - 1),
  eps = 0.1
)

base_out <- file.path("output", "simulations", "CPPcov", "fMA1")

run_one <- function(varType_tag) {
  dataInfo <- list(
    N = 200,
    type = "fMA1", varType = varType_tag, basisType = "fourier",
    kappa = 0.7,
    s.star = 0.5,
    factor = NULL,
    Psi = diag(21)
  )

  # N = 200
  out_dir <- file.path(base_out, paste0("varType_", varType_tag), "N200"); ensure_dir(out_dir)
  set.seed(101); writeStatistics(out_dir, factorVector, dataInfo, simInfo); writeStatistics_hong(out_dir, factorVector, dataInfo, simInfo)

  # N = 400
  dataInfo$N <- 400
  out_dir <- file.path(base_out, paste0("varType_", varType_tag), "N400"); ensure_dir(out_dir)
  set.seed(100); writeStatistics(out_dir, factorVector, dataInfo, simInfo); writeStatistics_hong(out_dir, factorVector, dataInfo, simInfo)

  # N = 1000
  dataInfo$N <- 1000
  out_dir <- file.path(base_out, paste0("varType_", varType_tag), "N1000"); ensure_dir(out_dir)
  set.seed(99); writeStatistics(out_dir, factorVector, dataInfo, simInfo); writeStatistics_hong(out_dir, factorVector, dataInfo, simInfo)
}

run_one("A")
run_one("B")

message("CPPcov simulations finished. Raw tables written to: ", normalizePath(base_out))
