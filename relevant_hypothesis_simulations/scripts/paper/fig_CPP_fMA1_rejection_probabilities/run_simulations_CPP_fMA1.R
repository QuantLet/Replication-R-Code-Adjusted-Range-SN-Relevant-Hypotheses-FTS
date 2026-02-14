rm(list = ls())

# Change-point simulations with fMA(1) errors for Figure fig:CPPfMA1_rejection_probabilities
# Assumes working directory is the repository root.

library(fda)
library(parallel)

source(file.path("R", "genData.R"))
source(file.path("R", "CPPmethods.R"))

ensure_dir <- function(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)

aVector <- 0:25/50

simInfo <- list(
  nsim = 1000,
  cores = max(1, detectCores() - 1),
  eps = 0.05
)

base_out <- file.path("output", "simulations", "CPP", "fMA1")

run_one <- function(tag, factor_val) {
  dataInfo <- list(
    n = 200,
    type = "fMA1", varType = "A", basisType = "bspline",
    kappa = 0.7,
    s.star = 0.5,
    factor = factor_val,
    muInfo = list(a = NULL, type = "horv")
  )

  # N = 200
  out_dir <- file.path(base_out, tag, "n200"); ensure_dir(out_dir)
  set.seed(193)
  writeStatistics(out_dir, aVector, dataInfo, simInfo)
  writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)

  # N = 500
  dataInfo$n <- 500
  out_dir <- file.path(base_out, tag, "n500"); ensure_dir(out_dir)
  set.seed(830)
  writeStatistics(out_dir, aVector, dataInfo, simInfo)
  writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)
}

# c = 1  (factor = sqrt(1))
run_one("c1", 1)

# c = 3  (factor = sqrt(3))
run_one("c3", sqrt(3))

message("CPP (fMA1) simulations finished. Raw tables written to: ", normalizePath(base_out))
