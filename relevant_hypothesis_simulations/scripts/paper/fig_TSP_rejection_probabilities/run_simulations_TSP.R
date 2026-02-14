rm(list = ls())

# Two-sample simulations for Figure fig:TSPrejection_probabilities
# Assumes working directory is the repository root.

library(fda)
library(parallel)

source(file.path("R", "genData.R"))
source(file.path("R", "TSPmethods.R"))

ensure_dir <- function(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)

base_out <- file.path("output", "simulations", "TSP")

simInfo <- list(
  nsim = 1000,
  cores = max(1, detectCores() - 1),
  normalizerType = "V_n"
)

# Helper to run the four sample-size configurations
run_sizes <- function(base_path, dataInfo, aVector) {
  # n = m = 50
  dataInfo$n <- 50; dataInfo$m <- 50
  out_dir <- file.path(base_path, "n50m50"); ensure_dir(out_dir)
  set.seed(193); writeStatistics(out_dir, aVector, dataInfo, simInfo); writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)

  # n = 50, m = 100
  dataInfo$n <- 50; dataInfo$m <- 100
  out_dir <- file.path(base_path, "n50m100"); ensure_dir(out_dir)
  set.seed(830); writeStatistics(out_dir, aVector, dataInfo, simInfo); writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)

  # n = m = 100
  dataInfo$n <- 100; dataInfo$m <- 100
  out_dir <- file.path(base_path, "n100m100"); ensure_dir(out_dir)
  set.seed(515); writeStatistics(out_dir, aVector, dataInfo, simInfo); writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)

  # n = 100, m = 200
  dataInfo$n <- 100; dataInfo$m <- 200
  out_dir <- file.path(base_path, "n100m200"); ensure_dir(out_dir)
  set.seed(515); writeStatistics(out_dir, aVector, dataInfo, simInfo); writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)
}

# -----------------------------
# Scenario 1: fIID errors
# -----------------------------
scenario <- "fIID"
aVector <- 0:50/100

dataInfo <- list(
  n = 50, m = 50,
  type = "fIID", varType = "A", basisType = "bspline",
  muInfo = list(a = NULL, type = "horv")
)

run_sizes(file.path(base_out, scenario), dataInfo, aVector)

# -----------------------------
# Scenario 2: fMA(1) errors (independent samples)
# -----------------------------
scenario <- "fMA1_independent"
aVector <- 0:50/100

dataInfo <- list(
  n = 50, m = 50,
  type = "fMA1", varType = "A", basisType = "bspline",
  kappa = 0.7,
  muInfo = list(a = NULL, type = "horv")
)

run_sizes(file.path(base_out, scenario), dataInfo, aVector)

# -----------------------------
# Scenario 3: Brownian bridge errors
# -----------------------------
scenario <- "BB"
aVector <- 0:(14-1)/10   # matches the original script

dataInfo <- list(
  n = 50, m = 50,
  type = "BB", basisType = "bspline", nArgvals = 300,
  muInfo = list(a = NULL, type = "horv")
)

run_sizes(file.path(base_out, scenario), dataInfo, aVector)

message("TSP simulations finished. Raw tables written to: ", normalizePath(base_out))
