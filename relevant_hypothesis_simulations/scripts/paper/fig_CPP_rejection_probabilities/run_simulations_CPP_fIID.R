rm(list = ls())

# Change-point simulations (IID errors) for Figure fig:CPP_rejection_probabilities
# Assumes working directory is the repository root.

library(fda)
library(parallel)

source(file.path("R", "genData.R"))
source(file.path("R", "CPPmethods.R"))

ensure_dir <- function(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# Parameter grid for the change size a
aVector <- 0:25/50   # a = 0, 0.02, ..., 0.5

simInfo <- list(
  nsim = 1000,
  cores = max(1, detectCores() - 1),
  eps = 0.05
)

dataInfo <- list(
  n = 200,
  type = "fIID", varType = "A", basisType = "bspline",
  s.star = 0.5,
  factor = 1,
  muInfo = list(a = NULL, type = "horv")
)

base_out <- file.path("output", "simulations", "CPP", "fIID")

# N = 200
out_dir <- file.path(base_out, "n200"); ensure_dir(out_dir)
set.seed(193)
writeStatistics(out_dir, aVector, dataInfo, simInfo)
writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)

# N = 500
dataInfo$n <- 500
out_dir <- file.path(base_out, "n500"); ensure_dir(out_dir)
set.seed(830)
writeStatistics(out_dir, aVector, dataInfo, simInfo)
writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)

message("CPP (IID) simulations finished. Raw tables written to: ", normalizePath(base_out))
