rm(list = ls())

# One-sample simulations for Figure fig:OSPrejection_probabilities
# Assumes working directory is the repository root.

library(fda)
library(parallel)

source(file.path("R", "genData.R"))
source(file.path("R", "OSPmethods.R"))

# -----------------------------
# Simulation design (matches the manuscript text)
# -----------------------------
Delta <- 0.02
deltaVector <- 10:40/1000
aVector <- sqrt(2 * deltaVector)

simInfo <- list(
  nsim = 1000,
  cores = max(1, detectCores() - 1)
)

# Helper to create output dirs
ensure_dir <- function(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)

base_out <- file.path("output", "simulations", "OSP")

# -----------------------------
# Scenario 1: fIID errors
# -----------------------------
scenario <- "fIID"
dataInfo <- list(
  n = 25, type = "fIID", varType = "A", basisType = "bspline",
  muInfo = list(a = NULL, type = "sin")
)

for (n_val in c(25, 50, 100)) {
  dataInfo$n <- n_val
  out_dir <- file.path(base_out, scenario, paste0("n", n_val))
  ensure_dir(out_dir)

  # Seeds chosen to match the original project scripts
  seed <- switch(as.character(n_val), "25" = 161, "50" = 418, "100" = 68)
  set.seed(seed)
  writeStatistics(out_dir, aVector, dataInfo, simInfo)
  writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)
}

# -----------------------------
# Scenario 2: fMA(1) errors
# -----------------------------
scenario <- "fMA1"
dataInfo <- list(
  n = 25, type = "fMA1", varType = "A", basisType = "bspline",
  kappa = 0.7,
  muInfo = list(a = NULL, type = "sin")
)

for (n_val in c(25, 50, 100)) {
  dataInfo$n <- n_val
  out_dir <- file.path(base_out, scenario, paste0("n", n_val))
  ensure_dir(out_dir)

  seed <- switch(as.character(n_val), "25" = 161, "50" = 418, "100" = 68)
  set.seed(seed)
  writeStatistics(out_dir, aVector, dataInfo, simInfo)
  writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)
}

# -----------------------------
# Scenario 3: Brownian bridge errors
# -----------------------------
scenario <- "BB"
dataInfo <- list(
  n = 25, type = "BB", basisType = "bspline", nArgvals = 300,
  muInfo = list(a = NULL, type = "sin")
)

for (n_val in c(25, 50, 100)) {
  dataInfo$n <- n_val
  out_dir <- file.path(base_out, scenario, paste0("n", n_val))
  ensure_dir(out_dir)

  seed <- switch(as.character(n_val), "25" = 161, "50" = 418, "100" = 68)
  set.seed(seed)
  writeStatistics(out_dir, aVector, dataInfo, simInfo)
  writeStatistics_hong(out_dir, aVector, dataInfo, simInfo)
}

message("OSP simulations finished. Raw tables written to: ", normalizePath(base_out))
