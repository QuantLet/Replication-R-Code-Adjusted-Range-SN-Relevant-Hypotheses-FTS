# R code for relevant functional mean changes

This repository contains the R programs for the paper *Relevant Mean Changes in Functional Time Series under Weak Directional Moments* by Zhuo Lin, Jiajing Sun, Wolfgang Karl Härdle, and Meiting Zhu.

The code implements adjusted-range and quadratic self-normalization and a same-block scalar Wald procedure. It covers the unknown-break simulations, nonuniform variance profiles, changing persistence, directional and orthogonal heavy tails, exact finite-grid Gaussian references, variance diagnostics, and the electricity-demand application. The empirical program downloads the public National Energy System Operator data and verifies the fixed file checksums before analysis.

## Requirements

Use R 4.2 or later and install:

```r
install.packages(c('digest','matrixStats','qrng','spacefillr','ggplot2','patchwork'))
```

## Run

Run commands from the repository root.

```sh
Rscript run_simulations.R
Rscript run_empirical.R
Rscript make_figures.R
Rscript validate.R
```

`run_simulations.R` reproduces the primary Monte Carlo study with 2,000 data samples per configuration and is computationally intensive. `run_empirical.R` downloads the public source files into `data/raw/` and writes calculated results to `results/`. `make_figures.R` writes the main simulation and empirical plots to `figures/`. Generated data, results, figures, and validation output are intentionally excluded from version control.

`run_all.R` executes all four steps in order.

## Code map

- `R/core.R`: functional statistics, split estimation, block profiles, Gaussian references, and bootstrap kernels.
- `R/simulation_helpers.R`: data-generating processes and common simulation utilities.
- `R/primary_simulations.R`: unknown-break, observed-loading, and ordinary-mean experiments.
- `R/size_adjusted_power.R`: independent size adjustment and paired uncertainty intervals.
- `R/variance_diagnostics.R`: block-variance bias diagnostics.
- `R/data_preparation.R`: public data download, checksum verification, and curve construction.
- `R/empirical_analysis.R`: level/shape decomposition and calendar sensitivity analysis.
- `R/figures.R`: color-blind-friendly figures with line-type distinctions.
- `R/validation.R`: numerical checks for segment rounding, reference construction, degeneracy handling, and scale behavior.
- `R/limiting_power.R`: Brownian limiting-power calculations.

The data source is the [NESO Historic Demand Data portal](https://www.neso.energy/data-portal/historic-demand-data).
