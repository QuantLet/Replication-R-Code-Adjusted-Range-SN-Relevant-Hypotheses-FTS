# Figure: Two-sample covariance-operator rejection probabilities (`fig:cov_tsp_rejection_probabilities`)

This folder reproduces Figure `fig:cov_tsp_rejection_probabilities` from the manuscript
(two-sample covariance operator comparison under fMA(1) errors, variance types A and B).

Panels:

- `TSPcov_fMA1_varTypeA_hong`, `TSPcov_fMA1_varTypeA`
- `TSPcov_fMA1_varTypeB_hong`, `TSPcov_fMA1_varTypeB`

## Files

- `run_simulations_TSPcov_fMA1.R`  
  Generates raw simulation tables under `output/simulations/TSPcov/fMA1/…`

- `make_figure_TSPcov_fMA1.R`  
  Reads the raw tables and writes figure panels under `output/figures/Covariance/`

## Workflow

```r
source("scripts/paper/fig_cov_tsp_rejection_probabilities/run_simulations_TSPcov_fMA1.R")
source("scripts/paper/fig_cov_tsp_rejection_probabilities/make_figure_TSPcov_fMA1.R")
```
