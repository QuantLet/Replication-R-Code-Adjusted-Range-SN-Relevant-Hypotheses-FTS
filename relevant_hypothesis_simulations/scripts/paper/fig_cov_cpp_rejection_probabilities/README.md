# Figure: Change-point covariance-operator rejection probabilities (`fig:cov_cpp_rejection_probabilities`)

This folder reproduces Figure `fig:cov_cpp_rejection_probabilities` from the manuscript
(change-point in covariance operator under fMA(1) errors, variance types A and B).

Panels:

- `CPPcov_fMA1_varTypeA_hong`, `CPPcov_fMA1_varTypeA`
- `CPPcov_fMA1_varTypeB_hong`, `CPPcov_fMA1_varTypeB`

## Files

- `run_simulations_CPPcov_fMA1.R`  
  Generates raw simulation tables under `output/simulations/CPPcov/fMA1/…`

- `make_figure_CPPcov_fMA1.R`  
  Reads the raw tables and writes figure panels under `output/figures/Covariance/`

## Workflow

```r
source("scripts/paper/fig_cov_cpp_rejection_probabilities/run_simulations_CPPcov_fMA1.R")
source("scripts/paper/fig_cov_cpp_rejection_probabilities/make_figure_CPPcov_fMA1.R")
```
