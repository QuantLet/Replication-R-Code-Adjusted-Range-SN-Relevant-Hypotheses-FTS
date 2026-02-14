# Figure: Change-point rejection probabilities (`fig:CPP_rejection_probabilities`)

This folder reproduces Figure `fig:CPP_rejection_probabilities` from the manuscript
(change-point in the mean function with IID errors).

It produces six panels:

- `a01power_hong`, `a01power`  (Δ = 0.1^2 / 30)
- `a02power_hong`, `a02power`  (Δ = 0.2^2 / 30)
- `a03power_hong`, `a03power`  (Δ = 0.3^2 / 30)

## Files

- `run_simulations_CPP_fIID.R`  
  Generates raw simulation tables under `output/simulations/CPP/fIID/…`

- `make_figure_CPP_power_curves.R`  
  Reads the raw tables and writes figure panels under `output/figures/CPP/`

## Workflow

```r
source("scripts/paper/fig_CPP_rejection_probabilities/run_simulations_CPP_fIID.R")
source("scripts/paper/fig_CPP_rejection_probabilities/make_figure_CPP_power_curves.R")
```
