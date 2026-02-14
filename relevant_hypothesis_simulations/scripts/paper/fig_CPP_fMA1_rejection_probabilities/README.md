# Figure: Change-point power curves with fMA(1) errors (`fig:CPPfMA1_rejection_probabilities`)

This folder reproduces Figure `fig:CPPfMA1_rejection_probabilities` from the manuscript.

Panels:

- `CPPfMA1_c1_hong`, `CPPfMA1_c1`
- `CPPfMA1_c3_hong`, `CPPfMA1_c3`

## Files

- `run_simulations_CPP_fMA1.R`  
  Generates raw simulation tables under `output/simulations/CPP/fMA1/…`

- `make_figure_CPP_fMA1_power_curves.R`  
  Reads the raw tables and writes figure panels under `output/figures/CPP/`

## Workflow

```r
source("scripts/paper/fig_CPP_fMA1_rejection_probabilities/run_simulations_CPP_fMA1.R")
source("scripts/paper/fig_CPP_fMA1_rejection_probabilities/make_figure_CPP_fMA1_power_curves.R")
```
