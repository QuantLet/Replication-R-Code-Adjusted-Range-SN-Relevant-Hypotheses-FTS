# Figure: Two-sample rejection probabilities (`fig:TSPrejection_probabilities`)

This folder reproduces Figure `fig:TSPrejection_probabilities` from the manuscript.

It produces six panels (as separate PNGs), comparing:

- **Adjusted-range** self-normalization (files ending in `_hong`)
- **Quadratic** self-normalization (baseline; Dette et al., 2020)

## Files

- `run_simulations_TSP.R`  
  Generates raw simulation tables under `output/simulations/TSP/…`

- `make_figure_TSP.R`  
  Reads the raw tables and writes figure panels under `output/figures/TSP/`

## Workflow

```r
source("scripts/paper/fig_TSP_rejection_probabilities/run_simulations_TSP.R")
source("scripts/paper/fig_TSP_rejection_probabilities/make_figure_TSP.R")
```
