# Figure: One-sample rejection probabilities (`fig:OSPrejection_probabilities`)

This folder reproduces Figure `fig:OSPrejection_probabilities` from the manuscript.

It produces six panels (as separate PNGs), comparing:

- **Adjusted-range** self-normalization (files ending in `_hong`)
- **Quadratic** self-normalization (baseline; Dette et al., 2020)

## Files

- `run_simulations_OSP.R`  
  Generates raw simulation tables under `output/simulations/OSP/…`

- `make_figure_OSP.R`  
  Reads the raw tables and writes figure panels under `output/figures/OSP/`

## Workflow

```r
source("scripts/paper/fig_OSP_rejection_probabilities/run_simulations_OSP.R")
source("scripts/paper/fig_OSP_rejection_probabilities/make_figure_OSP.R")
```
