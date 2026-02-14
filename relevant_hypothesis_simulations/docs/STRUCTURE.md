# Repository structure

This repository is split into three layers:

1. **Core R code (`R/`)**  
   Reusable functions implementing:
   - data generating processes used in simulations,
   - test statistics and self-normalizers,
   - Monte‑Carlo critical value utilities.

2. **Paper scripts (`scripts/paper/`)**  
   Reproducibility scripts grouped by manuscript figure.
   Each figure folder contains:
   - a simulation script (generates raw result tables under `output/simulations/…`)
   - a plotting script (generates figure files under `output/figures/…`)
   - a short README describing what the folder corresponds to in the manuscript.

## Output policy

This archive intentionally includes **no outputs**:

- No precomputed simulation tables (`*.csv`)
- No pre-rendered figures (`*.png`, `*.pdf`, …)
- No workspace artifacts (`*.RData`, `*.Rhistory`)

All scripts write outputs to:

- `output/simulations/…` (raw simulation tables)
- `output/figures/…` (figures for the manuscript)

