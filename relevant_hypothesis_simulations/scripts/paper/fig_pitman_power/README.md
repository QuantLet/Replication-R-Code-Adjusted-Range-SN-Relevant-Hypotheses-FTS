# Figure: Pitman local power curves (`fig:pitman-power`)

This folder contains code to reproduce the *Pitman local power* comparison figure in the manuscript.

- `make_figure.R` simulates Brownian-motion based limit objects and plots
  the local power curves for:
  - the adjusted-range (range-based) normalizer (proposed), and
  - the quadratic normalizer (baseline; Dette et al., 2020).

Outputs are written to `output/figures/`.
