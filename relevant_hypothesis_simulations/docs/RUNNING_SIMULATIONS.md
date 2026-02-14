# Running simulations and reproducing figures

All manuscript-facing scripts are under `scripts/paper/`.

Each figure folder usually contains two steps:

1. **Run simulations**  
   This writes raw tables (`shao_*.csv`, `hong_*.csv`) under `output/simulations/...`.

2. **Make figures**  
   This reads the raw tables and writes figures under `output/figures/...`.

## Working directory

All scripts assume your working directory is the repository root (the directory containing `README.md`).

In RStudio you can do:

```r
setwd(".../relevant_hypothesis_clean")
```

(or create an `.Rproj` file and open it).

## Output folders

Scripts create (if missing):

- `output/simulations/…`
- `output/figures/…`

To remove all generated files:

```r
unlink("output", recursive = TRUE)
```

## Typical workflow (example: one-sample figure)

```r
# Step 1: run simulations
source("scripts/paper/fig_OSP_rejection_probabilities/run_simulations_OSP.R")

# Step 2: build the figure panels
source("scripts/paper/fig_OSP_rejection_probabilities/make_figure_OSP.R")
```

Repeat similarly for other figure folders.
