 

This folder is a **clean, code-only** reorganization of the original project archive `relevant_hypothesis-20250423`.

It contains:

- The **core R implementations** of
  - the original (quadratic) self-normalization approach from Dette et al. (2020), and
  - the **adjusted-range (range-based) self-normalization** modifications used in the attached manuscript.
- **Paper-mapped scripts** grouped by figure (simulation drivers + plotting scripts).
- A `legacy/` snapshot that keeps the original working tree structure **without any outputs** (no precomputed `.csv` results, no `.png` figures, no `.RData`, etc.).

> **Important:** this repository intentionally contains **no simulation outputs** and **no pre-rendered figures**.
> Running the simulation/figure scripts will create outputs under `output/`, which is not included in this archive.

## Quick start

1. Open R/RStudio and set your working directory to the repository root (the folder that contains this README).

2. Install required packages (at minimum):

```r
install.packages(c("fda", "parallel"))
```

Some scripts also use optional packages (see the per-figure READMEs).

3. Reproduce simulations + figures from the manuscript:

- Scripts are under `scripts/paper/…`
- Each figure folder contains:
  - `run_simulations*.R` (writes raw simulation tables to `output/simulations/…`)
  - `make_figure*.R` (reads from `output/simulations/…` and writes figures to `output/figures/…`)

## Folder layout

- `R/`  
  Core functions (data generators, test statistics, self-normalizers, helpers).  
  These files are the “library” you source from scripts.

- `scripts/paper/`  
  **Main entry points** organized by manuscript figures (recommended).

- `legacy/`  
  Original code tree (kept for reference), cleaned to be code-only.

- `docs/`  
  High-level documentation: structure, mapping to the manuscript, and notes.

## Notes on reproducibility

- Many simulation scripts are computationally expensive (1000 replications × multiple parameter grids).
- Output paths are created automatically under `output/`.

See `docs/PAPER_MAPPING.md` for the detailed mapping from manuscript figures to scripts.
