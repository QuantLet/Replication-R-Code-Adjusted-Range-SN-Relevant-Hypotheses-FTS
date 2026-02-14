# Paper reproduction scripts

This folder contains the manuscript-facing entry points, organized **by figure**.

Each subfolder typically contains:

- `run_simulations_*.R` — generates raw simulation result tables under `output/simulations/…`
- `make_figure_*.R` — reads those tables and writes the figure panel PNGs under `output/figures/…`
- `README.md` — explains which manuscript figure the folder corresponds to

See `docs/PAPER_MAPPING.md` for the full mapping.
