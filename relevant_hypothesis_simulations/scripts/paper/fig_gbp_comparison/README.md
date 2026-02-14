# Figure: GBP/USD before vs after Brexit (`fig:gbp-comparison`)

This folder contains code to reproduce the two subfigures in Figure `fig:gbp-comparison`:

- `gbp-usd-before`
- `gbp-usd-after`

The original project used a proprietary 1‑minute GBP/USD dataset (Refinitiv / Eikon).
This clean repository ships **no data**.

## What you need

Provide a CSV file with at least:

- `DateTime` (parseable timestamp, e.g. `YYYY-mm-dd HH:MM:SS`)
- `Close` (GBP/USD close price)

Optionally you can also provide `Date` and `DayTime` (0–1), but the script can compute them from `DateTime`.

## Run

1. Edit `make_figure_gbp_usd.R` and set `csv_path`.
2. Run:

```r
source("scripts/paper/fig_gbp_comparison/make_figure_gbp_usd.R")
```

Figures are written to `output/figures/`:
- `gbp-usd-before.png`
- `gbp-usd-after.png`
