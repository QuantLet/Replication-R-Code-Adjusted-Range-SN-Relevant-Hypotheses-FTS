# Figure: Temperature mean functions (`fig:emp-mean-functions`)

This folder reproduces Figure `fig:emp-mean-functions` from the manuscript
(mean annual temperature curves for two American weather stations).

The script downloads daily **minimum** temperature data from NOAA GHCN-D via the `rnoaa` package,
projects each year onto a Fourier basis (49 basis functions), and plots the mean functions.

## Files

- `make_figure_temperature_means.R`

## Notes

- No data are shipped with this clean repository.
- Data download can take time because the time span is long (1940s–2024).

## Run

```r
source("scripts/paper/fig_emp_mean_functions/make_figure_temperature_means.R")
```
