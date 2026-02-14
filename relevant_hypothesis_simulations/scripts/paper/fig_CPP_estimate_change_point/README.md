# Figure: Change-point estimator histograms (`fig:CPP_estimate_change_point`)

This folder reproduces Figure `fig:CPP_estimate_change_point` from the manuscript.

The change-point estimator is the same in both methods (it depends on the argmax of a norm),
so the histograms can be produced from either the `shao_*.csv` or `hong_*.csv` tables.
This script uses the baseline (`shao_*.csv`) tables.

## Prerequisite

Run the IID change-point simulations first:

```r
source("scripts/paper/fig_CPP_rejection_probabilities/run_simulations_CPP_fIID.R")
```

## Then build the histograms

```r
source("scripts/paper/fig_CPP_estimate_change_point/make_figure_CPP_cp_histograms.R")
```
