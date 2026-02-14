# Mapping from the manuscript to code

This file maps **figures/simulations in `manuscript.tex`** to the corresponding scripts in this repository.

> The manuscript uses `\includegraphics{figs/<name>}`.
> The scripts here write figures to `output/figures/…`.
> If you want LaTeX compilation to work unchanged, copy the generated figure files into your `figs/` folder (or change the LaTeX paths).

---

## Figure `\label{fig:pitman-power}`

- Script: `scripts/paper/fig_pitman_power/make_figure.R`
- Output (when run): `output/figures/pitman_power_curves_m20.png`

---

## Figure `\label{fig:OSPrejection_probabilities}`

Six subfigures in the manuscript:

- `OSPfIID_hong`, `OSPfIID`
- `OSPfMA1_hong`, `OSPfMA1`
- `OSPBB_hong`, `OSPBB`

Scripts:

- Simulations:
  - `scripts/paper/fig_OSP_rejection_probabilities/run_simulations_OSP.R`
- Plotting:
  - `scripts/paper/fig_OSP_rejection_probabilities/make_figure_OSP.R`

Outputs (when run):

- Raw tables: `output/simulations/OSP/<scenario>/…`
- Figures: `output/figures/OSP/*.png`

---

## Figure `\label{fig:TSPrejection_probabilities}`

Six subfigures in the manuscript:

- `TSPfIID_hong`, `TSPfIID`
- `TSPfMA1_independent_hong`, `TSPfMA1_independent`
- `TSPBB_hong`, `TSPBB`

Scripts:

- Simulations:
  - `scripts/paper/fig_TSP_rejection_probabilities/run_simulations_TSP.R`
- Plotting:
  - `scripts/paper/fig_TSP_rejection_probabilities/make_figure_TSP.R`

---

## Figure `\label{fig:CPP_rejection_probabilities}`

Six subfigures in the manuscript:

- `a01power_hong`, `a01power`
- `a02power_hong`, `a02power`
- `a03power_hong`, `a03power`

Scripts:

- Simulations (IID errors):
  - `scripts/paper/fig_CPP_rejection_probabilities/run_simulations_CPP_fIID.R`
- Plotting:
  - `scripts/paper/fig_CPP_rejection_probabilities/make_figure_CPP_power_curves.R`

---

## Figure `\label{fig:CPP_estimate_change_point}`

Three histograms in the manuscript:

- `a01CPest`, `a02CPest`, `a03CPest`

Scripts:

- Uses the same simulations as Figure `fig:CPP_rejection_probabilities`
- Plotting:
  - `scripts/paper/fig_CPP_estimate_change_point/make_figure_CPP_cp_histograms.R`

---

## Figure `\label{fig:CPPfMA1_rejection_probabilities}`

Four subfigures in the manuscript:

- `CPPfMA1_c1_hong`, `CPPfMA1_c1`
- `CPPfMA1_c3_hong`, `CPPfMA1_c3`

Scripts:

- Simulations (fMA(1) errors):
  - `scripts/paper/fig_CPP_fMA1_rejection_probabilities/run_simulations_CPP_fMA1.R`
- Plotting:
  - `scripts/paper/fig_CPP_fMA1_rejection_probabilities/make_figure_CPP_fMA1_power_curves.R`

---

## Figure `\label{fig:cov_tsp_rejection_probabilities}`

Four subfigures in the manuscript:

- `TSPcov_fMA1_varTypeA_hong`, `TSPcov_fMA1_varTypeA`
- `TSPcov_fMA1_varTypeB_hong`, `TSPcov_fMA1_varTypeB`

Scripts:

- Simulations:
  - `scripts/paper/fig_cov_tsp_rejection_probabilities/run_simulations_TSPcov_fMA1.R`
- Plotting:
  - `scripts/paper/fig_cov_tsp_rejection_probabilities/make_figure_TSPcov_fMA1.R`

---

## Figure `\label{fig:cov_cpp_rejection_probabilities}`

Four subfigures in the manuscript:

- `CPPcov_fMA1_varTypeA_hong`, `CPPcov_fMA1_varTypeA`
- `CPPcov_fMA1_varTypeB_hong`, `CPPcov_fMA1_varTypeB`

Scripts:

- Simulations:
  - `scripts/paper/fig_cov_cpp_rejection_probabilities/run_simulations_CPPcov_fMA1.R`
- Plotting:
  - `scripts/paper/fig_cov_cpp_rejection_probabilities/make_figure_CPPcov_fMA1.R`

---

## Figure `\label{fig:emp-mean-functions}`

- Script: `scripts/paper/fig_emp_mean_functions/make_figure_temperature_means.R`
- Note: this script downloads data from NOAA GHCN via the `rnoaa` package (no data files are shipped here).

---

## Figure `\label{fig:gbp-comparison}`

- Script: `scripts/paper/fig_gbp_comparison/make_figure_gbp_usd.R`
- Note: the original project contained a local CSV. This clean repo ships **no data**.
  You must supply your own GBP/USD series (see script header).
