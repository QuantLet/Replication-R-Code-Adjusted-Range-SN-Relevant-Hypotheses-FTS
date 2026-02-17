# INFORMATION.md — BTC options IV smiles (relevant change-point empirical analysis)

This folder reproduces the **BTC options empirical analysis** in the paper:

1. **Construct daily 30-day constant-maturity BTC IV smiles** from executed Deribit option trades.
2. Apply the **rolling-window relevant mean change-point test** (quadratic SN vs adjusted-range SN).
3. Summarize **break dates**, **break direction** (up/down), and the **implied relevance boundary** (in RMS IV “vol points”).
4. Use the same rolling monitor as an **event detector** to drive a simple **break-direction, delta-hedged straddle** illustration.

> Units: Deribit trade IV is typically in **percent**; the code converts to **decimal IV** (`iv_dec = iv/100`).
> “Vol points” in the manuscript correspond to `100 × (decimal IV)`.

---

## Project structure

- `run_empirical_rolling.R`  
  Builds smiles (if needed) and runs the **rolling relevant change-point monitor**.

- `run_strategy_break_direction_straddle.R`  
  Uses the rolling monitor as an event detector and backtests a **break-direction delta-hedged straddle** (paper illustration).

- `run_strategy_relevant_trading.R` *(optional / research)*  
  A regime-gated **VRP straddle** strategy (not required for the paper’s main trading illustration).

- `R/`
  - `packages.R` — installs/loads required R packages
  - `parse_deribit_trades.R` — reads Deribit-style trade CSV/ZIP/dir and parses `instrument_name`
  - `binance_spot.R` — fetches daily BTCUSDT close from Binance (with caching; has a fallback)
  - `build_iv_smile.R` — constructs daily constant-maturity smiles on a fixed log-moneyness grid
  - `relevant_cpp_tests.R` — rolling-window relevant CP test (quadratic SN and adjusted-range SN) + pivotal quantile simulation/cache
  - `strategy_relevant_trading.R` — helpers for features, detector, and delta-hedged straddle backtest

---

## Data requirements

### 1) Deribit executed trades (required)

Input can be:
- a single CSV file,
- a directory containing multiple CSV/ZIP files,
- a ZIP file containing one or more trade CSVs,
- or a semicolon/comma-separated list of such paths.

**Expected columns (minimum):**
- `t` — millisecond timestamp
- `trade` — JSON string containing (at least)  
  `instrument_name`, `iv`, `amount`, `index_price`

`instrument_name` should look like: `BTC-30JAN26-80000-P`.

Put your file(s) under `data/` (recommended), e.g.
- `data/options_trades.csv`
- `data/part-00000-...csv.zip`
- `data/trades/` (directory)

### 2) BTC spot close (default via Binance; cached)

The code fetches daily **BTCUSDT** close via a public Binance endpoint and caches to:
- `data/btc_spot_binance.csv`

If Binance is blocked/unavailable, the code falls back to daily **median** `index_price` from Deribit trades (where applicable).

---

## Empirical pipeline (matches the manuscript)

### Step A — Build daily constant-maturity smiles

For each day `t` with sufficient liquidity:

1. **Select constant maturity**  
   Keep trades with time-to-maturity in  
   `tau_target_days ± tau_window_days` (default 30 ± 7).

2. **Compute log-moneyness**  
   `u = log(K / S_t)` where `S_t` is daily spot close.

3. **OTM filter (default)**  
   Use puts for `u<0` and calls for `u>0` to reduce microstructure noise.

4. **Bin + amount-weighted IV**  
   Bin trades to a fixed `u_grid` (default `[-0.35, 0.35]` with `61` points).  
   Within each bin: amount-weighted mean IV.

5. **Smooth within the day**  
   Fit a weighted smoothing spline across `u`, then evaluate on the fixed grid.

6. **Liquidity filters**  
   Drop days with fewer than `MIN_TRADES_DAY` trades (default 50) or
   fewer than `MIN_BINS_DAY` non-empty bins (default 10).

The resulting functional time series is the matrix `curves_mat` with rows = days and columns = grid points.

Cached output (if created):
- `output/tables/btc_daily_iv_smiles.csv`  
  (`date` + grid columns named `u_m0.350`, `u_p0.012`, etc.)

---

### Step B — Rolling-window relevant mean change-point test

On each rolling window of length `W` curves (paper: `W=30`), the code:

1. Estimates a candidate break location `k_hat` via a discretized functional objective.
2. Constructs the subsample path (lambda-grid) and computes:
   - quadratic SN normalizer
   - adjusted-range SN normalizer
3. Applies the **relevant** decision rule:

\[
\widehat{\mathbb{D}} > \Delta + q_{1-\alpha}\,\widehat{V},
\qquad
\Delta = \Delta_{\mathrm{rms}}^2,
\]

where `Delta_rms` is in **decimal IV units** (paper: `0.01` = 1 vol point RMS).

It also reports the **implied relevance boundary** per window (largest \(\Delta\) rejected):
\[
\widehat{\Delta}_\alpha = \max\{0,\ \widehat{\mathbb{D}} - q_{1-\alpha}\widehat{V}\},
\]
and converts to RMS vol points as `100*sqrt(Delta_hat)`.

---

### Step C — Break-direction delta-hedged straddle illustration (paper)

Using the rolling monitor as an **event detector**:

- When a window rejects, label direction as:
  - `up` if the average post-break mean smile is above the pre-break mean smile,
  - `down` otherwise.

Trading rule:
- enter a near-ATM constant-maturity **straddle** (1 call + 1 put) for `hold_days` days
- **long** straddle if direction is `up` (IV mark-up),
- **short** straddle if direction is `down` (IV mark-down),
- delta-hedge at entry (daily steps treat each smile date as a “trading day roll”),
- apply simple transaction costs (`tc_opt`, `tc_spot`).

Outputs include the cumulative P&L (USD per 1 straddle unit).

---

## How to run

### 0) One-time: open R in the project root
Make sure your working directory is the folder that contains `run_empirical_rolling.R`, `R/`, etc.

> **Note:** If `run_empirical_rolling.R` contains a hard-coded `setwd("...")`, delete/comment it or change it to your local project root.

### 1) Run the rolling empirical analysis (paper settings)

In R:

```r
Sys.setenv(OPTIONS_CSV  = "data/options_trades.csv")  # change to your path
Sys.setenv(SAMPLE_START = "2025-08-01")
Sys.setenv(SAMPLE_END   = "2026-01-30")

Sys.setenv(ROLL_WINDOW   = "30")
Sys.setenv(EPS_TRIM      = "0.10")    # trimming in CP estimation (paper)
Sys.setenv(ROLL_ALPHA    = "0.10")    # 10% level
Sys.setenv(ROLL_DELTA_RMS = "0.01")   # 1 vol point RMS

source("run_empirical_rolling.R")
```

Main outputs:
- `output/tables/btc_daily_iv_smiles.csv`
- `output/tables/btc_iv_rolling_cpp_W30_alpha0.10.csv`
- `output/tables/btc_iv_rolling_cpp_hits_W30_alpha0.10.csv` (frequency of cp_date among rejected windows)
- `output/figures/btc_iv_rolling_boundary_ar_W30_alpha0.10.png`

**How to recover the reported break dates** (e.g., 2025-10-08, 2025-12-03):  
look at `btc_iv_rolling_cpp_hits_*.csv` and take the most frequent `cp_date` values among rejected windows.

### 2) Optional: VRP regime-gated strategy (research)

```r
source("run_strategy_relevant_trading.R")
```

Outputs:
- `output/strategy/relevant_vrp_positions.csv`
- `output/strategy/relevant_vrp_pnl.csv`
- `output/strategy/relevant_vrp_cum_pnl.png`
- `output/strategy/relevant_vrp_summary.txt`

---

## Key environment variables

### Inputs / caching
- `OPTIONS_CSV` — Deribit trades CSV/ZIP/dir (or list separated by `;` or `,`)
- `SMILES_CSV` — cached smiles (default `output/tables/btc_daily_iv_smiles.csv`)
- `BINANCE_CACHE` — cached spot file (default `data/btc_spot_binance.csv`)
- `FORCE_REBUILD_SMILES` — `1` to rebuild smiles even if `SMILES_CSV` exists

### Sample window (recommended for the manuscript)
- `SAMPLE_START`, `SAMPLE_END` — e.g. `2025-08-01`, `2026-01-30`
- `SAMPLE_DAYS` — alternative: length counting back from `SAMPLE_END`

### Smile construction
- `TAU_TARGET_DAYS` (default 30), `TAU_WINDOW_DAYS` (default 7)
- `U_MIN`, `U_MAX` (default -0.35, 0.35), `U_GRID_N` (default 61)
- `MIN_TRADES_DAY` (default 50), `MIN_BINS_DAY` (default 10)

### Rolling relevant CP test (paper)
- `ROLL_WINDOW` (paper 30), `ROLL_STEP` (default 1)
- `EPS_TRIM` (paper 0.10)
- `ROLL_ALPHA` (paper 0.10)
- `ROLL_DELTA_RMS` (paper 0.01 = 1 vol point RMS)
- `ROLL_REQUIRE_BOTH` (0 = either normalizer rejects; 1 = both must reject)

### Pivotal quantiles cache (Brownian functional simulation)
- `Q_CACHE` (default `data/qvalues_cpp.rds`)
- `QSIM_N` (default 5000), `QSIM_GRID_N` (default 1000)

---

## Reproducibility notes

- Internet access is only needed for Binance spot fetch; cached data avoids repeated calls.
- Quantile simulation is cached in `data/qvalues_cpp.rds`; delete it if you want to re-simulate.
- The scripts are **research code** designed to match the manuscript’s empirical section; they do not model margin/funding/intraday hedging, and the P&L is “model-implied” per the paper’s description.
