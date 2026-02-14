# BTC options IV smiles – relevant change-point empirical script (R)

This folder is a **self-contained** R script that:

1. Reads trade-level BTC options data (Deribit-style CSV with a JSON column).
2. Fetches BTC spot prices from the **public Binance Spot REST API** (`/api/v3/klines`).
3. Constructs a **daily constant-maturity implied-volatility (IV) smile** as a function of log-moneyness.
4. Runs a **relevant mean change-point test** in functional time series form:
   - Quadratic self-normalization (Dette et al.-style benchmark)
   - Adjusted-range self-normalization (our proposed normalizer)

No FPCA/KL truncation is used: the functional objects are kept on a fixed grid.

---

## Quick start

1. copy your option csv to `data/options_trades.csv`).
2. In R, set your working directory to the folder where this README.md is.
3. Run:

```r
source("run_empirical.R")
```

Outputs are written to:

- `output/figures/`
- `output/tables/`

---

## Environment variables

- `OPTIONS_CSV` – full path to options CSV
- `BINANCE_CACHE` – where to cache Binance daily data
- `TAU_TARGET_DAYS` – constant maturity target (default 30)
- `TAU_WINDOW_DAYS` – selection window (default ±7)
- `EPS_TRIM` – trimming for change-point estimation (default 0.05)
- `U_MIN`, `U_MAX` – log-moneyness range (default -0.35, 0.35)
- `U_GRID_N` – number of grid points (default 61)
- `MIN_TRADES_DAY`, `MIN_BINS_DAY` – minimum liquidity thresholds
- `QSIM_N`, `QSIM_GRID_N` – Monte Carlo settings for pivotal quantiles (cached in `data/qvalues_cpp.rds`)

Example:

```r
Sys.setenv(OPTIONS_CSV = "/path/to/part-00000-...csv")
Sys.setenv(TAU_TARGET_DAYS = "30")
Sys.setenv(TAU_WINDOW_DAYS = "10")
source("run_empirical.R")
```

---

## Main outputs

**Tables**

- `output/tables/btc_iv_boundary_summary.csv`
  - estimated break date
  - statistic \(\widehat{\mathbb{D}}_N^{cp}\)
  - normalizers (quadratic vs adjusted-range)
  - implied relevance boundaries \(\widehat\Delta_{\alpha}\)
  - RMS boundaries in **vol points** (e.g., 3 vol points = 0.03)

- `output/tables/btc_iv_boundary.csv`
  - grid of \(\Delta\) values with reject/accept indicators for each \(\alpha\)

- `output/tables/btc_daily_iv_smiles.csv`
  - the **functional time series** itself: one row per day, columns are IV values on the fixed log-moneyness grid

**Figures**

- `output/figures/btc_iv_mean_pre.png`
- `output/figures/btc_iv_mean_post.png`
- `output/figures/btc_iv_means_overlay.png`
- `output/figures/btc_spot_with_cp.png`

---

## (Research) Trading strategy driven by relevant change

Because our empirical section is about **relevant change**, this project also includes
a daily-executable **research backtest** whose *gate* is the same relevant change-point
toolbox (quadratic SN and adjusted-range SN).

The code is intentionally simple and uses only daily data. It is designed for the
paper/backtest and is **not** production trading infrastructure.

### High-level idea

1. Treat the daily IV-smile as a functional time series `X_t(u)` on the fixed grid.
2. Run an **online rolling-window relevant CP test** on `X_t(u)` (no look-ahead).
3. If a relevant change is detected, the strategy enters a short **cooldown** period
   (stays flat) to avoid trading through structural breaks.
4. In stable periods (no recent relevant change), it trades a **delta-hedged near-ATM
   straddle** using a simple volatility-risk-premium (VRP) filter:
   - VRP proxy: implied ATM variance minus an EWMA forecast variance from spot returns.
   - Short vol when VRP is high and skew is not extreme.
   - Optionally long vol when VRP is very low.

### Run it

```r
source("run_empirical.R")                 # optional (builds smiles + tests)
source("run_strategy_relevant_trading.R") # strategy + P&L log
```

Outputs are written to `output/strategy/`:

- `relevant_vrp_positions.csv`
- `relevant_vrp_pnl.csv`
- `relevant_vrp_cum_pnl.png`
- `relevant_vrp_summary.txt`

### Useful environment variables

- `SMILES_CSV`  path to an existing daily-smiles CSV (defaults to `output/tables/btc_daily_iv_smiles.csv`)
- `DET_WINDOW`  rolling window length for the relevant-change detector
- `DET_DELTA_RMS` relevance threshold in IV units (e.g., 0.05 = 5 vol points)
- `DET_ALPHA` and `DET_REQUIRE_BOTH` (whether to require both tests to reject)
- `COOLDOWN_DAYS` flat days after a detected relevant change
- `EWMA_LAMBDA` and `VRP_WINDOW` for the VRP layer
- `STRIKE_STEP`, `TC_OPT`, `TC_SPOT` for execution assumptions

