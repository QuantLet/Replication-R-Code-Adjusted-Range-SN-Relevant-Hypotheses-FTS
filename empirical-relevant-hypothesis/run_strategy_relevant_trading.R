# BTC options IV smiles – trading strategy driven by *relevant change*
#
# This script is meant to be run AFTER (or instead of) run_empirical.R.
# It:
#   1) loads/constructs the daily functional IV smiles X_t(u)
#   2) runs an ONLINE rolling-window relevant change detector
#   3) uses that detector to gate a delta-hedged straddle VRP strategy
#   4) writes a trade/P&L log and a cumulative P&L plot
#
# Run:
#   source("run_strategy_relevant_trading.R")
#
# IMPORTANT: Research code for your paper; not production trading advice.

source("R/packages.R")
source("R/binance_spot.R")
source("R/parse_deribit_trades.R")
source("R/build_iv_smile.R")
source("R/relevant_cpp_tests.R")
source("R/strategy_relevant_trading.R")

dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create("output/strategy", recursive = TRUE, showWarnings = FALSE)

# -----------------------
# Configuration (env vars)
# -----------------------

config <- list(
  # inputs
  options_csv   = Sys.getenv("OPTIONS_CSV", "data/options_trades.csv"),
  smiles_csv    = Sys.getenv("SMILES_CSV", "output/tables/btc_daily_iv_smiles.csv"),
  binance_cache = Sys.getenv("BINANCE_CACHE", "data/btc_spot_binance.csv"),

  # smile construction (only used if smiles_csv doesn't exist)
  tau_target_days = as.numeric(Sys.getenv("TAU_TARGET_DAYS", 30)),
  tau_window_days = as.numeric(Sys.getenv("TAU_WINDOW_DAYS", 7)),
  u_min           = as.numeric(Sys.getenv("U_MIN", -0.35)),
  u_max           = as.numeric(Sys.getenv("U_MAX", 0.35)),
  u_grid_n        = as.integer(Sys.getenv("U_GRID_N", 61)),
  min_trades_day  = as.integer(Sys.getenv("MIN_TRADES_DAY", 50)),
  min_bins_day    = as.integer(Sys.getenv("MIN_BINS_DAY", 10)),

  # relevant-change detector
  det_window    = as.integer(Sys.getenv("DET_WINDOW", 20)),
  det_eps_trim  = as.numeric(Sys.getenv("DET_EPS_TRIM", 0.10)),
  det_Delta_rms = as.numeric(Sys.getenv("DET_DELTA_RMS", 0.05)), # 5 vol points
  det_alpha     = as.numeric(Sys.getenv("DET_ALPHA", 0.05)),
  det_require_both = as.logical(as.integer(Sys.getenv("DET_REQUIRE_BOTH", 1))),
  q_cache       = Sys.getenv("Q_CACHE", "data/qvalues_cpp.rds"),

  # VRP strategy layer
  ewma_lambda   = as.numeric(Sys.getenv("EWMA_LAMBDA", 0.94)),
  vrp_window    = as.integer(Sys.getenv("VRP_WINDOW", 60)),
  vrp_q_high    = as.numeric(Sys.getenv("VRP_Q_HIGH", 0.75)),
  vrp_q_low     = as.numeric(Sys.getenv("VRP_Q_LOW", 0.25)),
  skew_filter_q = as.numeric(Sys.getenv("SKEW_FILTER_Q", 0.80)),
  cooldown_days = as.integer(Sys.getenv("COOLDOWN_DAYS", 3)),

  # If 1, allow long straddles when VRP is very low. Default 0 (short-vol only).
  allow_long = as.logical(as.integer(Sys.getenv("ALLOW_LONG", 0))),

  # execution/backtest layer
  strike_step = as.numeric(Sys.getenv("STRIKE_STEP", 1000)),
  tc_opt      = as.numeric(Sys.getenv("TC_OPT", 0.0010)),
  tc_spot     = as.numeric(Sys.getenv("TC_SPOT", 0.0002)),
  day_count   = as.numeric(Sys.getenv("DAY_COUNT", 365))
)


# -------------------------
# Load or construct smiles
# -------------------------

if (file.exists(config$smiles_csv)) {
  message("Loading smiles from: ", config$smiles_csv)
  sm <- read_smiles_csv(config$smiles_csv)
  dates <- sm$dates
  u_grid <- sm$u_grid
  curves_mat <- sm$curves_mat
} else {
  message("Smiles CSV not found; constructing from options trades...")
  trades <- load_deribit_trade_csv(config$options_csv)
  date_range <- range(trades$date)
  spot <- tryCatch(
    get_binance_daily_close(
      start_date = date_range[1] - 2,
      end_date   = date_range[2] + 2,
      cache_path = config$binance_cache,
      symbol     = "BTCUSDT",
      interval   = "1d"
    ),
    error = function(e) {
      message("Binance fetch failed: ", conditionMessage(e))
      message("Falling back to daily median index_price from the options trades.")
      trades[, .(date, close = median(index_price, na.rm = TRUE))]
    }
  )
  u_grid <- seq(config$u_min, config$u_max, length.out = config$u_grid_n)
  smiles <- build_daily_iv_smiles(
    trades          = trades,
    spot            = spot,
    u_grid          = u_grid,
    tau_target_days = config$tau_target_days,
    tau_window_days = config$tau_window_days,
    min_trades_day  = config$min_trades_day,
    min_bins_day    = config$min_bins_day
  )
  curves_mat <- smiles$curves_mat
  dates <- smiles$dates
  u_grid <- smiles$u_grid
}

message(sprintf("Working with %d daily curves on a %d-point grid.", nrow(curves_mat), ncol(curves_mat)))


# -----------------
# Spot + features
# -----------------

date_range <- range(as.Date(dates))
spot <- tryCatch(
  get_binance_daily_close(
    start_date = date_range[1] - 60,
    end_date   = date_range[2] + 2,
    cache_path = config$binance_cache,
    symbol     = "BTCUSDT",
    interval   = "1d"
  ),
  error = function(e) {
    message("Binance fetch failed: ", conditionMessage(e))
    message("Falling back to NA spot; strategy needs spot to run.")
    data.table::data.table(date = sort(unique(as.Date(dates))), close = NA_real_)
  }
)

feat <- smile_features(curves_mat = curves_mat, dates = dates, u_grid = u_grid)
fcst <- ewma_vol_forecast(spot_dt = spot, lambda = config$ewma_lambda, day_count = config$day_count)


# ----------------------------
# Online relevant-change gate
# ----------------------------

det <- rolling_relevant_detector(
  curves_mat = curves_mat,
  dates      = dates,
  window     = min(config$det_window, nrow(curves_mat)),
  eps_trim   = config$det_eps_trim,
  Delta_rms  = config$det_Delta_rms,
  alpha      = config$det_alpha,
  require_both = config$det_require_both,
  q_cache    = config$q_cache
)


# ----------------------------
# Positioning + P&L backtest
# ----------------------------

pos_dt <- build_positions_relevant_vrp(
  features      = feat,
  spot_fcst     = fcst,
  detector      = det,
  tau_days      = config$tau_target_days,
  vrp_window    = config$vrp_window,
  q_high        = config$vrp_q_high,
  q_low         = config$vrp_q_low,
  skew_filter_q = config$skew_filter_q,
  cooldown_days = config$cooldown_days,
  allow_long    = config$allow_long
)


# Backtest helper: backtest uses spot_dt for close; if pos_dt already has a close column,
# the backtest function can handle it, but we keep this defensive to avoid merge suffix surprises.
pos_dt_bt <- data.table::copy(pos_dt)
if ("close" %in% names(pos_dt_bt)) pos_dt_bt[, close := NULL]

pnl <- backtest_delta_hedged_straddle(
  pos_dt     = pos_dt_bt,
  curves_mat = curves_mat,
  dates      = dates,
  u_grid     = u_grid,
  spot_dt    = spot,
  tau_days   = config$tau_target_days,
  strike_step = config$strike_step,
  tc_opt     = config$tc_opt,
  tc_spot    = config$tc_spot,
  day_count  = config$day_count
)

data.table::fwrite(pos_dt, file = "output/strategy/relevant_vrp_positions.csv")
data.table::fwrite(pnl,    file = "output/strategy/relevant_vrp_pnl.csv")
plot_cum_pnl(pnl, out_file = "output/strategy/relevant_vrp_cum_pnl.png")


# Simple summary
summary <- pnl[, .(
  n_days = .N,
  n_trade_days = sum(pos != 0, na.rm = TRUE),
  pnl_total = sum(pnl, na.rm = TRUE),
  pnl_mean = mean(pnl, na.rm = TRUE),
  pnl_sd   = stats::sd(pnl, na.rm = TRUE),
  sharpe_daily = data.table::fifelse(stats::sd(pnl, na.rm = TRUE) > 0, mean(pnl, na.rm = TRUE) / stats::sd(pnl, na.rm = TRUE), NA_real_),
  sharpe_ann_365 = data.table::fifelse(stats::sd(pnl, na.rm = TRUE) > 0, (mean(pnl, na.rm = TRUE) / stats::sd(pnl, na.rm = TRUE)) * sqrt(config$day_count), NA_real_)
)]

writeLines(capture.output(print(summary)), con = "output/strategy/relevant_vrp_summary.txt")

message("\nDone.")
message("- Positions: output/strategy/relevant_vrp_positions.csv")
message("- P&L log  : output/strategy/relevant_vrp_pnl.csv")
message("- Plot     : output/strategy/relevant_vrp_cum_pnl.png")
message("- Summary  : output/strategy/relevant_vrp_summary.txt")
print(summary)