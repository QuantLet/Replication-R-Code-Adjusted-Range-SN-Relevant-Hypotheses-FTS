# Profitable strategy driven by the relevant break test:
# trade delta-hedged straddle ONLY after a relevant rejection,
# long if direction=="up", short if direction=="down".

source("R/packages.R")
source("R/binance_spot.R")
source("R/relevant_cpp_tests.R")
source("R/strategy_relevant_trading.R")  # contains read_smiles_csv, smile_features, backtest, etc.

dir.create("output/strategy_break", recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load smiles (6-month file)
# -------------------------
smiles_csv <- "output/tables/btc_daily_iv_smiles.csv"
sm <- read_smiles_csv(smiles_csv)
dates <- sm$dates
u_grid <- sm$u_grid
curves_mat <- sm$curves_mat

# -------------------------
# Load spot (use cached CSV)
# -------------------------
spot <- data.table::fread("data/btc_spot_binance.csv")
spot[, date := as.Date(date)]
spot <- spot[, .(date, close)]
spot <- spot[date >= (min(dates) - 60) & date <= (max(dates) + 2)]
spot <- spot[order(date)]

# -------------------------
# Features (for bookkeeping)
# -------------------------
feat <- smile_features(curves_mat = curves_mat, dates = dates, u_grid = u_grid)

# -------------------------
# Rolling relevant detector
# -------------------------
# Key settings you asked for:
# - rolling window W = 30
# - alpha = 0.10 (90% quantile)
# - Delta_rms = 1 vol point = 0.01 in decimal IV units
# - require_both = FALSE (more rejections; uses both shells)
det <- rolling_relevant_detector(
  curves_mat = curves_mat,
  dates      = dates,
  window     = 30,
  eps_trim   = 0.10,
  Delta_rms  = 0.01,
  alpha      = 0.10,
  require_both = FALSE,
  q_cache    = "data/qvalues_cpp.rds"
)

# -------------------------
# Positions: break-direction straddle
# -------------------------
pos_dt <- build_positions_break_direction_straddle(
  features  = feat,
  detector  = det,
  hold_days = 3
)

# -------------------------
# Backtest: delta-hedged straddle
# -------------------------
# NOTE: pos_dt does not contain close, backtest merges from spot
pnl <- backtest_delta_hedged_straddle(
  pos_dt     = pos_dt[, .(date, pos)],
  curves_mat = curves_mat,
  dates      = dates,
  u_grid     = u_grid,
  spot_dt    = spot,
  tau_days   = 30,
  strike_step = 1000,
  tc_opt     = 0.0010,
  tc_spot    = 0.0002,
  day_count  = 365
)

# -------------------------
# Outputs
# -------------------------
data.table::fwrite(det,    "output/strategy_break/relevant_break_detector.csv")
data.table::fwrite(pos_dt, "output/strategy_break/relevant_break_positions.csv")
data.table::fwrite(pnl,    "output/strategy_break/relevant_break_pnl.csv")

# plot (reuse your helper)
p <- ggplot2::ggplot(pnl, ggplot2::aes(x = date, y = cum_pnl)) +
  ggplot2::geom_line() +
  ggplot2::labs(
    title = "Cumulative P&L (break-direction straddle, relevant-test driven)",
    x = "Date",
    y = "Cumulative P&L (USD per 1 straddle unit)"
  ) +
  ggplot2::theme_minimal()
ggplot2::ggsave("output/strategy_break/relevant_break_cum_pnl.png", p, width = 9, height = 4.5)

# summary
summary <- pnl[, .(
  n_days = .N,
  n_trade_days = sum(pos != 0, na.rm = TRUE),
  pnl_total = sum(pnl, na.rm = TRUE),
  pnl_mean = mean(pnl, na.rm = TRUE),
  pnl_sd   = stats::sd(pnl, na.rm = TRUE),
  sharpe_daily = data.table::fifelse(stats::sd(pnl, na.rm = TRUE) > 0,
                                     mean(pnl, na.rm = TRUE) / stats::sd(pnl, na.rm = TRUE),
                                     NA_real_),
  sharpe_ann_365 = data.table::fifelse(stats::sd(pnl, na.rm = TRUE) > 0,
                                       (mean(pnl, na.rm = TRUE) / stats::sd(pnl, na.rm = TRUE)) * sqrt(365),
                                       NA_real_)
)]
writeLines(capture.output(print(summary)), con = "output/strategy_break/relevant_break_summary.txt")
print(summary)

message("Done. Outputs in output/strategy_break/")
