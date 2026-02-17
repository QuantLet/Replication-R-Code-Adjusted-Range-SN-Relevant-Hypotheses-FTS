# Construct daily constant-maturity IV smiles from trade-level options data.
#
# Output is a daily functional time series represented on a fixed moneyness grid.
# No FPCA/KL truncation is used; we keep the full grid as the functional object.

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

build_daily_iv_smiles <- function(trades,
                                 spot,
                                 u_grid,
                                 tau_target_days = 30,
                                 tau_window_days = 7,
                                 min_trades_day = 50,
                                 min_bins_day = 10,
                                 otm_only = TRUE,
                                 spot_source = c("binance", "index_price")) {
  spot_source <- match.arg(spot_source)
  stopifnot(is.data.table(trades), is.data.table(spot))

  # Merge daily spot close onto trades
  spot_day <- spot[, .(date, spot_close = close)]
  data.table::setkey(spot_day, date)
  data.table::setkey(trades, date)
  dt <- spot_day[trades]

  # Fallback spot from index_price if Binance missing (or if requested)
  if (spot_source == "index_price") {
    dt[, spot_close := data.table::fifelse(is.finite(index_price), index_price, spot_close)]
  } else {
    dt[, spot_close := data.table::fifelse(is.finite(spot_close), spot_close, index_price)]
  }

  # time-to-maturity in days (calendar)
  dt[, tau_days := as.numeric(expiry_date - date)]
  dt <- dt[is.finite(tau_days) & tau_days > 0]
  dt <- dt[tau_days >= (tau_target_days - tau_window_days) & tau_days <= (tau_target_days + tau_window_days)]

  # log-moneyness u = log(K/S)
  dt <- dt[is.finite(spot_close) & spot_close > 0]
  dt[, u := log(strike / spot_close)]
  dt <- dt[is.finite(u) & u >= min(u_grid) & u <= max(u_grid)]

  # OTM filter: puts for u<0, calls for u>0
  if (otm_only) {
    dt <- dt[(cp_flag == "P" & u < 0) | (cp_flag == "C" & u > 0)]
  }

  # Bin setup: bins centered on u_grid; breaks are midpoints.
  mids <- u_grid
  breaks <- c(-Inf, (mids[-1] + mids[-length(mids)]) / 2, Inf)

  # For each day, compute weighted mean IV per bin and smooth.
  days <- sort(unique(dt$date))
  p <- length(u_grid)
  curves <- matrix(NA_real_, nrow = length(days), ncol = p)
  keep_day <- rep(FALSE, length(days))

  for (i in seq_along(days)) {
    d <- days[i]
    sub <- dt[date == d]
    if (nrow(sub) < min_trades_day) next

    sub[, bin := cut(u, breaks = breaks, labels = FALSE, include.lowest = TRUE)]
    bybin <- sub[, .(
      u_mid = mids[bin[1]],
      iv_bar = weighted_mean(iv_dec, amount),
      w_sum = sum(amount, na.rm = TRUE)
    ), by = bin]
    bybin <- bybin[is.finite(iv_bar) & is.finite(u_mid)]
    if (nrow(bybin) < min_bins_day) next

    # Smoothing spline (weighted)
    # Use only unique x (smooth.spline requires strictly increasing x)
    bybin <- bybin[order(u_mid)]
    bybin <- bybin[!duplicated(u_mid)]
    if (nrow(bybin) < min_bins_day) next

    fit <- tryCatch(
      stats::smooth.spline(x = bybin$u_mid, y = bybin$iv_bar, w = bybin$w_sum),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    pred <- stats::predict(fit, x = u_grid)$y
    if (any(!is.finite(pred))) {
      # simple fallback: linear interpolation across available bin means
      pred <- stats::approx(bybin$u_mid, bybin$iv_bar, xout = u_grid, rule = 2)$y
    }
    curves[i, ] <- pred
    keep_day[i] <- TRUE
  }

  days_kept <- days[keep_day]
  curves_kept <- curves[keep_day, , drop = FALSE]

  if (length(days_kept) < 10) {
    warning("Very few daily curves constructed (", length(days_kept), "). Consider relaxing filters.")
  }

  # Return as N x P matrix (N days)
  list(
    dates = days_kept,
    u_grid = u_grid,
    curves_mat = curves_kept
  )
}
