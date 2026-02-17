# ----------------------------------------------------------------------------
# Executable (research) trading strategy driven by *relevant change* detection
# ----------------------------------------------------------------------------
#
# Goal
# ----
# Turn the relevant-change toolbox used in the empirical section into an
# *implementable* daily strategy (with a clean, no-lookahead backtest).
#
# Inputs
# ------
# - Functional time series of IV smiles X_t(u) on a fixed u-grid.
# - Daily BTC spot close S_t.
#
# Strategy idea (high-level)
# --------------------------
# 1) Run an *online / rolling-window* relevant change-point test on X_t(u).
#    We only use information available up to day t.
# 2) Use detected relevant changes to define regimes and to impose a cooldown
#    (risk-off) period immediately after a structural break.
# 3) Within stable regimes, trade a delta-hedged, near-ATM constant-maturity
#    straddle based on a volatility-risk-premium (VRP) filter.
#
# IMPORTANT
# ---------
# This is research/backtest code for your paper. It is NOT production trading
# infrastructure and ignores many real-world frictions (margin, funding,
# microstructure, intraday hedging, etc.).


# ---- Helpers: parse smiles CSV into (dates, u_grid, curves_mat) ----

parse_u_from_colnames <- function(cols) {
  # cols like u_m0.350, u_p0.012
  stopifnot(is.character(cols))
  u <- rep(NA_real_, length(cols))
  is_m <- grepl("^u_m", cols)
  is_p <- grepl("^u_p", cols)
  u[is_m] <- -as.numeric(sub("^u_m", "", cols[is_m]))
  u[is_p] <-  as.numeric(sub("^u_p", "", cols[is_p]))
  u
}

read_smiles_csv <- function(path) {
  dt <- data.table::fread(path)
  if (!("date" %in% names(dt))) stop("smiles csv must contain a 'date' column")
  dt[, date := as.Date(date)]
  cols <- setdiff(names(dt), "date")
  u <- parse_u_from_colnames(cols)
  if (any(!is.finite(u))) stop("Could not parse u-grid from column names")
  # ensure columns are ordered by u ascending
  o <- order(u)
  u <- u[o]
  cols <- cols[o]
  curves <- as.matrix(dt[, ..cols])
  list(dates = dt$date, u_grid = u, curves_mat = curves)
}


# ---- Smile feature extraction (keeps functional nature; scalars are summaries) ----

interp_smile <- function(u_grid, iv_vec, u) {
  # linear interpolation (rule=2 extends flat outside range)
  stats::approx(x = u_grid, y = iv_vec, xout = u, rule = 2)$y
}

smile_features <- function(curves_mat, dates, u_grid) {
  stopifnot(is.matrix(curves_mat))
  N <- nrow(curves_mat)
  P <- ncol(curves_mat)
  if (length(dates) != N) stop("dates length must match nrow(curves_mat)")
  if (length(u_grid) != P) stop("u_grid length must match ncol(curves_mat)")

  # trapezoid weights over u-grid for a simple L2 proxy of level
  # (not used in tests; tests use their own trap_weights on [0,1])
  du <- diff(u_grid)
  w <- c(du[1], du[-1] + du[-length(du)], du[length(du)]) / 2

  out <- data.table::data.table(date = as.Date(dates))
  # ATM and skew at fixed u
  out[, iv_atm := vapply(seq_len(N), function(i) interp_smile(u_grid, curves_mat[i, ], 0), numeric(1))]
  out[, skew_020 := vapply(seq_len(N), function(i) {
    interp_smile(u_grid, curves_mat[i, ], -0.2) - interp_smile(u_grid, curves_mat[i, ], 0.2)
  }, numeric(1))]
  out[, level_l2 := vapply(seq_len(N), function(i) sum(w * curves_mat[i, ]), numeric(1))]
  out
}


# ---- Spot-based volatility forecast (EWMA on daily log returns) ----

ewma_vol_forecast <- function(spot_dt, lambda = 0.94, day_count = 365) {
  stopifnot(is.data.table(spot_dt))
  if (!all(c("date", "close") %in% names(spot_dt))) stop("spot_dt must have date, close")
  dt <- data.table::copy(spot_dt)
  dt[, date := as.Date(date)]
  dt <- dt[order(date)]
  dt[, r := c(NA_real_, diff(log(close)))]

  # initialize with sample variance of first few returns
  v <- rep(NA_real_, nrow(dt))
  r2 <- dt$r^2
  init_idx <- which(is.finite(r2))[1:min(10, sum(is.finite(r2)))]
  v0 <- if (length(init_idx) >= 2) mean(r2[init_idx], na.rm = TRUE) else NA_real_

  for (i in seq_len(nrow(dt))) {
    if (!is.finite(r2[i])) {
      v[i] <- if (i == 1) v0 else v[i - 1]
    } else if (i == 1) {
      v[i] <- v0
    } else {
      v[i] <- lambda * v[i - 1] + (1 - lambda) * r2[i]
    }
  }
  dt[, var_ewma := v]
  dt[, vol_ewma_ann := sqrt(pmax(var_ewma, 0) * day_count)]
  dt[, .(date, close, vol_ewma_ann)]
}


# ---- Rolling-window relevant change detector (online, no look-ahead) ----

rolling_relevant_detector <- function(curves_mat,
                                      dates,
                                      window = 60,
                                      eps_trim = 0.10,
                                      Delta_rms = 0.05,
                                      alpha = 0.05,
                                      require_both = TRUE,
                                      q_cache = "data/qvalues_cpp.rds") {
  # Delta_rms in *decimal* IV units (e.g., 0.05 = 5 vol points)
  stopifnot(is.matrix(curves_mat))
  N <- nrow(curves_mat)
  if (length(dates) != N) stop("dates length must equal nrow(curves_mat)")
  if (window < 10) stop("window too small")
  if (window > N) stop("window > sample size")

  # Pivotal quantiles (cached)
  qs <- get_cpp_q_values(alpha = alpha, cache_path = q_cache)
  q_quad <- as.numeric(qs$quadratic[as.character(1 - alpha)])
  q_ar   <- as.numeric(qs$adjusted_range[as.character(1 - alpha)])
  Delta <- Delta_rms^2

  out <- data.table::data.table(
    date_end = as.Date(dates),
    reject_quad = FALSE,
    reject_ar = FALSE,
    reject = FALSE,
    cp_date = as.Date(NA),
    direction = NA_character_,
    D_hat = NA_real_,
    norm_quad = NA_real_,
    norm_ar = NA_real_
  )

  for (t in window:N) {
    idx <- (t - window + 1):t
    cm  <- curves_mat[idx, , drop = FALSE]
    dd  <- as.Date(dates[idx])

    rq <- cpp_relevant_test(cm, dd, eps_trim = eps_trim, normalizer = "quadratic")
    rr <- cpp_relevant_test(cm, dd, eps_trim = eps_trim, normalizer = "adjusted_range")

    rej_q <- is.finite(rq$statistic) && is.finite(rq$normalizer) && (rq$statistic > Delta + q_quad * rq$normalizer)
    rej_r <- is.finite(rr$statistic) && is.finite(rr$normalizer) && (rr$statistic > Delta + q_ar   * rr$normalizer)
    rej <- if (require_both) (rej_q && rej_r) else (rej_q || rej_r)

    out[t, `:=`(
      reject_quad = rej_q,
      reject_ar   = rej_r,
      reject      = rej,
      cp_date     = as.Date(rq$cp_date),
      D_hat       = rq$statistic,
      norm_quad   = rq$normalizer,
      norm_ar     = rr$normalizer
    )]

    # direction: compare mean smile in "post" vs "pre" according to estimated k_hat
    if (rej) {
      k <- rq$k_hat
      pre_mean  <- colMeans(cm[1:k, , drop = FALSE], na.rm = TRUE)
      post_mean <- colMeans(cm[(k + 1):nrow(cm), , drop = FALSE], na.rm = TRUE)
      # use ATM difference sign as a simple direction label
      # (still derived from the functional object)
      u0 <- which.min(abs(seq_len(ncol(cm)) - ceiling(ncol(cm) / 2)))
      # safer: compute direction by average level shift
      lvl_shift <- mean(post_mean - pre_mean, na.rm = TRUE)
      out[t, direction := ifelse(lvl_shift >= 0, "up", "down")]
    }
  }

  out
}


# ---- Strategy: regime-gated VRP straddle (daily, delta-hedged) ----

bs_call_put <- function(S, K, tau, sigma, r = 0) {
  # Black-Scholes in USD (underlying is 1 BTC, price in USD)
  # r assumed continuously compounded, q=0.
  if (!all(is.finite(c(S, K, tau, sigma))) || S <= 0 || K <= 0 || tau <= 0 || sigma <= 0) {
    return(list(call = NA_real_, put = NA_real_, delta_call = NA_real_, delta_put = NA_real_))
  }
  sig_sqrt <- sigma * sqrt(tau)
  d1 <- (log(S / K) + (r + 0.5 * sigma^2) * tau) / sig_sqrt
  d2 <- d1 - sig_sqrt
  Nd1 <- stats::pnorm(d1)
  Nd2 <- stats::pnorm(d2)
  disc <- exp(-r * tau)
  call <- S * Nd1 - K * disc * Nd2
  put  <- K * disc * stats::pnorm(-d2) - S * stats::pnorm(-d1)
  list(call = call, put = put, delta_call = Nd1, delta_put = Nd1 - 1)
}


build_positions_relevant_vrp <- function(features,
                                        spot_fcst,
                                        detector,
                                        tau_days = 30,
                                        vrp_window = 60,
                                        q_high = 0.75,
                                        q_low  = 0.25,
                                        skew_filter_q = 0.80,
                                        cooldown_days = 3,
                                        allow_long = FALSE) {
  stopifnot(is.data.table(features), is.data.table(spot_fcst), is.data.table(detector))
  dt <- merge(features, spot_fcst, by = "date", all.x = TRUE)
  det <- detector[, .(date = date_end, reject, direction, cp_date)]
  dt <- merge(dt, det, by = "date", all.x = TRUE)
  dt <- dt[order(date)]

  # VRP proxy: implied ATM variance minus EWMA forecast variance
  dt[, vrp := iv_atm^2 - vol_ewma_ann^2]

  # rolling thresholds (no look-ahead)
  dt[, vrp_hi := NA_real_]
  dt[, vrp_lo := NA_real_]
  dt[, skew_cap := NA_real_]
  for (i in seq_len(nrow(dt))) {
    lo <- max(1, i - vrp_window)
    hist <- dt[lo:(i - 1)]
    if (nrow(hist) >= 10) {
      dt[i, vrp_hi := stats::quantile(hist$vrp, probs = q_high, na.rm = TRUE, names = FALSE)]
      dt[i, vrp_lo := stats::quantile(hist$vrp, probs = q_low,  na.rm = TRUE, names = FALSE)]
      dt[i, skew_cap := stats::quantile(hist$skew_020, probs = skew_filter_q, na.rm = TRUE, names = FALSE)]
    }
  }

  # cooldown flag: if a relevant change was detected within the last cooldown_days, stay flat
  dt[, cool := FALSE]
  if (cooldown_days > 0) {
    idx_rej <- which(dt$reject %in% TRUE)
    for (j in idx_rej) {
      lo <- j
      hi <- min(nrow(dt), j + cooldown_days)
      dt[lo:hi, cool := TRUE]
    }
  }

  # Position rule (daily rolled straddle):
  # - Only trade when not in cooldown AND thresholds exist.
  # - Default: SHORT-VOL ONLY (sell straddles when VRP is high and skew not extreme).
  #   Empirically, the long-straddle leg can be dominated by theta when hedging is only daily.
  # - Set allow_long = TRUE if you also want to buy straddles when VRP is very low.
  dt[, pos := 0L]

  # IMPORTANT: we are outside data.table's j here, so column names are NOT in scope.
  # Use dt$col (or compute inside dt[, ...]) to avoid errors like: object 'cool' not found.
  eligible <- !dt$cool & is.finite(dt$vrp_hi) & is.finite(dt$vrp_lo) & is.finite(dt$skew_cap) & is.finite(dt$vrp)

  # Short straddle when VRP is high and skew not extreme
  dt[eligible & (vrp > vrp_hi) & (skew_020 <= skew_cap), pos := -1L]

  # Optional: long straddle when VRP is very low
  if (isTRUE(allow_long)) {
    dt[eligible & (vrp < vrp_lo), pos := 1L]
  }

  dt
}


backtest_delta_hedged_straddle <- function(pos_dt,
                                           curves_mat,
                                           dates,
                                           u_grid,
                                           spot_dt,
                                           tau_days = 30,
                                           strike_step = 1000,
                                           tc_opt = 0.0010,
                                           tc_spot = 0.0002,
                                           day_count = 365) {
  # pos_dt must contain columns: date, pos, iv_atm, etc.
  stopifnot(is.data.table(pos_dt))
  stopifnot(is.matrix(curves_mat))
  stopifnot(length(dates) == nrow(curves_mat))
  stopifnot(length(u_grid) == ncol(curves_mat))

  spot <- data.table::copy(spot_dt)
  spot[, date := as.Date(date)]
  data.table::setkey(spot, date)
  # Merge spot close into pos_dt; avoid name clashes if pos_dt already contains 'close'
  dt <- merge(pos_dt, spot[, .(date, close)], by = "date", all.x = TRUE, suffixes = c("", "_spot"))
  if ("close_spot" %in% names(dt)) {
    # if pos_dt close is missing for a day, fill from spot close
    dt[!is.finite(close) & is.finite(close_spot), close := close_spot]
    dt[, close_spot := NULL]
  }
  dt <- dt[order(date)]

  # Map dates to row index in curves_mat
  idx_map <- data.table::data.table(date = as.Date(dates), row = seq_along(dates))
  dt <- merge(dt, idx_map, by = "date", all.x = TRUE)
  dt <- dt[is.finite(row) & is.finite(close) & close > 0]

  # One-day holding period backtest (daily roll)
  res <- dt[1:(.N - 1), .(date = date,
                          date_next = dt$date[.I + 1],
                          pos = pos,
                          S = close,
                          S_next = dt$close[.I + 1],
                          row = row,
                          row_next = dt$row[.I + 1])]

  res[, `:=`(
    K = NA_real_,
    u_entry = NA_real_,
    u_exit = NA_real_,
    sigma_entry = NA_real_,
    sigma_exit = NA_real_,
    price_entry = NA_real_,
    price_exit = NA_real_,
    delta_straddle = NA_real_,
    hedge_units = NA_real_,
    pnl_opt = 0,
    pnl_hedge = 0,
    costs = 0,
    pnl = 0
  )]

  tau0 <- tau_days / day_count
  tau1 <- (tau_days - 1) / day_count
  if (tau1 <= 0) stop("tau_days must be >= 2 for 1-day hold")

  for (i in seq_len(nrow(res))) {
    if (!is.finite(res$pos[i]) || res$pos[i] == 0L) next

    S  <- res$S[i]
    S1 <- res$S_next[i]
    if (!is.finite(S) || !is.finite(S1) || S <= 0 || S1 <= 0) next

    K <- round(S / strike_step) * strike_step
    if (!is.finite(K) || K <= 0) next

    u0 <- log(K / S)
    u1 <- log(K / S1)

    iv0 <- interp_smile(u_grid, curves_mat[res$row[i], ], u0)
    iv1 <- interp_smile(u_grid, curves_mat[res$row_next[i], ], u1)
    if (!is.finite(iv0) || !is.finite(iv1) || iv0 <= 0 || iv1 <= 0) next

    pr0 <- bs_call_put(S = S,  K = K, tau = tau0, sigma = iv0)
    pr1 <- bs_call_put(S = S1, K = K, tau = tau1, sigma = iv1)
    if (!is.finite(pr0$call) || !is.finite(pr0$put) || !is.finite(pr1$call) || !is.finite(pr1$put)) next

    str0 <- pr0$call + pr0$put
    str1 <- pr1$call + pr1$put

    # straddle delta at entry
    d_str <- pr0$delta_call + pr0$delta_put  # 2N(d1)-1
    hedge <- -as.numeric(res$pos[i]) * d_str

    pnl_opt <- as.numeric(res$pos[i]) * (str1 - str0)
    pnl_hedge <- hedge * (S1 - S)

    # simple proportional transaction costs
    cost_opt <- abs(res$pos[i]) * (str0 + str1) * tc_opt
    cost_spot <- (abs(hedge) * S + abs(hedge) * S1) * tc_spot
    cost <- cost_opt + cost_spot

    # IMPORTANT: data.table scoping
    # Avoid RHS names that collide with column names (e.g., K, pnl_opt).
    K_val         <- K
    pnl_opt_val   <- pnl_opt
    pnl_hedge_val <- pnl_hedge
    cost_val      <- cost

    res[i, `:=`(
      K = K_val,
      u_entry = u0,
      u_exit = u1,
      sigma_entry = iv0,
      sigma_exit = iv1,
      price_entry = str0,
      price_exit = str1,
      delta_straddle = d_str,
      hedge_units = hedge,
      pnl_opt = pnl_opt_val,
      pnl_hedge = pnl_hedge_val,
      costs = cost_val,
      pnl = pnl_opt_val + pnl_hedge_val - cost_val
    )]
  }

  res[, cum_pnl := cumsum(data.table::fifelse(is.finite(pnl), pnl, 0))]
  res
}


plot_cum_pnl <- function(pnl_dt, out_file) {
  stopifnot(is.data.table(pnl_dt))
  p <- ggplot2::ggplot(pnl_dt, ggplot2::aes(x = date, y = cum_pnl)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Cumulative P&L (delta-hedged straddle, regime-gated)",
                  x = "Date", y = "Cumulative P&L (USD per 1 straddle unit)") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(out_file, p, width = 9, height = 4.5)
}

# ---- Strategy: trade ONLY on relevant breaks, directionally in vol ----
# If direction == "up"  -> long straddle for hold_days
# If direction == "down"-> short straddle for hold_days
# Else flat.

build_positions_break_direction_straddle <- function(features,
                                                     detector,
                                                     hold_days = 3L) {
  stopifnot(data.table::is.data.table(features))
  stopifnot(data.table::is.data.table(detector))
  if (hold_days < 1) stop("hold_days must be >= 1")
  
  dt <- data.table::copy(features)
  det <- detector[, .(date = as.Date(date_end),
                      reject,
                      direction,
                      cp_date)]
  dt <- merge(dt, det, by = "date", all.x = TRUE)
  dt <- dt[order(date)]
  dt[is.na(reject), reject := FALSE]
  
  dt[, pos := 0L]
  
  idx_rej <- which(dt$reject %in% TRUE)
  if (length(idx_rej) == 0) return(dt)
  
  for (j in idx_rej) {
    dirj <- dt$direction[j]
    if (!is.character(dirj) || is.na(dirj) || !nzchar(dirj)) next
    
    sign <- if (dirj == "up") 1L else -1L
    hi <- min(nrow(dt), j + as.integer(hold_days) - 1L)
    
    # later signals overwrite earlier ones (intentional)
    dt[j:hi, pos := sign]
  }
  
  dt
}
