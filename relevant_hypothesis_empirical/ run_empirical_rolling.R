# run_empirical_rolling.R
# ------------------------------------------------------------
# Rolling-window relevant change-point tests for BTC IV smiles
# - Rolling window (default 30 observations)
# - Uses 90% quantile (alpha = 0.10) for more rejections
# - Writes per-window CP + relevance decisions + implied boundary
# ------------------------------------------------------------
#
# Run from project root:
#   source("run_empirical_rolling.R")
#
# ENV VARS (optional):
#   OPTIONS_CSV          path to Deribit trades CSV/ZIP/dir (default: data/options_trades.csv)
#   SMILES_CSV           path to cached daily smiles CSV (default: output/tables/btc_daily_iv_smiles.csv)
#   FORCE_REBUILD_SMILES 1 to rebuild smiles even if SMILES_CSV exists (default: 0)
#   SAMPLE_DAYS          length of sample in days counting back from SAMPLE_END (default: 183 ~ 6 months)
#   SAMPLE_START         YYYY-MM-DD (optional override)
#   SAMPLE_END           YYYY-MM-DD (optional override)
#
#   ROLL_WINDOW          rolling window length in #curves (default: 30)
#   ROLL_STEP            step size (default: 1)
#   ROLL_ALPHA           significance level (default: 0.10 -> 90% quantile)
#   ROLL_DELTA_RMS        relevance threshold in decimal IV units (default: 0.01 = 1 vol point RMS)
#   ROLL_REQUIRE_BOTH    1 require both SN normalizers to reject; 0 = either rejects (default: 0)
#
# OUTPUT:
#   output/tables/btc_iv_rolling_cpp_*.csv
#   output/figures/btc_iv_rolling_boundary_ar_*.png
setwd('/Users/sunjiajing/Desktop/R2026/btc_iv_empirical_relevant_tests 8')
source("R/packages.R")
source("R/binance_spot.R")
source("R/parse_deribit_trades.R")
source("R/build_iv_smile.R")
source("R/relevant_cpp_tests.R")

dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ---- Helpers to read cached smiles CSV (same format your project writes) ----
parse_u_from_colnames <- function(cols) {
  u <- rep(NA_real_, length(cols))
  is_m <- grepl("^u_m", cols)
  is_p <- grepl("^u_p", cols)
  u[is_m] <- -as.numeric(sub("^u_m", "", cols[is_m]))
  u[is_p] <-  as.numeric(sub("^u_p", "", cols[is_p]))
  u
}

read_smiles_csv_local <- function(path) {
  dt <- data.table::fread(path)
  if (!("date" %in% names(dt))) stop("smiles csv must contain a 'date' column")
  dt[, date := as.Date(date)]
  cols <- setdiff(names(dt), "date")
  u <- parse_u_from_colnames(cols)
  if (any(!is.finite(u))) stop("Could not parse u-grid from column names in smiles CSV")
  o <- order(u)
  u <- u[o]
  cols <- cols[o]
  curves <- as.matrix(dt[, ..cols])
  list(dates = dt$date, u_grid = u, curves_mat = curves)
}

# ---- Config ----
config <- list(
  options_csv   = Sys.getenv("OPTIONS_CSV", "data/options_trades.csv"),
  smiles_csv    = Sys.getenv("SMILES_CSV",  "output/tables/btc_daily_iv_smiles.csv"),
  binance_cache = Sys.getenv("BINANCE_CACHE", "data/btc_spot_binance.csv"),
  
  # 6-month sample control
  sample_days  = as.integer(Sys.getenv("SAMPLE_DAYS", 183)),
  sample_start = Sys.getenv("SAMPLE_START", ""),
  sample_end   = Sys.getenv("SAMPLE_END", ""),
  
  # smile construction
  tau_target_days = as.numeric(Sys.getenv("TAU_TARGET_DAYS", 30)),
  tau_window_days = as.numeric(Sys.getenv("TAU_WINDOW_DAYS", 7)),
  eps_trim        = as.numeric(Sys.getenv("EPS_TRIM", 0.05)),
  u_min           = as.numeric(Sys.getenv("U_MIN", -0.35)),
  u_max           = as.numeric(Sys.getenv("U_MAX",  0.35)),
  u_grid_n        = as.integer(Sys.getenv("U_GRID_N", 61)),
  min_trades_day  = as.integer(Sys.getenv("MIN_TRADES_DAY", 50)),
  min_bins_day    = as.integer(Sys.getenv("MIN_BINS_DAY", 10)),
  
  # rolling test settings
  roll_window       = as.integer(Sys.getenv("ROLL_WINDOW", 30)),
  roll_step         = as.integer(Sys.getenv("ROLL_STEP", 1)),
  roll_alpha        = as.numeric(Sys.getenv("ROLL_ALPHA", 0.10)),     # 90% quantile
  roll_Delta_rms    = as.numeric(Sys.getenv("ROLL_DELTA_RMS", 0.01)), # 1 vol point RMS
  roll_require_both = as.logical(as.integer(Sys.getenv("ROLL_REQUIRE_BOTH", 0))),
  
  # pivotal quantiles cache (Brownian functional)
  q_cache      = Sys.getenv("Q_CACHE", "data/qvalues_cpp.rds"),
  qsim_n       = as.integer(Sys.getenv("QSIM_N", 5000)),
  qsim_grid_N  = as.integer(Sys.getenv("QSIM_GRID_N", 1000)),
  
  force_rebuild_smiles = as.logical(as.integer(Sys.getenv("FORCE_REBUILD_SMILES", 0)))
)

# ---- Load or build smiles (6 months) ----
if (file.exists(config$smiles_csv) && !config$force_rebuild_smiles) {
  message("Loading cached smiles: ", config$smiles_csv)
  sm <- read_smiles_csv_local(config$smiles_csv)
  dates <- sm$dates
  u_grid <- sm$u_grid
  curves_mat <- sm$curves_mat
} else {
  message("Building smiles from trades (this can take time on large CSVs)...")
  trades <- load_deribit_trade_csv(config$options_csv)
  
  # enforce 6-month window
  sample_end <- if (nzchar(config$sample_end)) as.Date(config$sample_end) else max(trades$date)
  sample_start <- if (nzchar(config$sample_start)) as.Date(config$sample_start) else (sample_end - config$sample_days)
  trades <- trades[date >= sample_start & date <= sample_end]
  message(sprintf("Using trades sample window: %s to %s", sample_start, sample_end))
  
  date_range <- range(trades$date)
  spot <- tryCatch(
    get_binance_daily_close(
      start_date = date_range[1] - 5,
      end_date   = date_range[2] + 5,
      cache_path = config$binance_cache,
      symbol     = "BTCUSDT",
      interval   = "1d"
    ),
    error = function(e) {
      message("Binance fetch failed: ", conditionMessage(e))
      message("Falling back to daily median index_price from options trades.")
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
  
  # cache smiles (same format as your existing pipeline)
  smile_dt <- data.table::as.data.table(curves_mat)
  u_cols <- ifelse(u_grid >= 0,
                   paste0("u_p", sprintf("%.3f", abs(u_grid))),
                   paste0("u_m", sprintf("%.3f", abs(u_grid))))
  data.table::setnames(smile_dt, u_cols)
  smile_dt[, date := dates]
  data.table::setcolorder(smile_dt, c("date", u_cols))
  dir.create(dirname(config$smiles_csv), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(smile_dt, file = config$smiles_csv)
  message("Saved smiles to: ", config$smiles_csv)
}

message(sprintf("Working with %d daily curves on %d u-grid points.", nrow(curves_mat), ncol(curves_mat)))

# ---- Rolling window CP + relevance test at alpha=0.10 ----
# Ensure q-values exist for the chosen alpha
qs <- get_cpp_q_values(
  alpha      = config$roll_alpha,
  cache_path = config$q_cache,
  nsim       = config$qsim_n,
  N          = config$qsim_grid_N
)

prob <- as.character(1 - config$roll_alpha)
q_quad <- as.numeric(qs$quadratic[prob])
q_ar   <- as.numeric(qs$adjusted_range[prob])
Delta  <- config$roll_Delta_rms^2

rolling_cpp <- function(curves_mat, dates, window, step, eps_trim,
                        Delta, Delta_rms, q_quad, q_ar, require_both) {
  stopifnot(is.matrix(curves_mat))
  N <- nrow(curves_mat)
  dates <- as.Date(dates)
  if (length(dates) != N) stop("dates length must match nrow(curves_mat)")
  if (window < 10 || window > N) stop("Bad rolling window length")
  
  ends <- seq(window, N, by = step)
  out <- data.table::data.table(
    window_id = seq_along(ends),
    date_start = as.Date(NA),
    date_end   = as.Date(NA),
    n_obs      = window,
    k_hat      = NA_integer_,
    cp_date    = as.Date(NA),
    direction  = NA_character_,
    
    D_hat = NA_real_,
    RMS_shift_vol_pts = NA_real_,
    norm_quad = NA_real_,
    norm_ar   = NA_real_,
    
    Delta = Delta,
    Delta_rms_vol_pts = 100 * Delta_rms,
    q_quad = q_quad,
    q_ar   = q_ar,
    
    Delta_boundary_quad = NA_real_,
    Delta_boundary_ar   = NA_real_,
    RMS_boundary_vol_pts_quad = NA_real_,
    RMS_boundary_vol_pts_ar   = NA_real_,
    
    reject_quad = FALSE,
    reject_ar   = FALSE,
    reject      = FALSE,
    break_class = NA_character_
  )
  
  for (i in seq_along(ends)) {
    t <- ends[i]
    idx <- (t - window + 1):t
    cm <- curves_mat[idx, , drop = FALSE]
    dd <- dates[idx]
    
    out[i, `:=`(date_start = dd[1], date_end = dd[length(dd)])]
    
    rq <- tryCatch(cpp_relevant_test(cm, dd, eps_trim = eps_trim, normalizer = "quadratic"),
                   error = function(e) NULL)
    rr <- tryCatch(cpp_relevant_test(cm, dd, eps_trim = eps_trim, normalizer = "adjusted_range"),
                   error = function(e) NULL)
    if (is.null(rq) || is.null(rr)) next
    
    D <- rq$statistic
    Vq <- rq$normalizer
    Vr <- rr$normalizer
    
    rej_q <- is.finite(D) && is.finite(Vq) && (D > Delta + q_quad * Vq)
    rej_r <- is.finite(D) && is.finite(Vr) && (D > Delta + q_ar   * Vr)
    rej   <- if (isTRUE(require_both)) (rej_q && rej_r) else (rej_q || rej_r)
    
    # implied relevance boundary (largest Delta rejected):
    bq <- if (is.finite(D) && is.finite(Vq)) max(0, D - q_quad * Vq) else NA_real_
    br <- if (is.finite(D) && is.finite(Vr)) max(0, D - q_ar   * Vr) else NA_real_
    
    # direction: sign of average level shift (post - pre)
    dir <- NA_character_
    if (is.finite(rq$k_hat) && rq$k_hat >= 1 && rq$k_hat < nrow(cm)) {
      k <- rq$k_hat
      pre_mean  <- colMeans(cm[1:k, , drop = FALSE], na.rm = TRUE)
      post_mean <- colMeans(cm[(k + 1):nrow(cm), , drop = FALSE], na.rm = TRUE)
      lvl_shift <- mean(post_mean - pre_mean, na.rm = TRUE)
      dir <- ifelse(is.finite(lvl_shift) && lvl_shift >= 0, "up", "down")
    }
    
    out[i, `:=`(
      k_hat = rq$k_hat,
      cp_date = as.Date(rq$cp_date),
      direction = dir,
      
      D_hat = D,
      RMS_shift_vol_pts = 100 * sqrt(max(D, 0)),
      norm_quad = Vq,
      norm_ar   = Vr,
      
      Delta_boundary_quad = bq,
      Delta_boundary_ar   = br,
      RMS_boundary_vol_pts_quad = if (is.finite(bq)) 100 * sqrt(bq) else NA_real_,
      RMS_boundary_vol_pts_ar   = if (is.finite(br)) 100 * sqrt(br) else NA_real_,
      
      reject_quad = rej_q,
      reject_ar   = rej_r,
      reject      = rej,
      break_class = ifelse(rej, "structural", "irrelevant")
    )]
  }
  
  out
}

roll <- rolling_cpp(
  curves_mat = curves_mat,
  dates      = dates,
  window     = config$roll_window,
  step       = config$roll_step,
  eps_trim   = config$eps_trim,
  Delta      = Delta,
  Delta_rms  = config$roll_Delta_rms,
  q_quad     = q_quad,
  q_ar       = q_ar,
  require_both = config$roll_require_both
)

out_file <- sprintf("output/tables/btc_iv_rolling_cpp_W%d_alpha%.2f.csv",
                    config$roll_window, config$roll_alpha)
data.table::fwrite(roll, file = out_file)
message("Saved rolling results: ", out_file)

# Breakpoint frequency table (only windows where we reject)
hits <- roll[reject == TRUE & !is.na(cp_date),
             .N, by = .(cp_date, direction)][order(-N)]
hits_file <- sprintf("output/tables/btc_iv_rolling_cpp_hits_W%d_alpha%.2f.csv",
                     config$roll_window, config$roll_alpha)
data.table::fwrite(hits, file = hits_file)
message("Saved breakpoint hit counts: ", hits_file)

message("\nTop rolling-window break dates (when significant at alpha=0.10):")
print(head(hits, 20))
message("\nDone.")
