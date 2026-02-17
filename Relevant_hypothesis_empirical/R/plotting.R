# Plot helpers (ggplot2)

plot_mean_smiles <- function(u_grid,
                             curves_mat,
                             dates,
                             k_hat,
                             out_dir = "output/figures") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  N <- nrow(curves_mat)
  stopifnot(length(u_grid) == ncol(curves_mat))
  stopifnot(length(dates) == N)
  stopifnot(k_hat >= 2 && k_hat <= N - 2)

  pre_mean <- colMeans(curves_mat[1:k_hat, , drop = FALSE], na.rm = TRUE)
  post_mean <- colMeans(curves_mat[(k_hat + 1):N, , drop = FALSE], na.rm = TRUE)

  df_pre <- data.table::data.table(u = u_grid, iv = pre_mean)
  df_post <- data.table::data.table(u = u_grid, iv = post_mean)
  df_both <- data.table::rbindlist(list(
    data.table::data.table(u = u_grid, iv = pre_mean, regime = "Pre"),
    data.table::data.table(u = u_grid, iv = post_mean, regime = "Post")
  ))

  p1 <- ggplot2::ggplot(df_pre, ggplot2::aes(x = u, y = iv)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = sprintf("Mean IV smile (pre-break) – up to %s", as.character(dates[k_hat])),
      x = "log-moneyness u = log(K/S)",
      y = "Implied volatility (decimal)"
    ) +
    ggplot2::theme_minimal()

  p2 <- ggplot2::ggplot(df_post, ggplot2::aes(x = u, y = iv)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = sprintf("Mean IV smile (post-break) – after %s", as.character(dates[k_hat])),
      x = "log-moneyness u = log(K/S)",
      y = "Implied volatility (decimal)"
    ) +
    ggplot2::theme_minimal()

  p3 <- ggplot2::ggplot(df_both, ggplot2::aes(x = u, y = iv, linetype = regime)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = sprintf("Mean IV smiles (pre vs post), break at %s", as.character(dates[k_hat])),
      x = "log-moneyness u = log(K/S)",
      y = "Implied volatility (decimal)"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(filename = file.path(out_dir, "btc_iv_mean_pre.png"), plot = p1, width = 7, height = 4, dpi = 200)
  ggplot2::ggsave(filename = file.path(out_dir, "btc_iv_mean_post.png"), plot = p2, width = 7, height = 4, dpi = 200)
  ggplot2::ggsave(filename = file.path(out_dir, "btc_iv_means_overlay.png"), plot = p3, width = 7, height = 4, dpi = 200)
}


plot_btc_spot_with_cp <- function(spot,
                                 sample_dates,
                                 k_hat,
                                 out_file = "output/figures/btc_spot_with_cp.png") {
  stopifnot(is.data.table(spot), inherits(sample_dates, "Date"))
  cp_date <- sample_dates[k_hat]
  df <- spot[, .(date, close)]
  df <- df[date >= min(sample_dates) & date <= max(sample_dates)]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = date, y = close)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = as.numeric(cp_date), linetype = "dashed") +
    ggplot2::labs(
      title = sprintf("BTCUSDT (Binance) daily close with estimated break: %s", as.character(cp_date)),
      x = "Date",
      y = "Close price"
    ) +
    ggplot2::theme_minimal()

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = out_file, plot = p, width = 7, height = 4, dpi = 200)
}
