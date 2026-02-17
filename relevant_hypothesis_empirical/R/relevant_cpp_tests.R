# Relevant mean change-point tests for functional time series (discretized curves)
#
# Implements:
#   - Quadratic self-normalization (Dette et al. style)
#   - Adjusted-range self-normalization (your proposed normalizer)
#
# We work directly with the full curve evaluations on a fixed grid, i.e.
# no KL truncation / FPCA. The L^2 integrals are approximated by trapezoid weights.

trap_weights <- function(p) {
  # Trapezoidal-rule weights on [0,1] with p equally spaced points.
  stopifnot(p >= 2)
  h <- 1 / (p - 1)
  w <- rep(h, p)
  w[1] <- w[p] <- h / 2
  w
}

l2_sq <- function(x, w) {
  # x: numeric vector
  sum(w * x^2)
}

est_cp_matrix <- function(curves_mat, eps_trim = 0.05, w = NULL, min_seg = 2) {
  N <- nrow(curves_mat)
  P <- ncol(curves_mat)
  if (is.null(w)) w <- trap_weights(P)

  # 原来：k_min <- floor(N*eps_trim)+1
  #      k_max <- N - floor(N*eps_trim)
  # 改成：两侧至少 min_seg 个点
  k_min <- max(floor(N * eps_trim) + 1, min_seg)
  k_max <- min(N - floor(N * eps_trim), N - min_seg)

  if (k_max <= k_min) stop("Not enough observations for given eps_trim/min_seg")

  cs <- apply(curves_mat, 2, cumsum)
  total <- cs[N, ]

  best_k <- k_min
  best_val <- -Inf

  for (k in k_min:k_max) {
    mean_pre  <- cs[k, ] / k
    mean_post <- (total - cs[k, ]) / (N - k)
    diff <- mean_pre - mean_post
    f_hat <- (k / N) * (1 - k / N) * l2_sq(diff, w)
    if (is.finite(f_hat) && f_hat > best_val) {
      best_val <- f_hat
      best_k <- k
    }
  }
  list(k_hat = best_k, f_hat = best_val)
}

cpp_relevant_test <- function(curves_mat,
                             dates,
                             eps_trim = 0.05,
                             normalizer = c("quadratic", "adjusted_range"),
                             nPoints = 20) {
  normalizer <- match.arg(normalizer)
  stopifnot(is.matrix(curves_mat))
  N <- nrow(curves_mat)
  P <- ncol(curves_mat)
  if (length(dates) != N) stop("dates length must equal nrow(curves_mat)")

  w <- trap_weights(P)

  # 1) Estimate CP
  cp <- est_cp_matrix(curves_mat, eps_trim = eps_trim, w = w)
  k_hat <- cp$k_hat

  # 2) Build intValues over lambda grid
  n <- k_hat
  m <- N - k_hat
  if (n < 2 || m < 2) stop("Estimated change point leaves too short a segment")

  lambdas <- (1:nPoints) / nPoints
  intValues <- numeric(length(lambdas))

  # cumulative sums in each segment for fast partial means
  pre <- curves_mat[1:n, , drop = FALSE]
  post <- curves_mat[(n + 1):N, , drop = FALSE]
  cs_pre <- apply(pre, 2, cumsum)
  cs_post <- apply(post, 2, cumsum)
  total_pre <- cs_pre[n, ]
  total_post <- cs_post[m, ]

  for (i in seq_along(lambdas)) {
    lam <- lambdas[i]
    k1 <- floor(lam * n)
    k2 <- floor(lam * m)
    if (k1 < 1) k1 <- 1
    if (k2 < 1) k2 <- 1
    # IMPORTANT: scale by FULL segment lengths (n and m), not by partial counts (k1, k2)
    mean_pre  <- cs_pre[k1, ] / n
    mean_post <- cs_post[k2, ] / m
    diff <- mean_pre - mean_post
    intValues[i] <- l2_sq(diff, w)
    
   
  }

  statistic <- intValues[length(intValues)]
  lam_tilde <- (floor(lambdas * n) / n)
  integrand <- intValues - (lam_tilde^2) * statistic

  if (normalizer == "quadratic") {
    norm_val <- sqrt(sum(integrand^2) / (nPoints - 1))
  } else {
    norm_val <- max(integrand) - min(integrand)
  }

  list(
    k_hat = k_hat,
    cp_date = dates[k_hat],
    statistic = statistic,
    normalizer = norm_val,
    lambdas = lambdas,
    intValues = intValues,
    integrand = integrand
  )
}


# --- Pivotal quantiles for decision rule ---

simulate_BM_paths <- function(nsim = 5000, N = 1000) {
  # Returns matrix of size (N+1) x nsim of BM paths on grid 0,1/N,...,1
  inc <- matrix(rnorm(N * nsim), nrow = N, ncol = nsim)
  path <- apply(inc, 2, cumsum) / sqrt(N)
  rbind(rep(0, nsim), path)
}

calc_W_quadratic <- function(B, nPoints = 20) {
  # B: (N+1) x nsim matrix of BM paths
  N <- nrow(B) - 1
  nsim <- ncol(B)
  lambdas <- (1:nPoints) / nPoints
  idx <- floor(lambdas * N) + 1
  B1 <- B[N + 1, ]
  Bl <- B[idx, , drop = FALSE]
  lam_mat <- matrix(lambdas, nrow = nPoints, ncol = nsim)
  B1_mat <- matrix(rep(B1, each = nPoints), nrow = nPoints)
  values <- lam_mat * (Bl - lam_mat * B1_mat)
  res <- sqrt(colSums(values^2) / (nPoints - 1))
  B1 / res
}

calc_W_range <- function(B, nPoints = 20) {
  N <- nrow(B) - 1
  nsim <- ncol(B)
  lambdas <- (1:nPoints) / nPoints
  idx <- floor(lambdas * N) + 1
  B1 <- B[N + 1, ]
  Bl <- B[idx, , drop = FALSE]
  lam_mat <- matrix(lambdas, nrow = nPoints, ncol = nsim)
  B1_mat <- matrix(rep(B1, each = nPoints), nrow = nPoints)
  values <- lam_mat * (Bl - lam_mat * B1_mat)

  res <- apply(values, 2, function(v) diff(range(v)))
  B1 / res
}

get_cpp_q_values <- function(alpha = c(0.01, 0.05, 0.10),
                            cache_path = "data/qvalues_cpp.rds",
                            nsim = 5000,
                            N = 1000,
                            seed = 1) {
  # Returns a list with entries $quadratic and $adjusted_range, each a named
  # numeric vector of q_{1-alpha}.
  alpha <- sort(alpha)
  probs <- 1 - alpha

  if (file.exists(cache_path)) {
    obj <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(obj) &&
        is.list(obj) &&
        all(c("quadratic", "adjusted_range") %in% names(obj)) &&
        all(as.character(probs) %in% names(obj$quadratic)) &&
        all(as.character(probs) %in% names(obj$adjusted_range))) {
      return(obj)
    }
  }

  set.seed(seed)
  B <- simulate_BM_paths(nsim = nsim, N = N)
  Wq <- calc_W_quadratic(B)
  Wr <- calc_W_range(B)

  # Use names as the numeric probabilities ("0.99", "0.95", ...) to avoid the
  # default "99%" naming from quantile().
  q_quad <- as.numeric(stats::quantile(Wq, probs = probs, names = FALSE))
  q_rng  <- as.numeric(stats::quantile(Wr, probs = probs, names = FALSE))
  names(q_quad) <- as.character(probs)
  names(q_rng)  <- as.character(probs)

  out <- list(
    quadratic = q_quad,
    adjusted_range = q_rng,
    meta = list(nsim = nsim, N = N, seed = seed, nPoints = 20)
  )
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, cache_path)
  out
}


# --- Decision tables / boundaries ---

build_decision_table <- function(Delta_grid,
                                alpha,
                                res_quad,
                                res_rng,
                                qs_quad,
                                qs_rng) {
  alpha <- sort(alpha)
  probs <- as.character(1 - alpha)
  if (!all(probs %in% names(qs_quad))) stop("qs_quad names must include 1-alpha")
  if (!all(probs %in% names(qs_rng))) stop("qs_rng names must include 1-alpha")

  dt <- data.table::data.table(
    Delta = Delta_grid,
    RMS_vol_pts = 100 * sqrt(pmax(Delta_grid, 0))
  )

  for (a in alpha) {
    p <- as.character(1 - a)
    dt[[sprintf("quad_reject_alpha_%s", a)]] <- (res_quad$statistic > dt$Delta + as.numeric(qs_quad[p]) * res_quad$normalizer)
    dt[[sprintf("ar_reject_alpha_%s", a)]]   <- (res_rng$statistic  > dt$Delta + as.numeric(qs_rng[p])  * res_rng$normalizer)
  }

  dt
}

build_boundary_summary <- function(alpha,
                                  res_quad,
                                  res_rng,
                                  qs_quad,
                                  qs_rng) {
  alpha <- sort(alpha)
  probs <- as.character(1 - alpha)
  out <- data.table::data.table(
    alpha = alpha,
    q_quadratic = as.numeric(qs_quad[probs]),
    q_adjusted_range = as.numeric(qs_rng[probs]),
    k_hat = res_quad$k_hat,
    cp_date = as.character(res_quad$cp_date),
    statistic_D_hat = res_quad$statistic,
    normalizer_quad = res_quad$normalizer,
    normalizer_ar = res_rng$normalizer
  )

  out[, RMS_shift_vol_pts := 100 * sqrt(pmax(statistic_D_hat, 0))]

  out[, Delta_boundary_quad := pmax(0, statistic_D_hat - q_quadratic * normalizer_quad)]
  out[, Delta_boundary_ar   := pmax(0, statistic_D_hat - q_adjusted_range * normalizer_ar)]
  out[, RMS_boundary_vol_pts_quad := 100 * sqrt(Delta_boundary_quad)]
  out[, RMS_boundary_vol_pts_ar   := 100 * sqrt(Delta_boundary_ar)]
  out
}
