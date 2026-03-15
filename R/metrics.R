# ============================================================================
# Parameter Recovery Metrics
# ============================================================================

#' Compute rho trajectory recovery metrics
#'
#' Evaluates how well the estimated rho trajectory recovers the true values.
#'
#' @param rho_true Numeric vector of true rho values (length T-1).
#' @param rho_est Numeric vector of estimated rho (posterior mean).
#' @param rho_lower Numeric vector of lower CI bounds (2.5%).
#' @param rho_upper Numeric vector of upper CI bounds (97.5%).
#'
#' @return A named list with:
#'   - `bias`: mean bias
#'   - `relative_bias`: mean relative bias (%)
#'   - `coverage`: proportion of CIs containing true value
#'   - `interval_width`: mean CI width
#'   - `correlation`: Pearson correlation between true and estimated
#'   - `bias_start`, `bias_end`: bias at first/last time point
#'   - `coverage_start`, `coverage_end`: coverage at endpoints
#' @export
compute_rho_metrics <- function(rho_true, rho_est, rho_lower, rho_upper) {
  n <- length(rho_true)
  list(
    bias = mean(rho_est - rho_true),
    relative_bias = .relative_bias(rho_est, rho_true),
    coverage = mean(rho_true >= rho_lower & rho_true <= rho_upper),
    interval_width = mean(rho_upper - rho_lower),
    correlation = .safe_cor(rho_est, rho_true),
    bias_start = rho_est[1] - rho_true[1],
    bias_end = rho_est[n] - rho_true[n],
    coverage_start = as.numeric(rho_true[1] >= rho_lower[1] & rho_true[1] <= rho_upper[1]),
    coverage_end = as.numeric(rho_true[n] >= rho_lower[n] & rho_true[n] <= rho_upper[n])
  )
}


#' Compute scalar parameter recovery metrics
#'
#' @param true_value True parameter value.
#' @param est_mean Estimated posterior mean.
#' @param est_lower Lower CI bound (2.5%).
#' @param est_upper Upper CI bound (97.5%).
#'
#' @return A named list with `bias`, `relative_bias`, `covered`,
#'   `interval_width`.
#' @export
compute_param_metrics <- function(true_value, est_mean, est_lower, est_upper) {
  list(
    bias = est_mean - true_value,
    relative_bias = .relative_bias(est_mean, true_value),
    covered = as.numeric(true_value >= est_lower & true_value <= est_upper),
    interval_width = est_upper - est_lower
  )
}


#' Aggregate metrics across simulation replications
#'
#' @param metrics_list A list of metric lists (one per replication), each
#'   containing a `$rho` element as returned by [compute_rho_metrics()].
#'
#' @return A named list with:
#'   - `rho`: data frame of aggregated rho metrics (mean, SD, quantiles)
#'   - `n_reps`: number of replications
#' @export
aggregate_metrics <- function(metrics_list) {
  n_reps <- length(metrics_list)
  if (n_reps == 0) {
    cli_abort("{.arg metrics_list} must contain at least one replication.")
  }
  if (!all(vapply(metrics_list, function(x) is.list(x) && !is.null(x$rho), logical(1)))) {
    cli_abort("Each element of {.arg metrics_list} must be a list with a {.val rho} element.")
  }

  rho_df <- do.call(rbind, lapply(metrics_list, function(m) {
    data.frame(
      bias = m$rho$bias,
      relative_bias = m$rho$relative_bias,
      coverage = m$rho$coverage,
      interval_width = m$rho$interval_width,
      correlation = m$rho$correlation
    )
  }))

  summarise_col <- function(x) {
    c(mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      q2.5 = unname(quantile(x, 0.025, na.rm = TRUE)),
      q97.5 = unname(quantile(x, 0.975, na.rm = TRUE)))
  }

  rho_summ <- as.data.frame(t(sapply(rho_df, summarise_col)))
  rho_summ$metric <- rownames(rho_summ)
  rownames(rho_summ) <- NULL
  rho_summ <- rho_summ[, c("metric", "mean", "sd", "q2.5", "q97.5")]

  list(rho = rho_summ, n_reps = n_reps)
}


#' Internal: compute relative bias (%)
#' @noRd
.relative_bias <- function(estimated, true, epsilon = 1e-10) {
  denom <- ifelse(abs(true) < epsilon, epsilon, true)
  mean((estimated - true) / denom) * 100
}
