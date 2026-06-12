# ============================================================================
# Parameter Recovery Metrics
# ============================================================================

#' Compute rho trajectory recovery metrics
#'
#' Evaluates how well the estimated rho trajectory recovers the true values.
#'
#' @param rho_true Numeric vector of true rho values (length T-1).
#' @param rho_est Numeric vector of estimated rho (posterior mean).
#' @param rho_lower Numeric vector of lower interval bounds (2.5% quantiles).
#' @param rho_upper Numeric vector of upper interval bounds (97.5% quantiles).
#'
#' @return A named list with:
#'   - `bias`: mean bias
#'   - `relative_bias`: relative bias (%), the mean error normalized by the
#'     mean true value (by the mean absolute true value for mixed-sign
#'     trajectories; `NA` with a warning when all true values are near zero)
#'   - `coverage`: proportion of intervals containing true value
#'   - `interval_width`: mean interval width
#'   - `correlation`: Pearson correlation between true and estimated
#'   - `bias_start`, `bias_end`: bias at first/last time point
#'   - `coverage_start`, `coverage_end`: coverage at endpoints
#' @export
compute_rho_metrics <- function(rho_true, rho_est, rho_lower, rho_upper) {
  .validate_numeric_vector(rho_true, "rho_true")
  .validate_numeric_vector(rho_est, "rho_est")
  .validate_numeric_vector(rho_lower, "rho_lower")
  .validate_numeric_vector(rho_upper, "rho_upper")

  n <- length(rho_true)
  if (length(rho_est) != n || length(rho_lower) != n || length(rho_upper) != n) {
    cli_abort(c(
      "All rho metric vectors must have the same length.",
      "i" = "Expected {.val {n}} values for each of {.arg rho_est}, {.arg rho_lower}, and {.arg rho_upper}."
    ))
  }
  if (any(rho_lower > rho_upper)) {
    cli_abort(c(
      "{.arg rho_lower} must not exceed {.arg rho_upper}.",
      "i" = "Check that every credible interval is ordered lower <= upper."
    ))
  }

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
#' @param est_lower Lower interval bound (2.5% quantile).
#' @param est_upper Upper interval bound (97.5% quantile).
#'
#' @return A named list with `bias`, `relative_bias`, `covered`,
#'   `interval_width`.
#' @export
compute_param_metrics <- function(true_value, est_mean, est_lower, est_upper) {
  .validate_scalar_numeric(true_value, "true_value")
  .validate_scalar_numeric(est_mean, "est_mean")
  .validate_scalar_numeric(est_lower, "est_lower")
  .validate_scalar_numeric(est_upper, "est_upper")

  if (est_lower > est_upper) {
    cli_abort(c(
      "{.arg est_lower} must not exceed {.arg est_upper}.",
      "i" = "Check that the reported credible interval is ordered lower <= upper."
    ))
  }

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
#'   - `rho`: data frame of aggregated rho metrics (mean, SD, quantiles) for
#'     every field in the per-replication `rho` metric list
#'   - `n_reps`: number of replications
#' @export
aggregate_metrics <- function(metrics_list) {
  n_reps <- length(metrics_list)
  if (n_reps == 0) {
    cli_abort("{.arg metrics_list} must contain at least one replication.")
  }
  if (!is.list(metrics_list) || !all(vapply(metrics_list, function(x) is.list(x) && !is.null(x$rho), logical(1)))) {
    cli_abort("Each element of {.arg metrics_list} must be a list with a {.val rho} element.")
  }

  expected_metric_names <- c(
    "bias", "relative_bias", "coverage", "interval_width",
    "correlation", "bias_start", "bias_end",
    "coverage_start", "coverage_end"
  )
  metric_names <- names(metrics_list[[1]]$rho)
  if (is.null(metric_names) || any(metric_names == "")) {
    cli_abort("The {.val rho} element must be a named list of scalar metrics.")
  }
  if (!identical(metric_names, expected_metric_names)) {
    cli_abort("The {.val rho} element must contain the full set of compute_rho_metrics() outputs in the documented order.")
  }
  if (!all(vapply(metrics_list, function(x) identical(names(x$rho), expected_metric_names), logical(1)))) {
    cli_abort("All {.val rho} elements must contain the same metric names in the same order.")
  }

  rho_rows <- vapply(metrics_list, function(m) {
    vals <- m$rho[expected_metric_names]
    # NA is a documented value for relative_bias (null conditions) and for
    # correlation (degenerate series); the summaries below use na.rm = TRUE.
    ok <- vapply(
      vals,
      function(x) is.numeric(x) && length(x) == 1L && !is.nan(x) && !is.infinite(x),
      logical(1)
    )
    if (!all(ok)) {
      cli_abort("Each metric in {.val rho} must be a single finite or NA numeric value.")
    }
    unlist(vals, use.names = FALSE)
  }, numeric(length(expected_metric_names)))

  rho_df <- as.data.frame(t(rho_rows))
  names(rho_df) <- expected_metric_names

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
#'
#' Normalized form: mean error relative to the mean true value. Averaging
#' pointwise ratios instead would explode (or flip sign) whenever a single
#' true value is at or near zero -- e.g. a null-dependence rho trajectory --
#' silently dominating any aggregated summaries.
#'
#' When all true values share a sign, the denominator is `mean(true)`, which
#' preserves the conventional sign of relative bias for negative parameters
#' (attenuation toward zero is negative). For mixed-sign trajectories the
#' pointwise convention is ill-defined, so the mean error is normalized by the
#' mean absolute true value instead.
#' @noRd
.relative_bias <- function(estimated, true, epsilon = 1e-10) {
  denom <- mean(abs(true))
  if (denom < epsilon) {
    cli_warn(c(
      "Relative bias is undefined when all true values are (near) zero; returning {.val NA}.",
      "i" = "Use the absolute {.val bias} metric for null conditions."
    ))
    return(NA_real_)
  }
  if (all(true >= 0) || all(true <= 0)) {
    return(mean(estimated - true) / mean(true) * 100)
  }
  mean(estimated - true) / denom * 100
}

#' Internal: validate a numeric vector for metric calculations
#' @noRd
.validate_numeric_vector <- function(x, arg_name) {
  if (!is.numeric(x)) {
    cli_abort("{.arg {arg_name}} must be a numeric vector.")
  }
  if (length(x) == 0) {
    cli_abort("{.arg {arg_name}} must not be empty.")
  }
  if (any(!is.finite(x))) {
    cli_abort("{.arg {arg_name}} must contain only finite values.")
  }
}

#' Internal: validate a single finite numeric value for metric calculations
#' @noRd
.validate_scalar_numeric <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    cli_abort("{.arg {arg_name}} must be a single finite numeric value.")
  }
}
