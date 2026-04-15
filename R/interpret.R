# ============================================================================
# Clinical Interpretation Helper
# ============================================================================

# Default thresholds for classifying correlation strength (Cohen's conventions
# adapted for copula correlations). These are used by interpret_rho_trajectory()
# and can be overridden via the strength_breaks parameter.
.default_strength_breaks <- c(strong = 0.7, moderate = 0.4, weak = 0.2)

# Default thresholds for classifying the magnitude of change in a trajectory.
.default_magnitude_breaks <- c(large = 0.5, moderate = 0.3, small = 0.1)

#' Interpret a rho trajectory in clinical terms
#'
#' Generates a human-readable interpretation of the estimated rho trajectory,
#' describing the overall trend, magnitude of change, and key features.
#'
#' @param object A fitted model object (`dcvar_fit`, `dcvar_hmm_fit`, or
#'   `dcvar_constant_fit`).
#' @param threshold Minimum absolute change in posterior-mean rho to be
#'   considered "meaningful" (default: 0.1).
#' @param strength_breaks Named numeric vector of thresholds for classifying
#'   correlation strength (default: `c(strong = 0.7, moderate = 0.4,
#'   weak = 0.2)`). Values above the highest threshold are "strong", etc.
#' @param magnitude_breaks Named numeric vector of thresholds for classifying
#'   the magnitude of trajectory range (default: `c(large = 0.5, moderate = 0.3,
#'   small = 0.1)`).
#' @param fluctuation_threshold Proportion of sign changes in first differences
#'   to flag "substantial fluctuation" (default: 0.3).
#' @param ... Additional arguments (unused).
#'
#' @return A character string with the interpretation (invisibly). The
#'   interpretation is also printed to the console.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- simulate_dcvar(T = 100, rho_trajectory = rho_decreasing(100))
#' fit <- dcvar(sim$Y_df, vars = c("y1", "y2"))
#' interpret_rho_trajectory(fit)
#' }
interpret_rho_trajectory <- function(object, threshold = 0.1,
                                     strength_breaks = .default_strength_breaks,
                                     magnitude_breaks = .default_magnitude_breaks,
                                     fluctuation_threshold = 0.3,
                                     ...) {
  rho_df <- rho_trajectory(object)

  .classify <- function(value, breaks) {
    # breaks must be sorted descending (e.g., c(strong = 0.7, moderate = 0.4, weak = 0.2))
    for (i in seq_along(breaks)) {
      if (value > breaks[i]) return(names(breaks)[i])
    }
    "negligible"
  }

  if (object$model == "constant") {
    rho_val <- rho_df$mean[1]
    strength <- .classify(abs(rho_val), strength_breaks)
    direction <- if (rho_val > 0) "positive" else "negative"

    msg <- sprintf(
      "Constant copula model: The estimated correlation is %.3f (%s %s coupling).\nThis suggests a stable, time-invariant dependence structure between the variables.",
      rho_val, strength, direction
    )
    cat(msg, "\n")
    return(invisible(msg))
  }

  rho_mean <- rho_df$mean
  n <- length(rho_mean)
  rho_start <- rho_mean[1]
  rho_end <- rho_mean[n]
  delta <- rho_end - rho_start
  rho_min <- min(rho_mean)
  rho_max <- max(rho_mean)
  rho_range <- rho_max - rho_min

  # Classify overall trend
  if (abs(delta) < threshold) {
    trend <- "relatively stable"
  } else if (delta < -threshold) {
    trend <- "decreasing (decoupling)"
  } else {
    trend <- "increasing (strengthening coupling)"
  }

  # Classify magnitude
  magnitude <- .classify(rho_range, magnitude_breaks)

  # Detect non-monotonicity
  diffs <- diff(rho_mean)
  n_sign_changes <- sum(diff(sign(diffs)) != 0)
  nonmono <- if (n_sign_changes > (n - 2) * fluctuation_threshold) {
    " The trajectory shows substantial fluctuation, suggesting non-monotonic dynamics."
  } else {
    ""
  }

  model_label <- switch(object$model,
    dcvar = "DC-VAR",
    hmm = "HMM Copula",
    multilevel = "Multilevel Copula VAR",
    sem = "SEM Copula VAR",
    object$model
  )

  msg <- sprintf(
    paste0(
      "%s model: The dependence trajectory is %s (rho: %.3f -> %.3f, delta = %.3f).\n",
      "Range: [%.3f, %.3f] (%s variation).%s"
    ),
    model_label, trend, rho_start, rho_end, delta,
    rho_min, rho_max, magnitude, nonmono
  )

  if (object$model == "hmm") {
    states <- hmm_states(object)
    state_rhos <- states$rho_state$mean
    msg <- paste0(msg, sprintf(
      "\nHMM states: %s (rho = %s).",
      paste(paste0("State ", seq_along(state_rhos)), collapse = ", "),
      paste(sprintf("%.3f", state_rhos), collapse = ", ")
    ))
  }

  cat(msg, "\n")
  invisible(msg)
}
