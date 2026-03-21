# ============================================================================
# PIT (Probability Integral Transform) Diagnostics
# ============================================================================

utils::globalVariables("count")

#' Extract PIT values from a fitted model
#'
#' Computes approximate Probability Integral Transform values using posterior
#' mean residuals and posterior mean margin parameters. Large departures from
#' uniformity can indicate model misfit, but these are not exact posterior
#' predictive PIT values.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `time`, `variable`, `pit`.
#'
#' @details
#' PIT values are computed from posterior mean residuals and posterior mean
#' margin parameters. Treat them as a fast plug-in diagnostic rather than an
#' exact posterior predictive transform that integrates over full posterior
#' uncertainty.
#' @export
pit_values <- function(object, ...) {
  UseMethod("pit_values")
}

#' @rdname pit_values
#' @export
pit_values.default <- function(object, ...) {
  cli_abort("{.fun pit_values} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' @rdname pit_values
#' @export
pit_values.dcvar_model_fit <- function(object, ...) {
  margins <- object$margins %||% "normal"

  # Extract posterior mean residuals and parameters
  eps_draws <- posterior::as_draws_matrix(.fit_draws(object$fit, "eps", backend = object$backend))
  T_eff <- object$stan_data$T - 1
  D <- object$stan_data$D

  # Posterior mean residuals
  eps_mean <- matrix(NA_real_, T_eff, D)
  for (d in seq_len(D)) {
    cols <- paste0("eps[", seq_len(T_eff), ",", d, "]")
    eps_mean[, d] <- colMeans(eps_draws[, cols, drop = FALSE])
  }

  # Compute PIT values based on margin type
  pit_mat <- .pit_compute(object, eps_mean, margins)

  # Convert to long format
  var_names <- object$vars
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)
  data.frame(
    time = rep(time_values, D),
    variable = rep(var_names, each = T_eff),
    pit = as.vector(pit_mat)
  )
}


#' Internal: compute PIT values for given residuals and margin type
#' @noRd
.pit_compute <- function(object, eps_mean, margins) {
  T_eff <- nrow(eps_mean)
  D <- ncol(eps_mean)

  switch(margins,
    normal = {
      summ <- .fit_summary(object$fit, backend = object$backend)
      sigma_rows <- grep("^sigma_eps\\[", summ$variable)
      sigma_mean <- summ$mean[sigma_rows]
      pit_mat <- matrix(NA_real_, T_eff, D)
      for (d in seq_len(D)) {
        pit_mat[, d] <- stats::pnorm(eps_mean[, d] / sigma_mean[d])
      }
      pit_mat
    },
    exponential = {
      .pit_exponential(object, eps_mean)
    },
    skew_normal = {
      .pit_skew_normal(object, eps_mean)
    },
    gamma = {
      .pit_gamma(object, eps_mean)
    },
    cli_abort("Unknown margin type: {.val {margins}}")
  )
}


#' Internal: PIT for exponential margins
#' @noRd
.pit_exponential <- function(object, eps_mean) {
  T_eff <- nrow(eps_mean)
  D <- ncol(eps_mean)
  summ <- .fit_summary(object$fit, backend = object$backend)
  skew_dir <- object$skew_direction

  sigma_rows <- grep("^sigma_exp\\[", summ$variable)
  sigma_exp_mean <- summ$mean[sigma_rows]
  rate_exp <- 1.0 / sigma_exp_mean

  pit_mat <- matrix(NA_real_, T_eff, D)
  for (d in seq_len(D)) {
    x_shifted <- sigma_exp_mean[d] + skew_dir[d] * eps_mean[, d]
    u <- stats::pexp(x_shifted, rate = rate_exp[d])
    if (skew_dir[d] < 0) u <- 1.0 - u
    pit_mat[, d] <- u
  }
  pit_mat
}


#' Internal: PIT for skew-normal margins
#' @noRd
.pit_skew_normal <- function(object, eps_mean) {
  if (!requireNamespace("sn", quietly = TRUE)) {
    cli_abort("Package {.pkg sn} is required for skew-normal PIT values.")
  }
  T_eff <- nrow(eps_mean)
  D <- ncol(eps_mean)
  summ <- .fit_summary(object$fit, backend = object$backend)

  omega_rows <- grep("^omega\\[", summ$variable)
  delta_rows <- grep("^delta\\[", summ$variable)
  omega_mean <- summ$mean[omega_rows]
  delta_mean <- summ$mean[delta_rows]

  alpha_mean <- delta_mean / sqrt(1 - delta_mean^2)
  xi_mean <- -omega_mean * delta_mean * sqrt(2 / pi)

  pit_mat <- matrix(NA_real_, T_eff, D)
  for (d in seq_len(D)) {
    pit_mat[, d] <- sn::psn(eps_mean[, d],
      xi = xi_mean[d], omega = omega_mean[d], alpha = alpha_mean[d]
    )
  }
  pit_mat
}


#' Internal: PIT for gamma margins
#' @noRd
.pit_gamma <- function(object, eps_mean) {
  T_eff <- nrow(eps_mean)
  D <- ncol(eps_mean)
  summ <- .fit_summary(object$fit, backend = object$backend)
  skew_dir <- object$skew_direction

  sigma_rows <- grep("^sigma_gam\\[", summ$variable)
  shape_row <- grep("^shape_gam$", summ$variable)
  sigma_gam_mean <- summ$mean[sigma_rows]
  shape_gam_mean <- summ$mean[shape_row]

  sqrt_shape <- sqrt(shape_gam_mean)
  rate_gam <- sqrt_shape / sigma_gam_mean

  pit_mat <- matrix(NA_real_, T_eff, D)
  for (d in seq_len(D)) {
    mean_x <- sqrt_shape * sigma_gam_mean[d]
    x_shifted <- mean_x + skew_dir[d] * eps_mean[, d]
    u <- stats::pgamma(x_shifted, shape = shape_gam_mean, rate = rate_gam[d])
    if (skew_dir[d] < 0) u <- 1.0 - u
    pit_mat[, d] <- u
  }
  pit_mat
}


#' KS test for PIT uniformity
#'
#' Runs a Kolmogorov-Smirnov test per variable to assess whether PIT values
#' are approximately uniform. This is a heuristic check on the plug-in PIT
#' values returned by [pit_values()], not an exact posterior predictive test.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `variable`, `ks_statistic`, `p_value`, `n`.
#'
#' @details
#' This applies a Kolmogorov-Smirnov test to the approximate PIT values
#' returned by [pit_values()]. The result is a heuristic check and does not
#' account for serial dependence or full posterior uncertainty.
#' @export
pit_test <- function(object, ...) {
  UseMethod("pit_test")
}

#' @rdname pit_test
#' @export
pit_test.default <- function(object, ...) {
  cli_abort("{.fun pit_test} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' @rdname pit_test
#' @export
pit_test.dcvar_model_fit <- function(object, ...) {
  pit_df <- pit_values(object)

  results <- lapply(unique(pit_df$variable), function(v) {
    pit_v <- pit_df$pit[pit_df$variable == v]
    ks <- stats::ks.test(pit_v, "punif")
    data.frame(
      variable = v,
      ks_statistic = unname(ks$statistic),
      p_value = ks$p.value,
      n = length(pit_v)
    )
  })

  do.call(rbind, results)
}


#' @export
pit_values.dcvar_multilevel_fit <- function(object, ...) {
  cli_abort("PIT values are not yet supported for multilevel models.")
}

#' @export
pit_values.dcvar_sem_fit <- function(object, ...) {
  cli_abort("PIT values are not yet supported for SEM models.")
}

#' @export
pit_test.dcvar_multilevel_fit <- function(object, ...) {
  cli_abort("PIT tests are not yet supported for multilevel models.")
}

#' @export
pit_test.dcvar_sem_fit <- function(object, ...) {
  cli_abort("PIT tests are not yet supported for SEM models.")
}


#' Plot PIT histograms
#'
#' Creates faceted histograms of the approximate PIT values returned by
#' [pit_values()]. Under good model fit, these histograms should be roughly
#' uniform, but they remain plug-in diagnostics rather than exact posterior
#' predictive checks.
#'
#' @param object A fitted model object.
#' @param bins Number of histogram bins (default: 20).
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
plot_pit <- function(object, bins = 20, ...) {
  pit_df <- pit_values(object)

  ggplot2::ggplot(pit_df, ggplot2::aes(x = .data$pit)) +
    ggplot2::geom_histogram(
      bins = bins,
      binwidth = 1 / bins,
      boundary = 0,
      ggplot2::aes(y = ggplot2::after_stat(count / sum(count) * bins)),
      fill = "steelblue", alpha = 0.7, color = "white"
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::facet_wrap(~ variable) +
    ggplot2::labs(
      x = "PIT Value",
      y = "Density",
      title = "Approximate PIT Diagnostics",
      subtitle = "Red dashed line: uniform plug-in reference"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
