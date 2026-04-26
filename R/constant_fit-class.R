# ============================================================================
# S3 Class: dcvar_constant_fit
# ============================================================================

#' Construct a dcvar_constant_fit object
#' @noRd
new_dcvar_constant_fit <- function(fit, stan_data, vars, standardized,
                                   margins = "normal", copula = "gaussian",
                                   skew_direction = NULL,
                                   backend = "rstan", priors, meta) {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "constant",
      vars = vars,
      standardized = standardized,
      margins = margins,
      copula = copula,
      skew_direction = skew_direction,
      backend = backend,
      priors = priors,
      meta = meta
    ),
    class = c("dcvar_constant_fit", "dcvar_model_fit")
  )
}


#' S3 methods for dcvar_constant_fit objects
#'
#' @param x,object A `dcvar_constant_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcvar_constant_fit-methods
NULL

#' @describeIn dcvar_constant_fit-methods Print a concise overview of the
#'   constant copula fit.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_constant_fit <- function(x, ...) {
  .print_fit_header(x, "Constant Copula Model Fit")
  cat(sprintf("n_time = %d, D = %d\n", x$stan_data$n_time, x$stan_data$D))
  .print_fit_footer(x)

  if (identical(x$copula %||% "gaussian", "clayton")) {
    dep_df <- dependence_summary(x)
    cat(sprintf("Kendall's tau: %.3f [%.3f, %.3f]\n",
                dep_df$mean[1], dep_df$q2.5[1], dep_df$q97.5[1]))
  } else {
    rho_df <- rho_trajectory(x)
    cat(sprintf("rho: %.3f [%.3f, %.3f]\n",
                rho_df$mean[1], rho_df$q2.5[1], rho_df$q97.5[1]))
  }
  invisible(x)
}


#' @describeIn dcvar_constant_fit-methods Produce a detailed summary including
#'   constant rho, VAR parameters, and diagnostics.
#' @param probs Numeric vector of quantile probabilities for the rho estimate
#'   (default: `c(0.025, 0.5, 0.975)`).
#' @return A `dcvar_constant_summary` object (a list).
#' @export
summary.dcvar_constant_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  copula <- object$copula %||% "gaussian"
  rho_df <- if (identical(copula, "gaussian")) rho_trajectory(object, probs = probs) else NULL
  dep_df <- dependence_summary(object, probs = probs)
  vp <- var_params(object)
  diag <- dcvar_diagnostics(object)

  out <- list(
    model = "constant",
    copula = copula,
    n_time = object$stan_data$n_time,
    D = object$stan_data$D,
    rho = rho_df,
    dependence = dep_df,
    var_params = vp,
    diagnostics = diag
  )
  class(out) <- "dcvar_constant_summary"
  out
}


#' Print a dcvar_constant_summary object
#'
#' @param x A `dcvar_constant_summary` object as returned by
#'   [summary.dcvar_constant_fit()].
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_constant_summary <- function(x, ...) {
  dep <- if (identical(x$copula %||% "gaussian", "gaussian")) x$rho else x$dependence
  quantile_cols <- grep("^q", names(dep), value = TRUE)
  lower_col <- quantile_cols[1]
  upper_col <- quantile_cols[length(quantile_cols)]

  cat("Constant Copula Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("n_time = %d, D = %d\n\n", x$n_time, x$D))

  dep_label <- if (identical(x$copula %||% "gaussian", "clayton")) {
    "Kendall's tau (constant)"
  } else {
    "Rho (constant)"
  }
  cat(sprintf("%s: %.3f", dep_label, dep$mean[1]))
  if (length(quantile_cols) >= 2) {
    cat(sprintf(" [%.3f, %.3f]", dep[[lower_col]][1], dep[[upper_col]][1]))
  }
  cat("\n\n")

  cat("VAR(1) Parameters:\n")
  cat("  mu:\n")
  print(x$var_params$mu[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  cat("\n  Phi:\n")
  print(x$var_params$Phi[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  .print_margin_params(x$var_params)

  cat("\nDiagnostics:\n")
  cat(sprintf("  Divergences: %d\n", x$diagnostics$n_divergent))
  cat(sprintf("  Max Rhat: %.3f\n", x$diagnostics$max_rhat))
  cat(sprintf("  Min ESS bulk: %.0f\n", x$diagnostics$min_ess_bulk))

  invisible(x)
}


#' @describeIn dcvar_constant_fit-methods Extract posterior means of model
#'   coefficients.
#' @return A named list with elements `mu`, `Phi`, `sigma_eps`, and `rho`.
#' @export
coef.dcvar_constant_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  result <- list(
    mu = .extract_required_coef(summ, "^mu\\[", "mu", "coef.dcvar_constant_fit()"),
    Phi = .extract_required_coef(summ, "^Phi\\[", "Phi", "coef.dcvar_constant_fit()")
  )
  # Margin-specific scale params before rho
  margins <- object$margins %||% "normal"
  result <- c(result, .extract_margin_coefs(summ, margins))
  if (identical(object$copula %||% "gaussian", "clayton")) {
    result$theta <- .extract_required_coef(summ, "^theta$", "theta", "coef.dcvar_constant_fit()")
  } else {
    result$rho <- .extract_required_coef(summ, "^rho$", "rho", "coef.dcvar_constant_fit()")
  }
  result
}


#' @describeIn dcvar_constant_fit-methods Dispatch to a plot type: `"rho"`,
#'   `"phi"`, `"diagnostics"`, `"ppc"`, or `"pit"`.
#' @param type Character; one of `"rho"`, `"phi"`, `"diagnostics"`, `"ppc"`,
#'   `"pit"`.
#' @return A ggplot object.
#' @export
plot.dcvar_constant_fit <- function(x, type = c("rho", "phi", "diagnostics", "ppc", "pit"), ...) {
  type <- match.arg(type)
  if (identical(type, "rho") && identical(x$copula %||% "gaussian", "clayton")) {
    cli_abort("Clayton fits do not have a Gaussian copula rho. Use {.fun dependence_summary} for Kendall's tau.")
  }
  switch(type,
    rho = plot_rho(x, ...),
    phi = plot_phi(x, ...),
    diagnostics = plot_diagnostics(x, ...),
    ppc = plot_ppc(x, ...),
    pit = plot_pit(x, ...)
  )
}
