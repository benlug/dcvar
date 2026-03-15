# ============================================================================
# S3 Class: dcVar_constant_fit
# ============================================================================

#' Construct a dcVar_constant_fit object
#' @noRd
new_dcVar_constant_fit <- function(fit, stan_data, vars, standardized,
                                   margins = "normal", skew_direction = NULL,
                                   priors, meta) {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "constant",
      vars = vars,
      standardized = standardized,
      margins = margins,
      skew_direction = skew_direction,
      priors = priors,
      meta = meta
    ),
    class = c("dcVar_constant_fit", "dcVar_model_fit")
  )
}


#' S3 methods for dcVar_constant_fit objects
#'
#' @param x,object A `dcVar_constant_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcVar_constant_fit-methods
NULL

#' @describeIn dcVar_constant_fit-methods Print a concise overview of the
#'   constant copula fit.
#' @return Invisibly returns `x`.
#' @export
print.dcVar_constant_fit <- function(x, ...) {
  .print_fit_header(x, "Constant Copula Model Fit")
  cat(sprintf("T = %d, D = %d\n", x$stan_data$T, x$stan_data$D))
  .print_fit_footer(x)

  rho_df <- rho_trajectory(x)
  cat(sprintf("rho: %.3f [%.3f, %.3f]\n",
              rho_df$mean[1], rho_df$q2.5[1], rho_df$q97.5[1]))
  invisible(x)
}


#' @describeIn dcVar_constant_fit-methods Produce a detailed summary including
#'   constant rho, VAR parameters, and diagnostics.
#' @param probs Numeric vector of quantile probabilities for the rho estimate
#'   (default: `c(0.025, 0.5, 0.975)`).
#' @return A `dcVar_constant_summary` object (a list).
#' @export
summary.dcVar_constant_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  rho_df <- rho_trajectory(object, probs = probs)
  vp <- var_params(object)
  diag <- dcVar_diagnostics(object)

  out <- list(
    model = "constant",
    T = object$stan_data$T,
    D = object$stan_data$D,
    rho = rho_df,
    var_params = vp,
    diagnostics = diag
  )
  class(out) <- "dcVar_constant_summary"
  out
}


#' Print a dcVar_constant_summary object
#'
#' @param x A `dcVar_constant_summary` object as returned by
#'   [summary.dcVar_constant_fit()].
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcVar_constant_summary <- function(x, ...) {
  quantile_cols <- grep("^q", names(x$rho), value = TRUE)
  lower_col <- quantile_cols[1]
  upper_col <- quantile_cols[length(quantile_cols)]

  cat("Constant Copula Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("T = %d, D = %d\n\n", x$T, x$D))

  cat(sprintf("Rho (constant): %.3f", x$rho$mean[1]))
  if (length(quantile_cols) >= 2) {
    cat(sprintf(" [%.3f, %.3f]", x$rho[[lower_col]][1], x$rho[[upper_col]][1]))
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


#' @describeIn dcVar_constant_fit-methods Extract posterior means of model
#'   coefficients.
#' @return A named list with elements `mu`, `Phi`, `sigma_eps`, and `rho`.
#' @export
coef.dcVar_constant_fit <- function(object, ...) {
  summ <- object$fit$summary()
  result <- list(
    mu = .extract_coef(summ, "^mu\\["),
    Phi = .extract_coef(summ, "^Phi\\[")
  )
  # Margin-specific scale params before rho
  margins <- object$margins %||% "normal"
  result <- c(result, .extract_margin_coefs(summ, margins))
  result$rho <- .extract_coef(summ, "^rho$")
  result
}


#' @describeIn dcVar_constant_fit-methods Dispatch to a plot type: `"rho"`,
#'   `"phi"`, `"diagnostics"`, `"ppc"`, or `"pit"`.
#' @param type Character; one of `"rho"`, `"phi"`, `"diagnostics"`, `"ppc"`,
#'   `"pit"`.
#' @return A ggplot object.
#' @export
plot.dcVar_constant_fit <- function(x, type = c("rho", "phi", "diagnostics", "ppc", "pit"), ...) {
  type <- match.arg(type)
  switch(type,
    rho = plot_rho(x, ...),
    phi = plot_phi(x, ...),
    diagnostics = plot_diagnostics(x, ...),
    ppc = plot_ppc(x, ...),
    pit = plot_pit(x, ...)
  )
}
