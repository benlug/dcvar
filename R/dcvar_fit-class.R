# ============================================================================
# S3 Class: dcvar_fit
# ============================================================================

#' Construct a dcvar_fit object
#' @noRd
new_dcvar_fit <- function(fit, stan_data, vars, standardized,
                          margins = "normal", skew_direction = NULL,
                          backend = "rstan", priors, meta) {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "dcvar",
      vars = vars,
      standardized = standardized,
      margins = margins,
      skew_direction = skew_direction,
      backend = backend,
      priors = priors,
      meta = meta
    ),
    class = c("dcvar_fit", "dcvar_model_fit")
  )
}


#' S3 methods for dcvar_fit objects
#'
#' @param x,object A `dcvar_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcvar_fit-methods
NULL

#' @describeIn dcvar_fit-methods Print a concise overview of the DC-VAR fit.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_fit <- function(x, ...) {
  .print_fit_header(x, "DC-VAR Model Fit")
  cat(sprintf("n_time = %d, D = %d\n", x$stan_data$n_time, x$stan_data$D))
  .print_fit_footer(x)

  rho_df <- rho_trajectory(x)
  cat(sprintf("rho range: [%.3f, %.3f] (posterior mean)\n",
              min(rho_df$mean), max(rho_df$mean)))
  invisible(x)
}


#' @describeIn dcvar_fit-methods Produce a detailed summary including rho
#'   trajectory, VAR parameters, and diagnostics.
#' @param probs Numeric vector of quantile probabilities for the rho trajectory
#'   (default: `c(0.025, 0.5, 0.975)`).
#' @return A `dcvar_summary` object (a list).
#' @export
summary.dcvar_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  rho_df <- rho_trajectory(object, probs = probs)
  vp <- var_params(object)
  diag <- dcvar_diagnostics(object)

  out <- list(
    model = "dcvar",
    n_time = object$stan_data$n_time,
    D = object$stan_data$D,
    rho_trajectory = rho_df,
    rho_summary = data.frame(
      statistic = c("start", "end", "delta", "min", "max", "mean"),
      value = c(
        rho_df$mean[1],
        rho_df$mean[nrow(rho_df)],
        rho_df$mean[nrow(rho_df)] - rho_df$mean[1],
        min(rho_df$mean),
        max(rho_df$mean),
        mean(rho_df$mean)
      )
    ),
    var_params = vp,
    diagnostics = diag
  )
  class(out) <- "dcvar_summary"
  out
}


#' Print a dcvar_summary object
#'
#' @param x A `dcvar_summary` object as returned by [summary.dcvar_fit()].
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_summary <- function(x, ...) {
  cat("DC-VAR Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("n_time = %d, D = %d\n\n", x$n_time, x$D))

  cat("Rho Trajectory:\n")
  print(x$rho_summary, row.names = FALSE)

  cat("\nVAR(1) Parameters:\n")
  cat("  mu:\n")
  print(x$var_params$mu[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  cat("\n  Phi:\n")
  print(x$var_params$Phi[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  .print_margin_params(x$var_params)

  if (!is.null(x$var_params$sigma_omega)) {
    cat("\n  sigma_omega:\n")
    print(x$var_params$sigma_omega[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }

  cat("\nDiagnostics:\n")
  cat(sprintf("  Divergences: %d\n", x$diagnostics$n_divergent))
  cat(sprintf("  Max Rhat: %.3f\n", x$diagnostics$max_rhat))
  cat(sprintf("  Min ESS bulk: %.0f\n", x$diagnostics$min_ess_bulk))

  invisible(x)
}


#' @describeIn dcvar_fit-methods Extract posterior means of model coefficients.
#' @return A named list with elements `mu`, `Phi`, `sigma_eps`, and
#'   `sigma_omega`.
#' @export
coef.dcvar_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  result <- list(
    mu = .extract_required_coef(summ, "^mu\\[", "mu", "coef.dcvar_fit()"),
    Phi = .extract_required_coef(summ, "^Phi\\[", "Phi", "coef.dcvar_fit()")
  )
  # Margin-specific scale params before sigma_omega
  margins <- object$margins %||% "normal"
  result <- c(result, .extract_margin_coefs(summ, margins))
  result$sigma_omega <- .extract_required_coef(summ, "^sigma_omega$", "sigma_omega", "coef.dcvar_fit()")
  result
}


#' @describeIn dcvar_fit-methods Dispatch to a plot type: `"rho"`, `"phi"`,
#'   `"diagnostics"`, `"ppc"`, or `"pit"`.
#' @param type Character; one of `"rho"`, `"phi"`, `"diagnostics"`, `"ppc"`,
#'   or `"pit"`.
#' @return A ggplot object.
#' @export
plot.dcvar_fit <- function(x, type = c("rho", "phi", "diagnostics", "ppc", "pit"), ...) {
  type <- match.arg(type)
  switch(type,
    rho = plot_rho(x, ...),
    phi = plot_phi(x, ...),
    diagnostics = plot_diagnostics(x, ...),
    ppc = plot_ppc(x, ...),
    pit = plot_pit(x, ...)
  )
}
