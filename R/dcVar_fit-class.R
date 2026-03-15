# ============================================================================
# S3 Class: dcVar_fit
# ============================================================================

#' Construct a dcVar_fit object
#' @noRd
new_dcVar_fit <- function(fit, stan_data, vars, standardized,
                          margins = "normal", skew_direction = NULL,
                          priors, meta) {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "dcVar",
      vars = vars,
      standardized = standardized,
      margins = margins,
      skew_direction = skew_direction,
      priors = priors,
      meta = meta
    ),
    class = c("dcVar_fit", "dcVar_model_fit")
  )
}


#' S3 methods for dcVar_fit objects
#'
#' @param x,object A `dcVar_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcVar_fit-methods
NULL

#' @describeIn dcVar_fit-methods Print a concise overview of the DC-VAR fit.
#' @return Invisibly returns `x`.
#' @export
print.dcVar_fit <- function(x, ...) {
  .print_fit_header(x, "DC-VAR Model Fit")
  cat(sprintf("T = %d, D = %d\n", x$stan_data$T, x$stan_data$D))
  .print_fit_footer(x)

  rho_df <- rho_trajectory(x)
  cat(sprintf("rho range: [%.3f, %.3f] (posterior mean)\n",
              min(rho_df$mean), max(rho_df$mean)))
  invisible(x)
}


#' @describeIn dcVar_fit-methods Produce a detailed summary including rho
#'   trajectory, VAR parameters, and diagnostics.
#' @param probs Numeric vector of quantile probabilities for the rho trajectory
#'   (default: `c(0.025, 0.5, 0.975)`).
#' @return A `dcVar_summary` object (a list).
#' @export
summary.dcVar_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  rho_df <- rho_trajectory(object, probs = probs)
  vp <- var_params(object)
  diag <- dcVar_diagnostics(object)

  out <- list(
    model = "dcVar",
    T = object$stan_data$T,
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
  class(out) <- "dcVar_summary"
  out
}


#' Print a dcVar_summary object
#'
#' @param x A `dcVar_summary` object as returned by [summary.dcVar_fit()].
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcVar_summary <- function(x, ...) {
  cat("DC-VAR Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("T = %d, D = %d\n\n", x$T, x$D))

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


#' @describeIn dcVar_fit-methods Extract posterior means of model coefficients.
#' @return A named list with elements `mu`, `Phi`, `sigma_eps`, and
#'   `sigma_omega`.
#' @export
coef.dcVar_fit <- function(object, ...) {
  summ <- object$fit$summary()
  result <- list(
    mu = .extract_coef(summ, "^mu\\["),
    Phi = .extract_coef(summ, "^Phi\\[")
  )
  # Margin-specific scale params before sigma_omega
  margins <- object$margins %||% "normal"
  result <- c(result, .extract_margin_coefs(summ, margins))
  result$sigma_omega <- .extract_coef(summ, "^sigma_omega$")
  result
}


#' @describeIn dcVar_fit-methods Dispatch to a plot type: `"rho"`, `"phi"`,
#'   `"diagnostics"`, `"ppc"`, or `"pit"`.
#' @param type Character; one of `"rho"`, `"phi"`, `"diagnostics"`, `"ppc"`,
#'   or `"pit"`.
#' @return A ggplot object.
#' @export
plot.dcVar_fit <- function(x, type = c("rho", "phi", "diagnostics", "ppc", "pit"), ...) {
  type <- match.arg(type)
  switch(type,
    rho = plot_rho(x, ...),
    phi = plot_phi(x, ...),
    diagnostics = plot_diagnostics(x, ...),
    ppc = plot_ppc(x, ...),
    pit = plot_pit(x, ...)
  )
}
