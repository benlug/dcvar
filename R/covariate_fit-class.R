# ============================================================================
# S3 Class: dcvar_covariate_fit
# ============================================================================

#' Construct a dcvar_covariate_fit object
#' @noRd
new_dcvar_covariate_fit <- function(fit, stan_data, vars, covariates,
                                    standardized, standardized_covariates,
                                    drift, zero_init_eta,
                                    backend = "rstan", priors, meta) {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = if (drift) "dcvar_covariate" else "dcvar_covariate_nodrift",
      vars = vars,
      covariates = covariates,
      standardized = standardized,
      standardized_covariates = standardized_covariates,
      margins = "normal",
      drift = drift,
      zero_init_eta = zero_init_eta,
      backend = backend,
      priors = priors,
      meta = meta
    ),
    class = c("dcvar_covariate_fit", "dcvar_model_fit")
  )
}


#' S3 methods for covariate DC-VAR fits
#'
#' @param x,object A `dcvar_covariate_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcvar_covariate_fit-methods
NULL

#' @describeIn dcvar_covariate_fit-methods Print a concise overview of the
#'   covariate DC-VAR fit.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_covariate_fit <- function(x, ...) {
  .print_fit_header(x, "Covariate DC-VAR Model Fit")
  cat(sprintf("Covariates: %s\n", paste(x$covariates, collapse = ", ")))
  cat(sprintf("n_time = %d, D = %d, P = %d\n",
              x$stan_data$n_time, x$stan_data$D, x$stan_data$P))
  cat(sprintf("Residual drift: %s\n", if (x$drift) "yes" else "no"))
  .print_fit_footer(x)

  rho_df <- rho_trajectory(x)
  cat(sprintf("rho range: [%.3f, %.3f] (posterior mean)\n",
              min(rho_df$mean), max(rho_df$mean)))
  invisible(x)
}


#' @describeIn dcvar_covariate_fit-methods Produce a detailed summary including
#'   rho trajectory, VAR parameters, covariate effects, and diagnostics.
#' @param probs Numeric vector of quantile probabilities (default:
#'   `c(0.025, 0.5, 0.975)`).
#' @return A `dcvar_covariate_summary` object (a list).
#' @export
summary.dcvar_covariate_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  rho_df <- rho_trajectory(object, probs = probs)
  vp <- var_params(object)
  effects <- covariate_effects(object, probs = probs)
  diag <- dcvar_diagnostics(object)

  out <- list(
    model = object$model,
    n_time = object$stan_data$n_time,
    D = object$stan_data$D,
    P = object$stan_data$P,
    drift = object$drift,
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
    covariate_effects = effects,
    diagnostics = diag
  )
  class(out) <- "dcvar_covariate_summary"
  out
}


#' Print a dcvar_covariate_summary object
#'
#' @param x A `dcvar_covariate_summary` object as returned by `summary()`.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_covariate_summary <- function(x, ...) {
  cat("Covariate DC-VAR Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("n_time = %d, D = %d, P = %d\n", x$n_time, x$D, x$P))
  cat(sprintf("Residual drift: %s\n\n", if (x$drift) "yes" else "no"))

  cat("Rho Trajectory:\n")
  print(x$rho_summary, row.names = FALSE)

  cat("\nCovariate Effects:\n")
  print(x$covariate_effects, row.names = FALSE)

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


#' @describeIn dcvar_covariate_fit-methods Extract posterior means of model
#'   coefficients.
#' @return A named list with elements `mu`, `Phi`, `sigma_eps`, `beta_0`,
#'   `beta`, and, when `drift = TRUE`, `sigma_omega`.
#' @export
coef.dcvar_covariate_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  beta <- .extract_required_coef(summ, "^beta\\[", "beta", "coef.dcvar_covariate_fit()")
  names(beta) <- object$covariates

  result <- list(
    mu = .extract_required_coef(summ, "^mu\\[", "mu", "coef.dcvar_covariate_fit()"),
    Phi = .extract_required_coef(summ, "^Phi\\[", "Phi", "coef.dcvar_covariate_fit()"),
    sigma_eps = .extract_required_coef(summ, "^sigma_eps\\[", "sigma_eps", "coef.dcvar_covariate_fit()"),
    beta_0 = .extract_required_coef(summ, "^beta_0$", "beta_0", "coef.dcvar_covariate_fit()"),
    beta = beta
  )

  if (isTRUE(object$drift)) {
    result$sigma_omega <- .extract_required_coef(
      summ, "^sigma_omega$", "sigma_omega", "coef.dcvar_covariate_fit()"
    )
  }

  result
}


#' @describeIn dcvar_covariate_fit-methods Dispatch to a plot type: `"rho"`,
#'   `"phi"`, `"diagnostics"`, `"ppc"`, or `"pit"`.
#' @param type Character; one of `"rho"`, `"phi"`, `"diagnostics"`, `"ppc"`,
#'   or `"pit"`.
#' @return A ggplot object.
#' @export
plot.dcvar_covariate_fit <- function(x, type = c("rho", "phi", "diagnostics", "ppc", "pit"), ...) {
  type <- match.arg(type)
  switch(type,
    rho = plot_rho(x, ...),
    phi = plot_phi(x, ...),
    diagnostics = plot_diagnostics(x, ...),
    ppc = plot_ppc(x, ...),
    pit = plot_pit(x, ...)
  )
}
