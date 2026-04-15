# ============================================================================
# S3 Class: dcvar_sem_fit
# ============================================================================

#' Construct a dcvar_sem_fit object
#' @noRd
new_dcvar_sem_fit <- function(fit, stan_data, vars, J, lambda, sigma_e,
                               indicators, margins = "normal",
                               skew_direction = NULL, backend = "rstan",
                               priors, meta) {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "sem",
      vars = vars,
      standardized = FALSE,
      margins = margins,
      J = J,
      lambda = lambda,
      sigma_e = sigma_e,
      indicators = indicators,
      skew_direction = skew_direction,
      backend = backend,
      priors = priors,
      meta = meta
    ),
    class = c("dcvar_sem_fit", "dcvar_model_fit")
  )
}


#' S3 methods for dcvar_sem_fit objects
#'
#' @param x,object A `dcvar_sem_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcvar_sem_fit-methods
NULL

#' @describeIn dcvar_sem_fit-methods Print a concise overview.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_sem_fit <- function(x, ...) {
  .print_fit_header(x, "SEM Copula VAR Model Fit")
  cat(sprintf("T = %d, J = %d indicators per latent\n", x$stan_data$T, x$J))
  .print_fit_footer(x)

  cat(sprintf("rho: %.3f\n", coef(x)$rho[[1]]))
  invisible(x)
}


#' @describeIn dcvar_sem_fit-methods Produce a detailed summary.
#' @return A `dcvar_sem_summary` object (a list).
#' @export
summary.dcvar_sem_fit <- function(object, ...) {
  vp <- var_params(object)
  diag <- dcvar_diagnostics(object)

  out <- list(
    model = "sem",
    T = object$stan_data$T,
    J = object$J,
    lambda = object$lambda,
    sigma_e = object$sigma_e,
    var_params = vp,
    diagnostics = diag
  )
  class(out) <- "dcvar_sem_summary"
  out
}


#' Print a dcvar_sem_summary object
#' @param x A `dcvar_sem_summary` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_sem_summary <- function(x, ...) {
  cat("SEM Copula VAR Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("T = %d, J = %d indicators per latent\n", x$T, x$J))
  cat(sprintf("Fixed lambda: %s\n", paste(round(x$lambda, 3), collapse = ", ")))
  cat(sprintf("Fixed sigma_e: %.3f\n\n", x$sigma_e))

  cat("Latent VAR Parameters:\n")
  if (!is.null(x$var_params$mu)) {
    cat("  mu:\n")
    print(x$var_params$mu[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  if (!is.null(x$var_params$Phi)) {
    cat("\n  Phi:\n")
    print(x$var_params$Phi[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  if (!is.null(x$var_params$sigma)) {
    cat("\n  sigma:\n")
    print(x$var_params$sigma[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  } else {
    .print_margin_params(x$var_params)
  }
  if (!is.null(x$var_params$rho)) {
    cat("\n  rho:\n")
    print(x$var_params$rho[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }

  cat("\nDiagnostics:\n")
  cat(sprintf("  Divergences: %d\n", x$diagnostics$n_divergent))
  cat(sprintf("  Max Rhat: %.3f\n", x$diagnostics$max_rhat))
  cat(sprintf("  Min ESS bulk: %.0f\n", x$diagnostics$min_ess_bulk))

  invisible(x)
}


#' @describeIn dcvar_sem_fit-methods Extract posterior means of
#'   latent VAR coefficients.
#' @return A named list of posterior means.
#' @export
coef.dcvar_sem_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  result <- list(
    mu = .extract_required_coef(summ, "^mu\\[", "mu", "coef.dcvar_sem_fit()"),
    Phi = .extract_required_coef(summ, "^Phi\\[", "Phi", "coef.dcvar_sem_fit()"),
    rho = .extract_required_coef(summ, "^rho$", "rho", "coef.dcvar_sem_fit()")
  )
  margins <- object$margins %||% "normal"
  margin_coefs <- if (identical(margins, "normal")) {
    list(sigma = .extract_required_coef(summ, "^sigma\\[", "sigma", "coef.dcvar_sem_fit()"))
  } else {
    .extract_margin_coefs(summ, margins)
  }
  c(result[1:2], margin_coefs, result[3])
}


#' @describeIn dcvar_sem_fit-methods Dispatch to a plot type.
#' @param type Character; one of `"latent_states"`, `"rho"`, `"diagnostics"`.
#' @return A ggplot object.
#' @export
plot.dcvar_sem_fit <- function(x,
                                type = c("latent_states", "rho", "diagnostics"),
                                ...) {
  type <- match.arg(type)
  switch(type,
    latent_states = plot_latent_states(x, ...),
    rho = plot_rho(x, ...),
    diagnostics = plot_diagnostics(x, ...)
  )
}
