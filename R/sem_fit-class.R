# ============================================================================
# S3 Class: dcVar_sem_fit
# ============================================================================

#' Construct a dcVar_sem_fit object
#' @noRd
new_dcVar_sem_fit <- function(fit, stan_data, vars, J, lambda, sigma_e,
                               indicators, priors, meta) {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "sem",
      vars = vars,
      standardized = FALSE,
      margins = "normal",
      J = J,
      lambda = lambda,
      sigma_e = sigma_e,
      indicators = indicators,
      priors = priors,
      meta = meta
    ),
    class = c("dcVar_sem_fit", "dcVar_model_fit")
  )
}


#' S3 methods for dcVar_sem_fit objects
#'
#' @param x,object A `dcVar_sem_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcVar_sem_fit-methods
NULL

#' @describeIn dcVar_sem_fit-methods Print a concise overview.
#' @return Invisibly returns `x`.
#' @export
print.dcVar_sem_fit <- function(x, ...) {
  .print_fit_header(x, "SEM Copula VAR Model Fit")
  cat(sprintf("T = %d, J = %d indicators per latent\n", x$stan_data$T, x$J))
  .print_fit_footer(x)

  summ <- x$fit$summary()
  rho_row <- grep("^rho$", summ$variable)
  if (length(rho_row) > 0) {
    cat(sprintf("rho: %.3f\n", summ$mean[rho_row]))
  }
  invisible(x)
}


#' @describeIn dcVar_sem_fit-methods Produce a detailed summary.
#' @return A `dcVar_sem_summary` object (a list).
#' @export
summary.dcVar_sem_fit <- function(object, ...) {
  vp <- var_params(object)
  diag <- dcVar_diagnostics(object)

  out <- list(
    model = "sem",
    T = object$stan_data$T,
    J = object$J,
    lambda = object$lambda,
    sigma_e = object$sigma_e,
    var_params = vp,
    diagnostics = diag
  )
  class(out) <- "dcVar_sem_summary"
  out
}


#' Print a dcVar_sem_summary object
#' @param x A `dcVar_sem_summary` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcVar_sem_summary <- function(x, ...) {
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


#' @describeIn dcVar_sem_fit-methods Extract posterior means of
#'   latent VAR coefficients.
#' @return A named list of posterior means.
#' @export
coef.dcVar_sem_fit <- function(object, ...) {
  summ <- object$fit$summary()
  list(
    mu = .extract_coef(summ, "^mu\\["),
    Phi = .extract_coef(summ, "^Phi\\["),
    sigma = .extract_coef(summ, "^sigma\\["),
    rho = .extract_coef(summ, "^rho$")
  )
}


#' @describeIn dcVar_sem_fit-methods Dispatch to a plot type.
#' @param type Character; one of `"latent_states"`, `"rho"`, `"diagnostics"`.
#' @return A ggplot object.
#' @export
plot.dcVar_sem_fit <- function(x,
                                type = c("latent_states", "rho", "diagnostics"),
                                ...) {
  type <- match.arg(type)
  switch(type,
    latent_states = plot_latent_states(x, ...),
    rho = plot_rho(x, ...),
    diagnostics = plot_diagnostics(x, ...)
  )
}
