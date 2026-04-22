# ============================================================================
# S3 Class: dcvar_multilevel_fit
# ============================================================================

#' Construct a dcvar_multilevel_fit object
#' @noRd
new_dcvar_multilevel_fit <- function(fit, stan_data, N, vars, centered,
                                     person_means, priors, meta,
                                     standardized = FALSE,
                                     margins = "normal",
                                     backend = "rstan") {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "multilevel",
      N = N,
      vars = vars,
      standardized = standardized,
      margins = margins,
      centered = centered,
      person_means = person_means,
      backend = backend,
      priors = priors,
      meta = meta
    ),
    class = c("dcvar_multilevel_fit", "dcvar_model_fit")
  )
}


#' S3 methods for dcvar_multilevel_fit objects
#'
#' @param x,object A `dcvar_multilevel_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcvar_multilevel_fit-methods
NULL

#' @describeIn dcvar_multilevel_fit-methods Print a concise overview.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_multilevel_fit <- function(x, ...) {
  .print_fit_header(x, "Multilevel Copula VAR Model Fit")
  cat(sprintf("N = %d units, n_time = %d per unit\n", x$N, x$stan_data$n_time))
  .print_fit_footer(x)

  cat(sprintf("rho (global): %.3f\n", coef(x)$rho[[1]]))
  invisible(x)
}


#' @describeIn dcvar_multilevel_fit-methods Produce a detailed summary.
#' @return A `dcvar_multilevel_summary` object (a list).
#' @export
summary.dcvar_multilevel_fit <- function(object, ...) {
  vp <- var_params(object)
  diag <- dcvar_diagnostics(object)
  re <- random_effects(object)

  out <- list(
    model = "multilevel",
    N = object$N,
    n_time = object$stan_data$n_time,
    var_params = vp,
    random_effects = re,
    diagnostics = diag
  )
  class(out) <- "dcvar_multilevel_summary"
  out
}


#' Print a dcvar_multilevel_summary object
#' @param x A `dcvar_multilevel_summary` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_multilevel_summary <- function(x, ...) {
  cat("Multilevel Copula VAR Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("N = %d units, n_time = %d per unit\n\n", x$N, x$n_time))

  cat("Population-Level Parameters:\n")
  if (!is.null(x$var_params$phi_bar)) {
    cat("  phi_bar:\n")
    print(x$var_params$phi_bar[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  if (!is.null(x$var_params$tau_phi)) {
    cat("\n  tau_phi:\n")
    print(x$var_params$tau_phi[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
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


#' @describeIn dcvar_multilevel_fit-methods Extract posterior means of
#'   population-level coefficients.
#' @details
#' Unlike single-level models (where `coef()` returns `$Phi`), the multilevel
#' model returns hierarchical parameters:
#' \describe{
#'   \item{`phi_bar`}{Population-mean VAR coefficients (analogous to `Phi`
#'     in single-level models, vectorised as phi11, phi12, phi21, phi22).}
#'   \item{`tau_phi`}{Between-unit SD of VAR coefficients.}
#'   \item{`sigma`}{Innovation SDs.}
#'   \item{`rho`}{Copula correlation (constant across units).}
#' }
#' Use [random_effects()] to obtain unit-specific VAR coefficients.
#' @return A named list of posterior means.
#' @export
coef.dcvar_multilevel_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  list(
    phi_bar = .extract_required_coef(summ, "^phi_bar\\[", "phi_bar", "coef.dcvar_multilevel_fit()"),
    tau_phi = .extract_required_coef(summ, "^tau_phi\\[", "tau_phi", "coef.dcvar_multilevel_fit()"),
    sigma = .extract_required_coef(summ, "^sigma\\[", "sigma", "coef.dcvar_multilevel_fit()"),
    rho = .extract_required_coef(summ, "^rho$", "rho", "coef.dcvar_multilevel_fit()")
  )
}


#' @describeIn dcvar_multilevel_fit-methods Dispatch to a plot type.
#' @param type Character; one of `"random_effects"`, `"diagnostics"`.
#' @return A ggplot object.
#' @export
plot.dcvar_multilevel_fit <- function(x,
                                      type = c("random_effects", "diagnostics"),
                                      ...) {
  type <- match.arg(type)
  switch(type,
    random_effects = plot_random_effects(x, ...),
    diagnostics = plot_diagnostics(x, ...)
  )
}
