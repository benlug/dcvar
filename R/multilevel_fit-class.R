# ============================================================================
# S3 Class: dcVar_multilevel_fit
# ============================================================================

#' Construct a dcVar_multilevel_fit object
#' @noRd
new_dcVar_multilevel_fit <- function(fit, stan_data, N, vars, centered,
                                     person_means, priors, meta,
                                     standardized = FALSE,
                                     margins = "normal") {
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
      priors = priors,
      meta = meta
    ),
    class = c("dcVar_multilevel_fit", "dcVar_model_fit")
  )
}


#' S3 methods for dcVar_multilevel_fit objects
#'
#' @param x,object A `dcVar_multilevel_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcVar_multilevel_fit-methods
NULL

#' @describeIn dcVar_multilevel_fit-methods Print a concise overview.
#' @return Invisibly returns `x`.
#' @export
print.dcVar_multilevel_fit <- function(x, ...) {
  .print_fit_header(x, "Multilevel Copula VAR Model Fit")
  cat(sprintf("N = %d units, T = %d per unit\n", x$N, x$stan_data$T))
  .print_fit_footer(x)

  summ <- x$fit$summary()
  rho_row <- grep("^rho$", summ$variable)
  if (length(rho_row) > 0) {
    cat(sprintf("rho (global): %.3f\n", summ$mean[rho_row]))
  }
  invisible(x)
}


#' @describeIn dcVar_multilevel_fit-methods Produce a detailed summary.
#' @return A `dcVar_multilevel_summary` object (a list).
#' @export
summary.dcVar_multilevel_fit <- function(object, ...) {
  vp <- var_params(object)
  diag <- dcVar_diagnostics(object)
  re <- random_effects(object)

  out <- list(
    model = "multilevel",
    N = object$N,
    T = object$stan_data$T,
    var_params = vp,
    random_effects = re,
    diagnostics = diag
  )
  class(out) <- "dcVar_multilevel_summary"
  out
}


#' Print a dcVar_multilevel_summary object
#' @param x A `dcVar_multilevel_summary` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcVar_multilevel_summary <- function(x, ...) {
  cat("Multilevel Copula VAR Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("N = %d units, T = %d per unit\n\n", x$N, x$T))

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


#' @describeIn dcVar_multilevel_fit-methods Extract posterior means of
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
coef.dcVar_multilevel_fit <- function(object, ...) {
  summ <- object$fit$summary()
  list(
    phi_bar = .extract_coef(summ, "^phi_bar\\["),
    tau_phi = .extract_coef(summ, "^tau_phi\\["),
    sigma = .extract_coef(summ, "^sigma\\["),
    rho = .extract_coef(summ, "^rho$")
  )
}


#' @describeIn dcVar_multilevel_fit-methods Dispatch to a plot type.
#' @param type Character; one of `"random_effects"`, `"diagnostics"`.
#' @return A ggplot object.
#' @export
plot.dcVar_multilevel_fit <- function(x,
                                      type = c("random_effects", "diagnostics"),
                                      ...) {
  type <- match.arg(type)
  switch(type,
    random_effects = plot_random_effects(x, ...),
    diagnostics = plot_diagnostics(x, ...)
  )
}
