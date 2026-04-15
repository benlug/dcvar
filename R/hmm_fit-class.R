# ============================================================================
# S3 Class: dcvar_hmm_fit
# ============================================================================

#' Construct a dcvar_hmm_fit object
#' @noRd
new_dcvar_hmm_fit <- function(fit, stan_data, K, vars, standardized,
                              margins = "normal", skew_direction = NULL,
                              backend = "rstan", priors, meta) {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "hmm",
      K = K,
      vars = vars,
      standardized = standardized,
      margins = margins,
      skew_direction = skew_direction,
      backend = backend,
      priors = priors,
      meta = meta
    ),
    class = c("dcvar_hmm_fit", "dcvar_model_fit")
  )
}


#' S3 methods for dcvar_hmm_fit objects
#'
#' @param x,object A `dcvar_hmm_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcvar_hmm_fit-methods
NULL

#' @describeIn dcvar_hmm_fit-methods Print a concise overview of the HMM fit.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_hmm_fit <- function(x, ...) {
  .print_fit_header(x, "HMM Copula Model Fit")
  cat(sprintf("T = %d, D = %d, K = %d states\n",
              x$stan_data$T, x$stan_data$D, x$K))
  .print_fit_footer(x)

  states <- hmm_states(x)
  cat("\nState-specific rho:\n")
  for (k in seq_along(states$rho_state$mean)) {
    cat(sprintf("  State %d: %.3f [%.3f, %.3f]\n",
                k, states$rho_state$mean[k],
                states$rho_state$lower[k], states$rho_state$upper[k]))
  }

  cat("\nTransition matrix diagonal: ")
  cat(sprintf("%.3f", diag(states$A)), sep = ", ")
  cat("\n")

  invisible(x)
}


#' @describeIn dcvar_hmm_fit-methods Produce a detailed summary including
#'   state information, VAR parameters, and diagnostics.
#' @param probs Numeric vector of quantile probabilities for the rho trajectory
#'   (default: `c(0.025, 0.5, 0.975)`).
#' @return A `dcvar_hmm_summary` object (a list).
#' @export
summary.dcvar_hmm_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  rho_df <- rho_trajectory(object, probs = probs)
  vp <- var_params(object)
  states <- hmm_states(object)
  diag <- dcvar_diagnostics(object)

  out <- list(
    model = "hmm",
    T = object$stan_data$T,
    D = object$stan_data$D,
    K = object$K,
    rho_trajectory = rho_df,
    states = states,
    var_params = vp,
    diagnostics = diag
  )
  class(out) <- "dcvar_hmm_summary"
  out
}


#' Print a dcvar_hmm_summary object
#'
#' @param x A `dcvar_hmm_summary` object as returned by
#'   [summary.dcvar_hmm_fit()].
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_hmm_summary <- function(x, ...) {
  cat("HMM Copula Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("T = %d, D = %d, K = %d states\n\n", x$T, x$D, x$K))

  cat("State-Specific Rho:\n")
  for (k in seq_along(x$states$rho_state$mean)) {
    cat(sprintf("  State %d: %.3f [%.3f, %.3f]\n",
                k, x$states$rho_state$mean[k],
                x$states$rho_state$lower[k], x$states$rho_state$upper[k]))
  }

  cat("\nTransition Matrix:\n")
  print(round(x$states$A, 3))

  cat("\nVAR(1) Parameters:\n")
  cat("  mu:\n")
  print(x$var_params$mu[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  cat("\n  Phi:\n")
  print(x$var_params$Phi[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)

  cat("\nDiagnostics:\n")
  cat(sprintf("  Divergences: %d\n", x$diagnostics$n_divergent))
  cat(sprintf("  Max Rhat: %.3f\n", x$diagnostics$max_rhat))
  cat(sprintf("  Min ESS bulk: %.0f\n", x$diagnostics$min_ess_bulk))

  invisible(x)
}


#' @describeIn dcvar_hmm_fit-methods Extract posterior means of model
#'   coefficients including state-specific rho values.
#' @return A named list with elements `mu`, `Phi`, `sigma_eps`, `z_rho`, and
#'   `rho_state`.
#' @export
coef.dcvar_hmm_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  result <- list(
    mu = .extract_required_coef(summ, "^mu\\[", "mu", "coef.dcvar_hmm_fit()"),
    Phi = .extract_required_coef(summ, "^Phi\\[", "Phi", "coef.dcvar_hmm_fit()")
  )
  # Margin-specific scale params before HMM-specific params
  margins <- object$margins %||% "normal"
  result <- c(result, .extract_margin_coefs(summ, margins))
  result$z_rho <- .extract_required_coef(summ, "^z_rho\\[", "z_rho", "coef.dcvar_hmm_fit()")
  result$rho_state <- .extract_required_coef(summ, "^rho_state\\[", "rho_state", "coef.dcvar_hmm_fit()")
  result
}


#' @describeIn dcvar_hmm_fit-methods Dispatch to a plot type: `"rho"`,
#'   `"states"`, `"transition"`, `"phi"`, `"diagnostics"`, `"ppc"`, or `"pit"`.
#' @param type Character; one of `"rho"`, `"states"`, `"transition"`, `"phi"`,
#'   `"diagnostics"`, `"ppc"`, or `"pit"`.
#' @return A ggplot object.
#' @export
plot.dcvar_hmm_fit <- function(x,
                               type = c("rho", "states", "transition", "phi", "diagnostics", "ppc", "pit"),
                               ...) {
  type <- match.arg(type)
  switch(type,
    rho = plot_rho(x, ...),
    states = plot_hmm_states(x, ...),
    transition = {
      states <- hmm_states(x)
      .plot_transition_matrix(states$A, x$K)
    },
    phi = plot_phi(x, ...),
    diagnostics = plot_diagnostics(x, ...),
    ppc = plot_ppc(x, ...),
    pit = plot_pit(x, ...)
  )
}
