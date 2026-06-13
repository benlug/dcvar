# ============================================================================
# LOO-CV and Model Comparison
# ============================================================================

#' Compute LOO-CV for a fitted model
#'
#' @param x A fitted model object.
#' @param ... Additional arguments passed to [loo::loo()].
#'
#' @return A `loo` object from the loo package.
#'
#' @details PSIS-LOO is available for Gaussian and Clayton single-level fits,
#'   covariate fits, multilevel fits, and naive SEM score fits. Indicator SEM
#'   fits are not supported because their stored `log_lik` conditions on
#'   latent innovations rather than being an observation-level predictive
#'   density. Multilevel `log_lik` values are conditional one-step-ahead
#'   densities given the unit-level random effects.
#'
#'   Note that the pointwise `log_lik` estimands differ across model families:
#'   HMM fits store the state-marginalized one-step-ahead predictive density,
#'   while dynamic (DC-VAR and covariate) fits condition on the smoothed
#'   latent `rho` trajectory, which is itself informed by the observation.
#'   Comparisons across these families can systematically favor the
#'   latent-conditioned models; see [dcvar_compare()].
#' @importFrom loo loo
#' @name loo.dcvar
NULL

.loo_dcvar <- function(x, ...) {
  log_lik <- .fit_draws(
    x$fit, "log_lik",
    format = "draws_array",
    backend = x$backend,
    required = .stan_output_group_pattern("log_lik"),
    required_type = "pattern",
    context = ".loo_dcvar()",
    output_type = "generated quantity"
  )
  r_eff <- loo::relative_eff(exp(log_lik))
  loo::loo(log_lik, r_eff = r_eff, ...)
}

.abort_unsupported_loo <- function(x, reason) {
  model_label <- switch(class(x)[[1]],
    dcvar_sem_fit = "SEM",
    x$model %||% class(x)[[1]]
  )

  cli_abort(c(
    "{.fun loo} is not supported for {.val {model_label}} fits.",
    "i" = reason
  ))
}

.is_supported_loo_fit <- function(x) {
  if (inherits(x, "dcvar_sem_fit")) {
    return(identical(x$method %||% "indicator", "naive"))
  }
  TRUE
}

#' @rdname loo.dcvar
#' @method loo dcvar_fit
#' @export
loo.dcvar_fit <- function(x, ...) .loo_dcvar(x, ...)

#' @rdname loo.dcvar
#' @method loo dcvar_covariate_fit
#' @export
loo.dcvar_covariate_fit <- function(x, ...) .loo_dcvar(x, ...)

#' @rdname loo.dcvar
#' @method loo dcvar_hmm_fit
#' @export
loo.dcvar_hmm_fit <- function(x, ...) .loo_dcvar(x, ...)

#' @rdname loo.dcvar
#' @method loo dcvar_constant_fit
#' @export
loo.dcvar_constant_fit <- function(x, ...) .loo_dcvar(x, ...)

#' @rdname loo.dcvar
#' @method loo dcvar_multilevel_fit
#' @export
loo.dcvar_multilevel_fit <- function(x, ...) .loo_dcvar(x, ...)

#' @rdname loo.dcvar
#' @method loo dcvar_sem_fit
#' @export
loo.dcvar_sem_fit <- function(x, ...) {
  if (.is_supported_loo_fit(x)) {
    return(.loo_dcvar(x, ...))
  }
  .abort_unsupported_loo(
    x,
    paste(
      "The stored {.code log_lik} conditions on latent innovations rather than using",
      "an observation-level predictive density integrated over the latent states."
    )
  )
}


#' Compare multiple fitted models using LOO-CV
#'
#' Convenience wrapper around [loo::loo_compare()] that accepts named
#' dcvar model fits.
#'
#' @param ... Named fitted model objects (e.g., `dcvar = fit1, hmm = fit2`).
#'
#' @return A `loo_compare` matrix.
#'
#' @details The pointwise `log_lik` estimands are not identical across model
#'   families: HMM fits marginalize the discrete states out of each one-step
#'   predictive density, while dynamic (DC-VAR and covariate) fits condition
#'   on the smoothed latent `rho` trajectory, which is informed by the
#'   observation itself. Comparing an HMM fit against a dynamic fit can
#'   therefore systematically favor the dynamic model; a warning is emitted
#'   for such comparisons and the resulting elpd differences should be
#'   interpreted with caution.
#'
#' @seealso [loo::loo_compare()] for details on the comparison method,
#'   [dcvar()], [dcvar_hmm()], [dcvar_constant()] for fitting models.
#' @export
#'
#' @examples
#' \donttest{
#' sim <- simulate_dcvar(
#'   n_time = 12,
#'   rho_trajectory = rho_decreasing(12),
#'   seed = 1
#' )
#' fit_dcvar <- dcvar(
#'   sim$Y_df,
#'   vars = c("y1", "y2"),
#'   chains = 1,
#'   iter_warmup = 10,
#'   iter_sampling = 10,
#'   refresh = 0,
#'   seed = 1
#' )
#' fit_constant <- dcvar_constant(
#'   sim$Y_df,
#'   vars = c("y1", "y2"),
#'   chains = 1,
#'   iter_warmup = 10,
#'   iter_sampling = 10,
#'   refresh = 0,
#'   seed = 1
#' )
#' dcvar_compare(dcvar = fit_dcvar, constant = fit_constant)
#' }
dcvar_compare <- function(...) {
  fits <- list(...)
  if (is.null(names(fits)) || any(names(fits) == "")) {
    unnamed <- which(names(fits) == "" | is.na(names(fits)))
    cli_abort(c(
      "All arguments must be named.",
      "i" = "Unnamed argument{?s} at position{?s}: {.val {unnamed}}.",
      "i" = "Example: {.code dcvar_compare(dcvar = fit1, hmm = fit2)}"
    ))
  }

  for (nm in names(fits)) {
    if (!inherits(fits[[nm]], "dcvar_model_fit")) {
      cli_abort("Argument {.val {nm}} is not a dcvar model fit object.")
    }
  }

  unsupported <- names(fits)[!vapply(fits, .is_supported_loo_fit, logical(1))]
  if (length(unsupported) > 0) {
    cli_abort(c(
      "{.fun dcvar_compare} does not support one or more supplied fits.",
      "i" = "Unsupported argument{?s}: {.val {unsupported}}.",
      "i" = paste(
        "Indicator SEM fits store {.code log_lik} conditional on latent innovations,",
        "which is not a valid pointwise predictive density for PSIS-LOO."
      )
    ))
  }

  has_hmm <- vapply(fits, inherits, logical(1), what = "dcvar_hmm_fit")
  has_latent_conditioned <- vapply(
    fits,
    function(f) inherits(f, "dcvar_fit") || inherits(f, "dcvar_covariate_fit"),
    logical(1)
  )
  if (any(has_hmm) && any(has_latent_conditioned)) {
    cli_warn(c(
      "Comparing HMM and dynamic (DC-VAR/covariate) fits mixes two {.code log_lik} estimands.",
      "i" = paste(
        "HMM fits store the state-marginalized one-step predictive density, while dynamic",
        "fits condition on the smoothed latent rho trajectory (informed by the observation itself)."
      ),
      "!" = "elpd differences can systematically favor the dynamic model; interpret with caution."
    ))
  }

  has_tv <- vapply(fits, inherits, logical(1), what = "dcvar_tv_fit")
  if (any(has_tv) && any(!has_tv)) {
    cli_warn(c(
      "Time-varying DC-VAR fits condition their {.code log_lik} on additional smoothed latent paths (Phi(t) and/or sigma(t)).",
      "!" = "elpd differences can systematically favor the more heavily latent-conditioned model; interpret comparisons across the model ladder with caution."
    ))
  }

  loo_list <- lapply(fits, loo)
  loo::loo_compare(loo_list)
}
