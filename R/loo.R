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
#'   covariate fits, exponential-margin multilevel fits, and naive SEM score
#'   fits. Indicator SEM fits and normal-margin multilevel fits are not
#'   supported because their stored `log_lik` quantities are not comparable
#'   pointwise predictive densities.
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
    dcvar_multilevel_fit = "multilevel",
    dcvar_sem_fit = "SEM",
    x$model %||% class(x)[[1]]
  )

  cli_abort(c(
    "{.fun loo} is not supported for {.val {model_label}} fits.",
    "i" = reason
  ))
}

.is_supported_loo_fit <- function(x) {
  if (inherits(x, "dcvar_multilevel_fit")) {
    return(identical(x$margins %||% "normal", "exponential"))
  }
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
loo.dcvar_multilevel_fit <- function(x, ...) {
  if (.is_supported_loo_fit(x)) {
    return(.loo_dcvar(x, ...))
  }
  .abort_unsupported_loo(
    x,
    paste(
      "The stored {.code log_lik} is per-unit and conditional on unit-level random effects,",
      "so it is not a valid pointwise predictive density for PSIS-LOO."
    )
  )
}

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
        "Unsupported models store {.code log_lik} targets that are not valid, comparable",
        "pointwise predictive densities for PSIS-LOO."
      )
    ))
  }

  loo_list <- lapply(fits, loo)
  loo::loo_compare(loo_list)
}
