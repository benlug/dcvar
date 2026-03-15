# ============================================================================
# LOO-CV and Model Comparison
# ============================================================================

#' Compute LOO-CV for a fitted model
#'
#' @param x A fitted model object.
#' @param ... Additional arguments passed to [loo::loo()].
#'
#' @return A `loo` object from the loo package.
#' @importFrom loo loo
#' @name loo.dcVar
NULL

.loo_dcVar <- function(x, ...) {
  log_lik <- x$fit$draws("log_lik", format = "draws_array")
  r_eff <- loo::relative_eff(exp(log_lik))
  loo::loo(log_lik, r_eff = r_eff, ...)
}

.abort_unsupported_loo <- function(x, reason) {
  model_label <- switch(class(x)[[1]],
    dcVar_multilevel_fit = "multilevel",
    dcVar_sem_fit = "SEM",
    x$model %||% class(x)[[1]]
  )

  cli_abort(c(
    "{.fun loo} is not supported for {.val {model_label}} fits.",
    "i" = reason
  ))
}

.is_supported_loo_fit <- function(x) {
  !inherits(x, c("dcVar_multilevel_fit", "dcVar_sem_fit"))
}

#' @rdname loo.dcVar
#' @method loo dcVar_fit
#' @export
loo.dcVar_fit <- function(x, ...) .loo_dcVar(x, ...)

#' @rdname loo.dcVar
#' @method loo dcVar_hmm_fit
#' @export
loo.dcVar_hmm_fit <- function(x, ...) .loo_dcVar(x, ...)

#' @rdname loo.dcVar
#' @method loo dcVar_constant_fit
#' @export
loo.dcVar_constant_fit <- function(x, ...) .loo_dcVar(x, ...)

#' @rdname loo.dcVar
#' @method loo dcVar_multilevel_fit
#' @export
loo.dcVar_multilevel_fit <- function(x, ...) {
  .abort_unsupported_loo(
    x,
    paste(
      "The stored {.code log_lik} is per-unit and conditional on unit-level random effects,",
      "so it is not a valid pointwise predictive density for PSIS-LOO."
    )
  )
}

#' @rdname loo.dcVar
#' @method loo dcVar_sem_fit
#' @export
loo.dcVar_sem_fit <- function(x, ...) {
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
#' dcVar model fits.
#'
#' @param ... Named fitted model objects (e.g., `dcVar = fit1, hmm = fit2`).
#'
#' @return A `loo_compare` matrix.
#'
#' @seealso [loo::loo_compare()] for details on the comparison method,
#'   [dcVar()], [dcVar_hmm()], [dcVar_constant()] for fitting models.
#' @export
#'
#' @examples
#' \dontrun{
#' dcVar_compare(dcVar = fit1, hmm = fit2, constant = fit3)
#' }
dcVar_compare <- function(...) {
  fits <- list(...)
  if (is.null(names(fits)) || any(names(fits) == "")) {
    unnamed <- which(names(fits) == "" | is.na(names(fits)))
    cli_abort(c(
      "All arguments must be named.",
      "i" = "Unnamed argument{?s} at position{?s}: {.val {unnamed}}.",
      "i" = "Example: {.code dcVar_compare(dcVar = fit1, hmm = fit2)}"
    ))
  }

  for (nm in names(fits)) {
    if (!inherits(fits[[nm]], "dcVar_model_fit")) {
      cli_abort("Argument {.val {nm}} is not a dcVar model fit object.")
    }
  }

  unsupported <- names(fits)[!vapply(fits, .is_supported_loo_fit, logical(1))]
  if (length(unsupported) > 0) {
    cli_abort(c(
      "{.fun dcVar_compare} does not support SEM or multilevel fits.",
      "i" = "Unsupported argument{?s}: {.val {unsupported}}.",
      "i" = paste(
        "These models store {.code log_lik} targets that are not valid, comparable",
        "pointwise predictive densities for PSIS-LOO."
      )
    ))
  }

  loo_list <- lapply(fits, loo)
  loo::loo_compare(loo_list)
}
