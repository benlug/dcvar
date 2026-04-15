# ============================================================================
# as.data.frame methods for fit objects
# ============================================================================

#' Convert a dcvar fit summary to a data frame
#'
#' Returns the full parameter summary as a tidy data frame with
#' correct 2.5%/97.5% quantiles.
#'
#' @param x A fitted model object (`dcvar_fit`, `dcvar_hmm_fit`, or
#'   `dcvar_constant_fit`).
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `variable`, `mean`, `sd`, `q2.5`,
#'   `q97.5`, `rhat`, `ess_bulk`, `ess_tail`.
#' @export
as.data.frame.dcvar_model_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  summ <- .fit_summary(
    x$fit, variables = NULL, backend = x$backend,
    mean = mean,
    sd = sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975)),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  as.data.frame(summ)
}
