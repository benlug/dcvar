#' dcvar: Gaussian-copula VAR models for bivariate time series
#'
#' The package's core scope is the Gaussian-copula family for bivariate VAR(1)
#' time series. [dcvar()], [dcvar_hmm()], and [dcvar_constant()] support
#' normal, exponential, skew-normal, and gamma margins.
#'
#' Experimental extensions via [dcvar_multilevel()] and [dcvar_sem()] currently
#' support normal margins only and provide a narrower post-estimation interface
#' than the core single-level models.
#'
#' `dcvar` currently implements Gaussian-copula models only. Clayton copulas
#' and the paper-specific exponential-indicator SEM variant are not part of the
#' package API.
#'
#' PIT diagnostics and LOO-CV comparison are available for the supported
#' single-level models. See [pit_values()], [pit_test()], [loo::loo()], and
#' [dcvar_compare()].
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats complete.cases cor quantile qnorm rnorm runif sd setNames
#' @importFrom utils head
#' @importFrom rlang .data
#' @importFrom cli cli_abort cli_warn cli_inform cli_alert_info cli_alert_success
## usethis namespace: end
NULL
