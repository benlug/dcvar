#' dcvar: Gaussian-copula VAR models for bivariate time series
#'
#' `dcvar` focuses on Gaussian-copula VAR(1) workflows for bivariate time
#' series. The core single-level models [dcvar()], [dcvar_hmm()], and
#' [dcvar_constant()] support normal, exponential, skew-normal, and gamma
#' margins.
#'
#' Experimental extensions via [dcvar_multilevel()] and [dcvar_sem()] provide a
#' narrower diagnostic interface than the core models. The multilevel model
#' currently supports normal margins only, while the SEM model supports normal
#' and exponential latent innovation margins.
#'
#' The package implements Gaussian-copula models only. Clayton copulas are not
#' part of the public API.
#'
#' Estimation uses Stan through `rstan` by default, with optional `cmdstanr`
#' backend support. PIT diagnostics are approximate posterior-mean plug-in
#' diagnostics. PSIS-LOO model comparison is available for the supported
#' single-level models.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats coef complete.cases cor quantile qnorm rnorm runif sd setNames
#' @importFrom utils head
#' @importFrom rlang .data
#' @importFrom cli cli_abort cli_warn cli_inform cli_alert_info cli_alert_success
## usethis namespace: end
NULL
