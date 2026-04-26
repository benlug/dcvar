#' dcvar: Gaussian-copula VAR models for bivariate time series
#'
#' `dcvar` focuses on copula VAR(1) workflows for bivariate time
#' series. The core single-level models [dcvar()], [dcvar_hmm()], and
#' [dcvar_constant()] support normal, exponential, skew-normal, and gamma
#' margins. [dcvar_constant()] additionally supports a Clayton copula baseline
#' with normal margins.
#'
#' Experimental extensions via [dcvar_multilevel()] and [dcvar_sem()] provide a
#' narrower diagnostic interface than the core models. The multilevel model
#' supports normal and exponential margins, while the SEM model supports normal
#' and exponential latent innovation margins plus a naive row-mean score method.
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
#' @importFrom rlang .data %||%
#' @importFrom cli cli_abort cli_warn cli_inform cli_alert_info cli_alert_success
## usethis namespace: end
NULL
