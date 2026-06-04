#' dcvar: Copula VAR models for bivariate time series
#'
#' `dcvar` focuses on copula VAR(1) workflows for bivariate time series. Most
#' bundled models use Gaussian copulas. The core single-level Gaussian-copula
#' models [dcvar()], [dcvar_hmm()], and [dcvar_constant()] support normal,
#' exponential, skew-normal, and gamma margins. [dcvar_constant()] additionally
#' supports a Clayton copula baseline with normal margins or a per-variable
#' (mixed) margin vector.
#'
#' Per-variable (mixed) margins let each variable take its own marginal family.
#' Across all model families a length-2 `margins` vector may combine any of
#' normal, exponential, skew-normal, and gamma per dimension (routing to a
#' generic mixed-margins Stan model); a single-family fit is restricted to the
#' families that family's specialised model supports.
#'
#' Experimental extensions via [dcvar_multilevel()] and [dcvar_sem()] provide a
#' narrower diagnostic interface than the core models. Single-family multilevel
#' and SEM fits support normal and exponential (latent innovation) margins, but
#' a per-variable (mixed) margin vector extends both to all four margin families
#' per dimension; the SEM model also provides a naive row-mean score method.
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
