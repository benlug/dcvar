#' Fit the constant copula model
#'
#' Fits a copula VAR(1) model with a single time-invariant correlation
#' parameter. This serves as a baseline for comparison with the DC-VAR
#' and HMM copula models.
#'
#' @inheritParams dcvar
#' @param margins Character string specifying the marginal distribution.
#'   One of `"normal"` (default), `"exponential"`, `"skew_normal"`, or `"gamma"`.
#' @param copula Character string specifying the copula family. One of
#'   `"gaussian"` (default) or `"clayton"`. Clayton is currently available
#'   only with normal margins.
#' @param skew_direction Integer vector of length D indicating skew direction
#'   for asymmetric margins. Each element must be `1` (right-skewed) or `-1`
#'   (left-skewed). Required for `"exponential"` and `"gamma"` margins.
#' @param prior_z_rho_sd Prior SD for rho on Fisher-z scale (default: 1.0).
#' @param adapt_delta Target acceptance rate (default: 0.999). The constant
#'   model uses a higher default than DC-VAR (0.99) because its simpler
#'   posterior geometry benefits from tighter step-size adaptation without
#'   significant cost, reducing occasional divergences near the rho boundary.
#'
#' @details `adapt_delta` defaults to 0.999 because the constant-rho model has
#'   a simpler correlation structure that benefits from tighter step-size
#'   adaptation without significant computational cost, reducing occasional
#'   divergences near the rho boundary.
#'
#' @param backend Character: `"auto"` (default, uses rstan), `"rstan"`, or
#'   `"cmdstanr"`. Can also be set globally via
#'   `options(dcvar.backend = "cmdstanr")`.
#' @param ... Additional backend-specific sampling arguments.
#'
#' @return A `dcvar_constant_fit` object.
#'
#' @seealso [dcvar()] for the time-varying model, [dcvar_hmm()] for the
#'   regime-switching model, [dcvar_compare()] for LOO-CV model comparison.
#' @export
#'
#' @examples
#' \donttest{
#' sim <- simulate_dcvar(
#'   n_time = 12,
#'   rho_trajectory = rho_constant(12, rho = 0.5),
#'   seed = 1
#' )
#' fit <- dcvar_constant(
#'   sim$Y_df,
#'   vars = c("y1", "y2"),
#'   chains = 1,
#'   iter_warmup = 10,
#'   iter_sampling = 10,
#'   refresh = 0,
#'   seed = 1
#' )
#' print(fit)
#' }
dcvar_constant <- function(data, vars,
                           time_var = "time",
                           standardize = TRUE,
                           margins = "normal",
                           copula = "gaussian",
                           skew_direction = NULL,
                           allow_gaps = FALSE,
                           prior_mu_sd = 2,
                           prior_phi_sd = 0.5,
                           prior_sigma_eps_rate = 1,
                           prior_z_rho_sd = 1.0,
                           chains = 4,
                           iter_warmup = 2000,
                           iter_sampling = 4000,
                           adapt_delta = 0.999,
                           max_treedepth = 12,
                           seed = NULL,
                           cores = NULL,
                           refresh = 500,
                           init = NULL,
                           stan_file = NULL,
                           backend = getOption("dcvar.backend", "auto"),
                           ...) {
  backend <- .resolve_backend(backend)
  .validate_sampling_args(chains, iter_warmup, iter_sampling,
                          adapt_delta, max_treedepth)
  .validate_margins(margins, skew_direction)
  .validate_copula(copula)
  if (identical(copula, "clayton") && !identical(margins, "normal")) {
    cli_abort("Clayton copula support in {.fun dcvar_constant} currently requires {.arg margins = 'normal'}.")
  }

  # Prepare data
  stan_data <- prepare_constant_data(
    data, vars, time_var, standardize, margins, skew_direction,
    prior_mu_sd, prior_phi_sd, prior_sigma_eps_rate,
    prior_z_rho_sd,
    allow_gaps = allow_gaps
  )
  if (identical(copula, "clayton")) {
    stan_data$z_rho_prior_sd <- NULL
  }
  attr(stan_data, "copula") <- copula

  margins_label <- if (margins == "normal") "" else paste0(" [", margins, "]")
  copula_label <- if (identical(copula, "gaussian")) "" else paste0(" [", copula, " copula]")
  cli_inform("Fitting constant copula model{margins_label}{copula_label} (n_time = {stan_data$n_time}, D = {stan_data$D})...")

  # Compile model
  model <- .compile_model("constant", margins = margins, copula = copula,
                          stan_file = stan_file, backend = backend)

  # Default init
  if (is.null(init)) {
    D <- stan_data$D
    init <- function() .init_constant_params(D, margins, copula = copula)
  }

  cores <- .normalize_cores(cores, chains)

  # Fit
  fit <- .sample_model(
    compiled_model = model,
    stan_data = stan_data,
    backend = backend,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    cores = cores,
    init = init,
    refresh = refresh,
    ...
  )

  .report_sampling_outcome(fit, "Constant copula", chains = chains, backend = backend)

  # Wrap in S3 class
  new_dcvar_constant_fit(
    fit = fit,
    stan_data = stan_data,
    vars = vars,
    standardized = standardize,
    margins = margins,
    copula = copula,
    skew_direction = skew_direction,
    backend = backend,
    priors = list(
      mu_sd = prior_mu_sd,
      phi_sd = prior_phi_sd,
      sigma_eps_rate = prior_sigma_eps_rate,
      z_rho_sd = prior_z_rho_sd,
      copula = copula
    ),
    meta = list(
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      seed = seed
    )
  )
}
