#' Fit an experimental SEM copula VAR(1) model
#'
#' Fits a copula VAR(1) model with a fixed measurement model (factor
#' loadings and measurement error SD are not estimated). Latent innovations
#' are treated as parameters, making this model computationally intensive
#' for large T.
#'
#' @param data A data frame with time series of indicator variables.
#' @param indicators A list of two character vectors, each naming J indicator
#'   columns per latent variable. For example:
#'   `list(PA = c("y1_1", "y1_2", "y1_3"), NA_ = c("y2_1", "y2_2", "y2_3"))`.
#' @param J Number of indicators per latent variable.
#' @param lambda Numeric vector of length J with fixed factor loadings.
#' @param sigma_e Fixed measurement error SD (scalar).
#' @param margins Character string specifying the latent innovation margin.
#'   One of `"normal"` (default) or `"exponential"`.
#' @param skew_direction Integer vector of length 2 indicating skew direction
#'   for exponential margins. Each element must be `1` (right-skewed) or `-1`
#'   (left-skewed). Required when `margins = "exponential"`.
#' @param time_var Name of the time column (default: `"time"`).
#' @param prior_mu_sd Prior SD for intercepts: `mu ~ normal(0, prior_mu_sd)`.
#' @param prior_phi_sd Prior SD for VAR coefficients:
#'   `phi ~ normal(0, prior_phi_sd)`.
#' @param prior_sigma_sd Prior SD for the lognormal prior on the latent
#'   innovation scale parameter. For normal margins this is applied to
#'   `sigma`; for exponential margins it is applied to `sigma_exp`.
#' @param prior_rho_sd Prior SD for rho_raw:
#'   `rho_raw ~ normal(0, prior_rho_sd)`, with `rho = 0.97 * tanh(rho_raw)`.
#' @param chains Number of MCMC chains.
#' @param iter_warmup Warmup iterations per chain.
#' @param iter_sampling Sampling iterations per chain.
#' @param adapt_delta Target acceptance rate.
#' @param max_treedepth Maximum tree depth.
#' @param seed Random seed.
#' @param cores Number of parallel chains.
#' @param refresh How often to print progress.
#' @param init Custom init function or `NULL`.
#' @param stan_file Custom Stan file path or `NULL`.
#' @param backend Character: `"auto"` (default, uses rstan), `"rstan"`, or
#'   `"cmdstanr"`. Can also be set globally via
#'   `options(dcvar.backend = "cmdstanr")`.
#' @param ... Additional backend-specific sampling arguments.
#'
#' @return A `dcvar_sem_fit` object.
#'
#' @details
#' **Experimental extension.** This SEM variant supports `fitted()` and
#' `predict()`, but PIT diagnostics and PSIS-LOO are not yet implemented.
#'
#' **Boundary constraints.** The SEM model constrains each VAR coefficient
#' (Phi) to the interval \eqn{[-0.99, 0.99]}, unlike other dcvar models
#' where Phi is unconstrained. Very strong autoregressive or cross-lag
#' dynamics near \eqn{\pm 1} cannot be captured by this variant.
#'
#' The copula correlation \eqn{\rho} is constrained to \eqn{(-0.97, 0.97)}
#' via `rho = 0.97 * tanh(rho_raw)` to avoid boundary singularity in the
#' Gaussian copula density. Extremely high correlations near \eqn{\pm 1}
#' are truncated.
#'
#' **Margins.** The SEM model currently supports normal and exponential latent
#' innovation margins. Exponential margins use the same shifted-exponential
#' parameterization as the single-level models and therefore require
#' `skew_direction`. Other non-normal margins are not yet available.
#'
#' **Post-estimation.** `fitted()` and `predict()` are available for both the
#' latent-state scale (`type = "link"`) and the observed-indicator scale
#' (`type = "response"`). Use [latent_states()] when you specifically need the
#' full posterior summaries of the latent trajectories.
#'
#' @note This model currently supports normal and exponential latent margins.
#'   For skew-normal or gamma margins, use [dcvar()], [dcvar_constant()], or
#'   [dcvar_hmm()].
#'
#' @seealso [latent_states()] for extracting estimated latent states,
#'   [simulate_dcvar_sem()] for data generation.
#' @export
dcvar_sem <- function(data, indicators, J, lambda, sigma_e,
                      margins = "normal",
                      skew_direction = NULL,
                      time_var = "time",
                      prior_mu_sd = 0.25,
                      prior_phi_sd = 0.5,
                      prior_sigma_sd = 0.5,
                      prior_rho_sd = 0.75,
                      chains = 4,
                      iter_warmup = 2000,
                      iter_sampling = 4000,
                      adapt_delta = 0.95,
                      max_treedepth = 13,
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
  .validate_sem_margins(margins, skew_direction)

  if (!is.numeric(J) || length(J) != 1 || J < 1 || J != as.integer(J)) {
    cli_abort("{.arg J} must be a positive integer, got {.val {J}}.")
  }

  # Prepare SEM data
  stan_data <- prepare_sem_data(
    data, indicators, J, lambda, sigma_e,
    margins = margins,
    skew_direction = skew_direction,
    time_var = time_var,
    prior_mu_sd = prior_mu_sd,
    prior_phi_sd = prior_phi_sd,
    prior_sigma_sd = prior_sigma_sd,
    prior_rho_sd = prior_rho_sd
  )

  T_obs <- stan_data$T
  vars <- attr(stan_data, "vars")

  cli_inform("Fitting SEM copula VAR model [{margins}] (T = {T_obs}, J = {J})...")

  # Compile model
  model <- .compile_model("sem", margins = margins, stan_file = stan_file, backend = backend)

  # Default init
  if (is.null(init)) {
    init <- function() .init_sem_params(T_obs, margins)
  }

  cores <- .normalize_cores(cores, chains)

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

  .report_sampling_outcome(fit, "SEM copula VAR", chains = chains, backend = backend)

  new_dcvar_sem_fit(
    fit = fit,
    stan_data = stan_data,
    vars = vars,
    J = J,
    lambda = lambda,
    sigma_e = sigma_e,
    indicators = attr(stan_data, "indicators"),
    margins = margins,
    skew_direction = attr(stan_data, "skew_direction"),
    backend = backend,
    priors = list(
      mu_sd = prior_mu_sd,
      phi_sd = prior_phi_sd,
      sigma_sd = prior_sigma_sd,
      rho_sd = prior_rho_sd
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
