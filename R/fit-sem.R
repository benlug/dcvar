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
#' @param time_var Name of the time column (default: `"time"`).
#' @param prior_mu_sd Prior SD for intercepts: `mu ~ normal(0, prior_mu_sd)`.
#' @param prior_phi_sd Prior SD for VAR coefficients:
#'   `phi ~ normal(0, prior_phi_sd)`.
#' @param prior_sigma_sd Prior SD for lognormal on innovation SDs:
#'   `sigma ~ lognormal(0, prior_sigma_sd)`.
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
#' @param ... Additional arguments passed to `cmdstanr::CmdStanModel$sample()`.
#'
#' @return A `dcvar_sem_fit` object.
#'
#' @details
#' **Experimental extension.** This SEM variant currently has a narrower
#' post-estimation interface than the core single-level models. `fitted()`,
#' `predict()`, PIT diagnostics, and PSIS-LOO are not yet implemented.
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
#' **Margins.** The SEM model currently supports only normal marginal
#' distributions. Non-normal margins (e.g., exponential, gamma) are not
#' available; use [dcvar()], [dcvar_constant()], or [dcvar_hmm()] instead.
#'
#' Use [latent_states()] to extract estimated latent trajectories.
#'
#' @note This model currently supports normal marginal distributions only.
#'   For non-normal margins, use [dcvar()], [dcvar_constant()], or [dcvar_hmm()].
#'
#' @seealso [latent_states()] for extracting estimated latent states,
#'   [simulate_dcvar_sem()] for data generation.
#' @export
dcvar_sem <- function(data, indicators, J, lambda, sigma_e,
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
                      ...) {
  .check_cmdstanr()
  .validate_sampling_args(chains, iter_warmup, iter_sampling,
                          adapt_delta, max_treedepth)

  if (!is.numeric(J) || length(J) != 1 || J < 1 || J != as.integer(J)) {
    cli_abort("{.arg J} must be a positive integer, got {.val {J}}.")
  }

  # Prepare SEM data
  stan_data <- prepare_sem_data(
    data, indicators, J, lambda, sigma_e, time_var,
    prior_mu_sd, prior_phi_sd, prior_sigma_sd, prior_rho_sd
  )

  T_obs <- stan_data$T
  vars <- attr(stan_data, "vars")

  cli_inform("Fitting SEM copula VAR model (T = {T_obs}, J = {J})...")

  # Compile model
  model <- .compile_model("sem", stan_file = stan_file)

  # Default init
  if (is.null(init)) {
    init <- function() .init_sem_params(T_obs)
  }

  if (is.null(cores)) cores <- parallel::detectCores(logical = FALSE)

  fit <- model$sample(
    data = stan_data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    parallel_chains = cores,
    init = init,
    refresh = refresh,
    ...
  )

  .report_sampling_outcome(fit, "SEM copula VAR", chains = chains)

  new_dcvar_sem_fit(
    fit = fit,
    stan_data = stan_data,
    vars = vars,
    J = J,
    lambda = lambda,
    sigma_e = sigma_e,
    indicators = attr(stan_data, "indicators"),
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
