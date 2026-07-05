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
#' @param J Number of indicators per latent variable. If `method = "naive"`
#'   and `J` is `NULL`, it is inferred from `indicators`.
#' @param lambda Numeric vector of length J with fixed factor loadings.
#'   Required when `method = "indicator"`.
#' @param sigma_e Fixed measurement error SD (scalar). Required when
#'   `method = "indicator"`.
#' @param method Character string: `"indicator"` (default) fits the fixed
#'   measurement model; `"naive"` averages indicators within each latent block
#'   and fits the observed score VAR.
#' @param margins Latent-innovation marginal specification. A single string
#'   applies the same family to both latent variables. A length-2 character
#'   vector gives a per-variable (mixed) margin (for example
#'   `c("normal", "gamma")`). Normal, exponential, skew-normal, and gamma
#'   margins are supported; homogeneous skew-normal/gamma and per-variable specs
#'   route through the generic mixed-margin Stan models under the Gaussian
#'   copula. Applies to both the indicator and naive methods.
#' @param skew_direction Integer vector of length 2 of `1` (right-skewed) or
#'   `-1` (left-skewed). Required whenever any dimension uses an
#'   `"exponential"` or `"gamma"` margin; only those dimensions consult it.
#' @param time_var Name of the time column (default: `"time"`).
#' @param prior_mu_sd Prior SD for intercepts: `mu ~ normal(0, prior_mu_sd)`.
#' @param prior_phi_sd Prior SD for VAR coefficients:
#'   `phi ~ normal(0, prior_phi_sd)`.
#' @param prior_sigma_sd Prior SD for the lognormal prior on the latent
#'   innovation scale parameter. For normal margins this is applied to
#'   `sigma`; for exponential margins it is applied to `sigma_exp`.
#' @param prior_rho_sd Prior SD for rho_raw:
#'   `rho_raw ~ normal(0, prior_rho_sd)`, with `rho = 0.97 * tanh(rho_raw)`.
#'   For TV SEM fits this is the prior SD for the initial Fisher-z rho.
#' @param tv_phi Selects which latent VAR(1) coefficients evolve as random walks
#'   around their baseline. Either a logical scalar or a character selector as in
#'   [dcvar()]. Time-varying SEM currently requires `method = "indicator"`.
#' @param tv_sigma Logical; if `TRUE`, latent innovation scales evolve as
#'   log-scale random walks. Time-varying SEM currently requires
#'   `method = "indicator"`.
#' @param prior_sigma_omega_rate Prior mean for the TV SEM Fisher-z rho
#'   random-walk SD.
#' @param prior_tau_phi_rate Prior mean for Phi random-walk innovation SDs.
#' @param prior_tau_sigma_rate Prior mean for log-scale random-walk innovation
#'   SDs.
#' @param tv_sigma_k Soft-barrier sharpness for time-varying exponential/gamma
#'   scales.
#' @param chains Number of MCMC chains.
#' @param iter_warmup Warmup iterations per chain.
#' @param iter_sampling Sampling iterations per chain.
#' @param adapt_delta Target acceptance rate.
#' @param max_treedepth Maximum tree depth.
#' @param seed Random seed. When supplied, it seeds both the Stan sampler and
#'   the default per-chain initial values, so repeated fits are reproducible.
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
#' `predict()`. PSIS-LOO is available for naive score fits. PIT diagnostics are
#' not yet implemented.
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
#' **Margins.** SEM fits support normal, exponential, skew-normal, and gamma
#' latent innovation margins. Exponential and gamma margins use the same shifted
#' positive-support parameterization as the single-level models and therefore
#' require `skew_direction`. Homogeneous skew-normal/gamma and per-variable
#' margin specs route to the generic mixed-margins Stan model.
#'
#' **Post-estimation.** `fitted()` and `predict()` are available for both the
#' latent-state scale (`type = "link"`) and the observed-indicator scale
#' (`type = "response"`). Use [latent_states()] when you specifically need the
#' full posterior summaries of the latent trajectories.
#'
#' @seealso [latent_states()] for extracting estimated latent states,
#'   [simulate_dcvar_sem()] for data generation.
#' @export
dcvar_sem <- function(data, indicators, J = NULL, lambda = NULL, sigma_e = NULL,
                      margins = "normal",
                      skew_direction = NULL,
                      time_var = "time",
                      method = c("indicator", "naive"),
                      prior_mu_sd = 0.25,
                      prior_phi_sd = 0.5,
                      prior_sigma_sd = 0.5,
                      prior_rho_sd = 0.75,
                      tv_phi = FALSE,
                      tv_sigma = FALSE,
                      prior_sigma_omega_rate = 0.1,
                      prior_tau_phi_rate = 0.025,
                      prior_tau_sigma_rate = 0.05,
                      tv_sigma_k = 8,
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
  margins <- .normalize_margins_spec(margins)
  .validate_sem_margins(margins, skew_direction)
  method <- match.arg(method)
  phi_mask <- .resolve_phi_tv_mask(tv_phi)
  any_phi <- sum(phi_mask) > 0L
  .prep_validate_scalar_logical(tv_sigma, "tv_sigma")
  is_tv <- any_phi || tv_sigma

  if (is_tv && identical(method, "naive")) {
    cli_abort(c(
      "Time-varying SEM parameters are not implemented for {.code method = \"naive\"}.",
      "i" = "Use {.code method = \"indicator\"} with fixed measurement parameters, or fit {.fun dcvar} to score data."
    ))
  }
  if (!any_phi && !isTRUE(all.equal(prior_tau_phi_rate, 0.025))) {
    cli_warn("{.arg prior_tau_phi_rate} is ignored when no VAR coefficient is time-varying.")
  }
  if (!tv_sigma && !isTRUE(all.equal(prior_tau_sigma_rate, 0.05))) {
    cli_warn("{.arg prior_tau_sigma_rate} is ignored when {.code tv_sigma = FALSE}.")
  }
  if (!is_tv && !isTRUE(all.equal(prior_sigma_omega_rate, 0.1))) {
    cli_warn("{.arg prior_sigma_omega_rate} is ignored when SEM time-varying components are disabled.")
  }

  if (identical(method, "indicator") &&
      (!is.numeric(J) || length(J) != 1 || J < 1 || J != as.integer(J))) {
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
    prior_rho_sd = prior_rho_sd,
    tv_phi = phi_mask,
    tv_sigma = tv_sigma,
    prior_sigma_omega_rate = prior_sigma_omega_rate,
    prior_tau_phi_rate = prior_tau_phi_rate,
    prior_tau_sigma_rate = prior_tau_sigma_rate,
    tv_sigma_k = tv_sigma_k,
    method = method
  )

  n_time_obs <- stan_data$n_time
  vars <- attr(stan_data, "vars")
  J_fit <- attr(stan_data, "J") %||% J

  method_label <- if (identical(method, "naive")) "naive SEM" else if (is_tv) "TV SEM" else "SEM"
  margins_text <- paste(margins, collapse = ", ")
  cli_inform("Fitting {method_label} copula VAR model [{margins_text}] (n_time = {n_time_obs}, J = {J_fit})...")

  # Compile model
  model_type <- if (is_tv) {
    "sem_tv"
  } else if (identical(method, "naive")) {
    "sem_naive"
  } else {
    "sem"
  }
  model <- .compile_model(model_type, margins = margins, stan_file = stan_file, backend = backend)

  # Default init
  if (is.null(init)) {
    if (is_tv) {
      init <- function() .init_sem_tv_params(n_time_obs, margins,
                                             tv_phi = phi_mask,
                                             tv_sigma = tv_sigma)
    } else if (identical(method, "naive")) {
      y_scores <- stan_data$y
      init <- function() .init_sem_naive_params(y_scores, margins)
    } else {
      init <- function() .init_sem_params(n_time_obs, margins)
    }
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

  priors <- c(
    list(
      mu_sd = prior_mu_sd,
      phi_sd = prior_phi_sd,
      sigma_sd = prior_sigma_sd,
      rho_sd = prior_rho_sd
    ),
    if (is_tv) list(sigma_omega_rate = prior_sigma_omega_rate),
    if (any_phi) list(tau_phi_rate = prior_tau_phi_rate),
    if (tv_sigma) list(tau_sigma_rate = prior_tau_sigma_rate)
  )
  meta <- list(
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed
  )

  common_args <- list(
    fit = fit,
    stan_data = stan_data,
    vars = vars,
    J = J_fit,
    lambda = lambda,
    sigma_e = sigma_e,
    indicators = attr(stan_data, "indicators"),
    margins = margins,
    method = method,
    skew_direction = attr(stan_data, "skew_direction"),
    backend = backend,
    priors = priors,
    meta = meta
  )
  out <- if (is_tv) {
    do.call(new_dcvar_sem_tv_fit, c(common_args, list(
      tv_phi = any_phi,
      phi_tv_mask = phi_mask,
      tv_sigma = tv_sigma
    )))
  } else {
    do.call(new_dcvar_sem_fit, common_args)
  }

  .report_sampling_outcome(out$fit, if (is_tv) "TV SEM copula VAR" else if (identical(method, "naive")) "Naive SEM copula VAR" else "SEM copula VAR",
                           chains = chains, backend = backend, object = out)
  out
}
