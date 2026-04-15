#' Fit the DC-VAR model
#'
#' Fits a Dynamic Copula VAR(1) model with time-varying correlation following
#' a random walk on the Fisher-z scale. Uses non-centered parameterisation
#' for efficient HMC sampling.
#'
#' @param data A data frame with time series observations.
#' @param vars Character vector of two variable names to model.
#' @param time_var Name of the time column (default: `"time"`).
#' @param standardize Logical; whether to z-score variables (default: `TRUE`).
#' @param margins Character string specifying the marginal distribution.
#'   One of `"normal"` (default), `"exponential"`, `"skew_normal"`, or `"gamma"`.
#' @param skew_direction Integer vector of length D indicating skew direction
#'   for asymmetric margins. Each element must be `1` (right-skewed) or `-1`
#'   (left-skewed). Required for `"exponential"` and `"gamma"` margins.
#' @param allow_gaps Logical; if `FALSE` (default), interior missing values
#'   cause an error because they break VAR(1) time series adjacency. Set to
#'   `TRUE` to allow fitting with a warning instead.
#' @param prior_mu_sd Prior SD for intercepts: `mu ~ normal(0, prior_mu_sd)`.
#' @param prior_phi_sd Prior SD for VAR coefficients: `Phi ~ normal(0, prior_phi_sd)`.
#' @param prior_sigma_eps_rate Prior mean for innovation SDs (see
#'   [prepare_dcvar_data()]).
#' @param prior_sigma_omega_rate Prior mean for rho process SD (see
#'   [prepare_dcvar_data()]).
#' @param prior_rho_init_sd Prior SD for initial rho on Fisher-z scale.
#' @param chains Number of MCMC chains (default: 4).
#' @param iter_warmup Warmup iterations per chain (default: 2000).
#' @param iter_sampling Sampling iterations per chain (default: 4000).
#' @param adapt_delta Target acceptance rate (default: 0.99). The DC-VAR model
#'   uses a lower default than `dcvar_constant()` (0.999) because the
#'   non-centered parameterisation already handles posterior geometry well.
#' @param max_treedepth Maximum tree depth (default: 12).
#' @param seed Random seed.
#' @param cores Number of parallel chains. `NULL` uses all available cores.
#' @param refresh How often to print progress (default: 500). Set to 0 for
#'   silent operation.
#' @param init Custom init function or `NULL` for smart defaults.
#' @param stan_file Path to a custom Stan file, or `NULL` to use the bundled
#'   model.
#' @param backend Character: `"auto"` (default, uses rstan), `"rstan"`, or
#'   `"cmdstanr"`. Can also be set globally via
#'   `options(dcvar.backend = "cmdstanr")`.
#' @param ... Additional backend-specific sampling arguments.
#'
#' @return A `dcvar_fit` object.
#'
#' @seealso [dcvar_constant()] for the time-invariant baseline,
#'   [dcvar_hmm()] for the regime-switching model,
#'   [dcvar_compare()] for LOO-CV model comparison,
#'   [rho_trajectory()] and [plot_rho()] for inspecting results.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- simulate_dcvar(T = 100, rho_trajectory = rho_decreasing(100))
#' fit <- dcvar(sim$Y_df, vars = c("y1", "y2"))
#' print(fit)
#' summary(fit)
#' plot(fit)
#' }
dcvar <- function(data, vars, time_var = "time",
                  standardize = TRUE,
                  margins = "normal",
                  skew_direction = NULL,
                  allow_gaps = FALSE,
                  prior_mu_sd = 2,
                  prior_phi_sd = 0.5,
                  prior_sigma_eps_rate = 1,
                  prior_sigma_omega_rate = 0.1,
                  prior_rho_init_sd = 1,
                  chains = 4,
                  iter_warmup = 2000,
                  iter_sampling = 4000,
                  adapt_delta = 0.99,
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

  # Prepare data
  stan_data <- prepare_dcvar_data(
    data, vars, time_var, standardize, margins, skew_direction,
    prior_mu_sd, prior_phi_sd, prior_sigma_eps_rate,
    prior_sigma_omega_rate, prior_rho_init_sd,
    allow_gaps = allow_gaps
  )

  margins_label <- if (margins == "normal") "" else paste0(" [", margins, "]")
  cli_inform("Fitting DC-VAR model{margins_label} (T = {stan_data$T}, D = {stan_data$D})...")

  # Compile model
  model <- .compile_model("dcvar", margins = margins, stan_file = stan_file,
                          backend = backend)

  # Default init
  if (is.null(init)) {
    D <- stan_data$D
    T_obs <- stan_data$T
    init <- function() .init_dcvar_params(D, T_obs, margins)
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

  .report_sampling_outcome(fit, "DC-VAR", chains = chains, backend = backend)

  # Wrap in S3 class
  new_dcvar_fit(
    fit = fit,
    stan_data = stan_data,
    vars = vars,
    standardized = standardize,
    margins = margins,
    skew_direction = skew_direction,
    backend = backend,
    priors = list(
      mu_sd = prior_mu_sd,
      phi_sd = prior_phi_sd,
      sigma_eps_rate = prior_sigma_eps_rate,
      sigma_omega_rate = prior_sigma_omega_rate,
      rho_init_sd = prior_rho_init_sd
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
