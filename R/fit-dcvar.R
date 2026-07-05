#' Fit the DC-VAR model
#'
#' Fits a Dynamic Copula VAR(1) model with time-varying correlation following
#' a random walk on the Fisher-z scale. Uses non-centered parameterisation
#' for efficient HMC sampling.
#'
#' Beyond the always time-varying correlation, the VAR coefficients and the
#' residual scales can optionally vary over time as well (`tv_phi`,
#' `tv_sigma`): each enabled component evolves as a non-centered random walk
#' around its constant baseline, with innovation-SD priors that shrink toward
#' the constant-parameter model. With both flags off (the default) the fit is
#' exactly the classic DC-VAR. With any flag on, the fit uses a generic
#' time-varying Stan model and returns a `dcvar_tv_fit` object (a subclass of
#' `dcvar_fit`), with trajectory extractors [phi_trajectory()] and
#' [sigma_trajectory()].
#'
#' @section Time-varying scales and shifted margins:
#' `tv_sigma = TRUE` gives normal and skew-normal dimensions a multiplicative
#' log-scale random walk. Exponential and gamma dimensions use a **soft-barrier**
#' construction: the shifted variate `x = scale + skew * eps` (which has a hard
#' support boundary at 0) is replaced by `x = softplus_k(m_t + skew * eps)`
#' with a time-varying scale `m_t`. This keeps `x` positive, matches the exact
#' shifted margin in the interior, and rounds the boundary smoothly so the
#' scale can vary freely. `tv_sigma_k` controls the sharpness; larger values
#' track the exact exponential/gamma more closely at the cost of a stiffer
#' posterior geometry, smaller values are numerically gentler but introduce a
#' small residual-mean bias in the lower tail.
#'
#' @section Recommended workflow:
#' For typical series lengths (T of 100--300), fit the model ladder
#' `dcvar()` (constant Phi/sigma) -> `dcvar(tv_phi = TRUE)` ->
#' `dcvar(tv_phi = TRUE, tv_sigma = TRUE)` and compare with
#' [dcvar_compare()]. Up to seven latent random walks are estimated from a
#' single bivariate series; the tight innovation-SD priors are what keeps
#' this identifiable, so widen them deliberately.
#'
#' @param data A data frame with time series observations.
#' @param vars Character vector of two variable names to model.
#' @param time_var Name of the time column (default: `"time"`).
#' @param standardize Logical; whether to z-score variables (default: `TRUE`).
#' @param margins Marginal distribution specification. Either a single string
#'   applied to both variables, or a length-2 character vector giving a
#'   per-variable (mixed) margin, e.g. `c("normal", "exponential")`. Each entry
#'   is one of `"normal"` (default), `"exponential"`, `"skew_normal"`, or
#'   `"gamma"`. When the two entries differ the fit uses a generic
#'   mixed-margins Stan model under the Gaussian copula; identical entries reuse
#'   the specialised single-family model.
#' @param skew_direction Integer vector of length 2 indicating skew direction
#'   for asymmetric margins. Each element must be `1` (right-skewed) or `-1`
#'   (left-skewed). Required whenever any dimension uses an `"exponential"` or
#'   `"gamma"` margin; only those dimensions consult it.
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
#' @param tv_phi Selects which VAR(1) coefficients evolve as independent
#'   random walks around the constant baseline `Phi`. Either a logical
#'   (`TRUE` = all four, `FALSE` = none, the default) or a character selector:
#'   `"ar"` (the autoregressive effects phi11, phi22 -- e.g. for modelling
#'   changing emotional inertia / critical slowing down), `"cross"` (the
#'   cross-lagged effects phi12, phi21), or specific names from
#'   `c("phi11", "phi12", "phi21", "phi22")`.
#' @param tv_sigma Logical; if `TRUE`, the residual scales evolve as log-scale
#'   random walks around their constant baselines (default: `FALSE`). Applies
#'   to all margin families; see the section on shifted margins for how
#'   exponential and gamma dimensions are handled.
#' @param prior_tau_phi_rate Prior mean for the Phi random-walk innovation
#'   SDs (`tau_phi ~ exponential(1/prior_tau_phi_rate)`; default `0.025`).
#'   Used only when `tv_phi = TRUE`.
#' @param prior_tau_sigma_rate Prior mean for the log-scale random-walk
#'   innovation SDs (`tau_sigma ~ exponential(1/prior_tau_sigma_rate)`;
#'   default `0.05`). Used only when `tv_sigma = TRUE`.
#' @param tv_sigma_k Soft-barrier sharpness for time-varying exponential/gamma
#'   scales (default `8`). Larger values approximate the exact shifted margin
#'   more closely but stiffen the geometry. Used only when `tv_sigma = TRUE`
#'   and an exponential or gamma margin is present.
#' @param chains Number of MCMC chains (default: 4).
#' @param iter_warmup Warmup iterations per chain (default: 2000).
#' @param iter_sampling Sampling iterations per chain (default: 4000).
#' @param adapt_delta Target acceptance rate (default: 0.99). The DC-VAR model
#'   uses a lower default than `dcvar_constant()` (0.999) because the
#'   non-centered parameterisation already handles posterior geometry well.
#' @param max_treedepth Maximum tree depth (default: 12).
#' @param seed Random seed. When supplied, it seeds both the Stan sampler and
#'   the default per-chain initial values, so repeated fits are reproducible.
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
#' \donttest{
#' sim <- simulate_dcvar(
#'   n_time = 12,
#'   rho_trajectory = rho_decreasing(12),
#'   seed = 1
#' )
#' fit <- dcvar(
#'   sim$Y_df,
#'   vars = c("y1", "y2"),
#'   chains = 1,
#'   iter_warmup = 10,
#'   iter_sampling = 10,
#'   refresh = 0,
#'   seed = 1
#' )
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
                  tv_phi = FALSE,
                  tv_sigma = FALSE,
                  prior_tau_phi_rate = 0.025,
                  prior_tau_sigma_rate = 0.05,
                  tv_sigma_k = 8,
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
  margins <- .normalize_margins_spec(margins)
  .validate_margins(margins, skew_direction)
  phi_mask <- .resolve_phi_tv_mask(tv_phi)
  any_phi <- sum(phi_mask) > 0L
  .prep_validate_scalar_logical(tv_sigma, "tv_sigma")

  is_tv <- any_phi || tv_sigma
  model_type <- if (is_tv) "dcvar_tv" else "dcvar"
  margins_vec <- rep(margins, length.out = 2L)

  if (!any_phi && !isTRUE(all.equal(prior_tau_phi_rate, 0.025))) {
    cli_warn("{.arg prior_tau_phi_rate} is ignored when no VAR coefficient is time-varying.")
  }
  if (!tv_sigma && !isTRUE(all.equal(prior_tau_sigma_rate, 0.05))) {
    cli_warn("{.arg prior_tau_sigma_rate} is ignored when {.code tv_sigma = FALSE}.")
  }
  if (tv_sigma) {
    .prep_validate_positive_scalar(tv_sigma_k, "tv_sigma_k")
  }

  # Prepare data
  stan_data <- prepare_dcvar_data(
    data, vars, time_var, standardize, margins, skew_direction,
    prior_mu_sd, prior_phi_sd, prior_sigma_eps_rate,
    prior_sigma_omega_rate, prior_rho_init_sd,
    tv_phi = tv_phi, tv_sigma = tv_sigma,
    prior_tau_phi_rate = prior_tau_phi_rate,
    prior_tau_sigma_rate = prior_tau_sigma_rate,
    tv_sigma_k = tv_sigma_k,
    allow_gaps = allow_gaps
  )

  margins_label <- if (all(margins == "normal")) {
    ""
  } else {
    paste0(" [", paste(margins, collapse = ", "), "]")
  }
  if (is_tv) {
    phi_label <- if (any_phi) {
      if (all(phi_mask == 1L)) "Phi(t)" else paste0("Phi(t):", paste(names(phi_mask)[phi_mask == 1L], collapse = "/"))
    }
    components <- c("rho(t)", phi_label, if (tv_sigma) "sigma(t)")
    cli_inform("Fitting TV DC-VAR model{margins_label} [{paste(components, collapse = ', ')}] (n_time = {stan_data$n_time}, D = {stan_data$D})...")
  } else {
    cli_inform("Fitting DC-VAR model{margins_label} (n_time = {stan_data$n_time}, D = {stan_data$D})...")
  }

  uses_engine <- is_tv && .uses_dynamic_stan_file(stan_file, margins = margins)
  compile_stan_file <- if (uses_engine) {
    dcvar_stan_path("dcvar_dynamic", margins = margins)
  } else {
    stan_file
  }

  # Compile model
  model <- .compile_model(model_type, margins = margins, stan_file = compile_stan_file,
                          backend = backend)

  # The time-varying path is served by the unified dynamic engine
  # (dcvar_dynamic.stan), whose data block is a superset of the TV output of
  # prepare_dcvar_data(). Augment the prepared data with the engine's covariate
  # predictor / drift / fast-path fields just before sampling when the bundled
  # engine is selected. A user-supplied legacy TV Stan file keeps the legacy
  # TV block.
  if (uses_engine) {
    stan_data <- .as_dynamic_stan_data(stan_data)
  }

  # Default init
  if (is.null(init)) {
    D <- stan_data$D
    n_time_obs <- stan_data$n_time
    init <- if (is_tv) {
      # The TV path runs on the unified dynamic engine, so use its init (P = 0
      # omits beta, zero_init_eta = FALSE matches the TV omega_raw length). This
      # is also valid for a custom legacy dcvar_tv_mixed.stan file, whose
      # parameter block is the same minus the (omitted) covariate effects.
      function() .init_dcvar_dynamic_params(D, n_time_obs, margins,
                                            tv_phi = phi_mask, tv_sigma = tv_sigma,
                                            P = 0L, zero_init_eta = FALSE)
    } else {
      function() .init_dcvar_params(D, n_time_obs, margins)
    }
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

  priors <- c(
    list(
      mu_sd = prior_mu_sd,
      phi_sd = prior_phi_sd,
      sigma_eps_rate = prior_sigma_eps_rate,
      sigma_omega_rate = prior_sigma_omega_rate,
      rho_init_sd = prior_rho_init_sd
    ),
    # tau priors only act on enabled components; recording them otherwise
    # would suggest they were used.
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

  # Wrap in S3 class
  out <- if (is_tv) {
    new_dcvar_tv_fit(
      fit = fit,
      stan_data = stan_data,
      vars = vars,
      standardized = standardize,
      margins = margins,
      skew_direction = skew_direction,
      tv_phi = any_phi,
      phi_tv_mask = phi_mask,
      tv_sigma = tv_sigma,
      backend = backend,
      priors = priors,
      meta = meta
    )
  } else {
    new_dcvar_fit(
      fit = fit,
      stan_data = stan_data,
      vars = vars,
      standardized = standardize,
      margins = margins,
      skew_direction = skew_direction,
      backend = backend,
      priors = priors,
      meta = meta
    )
  }

  .report_sampling_outcome(out$fit, if (is_tv) "TV DC-VAR" else "DC-VAR",
                           chains = chains, backend = backend, object = out)
  out
}
