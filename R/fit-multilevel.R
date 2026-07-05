#' Fit an experimental multilevel copula VAR(1) model
#'
#' Fits a hierarchical copula VAR(1) model with unit-specific VAR coefficients
#' (random effects) and a global copula correlation. Uses non-centered
#' parameterization for the random Phi coefficients.
#'
#' @param data A data frame in long (panel) format with columns for unit ID,
#'   time, and two outcome variables.
#' @param vars Character vector of two variable names to model.
#' @param id_var Name of the unit/person ID column (default: `"id"`).
#' @param time_var Name of the time column (default: `"time"`).
#' @param center Logical; whether to person-mean center the data
#'   (default: `TRUE`). The bundled multilevel Stan model requires
#'   `center = TRUE`; set `center = FALSE` only with a custom
#'   `stan_file` that includes intercept terms.
#' @param margins Marginal distribution specification. A single string applies
#'   the same family to both variables. A length-2 character vector gives a
#'   per-variable (mixed) margin (for example `c("normal", "gamma")`). Normal,
#'   exponential, skew-normal, and gamma margins are supported; homogeneous
#'   skew-normal/gamma and per-variable specs route through the generic
#'   mixed-margin Stan model under the Gaussian copula.
#' @param skew_direction Integer vector of length 2 of `1`/`-1`. Required
#'   whenever any dimension uses an `"exponential"` or `"gamma"` margin; only
#'   those dimensions consult it.
#' @param prior_phi_bar_sd Prior SD for population-mean VAR coefficients.
#' @param prior_tau_phi_scale Prior scale for half-t(3) on tau_phi.
#' @param prior_sigma_sd Prior SD for the innovation scales. For normal
#'   margins this is the SD of a half-normal prior on `sigma`; for
#'   exponential/gamma margins it is the SD of the lognormal prior on the
#'   marginal scale.
#' @param prior_rho_sd Prior SD for normal on rho.
#' @param tv_phi Selects which population VAR(1) coefficients carry a shared
#'   time-varying drift around the unit-specific baselines. Either a logical
#'   scalar or a character selector as in [dcvar()].
#' @param tv_sigma Logical; if `TRUE`, residual scales evolve as shared
#'   log-scale random walks across time.
#' @param prior_tau_phi_rate Prior mean for the time-drift Phi random-walk
#'   innovation SDs.
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
#' @return A `dcvar_multilevel_fit` object.
#'
#' @details **Experimental extension.** This multilevel variant supports
#'   `fitted()` and `predict()`. PSIS-LOO is available for all multilevel
#'   fits; the stored `log_lik` values are conditional one-step-ahead
#'   densities given the unit-level random effects. PIT diagnostics are not yet implemented.
#'
#'   `adapt_delta` defaults to 0.90 and `max_treedepth` to 14 because the
#'   hierarchical structure with random effects benefits from deeper trees but
#'   does not require aggressive step-size adaptation.
#'
#' @note The bundled multilevel Stan program is defined for person-mean
#'   centered data and omits intercept terms. With the bundled model,
#'   `center = FALSE` is therefore not supported.
#'
#' @seealso [random_effects()] for extracting unit-specific coefficients,
#'   [simulate_dcvar_multilevel()] for data generation.
#' @export
dcvar_multilevel <- function(data, vars,
                             id_var = "id",
                             time_var = "time",
                             center = TRUE,
                             margins = "normal",
                             skew_direction = NULL,
                             prior_phi_bar_sd = 0.5,
                             prior_tau_phi_scale = 0.2,
                             prior_sigma_sd = 1,
                             prior_rho_sd = 0.5,
                             tv_phi = FALSE,
                             tv_sigma = FALSE,
                             prior_tau_phi_rate = 0.025,
                             prior_tau_sigma_rate = 0.05,
                             tv_sigma_k = 8,
                             chains = 4,
                             iter_warmup = 2000,
                             iter_sampling = 4000,
                             adapt_delta = 0.90,
                             max_treedepth = 14,
                             seed = NULL,
                             cores = NULL,
                             refresh = 500,
                             init = NULL,
                             stan_file = NULL,
                             backend = getOption("dcvar.backend", "auto"),
                             ...) {
  margins <- .normalize_margins_spec(margins)
  .validate_margins(margins, skew_direction)
  phi_mask <- .resolve_phi_tv_mask(tv_phi)
  any_phi <- sum(phi_mask) > 0L
  .prep_validate_scalar_logical(tv_sigma, "tv_sigma")
  is_tv <- any_phi || tv_sigma

  if (!any_phi && !isTRUE(all.equal(prior_tau_phi_rate, 0.025))) {
    cli_warn("{.arg prior_tau_phi_rate} is ignored when no VAR coefficient is time-varying.")
  }
  if (!tv_sigma && !isTRUE(all.equal(prior_tau_sigma_rate, 0.05))) {
    cli_warn("{.arg prior_tau_sigma_rate} is ignored when {.code tv_sigma = FALSE}.")
  }

  bundled_stan <- dcvar_stan_path(if (is_tv) "multilevel_tv" else "multilevel", margins = margins)
  uses_bundled_stan <- is.null(stan_file) || identical(
    normalizePath(stan_file, winslash = "/", mustWork = TRUE),
    normalizePath(bundled_stan, winslash = "/", mustWork = TRUE)
  )

  if (!isTRUE(center) && uses_bundled_stan) {
    cli_abort(c(
      "{.arg center = FALSE} is not supported by the bundled multilevel model.",
      "i" = "The bundled Stan program assumes person-mean centered data and omits intercept terms.",
      "i" = "Use {.code center = TRUE} or supply a custom {.arg stan_file} that implements intercepts."
    ))
  }

  backend <- .resolve_backend(backend)
  .validate_sampling_args(chains, iter_warmup, iter_sampling,
                          adapt_delta, max_treedepth)

  # Prepare panel data
  stan_data <- prepare_multilevel_data(
    data = data,
    vars = vars,
    id_var = id_var,
    time_var = time_var,
    center = center,
    prior_phi_bar_sd = prior_phi_bar_sd,
    prior_tau_phi_scale = prior_tau_phi_scale,
    prior_sigma_sd = prior_sigma_sd,
    prior_rho_sd = prior_rho_sd,
    tv_phi = phi_mask,
    tv_sigma = tv_sigma,
    prior_tau_phi_rate = prior_tau_phi_rate,
    prior_tau_sigma_rate = prior_tau_sigma_rate,
    tv_sigma_k = tv_sigma_k,
    margins = margins,
    skew_direction = skew_direction
  )

  N <- stan_data$N
  n_time_obs <- stan_data$n_time

  margins_label <- if (all(margins == "normal")) "" else paste0(" [", paste(margins, collapse = ", "), "]")
  model_label <- if (is_tv) "TV multilevel" else "multilevel"
  cli_inform("Fitting {model_label} copula VAR model{margins_label} (N = {N}, n_time = {n_time_obs})...")

  # Compile model
  model <- .compile_model(if (is_tv) "multilevel_tv" else "multilevel",
                          margins = margins, stan_file = stan_file, backend = backend)

  # Default init
  if (is.null(init)) {
    init <- if (is_tv) {
      function() .init_multilevel_tv_params(2, N, n_time_obs, margins = margins,
                                            tv_phi = phi_mask,
                                            tv_sigma = tv_sigma)
    } else {
      function() .init_multilevel_params(2, N, margins = margins)
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
      phi_bar_sd = prior_phi_bar_sd,
      tau_phi_scale = prior_tau_phi_scale,
      sigma_sd = prior_sigma_sd,
      rho_sd = prior_rho_sd
    ),
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
    N = N,
    vars = vars,
    centered = center,
    person_means = attr(stan_data, "person_means"),
    margins = margins,
    skew_direction = attr(stan_data, "skew_direction"),
    backend = backend,
    priors = priors,
    meta = meta
  )
  out <- if (is_tv) {
    do.call(new_dcvar_multilevel_tv_fit, c(common_args, list(
      tv_phi = any_phi,
      phi_tv_mask = phi_mask,
      tv_sigma = tv_sigma
    )))
  } else {
    do.call(new_dcvar_multilevel_fit, common_args)
  }

  .report_sampling_outcome(out$fit, if (is_tv) "TV multilevel copula VAR" else "Multilevel copula VAR",
                           chains = chains, backend = backend, object = out)
  out
}
