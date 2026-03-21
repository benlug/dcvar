#' Fit a multilevel copula VAR(1) model
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
#' @param prior_phi_bar_sd Prior SD for population-mean VAR coefficients.
#' @param prior_tau_phi_scale Prior scale for half-t(3) on tau_phi.
#' @param prior_sigma_sd Prior SD for half-normal on innovation SDs.
#' @param prior_rho_sd Prior SD for normal on rho.
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
#' @return A `dcvar_multilevel_fit` object.
#'
#' @details `adapt_delta` defaults to 0.90 and `max_treedepth` to 14 because
#'   the hierarchical structure with random effects benefits from deeper trees
#'   but does not require aggressive step-size adaptation.
#'
#' @note This model currently supports normal marginal distributions only.
#'   For non-normal margins, use [dcvar()], [dcvar_constant()], or [dcvar_hmm()].
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
                             prior_phi_bar_sd = 0.5,
                             prior_tau_phi_scale = 0.2,
                             prior_sigma_sd = 1,
                             prior_rho_sd = 0.5,
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
  bundled_stan <- dcvar_stan_path("multilevel")
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
    data, vars, id_var, time_var, center,
    prior_phi_bar_sd, prior_tau_phi_scale, prior_sigma_sd, prior_rho_sd
  )

  N <- stan_data$N
  T_obs <- stan_data$T

  cli_inform("Fitting multilevel copula VAR model (N = {N}, T = {T_obs})...")

  # Compile model
  model <- .compile_model("multilevel", stan_file = stan_file, backend = backend)

  # Default init
  if (is.null(init)) {
    init <- function() .init_multilevel_params(2, N)
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

  .report_sampling_outcome(fit, "Multilevel copula VAR", chains = chains, backend = backend)

  new_dcvar_multilevel_fit(
    fit = fit,
    stan_data = stan_data,
    N = N,
    vars = vars,
    centered = center,
    person_means = attr(stan_data, "person_means"),
    backend = backend,
    priors = list(
      phi_bar_sd = prior_phi_bar_sd,
      tau_phi_scale = prior_tau_phi_scale,
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
