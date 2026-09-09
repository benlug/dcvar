#' Fit a Bayesian causal SEM with ordinal indicators
#'
#' Fits one of four models from the SEM causal effects study. Each model has
#' two treatment groups and one factor with three ordinal indicators. The
#' models use normal structural errors and an ordinal probit measurement model.
#'
#' @inheritParams prepare_causal_sem_data
#' @param chains Number of MCMC chains.
#' @param iter_warmup Number of warmup iterations per chain.
#' @param iter_sampling Number of retained iterations per chain.
#' @param adapt_delta Target acceptance probability.
#' @param max_treedepth Maximum tree depth.
#' @param seed Optional integer seed.
#' @param cores Number of parallel chains. `NULL` uses the package default.
#' @param refresh Sampler progress interval. Use zero to suppress progress.
#' @param init Initial values accepted by the selected backend. `NULL` uses
#'   small random starts and ordered thresholds.
#' @param backend Sampling backend: `"auto"`, `"rstan"`, or `"cmdstanr"`.
#'   `"auto"` uses RStan.
#' @param ... Additional arguments for the selected sampler.
#'
#' @details
#' The first loading is one. The factor mean in the control group is zero.
#' Item residual standard deviations are one in the control group and free in
#' the treated group. Loadings and thresholds are shared across groups.
#' Endogenous factor variances follow from their regression and residual
#' variance. The factor scale does not require known simulation parameters.
#' The sampler conditions person factors on continuous outcomes when present.
#' This exact normal factorization reduces dependence during sampling. The
#' ordinal likelihood then updates these factors. The fit records the factor
#' parameterization in `sampling$factor_parameterization`.
#'
#' Treatment slopes and structural variances can differ across groups.
#' [causal_effects()] computes effects for each joint posterior draw. The
#' direct effect follows the paper's adjustment at the pooled mediator mean.
#' It is not a natural direct effect. A causal interpretation requires the
#' assumptions that identify the effects in the study design.
#'
#' `prior_only = TRUE` omits all observed likelihood terms. Predictive draws
#' still condition on supplied manifest predictors where the model uses them.
#' Use [simulate_dcvar_causal_sem()] with `parameter_prior = TRUE` to generate
#' a complete data set from the prior.
#'
#' @return A `dcvar_causal_sem_fit` object. It contains the backend fit, Stan
#'   data, priors, roles, group labels, scales, and sampling settings.
#' @seealso [simulate_dcvar_causal_sem()], [causal_effects()],
#'   [prepare_causal_sem_data()], [dcvar_diagnostics()]
#' @export
#' @examples
#' sim <- simulate_dcvar_causal_sem(n = 50, model = "latent_covariate", seed = 42)
#' prepared <- do.call(prepare_causal_sem_data, sim$args)
#' \dontrun{
#' fit <- do.call(dcvar_causal_sem, c(sim$args, list(seed = 42)))
#' causal_effects(fit)
#' }
dcvar_causal_sem <- function(data, model, treatment, control = 0, treated = 1,
                             outcome, covariate = NULL, mediator = NULL,
                             mediator_baseline = NULL, outcome_baseline = NULL,
                             indicators, target_weights = c(0.5, 0.5),
                             standardize = TRUE, priors = list(),
                             prior_only = FALSE, chains = 4,
                             iter_warmup = 1000, iter_sampling = 1000,
                             adapt_delta = 0.95, max_treedepth = 12,
                             seed = NULL, cores = NULL, refresh = 100,
                             init = NULL,
                             backend = getOption("dcvar.backend", "auto"), ...) {
  .validate_sampling_args(chains, iter_warmup, iter_sampling,
                          adapt_delta, max_treedepth)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      seed != floor(seed))) {
    cli_abort("{.arg seed} must be an integer from 0 to 2147483647.")
  }
  if (!is.numeric(refresh) || length(refresh) != 1L || !is.finite(refresh) ||
      refresh < 0 || refresh != floor(refresh)) {
    cli_abort("{.arg refresh} must be a nonnegative integer.")
  }
  prepared <- prepare_causal_sem_data(
    data = data, model = model, treatment = treatment, control = control,
    treated = treated, outcome = outcome, covariate = covariate,
    mediator = mediator, mediator_baseline = mediator_baseline,
    outcome_baseline = outcome_baseline, indicators = indicators,
    target_weights = target_weights, standardize = standardize,
    priors = priors, prior_only = prior_only
  )
  backend <- .resolve_backend(backend)
  cores <- .normalize_cores(cores, chains)
  stan_file <- .causal_sem_stan_path(prepared$meta$model)
  compiled <- .compile_model_backend(stan_file, backend, quiet = refresh == 0)
  if (is.null(init)) init <- .causal_sem_init(prepared$stan_data, prepared$meta$model)
  raw_fit <- .sample_model(
    compiled, prepared$stan_data, backend, chains, iter_warmup, iter_sampling,
    adapt_delta, max_treedepth, seed, cores, init, refresh, ...
  )
  out <- structure(list(
    fit = raw_fit, model = prepared$meta$model, backend = backend,
    stan_data = prepared$stan_data, priors = prepared$meta$priors,
    meta = prepared$meta,
    sampling = list(chains = chains, iter_warmup = iter_warmup,
                    iter_sampling = iter_sampling, adapt_delta = adapt_delta,
                    max_treedepth = max_treedepth, seed = seed, cores = cores,
                    factor_parameterization = if (prior_only || prepared$meta$model == "latent_outcome") {
                      "noncentered"
                    } else "conditional_noncentered")
  ), class = "dcvar_causal_sem_fit")
  diagnostics <- dcvar_diagnostics(out)
  if (isTRUE(diagnostics$n_divergent > 0) ||
      isTRUE(diagnostics$n_max_treedepth > 0) ||
      isTRUE(diagnostics$max_rhat >= 1.01) ||
      isTRUE(diagnostics$min_ess_bulk < 400) ||
      isTRUE(diagnostics$min_ess_tail < 400) ||
      isTRUE(diagnostics$high_effect_mcse) ||
      isTRUE(diagnostics$incomplete_diagnostics) ||
      isTRUE(diagnostics$low_ebfmi)) {
    cli_warn(c("Causal SEM sampling needs review.",
               "i" = "Inspect {.fun dcvar_diagnostics} before using the effects."))
  }
  out
}

#' Locate a causal SEM Stan model
#' @noRd
.causal_sem_stan_path <- function(model) {
  model <- .causal_sem_validate_model(model)
  system.file("stan", paste0("causal_sem_", model, ".stan"),
              package = "dcvar", mustWork = TRUE)
}

#' Set stable starts for the measurement model
#' @noRd
.causal_sem_init <- function(stan_data, model) {
  function() {
    out <- list(lambda_free = rep(1, 2), item_sd_treated = rep(1, 3),
                latent_raw = stats::rnorm(stan_data$N, 0, 0.1), sigma_y = rep(1, 2))
    for (j in seq_len(3L)) {
      out[[paste0("threshold_", j)]] <-
        stats::qnorm(seq_len(stan_data$K[j] - 1L) / stan_data$K[j])
    }
    if (model %in% c("latent_covariate", "latent_outcome")) {
      out$sigma_v <- rep(1, 2)
      out$beta <- rep(0, 2)
      if (model == "latent_covariate") {
        out$alpha <- rep(0, 2)
        out$mu_v_treated <- 0
      } else {
        out$alpha_treated <- 0
        out$mu_v <- rep(0, 2)
      }
    } else {
      out$sigma_q <- matrix(1, 2, 2)
      out$L_q <- array(0, c(2, 2, 2))
      for (g in seq_len(2L)) out$L_q[g, , ] <- diag(2)
      out$B <- out$D <- matrix(0, 2, 2)
      out$alpha_y <- out$d <- rep(0, 2)
      out$sigma_m <- rep(1, 2)
      if (model == "latent_mediator") {
        out$mu_q <- matrix(0, 2, 2)
        out$alpha_m_treated <- 0
      } else {
        out$mu_q1_treated <- 0
        out$mu_q2 <- rep(0, 2)
        out$alpha_m <- rep(0, 2)
      }
    }
    out
  }
}
