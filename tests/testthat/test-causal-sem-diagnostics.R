make_causal_diagnostic_fit <- function() {
  set.seed(71)
  n <- 2000L
  chains <- 4L
  beta <- matrix(stats::rnorm(n * chains), n, chains)
  # Each marginal mixes well. The joint dependence differs between chains.
  mu <- beta * rep(rep(c(1, -1), each = n), 2)
  small <- function() matrix(stats::rnorm(n * chains, sd = 0.001), n, chains)
  fields <- list(
    "alpha[1]" = small(), "alpha[2]" = small(),
    "beta[1]" = small(), "beta[2]" = beta,
    "mu_v[1]" = matrix(0, n, chains), "mu_v[2]" = mu,
    mu_v_treated = mu
  )
  a <- array(unlist(fields, use.names = FALSE), c(n, chains, length(fields)),
              dimnames = list(NULL, NULL, names(fields)))
  structure(list(
    fit = posterior::as_draws_array(a), model = "latent_covariate",
    backend = "rstan",
    meta = list(target_weights = c(0.5, 0.5), roles = list(outcome = "y"),
                scales = list(outcome = list(center = 0, scale = 1),
                              covariate = list(center = 0, scale = 1)))
  ), class = "dcvar_causal_sem_fit")
}

test_that("causal diagnostics detect failed dependence in effect draws", {
  object <- make_causal_diagnostic_fit()
  result <- dcvar_diagnostics(object)
  expect_lt(result$max_parameter_rhat, 1.01)
  expect_gt(result$min_parameter_ess_bulk, 400)
  expect_gt(result$min_parameter_ess_tail, 400)
  expect_gt(result$max_effect_rhat, 1.5)
  expect_lt(result$min_effect_ess_bulk, 20)
  expect_equal(result$max_rhat, result$max_effect_rhat)
  expect_equal(result$min_ess_bulk, result$min_effect_ess_bulk)
  expect_equal(result$min_ess_tail, result$min_effect_ess_tail)
  expect_true(result$high_effect_mcse)
  expect_gt(result$max_effect_mcse_ratio, 0.05)
  effects <- causal_effects(object)
  expect_equal(result$effects$mcse_ratio, effects$mcse_mean / effects$sd)
  expect_false(result$low_ebfmi)
  expect_true(all(is.na(result$ebfmi)))
})

test_that("causal diagnostics report low energy movement and missing energy", {
  object <- make_causal_diagnostic_fit()
  energy <- array(c(seq_len(100), rep(NA_real_, 100)), c(100, 2, 1),
                  dimnames = list(NULL, NULL, "energy__"))
  testthat::local_mocked_bindings(
    .fit_sampler_diagnostics = function(fit, backend) energy,
    .package = "dcvar"
  )
  result <- dcvar_diagnostics(object)
  expect_true(result$low_ebfmi)
  expect_equal(result$ebfmi[1], 1 / stats::var(seq_len(100)))
  expect_true(is.na(result$ebfmi[2]))
  energy[] <- NA_real_
  result <- dcvar_diagnostics(object)
  expect_false(result$low_ebfmi)
  expect_true(all(is.na(result$ebfmi)))
})

test_that("the causal fit warns when effect dependence fails", {
  object <- make_causal_diagnostic_fit()
  testthat::local_mocked_bindings(
    .causal_sem_stan_path = function(model) "unused.stan",
    .compile_model_backend = function(...) NULL,
    .sample_model = function(...) object$fit,
    .fit_diagnostic_summary = function(fit, backend) {
      list(num_divergent = 0L, num_max_treedepth = 0L)
    },
    .package = "dcvar"
  )
  sim <- simulate_dcvar_causal_sem(n = 20, model = "latent_covariate", seed = 71)
  expect_warning(
    do.call(dcvar_causal_sem, c(sim$args, list(backend = "rstan", refresh = 0))),
    "Causal SEM sampling needs review"
  )
})

test_that("the causal fit warns for energy and Monte Carlo flags", {
  testthat::local_mocked_bindings(
    .causal_sem_stan_path = function(model) "unused.stan",
    .compile_model_backend = function(...) NULL,
    .sample_model = function(...) NULL,
    dcvar_diagnostics = function(object, ...) diagnostics,
    .package = "dcvar"
  )
  sim <- simulate_dcvar_causal_sem(n = 20, model = "latent_covariate", seed = 71)
  diagnostics <- list(n_divergent = 0, n_max_treedepth = 0, max_rhat = 1,
                      min_ess_bulk = 1000, min_ess_tail = 1000,
                      high_effect_mcse = FALSE, low_ebfmi = TRUE)
  expect_warning(
    do.call(dcvar_causal_sem, c(sim$args, list(backend = "rstan", refresh = 0))),
    "Causal SEM sampling needs review"
  )
  diagnostics$low_ebfmi <- FALSE
  diagnostics$high_effect_mcse <- TRUE
  expect_warning(
    do.call(dcvar_causal_sem, c(sim$args, list(backend = "rstan", refresh = 0))),
    "Causal SEM sampling needs review"
  )
  diagnostics$high_effect_mcse <- FALSE
  expect_no_warning(
    do.call(dcvar_causal_sem, c(sim$args, list(backend = "rstan", refresh = 0)))
  )
})

test_that("the causal fit warns when a sampled parameter has no diagnostics", {
  object <- make_causal_diagnostic_fit()
  set.seed(176)
  mu <- matrix(stats::rnorm(8000), 2000, 4)
  object$fit[, , "mu_v[2]"] <- mu
  object$fit[, , "mu_v_treated"] <- mu
  constant <- array(1, c(2000, 4, 1),
    dimnames = list(NULL, NULL, "lambda_free[1]"))
  object$fit <- posterior::bind_draws(object$fit,
    posterior::as_draws_array(constant), along = "variable")
  testthat::local_mocked_bindings(
    .causal_sem_stan_path = function(model) "unused.stan",
    .compile_model_backend = function(...) NULL,
    .sample_model = function(...) object$fit,
    .fit_sampler_diagnostics = function(fit, backend) {
      array(stats::rnorm(8000), c(2000, 4, 1),
        dimnames = list(NULL, NULL, "energy__"))
    },
    .fit_diagnostic_summary = function(fit, backend) {
      list(num_divergent = 0L, num_max_treedepth = 0L)
    },
    .package = "dcvar"
  )
  result <- dcvar_diagnostics(object)
  expect_lt(result$max_rhat, 1.01)
  expect_gt(result$min_ess_bulk, 400)
  expect_gt(result$min_ess_tail, 400)
  expect_true(result$incomplete_diagnostics)
  sim <- simulate_dcvar_causal_sem(n = 20, model = "latent_covariate", seed = 71)
  expect_warning(
    do.call(dcvar_causal_sem, c(sim$args, list(backend = "rstan", refresh = 0))),
    "Causal SEM sampling needs review"
  )
})
