test_that("causal SEM simulation is reproducible and keeps factors separate", {
  models <- c("latent_covariate", "latent_outcome", "latent_mediator", "latent_mediator_baseline")
  for (model in models) {
    sim <- simulate_dcvar_causal_sem(n = 101, model = model, seed = 283)
    expect_identical(sim, simulate_dcvar_causal_sem(n = 101, model = model, seed = 283))
    expect_equal(nrow(sim$data), 101)
    expect_length(sim$latent, 101)
    expect_false(names(sim$args$indicators) %in% names(sim$data))
    expect_true(all(vapply(sim$data[paste0("u", 1:3)], is.ordered, logical(1))))
    expect_identical(sim$args$standardize, FALSE)
    expect_equal(sim$parameters$lambda[1], 1)
    expect_equal(sim$parameters$item_sd[1, ], rep(1, 3))
    expect_silent(do.call(prepare_causal_sem_data, sim$args))
  }
})

test_that("causal SEM simulation calculates independent regression effects", {
  sim <- simulate_dcvar_causal_sem(
    model = "latent_covariate", target_weights = c(0.2, 0.8), seed = 172,
    parameters = list(mu_v = c(0, 2), alpha = c(-1, 1), beta = c(0.2, 0.7))
  )
  expect_equal(sim$effects, c(ATE = 2.8, interaction = 0.5))
  zero <- simulate_dcvar_causal_sem(
    model = "latent_covariate", seed = 174,
    parameters = list(alpha = c(0.3, 0.3), beta = c(-0.2, -0.2))
  )
  expect_equal(zero$effects, c(ATE = 0, interaction = 0))
  outcome <- simulate_dcvar_causal_sem(
    model = "latent_outcome", seed = 172,
    parameters = list(mu_v = c(2, 3), beta = c(0.4, 0.6))
  )
  expect_equal(outcome$parameters$alpha[1], -0.8)
  expect_equal(outcome$effects, c(ATE = 1.8, interaction = 0.2))
})

test_that("mediation effects use every path and the weighted mediator mean", {
  model <- "latent_mediator_baseline"
  parameters <- list(
    mu_q = matrix(c(0, 1, 2, -1), 2, byrow = TRUE),
    alpha_m = c(1, 2), B = matrix(c(1, 2, 3, 4), 2, byrow = TRUE),
    alpha_y = c(-1, 1), D = matrix(c(0.5, 1, 1.5, 2), 2, byrow = TRUE), d = c(2, 3)
  )
  sim <- simulate_dcvar_causal_sem(model = model, parameters = parameters,
                                  target_weights = c(0.25, 0.75), seed = 923)
  # The target baseline mean is (1.5, -0.5). The target mediator mean is 3.75.
  # Total means are 2.25 and 15.75. Direct means are 6.75 and 13.5.
  expect_equal(sim$effects, c(ATE = 13.5, ADE = 6.75, AIE = 6.75))
  expect_equal(unname(sim$effects["ATE"]), unname(sum(sim$effects[c("ADE", "AIE")])))
  zero <- simulate_dcvar_causal_sem(model = "latent_mediator", seed = 921,
    parameters = list(B = matrix(0, 2, 2), D = matrix(0.2, 2, 2),
                      alpha_m = c(0, 0), alpha_y = c(1, 1), d = c(0.4, 0.4)))
  expect_equal(zero$effects, c(ATE = 0, ADE = 0, AIE = 0))
})

test_that("each mediation path changes the simulated outcome", {
  for (model in c("latent_mediator", "latent_mediator_baseline")) {
    reference <- simulate_dcvar_causal_sem(model = model, n = 120, seed = 211)
    for (name in c("B", "D")) {
      for (column in 1:2) {
        p <- reference$parameters
        p[[name]][2, column] <- p[[name]][2, column] + 1
        changed <- simulate_dcvar_causal_sem(model = model, n = 120, parameters = p, seed = 211)
        expect_false(isTRUE(all.equal(changed$data$y, reference$data$y)))
      }
    }
  }
})

test_that("causal SEM simulation preserves mixed category counts", {
  sim <- simulate_dcvar_causal_sem(seed = 671,
    thresholds = list(c(-0.5, 0.5), c(-1, 0, 1), c(-1.5, -0.5, 0.5, 1.5)))
  expect_identical(vapply(sim$data[paste0("u", 1:3)], nlevels, integer(1)),
                   c(u1 = 3L, u2 = 4L, u3 = 5L))
  expect_identical(do.call(prepare_causal_sem_data, sim$args)$stan_data$K, 3:5)
})

test_that("causal SEM simulation rejects invalid parameters and anchor conflicts", {
  expect_error(simulate_dcvar_causal_sem(n = 1), "at least one")
  expect_error(simulate_dcvar_causal_sem(n = c(10, NA)), "positive integer")
  expect_error(simulate_dcvar_causal_sem(parameters = list(unknown = 1)), "Unknown parameter")
  expect_error(simulate_dcvar_causal_sem(parameters = list(beta = NULL)), "finite numeric")
  expect_error(simulate_dcvar_causal_sem(parameters = list(lambda = c(0.8, 1, 1))), "first loading")
  expect_error(simulate_dcvar_causal_sem(parameters = list(item_sd = matrix(2, 2, 3))), "control item")
  expect_error(simulate_dcvar_causal_sem(parameters = list(mu_v = c(1, 1))), "factor mean")
  expect_error(simulate_dcvar_causal_sem(model = "latent_outcome",
    parameters = list(mu_v = c(1, 1), alpha = c(0, 0))), "factor mean")
  expect_error(simulate_dcvar_causal_sem(model = "latent_mediator",
    parameters = list(mu_q = matrix(1, 2, 2), alpha_m = c(0, 0))), "factor mean")
  expect_error(simulate_dcvar_causal_sem(model = "latent_mediator_baseline",
    parameters = list(mu_q = matrix(1, 2, 2))), "factor mean")
  expect_error(simulate_dcvar_causal_sem(thresholds = c(0, 0)), "strictly increasing")
  expect_error(simulate_dcvar_causal_sem(thresholds = list(NULL, c(-1, 1), c(-1, 1))),
               "finite numeric")
  expect_error(simulate_dcvar_causal_sem(parameters = list(sigma_y = c(0, 1))), "positive SDs")
  expect_error(simulate_dcvar_causal_sem(model = "latent_mediator",
    parameters = list(rho_q = c(-1, 0))), "strictly between")
})

test_that("causal SEM simulation accepts names on numeric parameters", {
  plain <- simulate_dcvar_causal_sem(parameters = list(mu_v = c(0, 1)), seed = 613)
  named <- simulate_dcvar_causal_sem(
    parameters = list(mu_v = c(control = 0, treated = 1)), seed = 613)
  expect_identical(named, plain)
})

test_that("prior draws obey all identification restrictions", {
  for (model in c("latent_covariate", "latent_outcome", "latent_mediator", "latent_mediator_baseline")) {
    set.seed(168)
    p <- .causal_sem_draw_prior(model, priors = list(loading_sd = 0.2), K = 3:5)
    expect_equal(p$lambda[1], 1)
    expect_equal(p$item_sd[1, ], rep(1, 3))
    expect_equal(vapply(p[paste0("threshold_", 1:3)], length, integer(1)),
                 c(threshold_1 = 2L, threshold_2 = 3L, threshold_3 = 4L))
    expect_true(all(vapply(p[paste0("threshold_", 1:3)], function(x) all(diff(x) > 0), logical(1))))
    if (model == "latent_covariate") expect_equal(p$mu_v[1], 0)
    if (model == "latent_outcome") expect_equal(p$alpha[1] + p$beta[1] * p$mu_v[1], 0)
    if (model == "latent_mediator") expect_equal(p$alpha_m[1] + sum(p$B[1, ] * p$mu_q[1, ]), 0)
    if (model == "latent_mediator_baseline") expect_equal(p$mu_q[1, 1], 0)
    expect_silent(simulate_dcvar_causal_sem(model = model, parameters = p, seed = 12))
  }
})

test_that("prior simulation draws once and retains constant indicators", {
  sim <- simulate_dcvar_causal_sem(parameter_prior = TRUE, seed = 641)
  expect_identical(sim, simulate_dcvar_causal_sem(parameter_prior = TRUE, seed = 641))
  expect_error(simulate_dcvar_causal_sem(parameter_prior = TRUE, thresholds = c(-1, 1)),
               "without explicit")
  constant <- simulate_dcvar_causal_sem(n = 20, thresholds = c(100, 101), seed = 417)
  expect_equal(length(unique(constant$data$u1)), 1)
  expect_equal(nrow(constant$data), 20)
  expect_error(do.call(prepare_causal_sem_data, constant$args), "constant")
})

test_that("causal SEM simulation validates seeds and finite parameter values", {
  for (seed in list(-1, 0.5, NA_real_, Inf, numeric(), c(1, 2), matrix(1), 2^31)) {
    expect_error(simulate_dcvar_causal_sem(seed = seed), "seed")
  }
  for (n in list(0, -1, 1.5, numeric(), c(2, 2, 2), Inf, matrix(2), 2^31,
                c(.Machine$integer.max, 1))) {
    expect_error(simulate_dcvar_causal_sem(n = n), "size")
  }
  for (value in c(NA_real_, NaN, Inf, -Inf)) {
    expect_error(simulate_dcvar_causal_sem(parameters = list(beta = c(1, value))), "finite")
    expect_error(simulate_dcvar_causal_sem(thresholds = c(-1, value)), "finite")
  }
  expect_error(.causal_sem_draw_prior("latent_covariate", K = c(3, 4, 1e100)), "integers")
  expect_error(.causal_sem_draw_prior("latent_covariate", K = c(3, 4, 4.5)), "integers")
  expect_error(.causal_sem_draw_prior("latent_covariate", K = matrix(3:5)), "integers")
  expect_error(simulate_dcvar_causal_sem(parameters = list(beta = c(0, 1), beta = c(1, 2))),
               "unique")
  expect_silent(simulate_dcvar_causal_sem(n = 20, seed = 0))
})

test_that("causal SEM simulation normalizes accepted target weights", {
  weights <- c(0.5, 0.500000005)
  normalized <- weights / sum(weights)
  approximate <- simulate_dcvar_causal_sem(seed = 15, target_weights = weights)
  reference <- simulate_dcvar_causal_sem(seed = 15, target_weights = normalized)
  expect_equal(approximate, reference, tolerance = 1e-14)
  expect_equal(sum(approximate$args$target_weights), 1, tolerance = 1e-14)
})
