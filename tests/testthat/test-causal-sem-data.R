test_that("all causal SEM models produce the documented Stan data", {
  models <- c("latent_covariate", "latent_outcome", "latent_mediator", "latent_mediator_baseline")
  for (model in models) {
    sim <- simulate_dcvar_causal_sem(n = c(60, 80), model = model, seed = 94)
    prepared <- do.call(prepare_causal_sem_data, sim$args)
    expect_equal(prepared$stan_data$N, 140L)
    expect_identical(prepared$stan_data$J, 3L)
    expect_identical(prepared$stan_data$K, rep(3L, 3))
    expect_identical(dim(prepared$stan_data$u), c(140L, 3L))
    expect_identical(prepared$stan_data$group, c(rep(1L, 60), rep(2L, 80)))
    expect_identical(prepared$meta$n_group, c(control = 60L, treated = 80L))
    expect_equal(prepared$meta$model, model)
    expect_equal(prepared$meta$scales$outcome, list(center = 0, scale = 1))
    expect_equal(prepared$stan_data$prior_beta_sd, 1)
    expect_identical(prepared$stan_data$prior_only, 0L)
    expect_true(all(vapply(prepared$stan_data[c("y", "z", "m", "q1", "q2")],
                           length, integer(1)) == 140L))
  }
})

test_that("causal SEM preparation retains global ordinal levels", {
  sim <- simulate_dcvar_causal_sem(n = 120, seed = 11)
  sim$args$data$u1 <- ordered(c(rep("low", 60), rep("high", 60)),
                            levels = c("low", "middle", "high", "unused"))
  prepared <- do.call(prepare_causal_sem_data, sim$args)
  expect_identical(prepared$stan_data$K[1], 4L)
  expect_identical(prepared$meta$levels$u1, c("low", "middle", "high", "unused"))
  expect_identical(prepared$stan_data$u[, 1], c(rep(1L, 60), rep(3L, 60)))
  expect_equal(prepared$stan_data$N, nrow(sim$args$data))
})

test_that("causal SEM standardization uses both groups and records each role", {
  sim <- simulate_dcvar_causal_sem(n = c(40, 80), model = "latent_mediator", seed = 412)
  sim$args$standardize <- TRUE
  original <- sim$args$data
  prepared <- do.call(prepare_causal_sem_data, sim$args)
  mapping <- c(outcome = "y", mediator_baseline = "q1", outcome_baseline = "q2")
  for (role in names(mapping)) {
    values <- original[[sim$roles[[role]]]]
    transform <- prepared$meta$scales[[role]]
    expect_equal(transform$center, mean(values))
    expect_equal(transform$scale, sd(values))
    expect_equal(prepared$stan_data[[mapping[[role]]]], as.numeric(scale(values)))
  }
  expect_identical(prepared$meta$scales$mediator, list(center = 0, scale = 1))
  expect_identical(sim$args$data, original)
})

test_that("causal SEM groups use explicit labels", {
  sim <- simulate_dcvar_causal_sem(seed = 74)
  sim$args$data$x <- factor(ifelse(sim$args$data$x == 0, "usual", "new"),
                            levels = c("new", "usual"))
  sim$args$control <- "usual"
  sim$args$treated <- "new"
  prepared <- do.call(prepare_causal_sem_data, sim$args)
  expect_equal(prepared$stan_data$group, rep(1:2, each = 200))
  expect_identical(prepared$meta$group_labels, list(control = "usual", treated = "new"))
  sim$args$data$x[1] <- NA
  expect_error(do.call(prepare_causal_sem_data, sim$args), "complete vector")
})

test_that("causal SEM roles reject ambiguous or unused assignments", {
  sim <- simulate_dcvar_causal_sem(seed = 12)
  check <- function(change, pattern) {
    args <- utils::modifyList(sim$args, change, keep.null = TRUE)
    expect_error(do.call(prepare_causal_sem_data, args), pattern)
  }
  check(list(covariate = NULL), "covariate")
  check(list(mediator = "m"), "not used")
  check(list(outcome = "x"), "distinct")
  check(list(indicators = list(y = c("u1", "u2", "u3"))), "xi")
  check(list(indicators = list(xi = c("u1", "u2", "x"))), "cannot also")
  check(list(indicators = list(xi = c("u1", "u1", "u3"))), "distinct")
  check(list(outcome = "missing"), "not found")
  latent_data <- sim$args$data
  latent_data$xi <- sim$latent
  check(list(data = latent_data), "must not also be a data column")
  check(list(control = 1), "must differ")
  check(list(treated = 2), "unknown group")
  missing_group <- sim$args$data
  missing_group$x[] <- 0
  check(list(data = missing_group), "Both treatment groups")
})

test_that("causal SEM data validation rejects missing and invalid values", {
  sim <- simulate_dcvar_causal_sem(seed = 19)
  check_data <- function(data, pattern) {
    args <- sim$args
    args$data <- data
    expect_error(do.call(prepare_causal_sem_data, args), pattern)
  }
  bad <- sim$data
  bad$y[1] <- NA_real_
  check_data(bad, "complete finite numeric")
  bad$y[1] <- Inf
  check_data(bad, "complete finite numeric")
  bad <- sim$data
  bad$u1 <- factor(bad$u1, ordered = FALSE)
  check_data(bad, "ordered factor")
  bad$u1 <- ordered(rep(1, nrow(bad)), levels = 1:3)
  check_data(bad, "constant")
  bad$u1 <- ordered(rep(c(1, 2), length.out = nrow(bad)), levels = 1:2)
  check_data(bad, "three to five")
  bad <- sim$data
  bad$u1[1] <- NA
  check_data(bad, "missing values")
  bad <- sim$data
  bad$unused <- NA_real_
  args <- sim$args
  args$data <- bad
  expect_equal(do.call(prepare_causal_sem_data, args)$stan_data$N, nrow(bad))
})

test_that("causal SEM validates priors, weights, flags, and constant scales", {
  sim <- simulate_dcvar_causal_sem(seed = 171)
  check <- function(change, pattern) {
    args <- utils::modifyList(sim$args, change, keep.null = TRUE)
    expect_error(do.call(prepare_causal_sem_data, args), pattern)
  }
  check(list(target_weights = c(0.2, 0.2)), "sum to one")
  check(list(target_weights = c(-0.1, 1.1)), "non-negative")
  check(list(standardize = NA), "TRUE or FALSE")
  check(list(prior_only = 1), "TRUE or FALSE")
  check(list(priors = list(beta_sd = 0)), "positive finite")
  check(list(priors = list(beta_sd = NULL)), "positive finite")
  check(list(priors = list(unknown = 1)), "Unknown prior")
  sim$args$data$y[] <- 4
  sim$args$standardize <- TRUE
  expect_error(do.call(prepare_causal_sem_data, sim$args), "positive finite SD")
  sim$args$standardize <- FALSE
  sim$args$prior_only <- TRUE
  sim$args$priors <- list(beta_sd = 0.25, loading_mean = -0.2)
  prepared <- do.call(prepare_causal_sem_data, sim$args)
  expect_identical(prepared$stan_data$prior_only, 1L)
  expect_equal(prepared$stan_data$prior_beta_sd, 0.25)
  expect_equal(prepared$stan_data$prior_loading_mean, -0.2)
  expect_equal(prepared$stan_data$y, rep(4, 400))
})

test_that("weights with roundoff pass preparation and effect extraction", {
  weights <- c(0.5, 0.500000005)
  sim <- simulate_dcvar_causal_sem(n = 120, seed = 194)
  sim$args$target_weights <- weights
  prepared <- do.call(prepare_causal_sem_data, sim$args)
  expect_equal(prepared$meta$target_weights, weights / sum(weights), tolerance = 1e-14)
  expect_equal(sum(prepared$meta$target_weights), 1, tolerance = 1e-14)
  expect_identical(prepared$stan_data$target_weights, prepared$meta$target_weights)
  p <- sim$parameters
  values <- c(p$alpha, p$beta, p$mu_v)
  parameter_names <- c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]", "mu_v[1]", "mu_v[2]")
  a <- array(rep(values, each = 8), dim = c(4, 2, 6),
             dimnames = list(NULL, NULL, parameter_names))
  fit <- structure(list(fit = posterior::as_draws_array(a), model = sim$args$model,
                         backend = "rstan", meta = prepared$meta),
                    class = "dcvar_causal_sem_fit")
  effects <- causal_effects(fit)
  explicit <- causal_effects(fit, target_weights = weights)
  expect_equal(effects, explicit)
  expect_equal(attr(effects, "target_weights"), weights / sum(weights), tolerance = 1e-14)
  expected <- p$alpha[2] - p$alpha[1] +
    (p$beta[2] - p$beta[1]) * sum(weights * p$mu_v) / sum(weights)
  expect_equal(effects$mean[1], expected, tolerance = 1e-14)
})

test_that("causal SEM priors and weights reject nonfinite and nonscalar input", {
  for (value in list(NA_real_, NaN, Inf, -Inf, numeric(), c(1, 2), matrix(1))) {
    expect_error(.causal_sem_priors(list(beta_sd = value)), "positive finite")
  }
  expect_error(.causal_sem_priors(list(beta_sd = 1, beta_sd = 2)), "unique")
  expect_error(.causal_sem_priors(list(1)), "unique")
  expect_equal(.causal_sem_priors(list(beta_sd = c(custom = 2)))$beta_sd, 2)
  for (weights in list(c(NA, 1), c(Inf, -Inf), matrix(c(0.5, 0.5)), c(0, 0))) {
    expect_error(.causal_sem_validate_weights(weights), "target_weights")
  }
  expect_error(.causal_sem_validate_weights(c(0.5, 0.5000001)), "sum to one")
})
