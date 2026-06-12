# Regression test: the seed argument must reproduce the full MCMC output,
# including the default chain initial values (which previously came from the
# unseeded global R RNG).

test_that("identical seeds reproduce identical draws", {
  skip_if_no_rstan()

  sim <- simulate_dcvar(n_time = 25, rho_trajectory = rho_constant(25, 0.4), seed = 9)
  fit_args <- list(
    sim$Y_df,
    vars = c("y1", "y2"),
    chains = 2,
    iter_warmup = 50,
    iter_sampling = 50,
    refresh = 0,
    seed = 123
  )
  fit1 <- suppressWarnings(do.call(dcvar_constant, fit_args))
  fit2 <- suppressWarnings(do.call(dcvar_constant, fit_args))

  d1 <- posterior::as_draws_matrix(draws(fit1))
  d2 <- posterior::as_draws_matrix(draws(fit2))
  expect_equal(unclass(d1), unclass(d2))
})

test_that("seeded default inits do not disturb the caller's RNG stream", {
  init <- function() list(x = rnorm(1))
  set.seed(77)
  expected <- rnorm(3)

  set.seed(77)
  inits <- .seeded_chain_inits(init, chains = 2, seed = 5)
  expect_length(inits, 2)
  expect_false(identical(inits[[1]]$x, inits[[2]]$x))
  expect_equal(rnorm(3), expected)

  # Same seed -> same inits
  inits2 <- .seeded_chain_inits(init, chains = 2, seed = 5)
  expect_identical(inits, inits2)
})
