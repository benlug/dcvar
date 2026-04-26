make_constant_clayton_stub_fit <- function() {
  variables <- c(
    "mu[1]", "mu[2]",
    "Phi[1,1]", "Phi[1,2]", "Phi[2,1]", "Phi[2,2]",
    "sigma_eps[1]", "sigma_eps[2]",
    "theta",
    paste0("log_lik[", 1:4, "]")
  )
  draws <- array(
    rnorm(40 * length(variables), 0, 0.1),
    dim = c(40, 1, length(variables)),
    dimnames = list(NULL, NULL, variables)
  )
  draws[, , "sigma_eps[1]"] <- runif(40, 0.8, 1.2)
  draws[, , "sigma_eps[2]"] <- runif(40, 0.8, 1.2)
  draws[, , "theta"] <- runif(40, 0.2, 2.5)
  draws[, , paste0("log_lik[", 1:4, "]")] <- -abs(draws[, , paste0("log_lik[", 1:4, "]")])

  stan_data <- structure(
    list(n_time = 5, D = 2),
    time_values = 101:105,
    copula = "clayton"
  )

  new_dcvar_constant_fit(
    fit = posterior::as_draws_array(draws),
    stan_data = stan_data,
    vars = c("y1", "y2"),
    standardized = TRUE,
    margins = "normal",
    copula = "clayton",
    backend = "rstan",
    priors = list(),
    meta = list(chains = 1, iter_warmup = 10, iter_sampling = 40)
  )
}

test_that("Clayton constant dispatch resolves the bundled Stan file", {
  expect_equal(.margin_stan_file("constant", "normal", copula = "clayton"),
               "constant_NCl.stan")
  expect_true(file.exists(dcvar_stan_path("constant", copula = "clayton")))
  expect_error(.margin_stan_file("constant", "exponential", copula = "clayton"),
               "normal margins")
})

test_that("Clayton constant stub fit extractors use theta and Kendall tau", {
  fit <- make_constant_clayton_stub_fit()

  co <- coef(fit)
  expect_named(co, c("mu", "Phi", "sigma_eps", "theta"))
  expect_false("rho" %in% names(co))

  vp <- var_params(fit)
  expect_true("theta" %in% names(vp))

  dep <- dependence_summary(fit, probs = c(0.025, 0.5, 0.975))
  expect_s3_class(dep, "data.frame")
  expect_equal(dep$time, 102:105)
  expect_true(all(dep$mean > 0 & dep$mean < 1))

  expect_error(rho_trajectory(fit), "dependence_summary")
})

test_that("Clayton constant summary, print, and loo work on stubs", {
  fit <- make_constant_clayton_stub_fit()

  expect_no_error(capture.output(print(fit)))
  s <- summary(fit)
  expect_s3_class(s, "dcvar_constant_summary")
  expect_equal(s$copula, "clayton")
  expect_s3_class(s$dependence, "data.frame")

  out <- suppressWarnings(loo::loo(fit))
  expect_s3_class(out, "loo")
})
