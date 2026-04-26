make_multilevel_exponential_stub_fit <- function() {
  N <- 2
  n_time <- 4
  variables <- c(
    paste0("phi_bar[", 1:4, "]"),
    paste0("tau_phi[", 1:4, "]"),
    paste0(
      "z_phi[",
      rep(seq_len(N), each = 4),
      ",",
      rep(seq_len(4), times = N),
      "]"
    ),
    "eta[1]", "eta[2]",
    paste0(
      "phi_unit[",
      rep(seq_len(N), each = 4),
      ",",
      rep(seq_len(4), times = N),
      "]"
    ),
    "sigma_exp[1]", "sigma_exp[2]",
    "rho",
    paste0(
      "log_lik[",
      rep(seq_len(N), each = n_time - 1),
      ",",
      rep(seq_len(n_time - 1), times = N),
      "]"
    )
  )
  draws <- array(
    rnorm(40 * length(variables), 0, 0.1),
    dim = c(40, 1, length(variables)),
    dimnames = list(NULL, NULL, variables)
  )
  draws[, , paste0("tau_phi[", 1:4, "]")] <- runif(40 * 4, 0.05, 0.2)
  draws[, , c("sigma_exp[1]", "sigma_exp[2]")] <- runif(40 * 2, 0.5, 1.5)
  draws[, , "rho"] <- runif(40, -0.4, 0.4)
  log_lik_vars <- grep("^log_lik\\[", variables, value = TRUE)
  draws[, , log_lik_vars] <- -abs(draws[, , log_lik_vars])

  stan_data <- structure(
    list(
      N = N,
      n_time = n_time,
      y = list(matrix(rnorm(n_time * 2), n_time, 2), matrix(rnorm(n_time * 2), n_time, 2)),
      skew_direction = c(1, -1)
    ),
    ids = c("a", "b"),
    time_values = 11:14,
    margins = "exponential",
    skew_direction = c(1, -1),
    person_means = matrix(0, N, 2)
  )

  new_dcvar_multilevel_fit(
    fit = posterior::as_draws_array(draws),
    stan_data = stan_data,
    N = N,
    vars = c("y1", "y2"),
    centered = TRUE,
    person_means = attr(stan_data, "person_means"),
    margins = "exponential",
    skew_direction = c(1, -1),
    backend = "rstan",
    priors = list(),
    meta = list(chains = 1, iter_warmup = 10, iter_sampling = 40)
  )
}

test_that("prepare_multilevel_data handles exponential margins", {
  df <- data.frame(
    id = rep(1:2, each = 4),
    time = rep(1:4, times = 2),
    y1 = rnorm(8),
    y2 = rnorm(8)
  )

  expect_error(
    prepare_multilevel_data(df, vars = c("y1", "y2"), margins = "exponential"),
    "skew_direction"
  )

  out <- prepare_multilevel_data(
    df,
    vars = c("y1", "y2"),
    margins = "exponential",
    skew_direction = c(1, -1)
  )
  expect_equal(out$skew_direction, c(1, -1))
  expect_equal(attr(out, "margins"), "exponential")
  expect_equal(.margin_stan_file("multilevel", "exponential"), "multilevel_EG.stan")
})

test_that("multilevel exponential stub fit extractors and loo work", {
  fit <- make_multilevel_exponential_stub_fit()

  co <- coef(fit)
  expect_named(co, c("phi_bar", "tau_phi", "sigma_exp", "rho"))
  expect_true(length(co$sigma_exp) == 2)

  vp <- var_params(fit)
  expect_named(vp, c("phi_bar", "tau_phi", "sigma_exp", "rho"))

  s <- summary(fit)
  expect_s3_class(s, "dcvar_multilevel_summary")
  expect_true("random_effects" %in% names(s))

  out <- suppressWarnings(loo::loo(fit))
  expect_s3_class(out, "loo")
})
