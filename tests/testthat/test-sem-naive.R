make_sem_naive_stub_fit <- function(margins = c("normal", "exponential")) {
  margins <- match.arg(margins)
  n_time <- 4
  variables <- c(
    "mu[1]", "mu[2]",
    "phi11", "phi12", "phi21", "phi22",
    "Phi[1,1]", "Phi[1,2]", "Phi[2,1]", "Phi[2,2]",
    "rho_raw", "rho",
    paste0("log_lik[", seq_len(n_time), "]")
  )
  if (identical(margins, "exponential")) {
    variables <- c(variables, "eta[1]", "eta[2]", "sigma_exp[1]", "sigma_exp[2]")
  } else {
    variables <- c(variables, "sigma[1]", "sigma[2]")
  }

  draws <- array(
    rnorm(40 * length(variables), 0, 0.1),
    dim = c(40, 1, length(variables)),
    dimnames = list(NULL, NULL, variables)
  )
  draws[, , "rho"] <- runif(40, -0.5, 0.5)
  draws[, , paste0("log_lik[", seq_len(n_time), "]")] <-
    -abs(draws[, , paste0("log_lik[", seq_len(n_time), "]")])
  if (identical(margins, "exponential")) {
    draws[, , c("sigma_exp[1]", "sigma_exp[2]")] <- runif(40 * 2, 0.5, 1.5)
  } else {
    draws[, , c("sigma[1]", "sigma[2]")] <- runif(40 * 2, 0.5, 1.5)
  }

  indicators <- list(
    latent1 = c("y1_1", "y1_2"),
    latent2 = c("y2_1", "y2_2")
  )
  stan_data <- structure(
    list(
      n_time = n_time,
      y = matrix(rnorm(n_time * 2), n_time, 2)
    ),
    time_values = 21:24,
    indicators = indicators,
    vars = names(indicators),
    margins = margins,
    method = "naive",
    J = 2
  )
  if (identical(margins, "exponential")) {
    stan_data$skew_direction <- c(1, -1)
    attr(stan_data, "skew_direction") <- c(1, -1)
  }

  new_dcvar_sem_fit(
    fit = posterior::as_draws_array(draws),
    stan_data = stan_data,
    vars = names(indicators),
    J = 2,
    lambda = NULL,
    sigma_e = NULL,
    indicators = indicators,
    margins = margins,
    method = "naive",
    skew_direction = attr(stan_data, "skew_direction"),
    backend = "rstan",
    priors = list(),
    meta = list(chains = 1, iter_warmup = 10, iter_sampling = 40)
  )
}

test_that("prepare_sem_data builds row-mean scores for naive SEM", {
  df <- data.frame(
    time = 1:4,
    y1_1 = c(1, 2, 3, 4),
    y1_2 = c(3, 4, 5, 6),
    y2_1 = c(2, 4, 6, 8),
    y2_2 = c(4, 6, 8, 10)
  )
  indicators <- list(
    latent1 = c("y1_1", "y1_2"),
    latent2 = c("y2_1", "y2_2")
  )

  out <- prepare_sem_data(
    df,
    indicators = indicators,
    J = 2,
    lambda = NULL,
    sigma_e = NULL,
    method = "naive"
  )

  expect_equal(out$y[, 1], rowMeans(df[, indicators[[1]]]))
  expect_equal(out$y[, 2], rowMeans(df[, indicators[[2]]]))
  expect_equal(attr(out, "method"), "naive")
  expect_equal(.margin_stan_file("sem_naive", "normal"), "sem_naive_NG.stan")
  expect_equal(.margin_stan_file("sem_naive", "exponential"), "sem_naive_EG.stan")
})

test_that("SEM naive normal stub fit reports sigma and supports loo", {
  fit <- make_sem_naive_stub_fit("normal")

  co <- coef(fit)
  expect_named(co, c("mu", "Phi", "sigma", "rho"))

  vp <- var_params(fit)
  expect_named(vp, c("mu", "Phi", "sigma", "rho"))

  s <- summary(fit)
  expect_s3_class(s, "dcvar_sem_summary")
  expect_equal(s$method, "naive")

  out <- suppressWarnings(loo::loo(fit))
  expect_s3_class(out, "loo")
  expect_error(latent_states(fit), "naive SEM")
})

test_that("SEM naive exponential stub fit reports sigma_exp and supports loo", {
  fit <- make_sem_naive_stub_fit("exponential")

  co <- coef(fit)
  expect_named(co, c("mu", "Phi", "sigma_exp", "rho"))

  vp <- var_params(fit)
  expect_named(vp, c("mu", "Phi", "sigma_exp", "rho"))

  out <- suppressWarnings(loo::loo(fit))
  expect_s3_class(out, "loo")
})
