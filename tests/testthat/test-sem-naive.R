make_sem_naive_stub_fit <- function(margins = c("normal", "exponential", "skew_normal", "gamma")) {
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
  } else if (identical(margins, "skew_normal")) {
    variables <- c(variables, "omega[1]", "omega[2]", "delta[1]", "delta[2]")
  } else if (identical(margins, "gamma")) {
    variables <- c(
      variables,
      "eta[1]", "eta[2]",
      "sigma_gam[1]", "sigma_gam[2]",
      "shape_gam[1]", "shape_gam[2]"
    )
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
  } else if (identical(margins, "skew_normal")) {
    draws[, , c("omega[1]", "omega[2]")] <- runif(40 * 2, 0.5, 1.5)
    draws[, , c("delta[1]", "delta[2]")] <- runif(40 * 2, -0.5, 0.5)
  } else if (identical(margins, "gamma")) {
    draws[, , c("sigma_gam[1]", "sigma_gam[2]")] <- runif(40 * 2, 0.5, 1.5)
    draws[, , c("shape_gam[1]", "shape_gam[2]")] <- runif(40 * 2, 0.5, 2.0)
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
  if (margins %in% c("exponential", "gamma")) {
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

test_that("SEM naive homogeneous skew_normal and gamma stubs report mixed-engine params", {
  skew_fit <- make_sem_naive_stub_fit("skew_normal")
  skew_co <- coef(skew_fit)
  skew_vp <- var_params(skew_fit)
  expect_named(skew_co, c("mu", "Phi", "omega", "delta", "rho"))
  expect_equal(names(skew_co$omega), c("omega[1]", "omega[2]"))
  expect_true(all(c("omega", "delta") %in% names(skew_vp)))

  gamma_fit <- make_sem_naive_stub_fit("gamma")
  gamma_co <- coef(gamma_fit)
  gamma_vp <- var_params(gamma_fit)
  expect_named(gamma_co, c("mu", "Phi", "sigma_gam", "shape_gam", "rho"))
  expect_equal(names(gamma_co$shape_gam), c("shape_gam[1]", "shape_gam[2]"))
  expect_true(all(c("sigma_gam", "shape_gam") %in% names(gamma_vp)))

  out <- suppressWarnings(loo::loo(gamma_fit))
  expect_s3_class(out, "loo")
})

test_that("SEM naive fits homogeneous skew_normal and gamma margins", {
  skip_if_no_rstan()
  skip_if_not_slow()

  for (margin in c("skew_normal", "gamma")) {
    J <- 2
    skew_direction <- if (identical(margin, "gamma")) c(1, 1) else NULL
    sim <- simulate_dcvar_sem(
      n_time = 80, J = J, lambda = rep(0.8, J), sigma_e = sqrt(0.2),
      margins = margin,
      skew_direction = skew_direction,
      skew_params = if (identical(margin, "skew_normal")) list(alpha = c(2, 2)) else list(shape = 2),
      rho = 0.5,
      seed = 700 + match(margin, c("skew_normal", "gamma"))
    )
    fit <- dcvar_sem(
      sim$data,
      indicators = list(latent1 = paste0("y1_", seq_len(J)), latent2 = paste0("y2_", seq_len(J))),
      J = J,
      method = "naive",
      margins = margin,
      skew_direction = skew_direction,
      chains = 2, iter_warmup = 250, iter_sampling = 250,
      adapt_delta = 0.999, max_treedepth = 13, refresh = 0, seed = 123
    )

    expect_equal(fit$method, "naive")
    expect_equal(fit$margins, margin)
    expect_equal(fit$stan_data$family, rep(unname(.family_codes[[margin]]), 2))
    expect_true(coef(fit)$rho > -0.2)
  }
})
