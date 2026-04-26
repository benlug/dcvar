test_that("prepare_dcvar_covariate_data aligns covariates after sorting and outcome filtering", {
  df <- data.frame(
    time = c(4, 1, 3, 2),
    y1 = c(40, 10, NA, 20),
    y2 = c(44, 11, 33, 22),
    phase = c(4, 1, 3, 2)
  )

  expect_warning(
    out <- prepare_dcvar_covariate_data(
      df,
      vars = c("y1", "y2"),
      covariates = "phase",
      standardize = FALSE,
      allow_gaps = TRUE
    ),
    "Removing 1 row"
  )

  expect_equal(out$n_time, 3)
  expect_equal(as.numeric(out$Y[, 1]), c(10, 20, 40))
  expect_equal(as.numeric(out$X[, 1]), c(1, 2, 4))
  expect_equal(attr(out, "time_values"), c(1, 2, 4))
  expect_equal(attr(out, "covariates"), "phase")
})

test_that("prepare_dcvar_covariate_data validates covariates", {
  df <- data.frame(
    time = 1:5,
    y1 = rnorm(5),
    y2 = rnorm(5),
    phase = c(0, 0, 1, 1, 1),
    label = letters[1:5],
    missing_x = c(0, 1, NA, 1, 0),
    constant_x = 1
  )

  expect_error(
    prepare_dcvar_covariate_data(df, vars = c("y1", "y2"), covariates = "missing_col"),
    "not found"
  )
  expect_error(
    prepare_dcvar_covariate_data(df, vars = c("y1", "y2"), covariates = "label"),
    "must be numeric"
  )
  expect_error(
    prepare_dcvar_covariate_data(df, vars = c("y1", "y2"), covariates = "missing_x"),
    "complete finite"
  )
  expect_error(
    prepare_dcvar_covariate_data(
      df,
      vars = c("y1", "y2"),
      covariates = "constant_x",
      standardize_covariates = TRUE
    ),
    "zero"
  )
})

test_that("prepare_dcvar_covariate_data allows binary phase indicators by default", {
  df <- data.frame(
    time = 1:6,
    y1 = rnorm(6),
    y2 = rnorm(6),
    phase = c(0, 0, 0, 1, 1, 1)
  )

  out <- prepare_dcvar_covariate_data(df, vars = c("y1", "y2"), covariates = "phase")

  expect_equal(as.numeric(out$X[, 1]), df$phase)
  expect_false(attr(out, "standardized_covariates"))
  expect_equal(out$P, 1)
})

test_that("dcvar_stan_path exposes covariate Stan models", {
  drift_path <- dcvar_stan_path("dcvar_covariate")
  nodrift_path <- dcvar_stan_path("dcvar_covariate_nodrift")

  expect_true(file.exists(drift_path))
  expect_true(file.exists(nodrift_path))
  expect_match(basename(drift_path), "^dcvar_covariate_ncp\\.stan$")
  expect_match(basename(nodrift_path), "^dcvar_covariate_nodrift\\.stan$")
  expect_error(dcvar_stan_path("dcvar_covariate", margins = "gamma"), "normal")
})

make_covariate_stub_fit <- function(drift = TRUE) {
  variables <- c(
    "mu[1]", "mu[2]",
    "Phi[1,1]", "Phi[2,1]", "Phi[1,2]", "Phi[2,2]",
    "sigma_eps[1]", "sigma_eps[2]",
    "beta_0", "beta[1]", "beta[2]",
    paste0("rho[", 1:4, "]"),
    paste0("log_lik[", 1:4, "]")
  )
  if (drift) {
    variables <- c(variables, "sigma_omega")
  }

  draws <- array(
    rnorm(40 * length(variables), 0, 0.1),
    dim = c(40, 1, length(variables)),
    dimnames = list(NULL, NULL, variables)
  )
  rho_cols <- match(paste0("rho[", 1:4, "]"), variables)
  draws[, , rho_cols] <- tanh(draws[, , rho_cols])
  log_lik_cols <- match(paste0("log_lik[", 1:4, "]"), variables)
  draws[, , log_lik_cols] <- -abs(draws[, , log_lik_cols])

  stan_data <- structure(
    list(n_time = 5, D = 2, P = 2),
    time_values = 101:105
  )

  new_dcvar_covariate_fit(
    fit = posterior::as_draws_array(draws),
    stan_data = stan_data,
    vars = c("y1", "y2"),
    covariates = c("phase", "load"),
    standardized = TRUE,
    standardized_covariates = FALSE,
    drift = drift,
    zero_init_eta = TRUE,
    backend = "rstan",
    priors = list(),
    meta = list(chains = 1, iter_warmup = 10, iter_sampling = 40)
  )
}

test_that("covariate fit extractors work for drift and no-drift objects", {
  fit <- make_covariate_stub_fit(drift = TRUE)
  no_drift_fit <- make_covariate_stub_fit(drift = FALSE)

  expect_s3_class(fit, "dcvar_covariate_fit")
  expect_equal(fit$model, "dcvar_covariate")
  expect_equal(no_drift_fit$model, "dcvar_covariate_nodrift")

  co <- coef(fit)
  expect_named(co, c("mu", "Phi", "sigma_eps", "beta_0", "beta", "sigma_omega"))
  expect_equal(names(co$beta), c("phase", "load"))

  co_no_drift <- coef(no_drift_fit)
  expect_named(co_no_drift, c("mu", "Phi", "sigma_eps", "beta_0", "beta"))

  rho_df <- rho_trajectory(fit)
  expect_equal(rho_df$time, 102:105)
  expect_true(all(rho_df$mean >= -1 & rho_df$mean <= 1))

  effects <- covariate_effects(fit)
  expect_true(all(c("(Intercept)", "phase", "load", "sigma_omega") %in% effects$term))

  effects_no_drift <- covariate_effects(no_drift_fit)
  expect_false("sigma_omega" %in% effects_no_drift$term)
})

test_that("loo works for covariate fits with log_lik", {
  fit <- make_covariate_stub_fit(drift = TRUE)

  out <- suppressWarnings(loo::loo(fit))

  expect_s3_class(out, "loo")
})
