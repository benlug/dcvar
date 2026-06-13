test_that("cmdstanr backend fits the core models", {
  skip_if_no_cmdstanr_backend()

  dc_fit <- get_dcvar_cmdstanr_fit()
  hmm_fit <- get_hmm_cmdstanr_fit()
  con_fit <- get_constant_cmdstanr_fit()

  expect_s3_class(dc_fit, "dcvar_fit")
  expect_s3_class(hmm_fit, "dcvar_hmm_fit")
  expect_s3_class(con_fit, "dcvar_constant_fit")

  expect_equal(dc_fit$backend, "cmdstanr")
  expect_equal(hmm_fit$backend, "cmdstanr")
  expect_equal(con_fit$backend, "cmdstanr")
})

test_that("cmdstanr fits support the standard extractors", {
  skip_if_no_cmdstanr_backend()

  dc_fit <- get_dcvar_cmdstanr_fit()
  hmm_fit <- get_hmm_cmdstanr_fit()
  con_fit <- get_constant_cmdstanr_fit()

  expect_s3_class(summary(dc_fit), "dcvar_summary")
  expect_named(coef(dc_fit), c("mu", "Phi", "sigma_eps", "sigma_omega"))
  expect_s3_class(rho_trajectory(dc_fit), "data.frame")
  expect_type(dcvar_diagnostics(dc_fit), "list")
  expect_s3_class(fitted(dc_fit), "data.frame")
  expect_s3_class(predict(dc_fit), "data.frame")
  expect_s3_class(as.data.frame(dc_fit), "data.frame")
  expect_s3_class(loo(dc_fit), "loo")

  expect_named(coef(hmm_fit), c("mu", "Phi", "sigma_eps", "z_rho", "rho_state"))
  expect_type(dcvar_diagnostics(hmm_fit), "list")
  expect_type(hmm_states(hmm_fit), "list")
  expect_s3_class(fitted(hmm_fit), "data.frame")
  expect_s3_class(predict(hmm_fit), "data.frame")
  expect_s3_class(loo(hmm_fit), "loo")

  expect_named(coef(con_fit), c("mu", "Phi", "sigma_eps", "rho"))
  expect_type(dcvar_diagnostics(con_fit), "list")
  expect_s3_class(fitted(con_fit), "data.frame")
  expect_s3_class(predict(con_fit), "data.frame")
  expect_s3_class(plot(con_fit, type = "rho"), "ggplot")
  expect_s3_class(loo(con_fit), "loo")
})

test_that("cmdstanr backend fits partial time-varying models (zero-sized init params omitted)", {
  skip_if_no_cmdstanr_backend()

  sim <- simulate_dcvar(
    n_time = 30,
    rho_trajectory = rho_constant(30, 0.4),
    phi_trajectory = list(
      rho_decreasing(30, 0.3, 0.05), rho_constant(30, 0.1),
      rho_constant(30, 0.1), rho_constant(30, 0.3)
    ),
    seed = 91
  )

  # tv_phi only: sigma_raw / tau_sigma are inactive (zero-sized) -> must be
  # omitted from the default init so cmdstanr does not abort on a shapeless
  # empty matrix.
  fit_phi <- dcvar(
    sim$Y_df, vars = c("y1", "y2"),
    tv_phi = "ar", tv_sigma = FALSE,
    backend = "cmdstanr",
    chains = 1, iter_warmup = 50, iter_sampling = 50, refresh = 0, seed = 7
  )
  expect_s3_class(fit_phi, "dcvar_tv_fit")
  expect_equal(fit_phi$backend, "cmdstanr")
  expect_s3_class(phi_trajectory(fit_phi), "data.frame")

  # tv_sigma only: phi_raw / tau_phi are inactive (zero-sized).
  sim2 <- simulate_dcvar(
    n_time = 30,
    rho_trajectory = rho_constant(30, 0.4),
    sigma_trajectory = cbind(seq(0.8, 1.4, length.out = 29), rep(1, 29)),
    seed = 92
  )
  fit_sig <- dcvar(
    sim2$Y_df, vars = c("y1", "y2"),
    tv_phi = FALSE, tv_sigma = TRUE,
    backend = "cmdstanr",
    chains = 1, iter_warmup = 50, iter_sampling = 50, refresh = 0, seed = 8
  )
  expect_s3_class(fit_sig, "dcvar_tv_fit")
  expect_equal(fit_sig$backend, "cmdstanr")
  expect_s3_class(sigma_trajectory(fit_sig), "data.frame")
})
