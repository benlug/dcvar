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
