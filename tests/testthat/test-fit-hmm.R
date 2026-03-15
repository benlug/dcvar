test_that("dcVar_hmm() returns a dcVar_hmm_fit object", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  expect_s3_class(fit, "dcVar_hmm_fit")
  expect_s3_class(fit, "dcVar_model_fit")
  expect_equal(fit$model, "hmm")
  expect_equal(fit$K, 2)
})

test_that("print.dcVar_hmm_fit runs without error", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  out <- capture.output(print(fit))
  expect_true(any(grepl("HMM", out)))
  expect_true(any(grepl("State", out)))
})

test_that("summary.dcVar_hmm_fit returns expected class", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  s <- summary(fit)
  expect_s3_class(s, "dcVar_hmm_summary")
  out <- capture.output(print(s))
  expect_true(any(grepl("HMM", out)))
})

test_that("coef.dcVar_hmm_fit returns named list", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  co <- coef(fit)
  expect_type(co, "list")
  expect_named(co, c("mu", "Phi", "sigma_eps", "z_rho", "rho_state"))
})

test_that("plot.dcVar_hmm_fit dispatches without error", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  expect_s3_class(plot(fit, type = "rho"), "ggplot")
  expect_s3_class(plot(fit, type = "phi"), "ggplot")
})

test_that("dcVar_hmm() fits exponential margins", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_exponential_fit()
  expect_s3_class(fit, "dcVar_hmm_fit")
  expect_equal(fit$margins, "exponential")
  expect_equal(fit$K, 2)
})
