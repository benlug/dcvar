test_that("dcVar_constant() returns a dcVar_constant_fit object", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  expect_s3_class(fit, "dcVar_constant_fit")
  expect_s3_class(fit, "dcVar_model_fit")
  expect_equal(fit$model, "constant")
})

test_that("print.dcVar_constant_fit runs without error", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  out <- capture.output(print(fit))
  expect_true(any(grepl("Constant", out)))
  expect_true(any(grepl("rho:", out)))
})

test_that("summary.dcVar_constant_fit returns expected class", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  s <- summary(fit)
  expect_s3_class(s, "dcVar_constant_summary")
  out <- capture.output(print(s))
  expect_true(any(grepl("Constant", out)))
})

test_that("summary.dcVar_constant_fit prints custom quantiles", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  s <- summary(fit, probs = c(0.1, 0.9))
  expect_no_error(capture.output(print(s)))
})

test_that("coef.dcVar_constant_fit returns named list", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  co <- coef(fit)
  expect_type(co, "list")
  expect_named(co, c("mu", "Phi", "sigma_eps", "rho"))
})

test_that("plot.dcVar_constant_fit dispatches without error", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  expect_s3_class(plot(fit, type = "phi"), "ggplot")
})

test_that("dcVar_constant() fits gamma margins", {
  skip_if_no_cmdstanr()

  fit <- get_constant_gamma_fit()
  expect_s3_class(fit, "dcVar_constant_fit")
  expect_equal(fit$margins, "gamma")
})

test_that("dcVar_constant() fits skew-normal margins", {
  skip_if_no_cmdstanr()
  skip_if_not_installed("sn")

  fit <- get_constant_skew_normal_fit()
  expect_s3_class(fit, "dcVar_constant_fit")
  expect_equal(fit$margins, "skew_normal")
})
