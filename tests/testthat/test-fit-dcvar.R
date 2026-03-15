test_that("dcVar() returns a dcVar_fit object", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  expect_s3_class(fit, "dcVar_fit")
  expect_s3_class(fit, "dcVar_model_fit")
  expect_equal(fit$model, "dcVar")
})

test_that("print.dcVar_fit runs without error", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  out <- capture.output(print(fit))
  expect_true(any(grepl("DC-VAR", out)))
  expect_true(any(grepl("rho range", out)))
})

test_that("summary.dcVar_fit returns expected class", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  s <- summary(fit)
  expect_s3_class(s, "dcVar_summary")
  out <- capture.output(print(s))
  expect_true(any(grepl("DC-VAR", out)))
})

test_that("coef.dcVar_fit returns named list", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  co <- coef(fit)
  expect_type(co, "list")
  expect_named(co, c("mu", "Phi", "sigma_eps", "sigma_omega"))
})

test_that("plot.dcVar_fit dispatches without error", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  expect_s3_class(plot(fit, type = "rho"), "ggplot")
  expect_s3_class(plot(fit, type = "phi"), "ggplot")
})

test_that("dcVar() fits exponential margins", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_exponential_fit()
  expect_s3_class(fit, "dcVar_fit")
  expect_equal(fit$margins, "exponential")
})
