test_that("dcvar() returns a dcvar_fit object", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  expect_s3_class(fit, "dcvar_fit")
  expect_s3_class(fit, "dcvar_model_fit")
  expect_equal(fit$model, "dcvar")
})

test_that("print.dcvar_fit runs without error", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  out <- capture.output(print(fit))
  expect_true(any(grepl("DC-VAR", out)))
  expect_true(any(grepl("rho range", out)))
})

test_that("summary.dcvar_fit returns expected class", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  s <- summary(fit)
  expect_s3_class(s, "dcvar_summary")
  out <- capture.output(print(s))
  expect_true(any(grepl("DC-VAR", out)))
})

test_that("coef.dcvar_fit returns named list", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  co <- coef(fit)
  expect_type(co, "list")
  expect_named(co, c("mu", "Phi", "sigma_eps", "sigma_omega"))
})

test_that("plot.dcvar_fit dispatches without error", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  expect_s3_class(plot(fit, type = "rho"), "ggplot")
  expect_s3_class(plot(fit, type = "phi"), "ggplot")
})

test_that("dcvar() fits exponential margins", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_exponential_fit()
  expect_s3_class(fit, "dcvar_fit")
  expect_equal(fit$margins, "exponential")
})
