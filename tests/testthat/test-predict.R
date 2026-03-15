test_that("fitted() returns correct structure for dcVar", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  fit_df <- fitted(fit)

  expect_s3_class(fit_df, "data.frame")
  expect_true(all(c("time", "y1", "y2") %in% names(fit_df)))
  expect_equal(nrow(fit_df), fit$stan_data$T - 1)
  expect_equal(fit_df$time, attr(fit$stan_data, "time_values")[-1])
})

test_that("fitted() returns correct structure for hmm", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  fit_df <- fitted(fit)

  expect_s3_class(fit_df, "data.frame")
  expect_equal(nrow(fit_df), fit$stan_data$T - 1)
})

test_that("fitted() returns correct structure for constant", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  fit_df <- fitted(fit)

  expect_s3_class(fit_df, "data.frame")
  expect_equal(nrow(fit_df), fit$stan_data$T - 1)
})

test_that("predict() returns correct structure", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  pred_df <- predict(fit)

  expect_s3_class(pred_df, "data.frame")
  expect_true(all(c("time", "variable", "mean", "lower", "upper") %in% names(pred_df)))
  # Two variables, each with T-1 rows
  expect_equal(nrow(pred_df), 2 * (fit$stan_data$T - 1))
  expect_true(all(pred_df$lower <= pred_df$mean))
  expect_true(all(pred_df$upper >= pred_df$mean))
})

test_that("predict() respects ci_level parameter", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  pred_95 <- predict(fit, ci_level = 0.95)
  pred_80 <- predict(fit, ci_level = 0.80)

  # 80% intervals should be narrower than 95%
  width_95 <- pred_95$upper - pred_95$lower
  width_80 <- pred_80$upper - pred_80$lower
  expect_true(all(width_80 < width_95))
})

test_that("predict() works for hmm and constant", {
  skip_if_no_cmdstanr()

  pred_hmm <- predict(get_hmm_fit())
  pred_con <- predict(get_constant_fit())

  expect_s3_class(pred_hmm, "data.frame")
  expect_s3_class(pred_con, "data.frame")
  expect_true(all(c("time", "variable", "mean", "lower", "upper") %in% names(pred_hmm)))
  expect_true(all(c("time", "variable", "mean", "lower", "upper") %in% names(pred_con)))
})

test_that("fitted() works and predict() errors for non-normal fits", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_exponential_fit()
  fit_df <- fitted(fit)

  expect_s3_class(fit_df, "data.frame")
  expect_error(
    predict(fit),
    "only supported for normal margins"
  )
})

test_that("fitted() type='response' unstandardizes correctly", {
  skip_if_no_cmdstanr()
  fit <- get_dcVar_fit()
  if (!isTRUE(fit$standardized)) skip("fit not standardized")

  fit_link <- fitted(fit, type = "link")
  fit_resp <- fitted(fit, type = "response")

  # Response should differ from link when standardized
  expect_false(identical(fit_link, fit_resp))

  # Response values should be on original data scale (larger magnitude)
  expect_true(any(abs(fit_resp[[fit$vars[1]]]) != abs(fit_link[[fit$vars[1]]])))
})

test_that("fitted() default type='link' is backward compatible", {
  skip_if_no_cmdstanr()
  fit <- get_dcVar_fit()

  # Default should be same as explicit "link"
  fit_default <- fitted(fit)
  fit_link <- fitted(fit, type = "link")
  expect_equal(fit_default, fit_link)
})

test_that("predict() type='response' unstandardizes correctly", {
  skip_if_no_cmdstanr()
  fit <- get_dcVar_fit()
  if (!isTRUE(fit$standardized)) skip("fit not standardized")

  pred_link <- predict(fit, type = "link")
  pred_resp <- predict(fit, type = "response")

  expect_false(identical(pred_link$mean, pred_resp$mean))
})

test_that("fitted() and predict() honor preserved time values", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  attr(fit$stan_data, "time_values") <- seq.Date(
    as.Date("2021-01-01"),
    by = "day",
    length.out = fit$stan_data$T
  )

  fit_df <- fitted(fit)
  pred_df <- predict(fit)

  expect_equal(fit_df$time, attr(fit$stan_data, "time_values")[-1])
  expect_equal(unique(pred_df$time), attr(fit$stan_data, "time_values")[-1])
})
