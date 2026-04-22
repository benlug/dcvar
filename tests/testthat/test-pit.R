# ============================================================================
# Tests for PIT diagnostics
# ============================================================================

test_that("pit_values returns correct structure for constant fit", {
  skip_if_no_rstan()
  fit <- get_constant_fit()
  pit_df <- pit_values(fit)

  expect_s3_class(pit_df, "data.frame")
  expect_named(pit_df, c("time", "variable", "pit"))
  expect_equal(nrow(pit_df), (fit$stan_data$n_time - 1) * fit$stan_data$D)
  expect_true(all(pit_df$pit >= 0 & pit_df$pit <= 1))
  expect_equal(unique(pit_df$time), attr(fit$stan_data, "time_values")[-1])
})

test_that("pit_values returns correct structure for dcvar fit", {
  skip_if_no_rstan()
  fit <- get_dcvar_fit()
  pit_df <- pit_values(fit)

  expect_s3_class(pit_df, "data.frame")
  expect_named(pit_df, c("time", "variable", "pit"))
  expect_true(all(pit_df$pit >= 0 & pit_df$pit <= 1))
  expect_equal(unique(pit_df$time), attr(fit$stan_data, "time_values")[-1])
})

test_that("pit_values supports gamma margins", {
  skip_if_no_rstan()
  fit <- get_constant_gamma_fit()
  pit_df <- pit_values(fit)

  expect_s3_class(pit_df, "data.frame")
  expect_named(pit_df, c("time", "variable", "pit"))
  expect_true(all(pit_df$pit >= 0 & pit_df$pit <= 1))
})

test_that("pit_values supports skew-normal margins", {
  skip_if_no_rstan()
  skip_if_not_installed("sn")
  fit <- get_constant_skew_normal_fit()
  pit_df <- pit_values(fit)

  expect_s3_class(pit_df, "data.frame")
  expect_named(pit_df, c("time", "variable", "pit"))
  expect_true(all(pit_df$pit >= 0 & pit_df$pit <= 1))
})

test_that("pit_test returns KS test results", {
  skip_if_no_rstan()
  fit <- get_constant_fit()
  ks_df <- pit_test(fit)

  expect_s3_class(ks_df, "data.frame")
  expect_named(ks_df, c("variable", "ks_statistic", "p_value", "n"))
  expect_equal(nrow(ks_df), 2)
  expect_true(all(ks_df$ks_statistic >= 0))
  expect_true(all(ks_df$p_value >= 0 & ks_df$p_value <= 1))
})

test_that("plot_pit returns a ggplot object", {
  skip_if_no_rstan()
  fit <- get_constant_fit()
  p <- plot_pit(fit)

  expect_s3_class(p, "ggplot")
})

test_that("pit plot dispatch works via plot(fit, type = 'pit')", {
  skip_if_no_rstan()
  fit <- get_constant_fit()
  p <- plot(fit, type = "pit")

  expect_s3_class(p, "ggplot")
})

test_that("pit_values.default errors for unsupported objects", {
  expect_error(pit_values(list()), class = "rlang_error")
})

test_that("pit_test.default errors for unsupported objects", {
  expect_error(pit_test(list()), class = "rlang_error")
})

test_that("pit_values honor preserved time values", {
  skip_if_no_rstan()

  fit <- get_constant_fit()
  attr(fit$stan_data, "time_values") <- seq.Date(
    as.Date("2022-01-01"),
    by = "week",
    length.out = fit$stan_data$n_time
  )
  pit_df <- pit_values(fit)

  expect_equal(unique(pit_df$time), attr(fit$stan_data, "time_values")[-1])
})
