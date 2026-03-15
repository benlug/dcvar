test_that("interpret_rho_trajectory() works for dcvar", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  msg <- expect_output(interpret_rho_trajectory(fit), "DC-VAR")
  expect_type(msg, "character")
  expect_true(nchar(msg) > 0)
})

test_that("interpret_rho_trajectory() works for hmm", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  msg <- expect_output(interpret_rho_trajectory(fit), "HMM")
  expect_type(msg, "character")
  expect_true(grepl("HMM states", msg))
})

test_that("interpret_rho_trajectory() works for constant", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  msg <- expect_output(interpret_rho_trajectory(fit), "Constant")
  expect_type(msg, "character")
  expect_true(grepl("time-invariant", msg))
})

test_that("interpret_rho_trajectory() respects threshold", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  msg1 <- capture.output(interpret_rho_trajectory(fit, threshold = 0.01))
  msg2 <- capture.output(interpret_rho_trajectory(fit, threshold = 10))
  # With very high threshold, everything should be "stable"
  expect_true(any(grepl("stable", msg2)))
})
