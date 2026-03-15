test_that("as.data.frame() works for dcVar_fit", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  df <- as.data.frame(fit)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("variable", "mean", "sd", "q2.5", "q97.5",
                     "rhat", "ess_bulk", "ess_tail") %in% names(df)))
  expect_true(nrow(df) > 0)
})

test_that("as.data.frame() works for dcVar_hmm_fit", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  df <- as.data.frame(fit)

  expect_s3_class(df, "data.frame")
  expect_true("variable" %in% names(df))
  expect_true(nrow(df) > 0)
})

test_that("as.data.frame() works for dcVar_constant_fit", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  df <- as.data.frame(fit)

  expect_s3_class(df, "data.frame")
  expect_true("variable" %in% names(df))
  expect_true(nrow(df) > 0)
})
