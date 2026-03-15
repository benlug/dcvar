test_that("loo() works for supported model types", {
  skip_if_no_cmdstanr()

  loo_dc <- loo(get_dcvar_fit())
  loo_hmm <- loo(get_hmm_fit())
  loo_con <- loo(get_constant_fit())

  expect_s3_class(loo_dc, "loo")
  expect_s3_class(loo_hmm, "loo")
  expect_s3_class(loo_con, "loo")
})

test_that("loo() rejects multilevel and SEM fits with informative errors", {
  skip_if_no_cmdstanr()

  expect_error(
    loo(get_multilevel_fit()),
    "not supported.*multilevel|multilevel.*not supported"
  )
  expect_error(
    loo(get_sem_fit()),
    "not supported.*SEM|SEM.*not supported"
  )
})

test_that("dcvar_compare() works with named fits", {
  skip_if_no_cmdstanr()

  result <- dcvar_compare(
    dcvar = get_dcvar_fit(),
    constant = get_constant_fit()
  )

  expect_true(is.matrix(result))
  expect_true(nrow(result) == 2)
  expect_true("elpd_diff" %in% colnames(result))
})

test_that("dcvar_compare() rejects unnamed arguments", {
  skip_if_no_cmdstanr()

  expect_error(
    dcvar_compare(get_dcvar_fit(), get_constant_fit()),
    "must be named"
  )
})

test_that("dcvar_compare() rejects non-fit objects", {
  skip_if_no_cmdstanr()

  expect_error(
    dcvar_compare(dcvar = get_dcvar_fit(), bad = "not a fit"),
    "not a dcvar model fit"
  )
})

test_that("dcvar_compare() rejects unsupported LOO targets", {
  skip_if_no_cmdstanr()

  expect_error(
    dcvar_compare(dcvar = get_dcvar_fit(), multilevel = get_multilevel_fit()),
    "does not support SEM or multilevel"
  )
  expect_error(
    dcvar_compare(dcvar = get_dcvar_fit(), sem = get_sem_fit()),
    "does not support SEM or multilevel"
  )
})
