test_that("loo() works for supported model types", {
  skip_if_no_rstan()

  loo_dc <- loo(get_dcvar_fit())
  loo_hmm <- loo(get_hmm_fit())
  loo_con <- loo(get_constant_fit())

  expect_s3_class(loo_dc, "loo")
  expect_s3_class(loo_hmm, "loo")
  expect_s3_class(loo_con, "loo")
})

test_that("loo() supports multilevel fits and rejects indicator SEM fits", {
  skip_if_no_rstan()

  # multilevel_copula_var.stan now stores per-observation log_lik like the EG
  # and mixed variants, so all multilevel fits are valid PSIS-LOO targets.
  loo_ml <- loo(get_multilevel_fit())
  expect_s3_class(loo_ml, "loo")

  expect_error(
    loo(get_sem_fit()),
    "not supported.*SEM|SEM.*not supported"
  )
})

test_that("dcvar_compare() works with named fits", {
  skip_if_no_rstan()

  result <- dcvar_compare(
    dcvar = get_dcvar_fit(),
    constant = get_constant_fit()
  )

  # loo_compare() returns a two-model comparison table. Whether that object
  # carries the base `matrix` class is loo-build-dependent (the dev build that
  # some CI runners resolve from r-universe does not always), so coerce before
  # asserting the structure rather than relying on is.matrix().
  comparison <- as.matrix(result)
  expect_equal(nrow(comparison), 2L)
  expect_true("elpd_diff" %in% colnames(comparison))
})

test_that("dcvar_compare() rejects unnamed arguments", {
  skip_if_no_rstan()

  expect_error(
    dcvar_compare(get_dcvar_fit(), get_constant_fit()),
    "must be named"
  )
})

test_that("dcvar_compare() rejects non-fit objects", {
  skip_if_no_rstan()

  expect_error(
    dcvar_compare(dcvar = get_dcvar_fit(), bad = "not a fit"),
    "not a dcvar model fit"
  )
})

test_that("dcvar_compare() rejects unsupported LOO targets", {
  skip_if_no_rstan()

  expect_error(
    dcvar_compare(dcvar = get_dcvar_fit(), sem = get_sem_fit()),
    "does not support one or more supplied fits"
  )
})

test_that("dcvar_compare() warns when mixing HMM and latent-conditioned fits", {
  skip_if_no_rstan()

  expect_warning(
    dcvar_compare(dcvar = get_dcvar_fit(), hmm = get_hmm_fit()),
    "mixes two"
  )
})

test_that("loo() fails clearly when custom Stan output omits log_lik", {
  fit <- structure(
    list(
      fit = posterior::as_draws_array(
        array(
          1:8,
          dim = c(4L, 1L, 2L),
          dimnames = list(NULL, NULL, c("mu[1]", "mu[2]"))
        )
      ),
      stan_data = list(n_time = 5, D = 2),
      model = "dcvar",
      vars = c("y1", "y2"),
      standardized = TRUE,
      margins = "normal",
      backend = "rstan",
      priors = list(),
      meta = list()
    ),
    class = c("dcvar_fit", "dcvar_model_fit")
  )

  expect_error(
    loo(fit),
    "Custom Stan files must preserve the expected parameter and generated-quantity names"
  )
})
