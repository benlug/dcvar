test_that("fitted() returns correct structure for dcvar", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
  fit_df <- fitted(fit)

  expect_s3_class(fit_df, "data.frame")
  expect_true(all(c("time", "y1", "y2") %in% names(fit_df)))
  expect_equal(nrow(fit_df), fit$stan_data$T - 1)
  expect_equal(fit_df$time, attr(fit$stan_data, "time_values")[-1])
})

test_that("fitted() returns correct structure for hmm", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  fit_df <- fitted(fit)

  expect_s3_class(fit_df, "data.frame")
  expect_equal(nrow(fit_df), fit$stan_data$T - 1)
})

test_that("fitted() returns correct structure for constant", {
  skip_if_no_rstan()

  fit <- get_constant_fit()
  fit_df <- fitted(fit)

  expect_s3_class(fit_df, "data.frame")
  expect_equal(nrow(fit_df), fit$stan_data$T - 1)
})

test_that("predict() returns correct structure", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
  pred_df <- predict(fit)

  expect_s3_class(pred_df, "data.frame")
  expect_true(all(c("time", "variable", "mean", "lower", "upper") %in% names(pred_df)))
  # Two variables, each with T-1 rows
  expect_equal(nrow(pred_df), 2 * (fit$stan_data$T - 1))
  expect_true(all(pred_df$lower <= pred_df$mean))
  expect_true(all(pred_df$upper >= pred_df$mean))
})

test_that("predict() respects ci_level parameter", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
  pred_95 <- predict(fit, ci_level = 0.95)
  pred_80 <- predict(fit, ci_level = 0.80)

  # 80% intervals should be narrower than 95%
  width_95 <- pred_95$upper - pred_95$lower
  width_80 <- pred_80$upper - pred_80$lower
  expect_true(all(width_80 < width_95))
})

test_that("predict() rejects invalid ci_level values", {
  fit <- structure(
    list(margins = "normal"),
    class = c("dcvar_fit", "dcvar_model_fit")
  )

  expect_error(predict(fit, ci_level = 0), "ci_level")
  expect_error(predict(fit, ci_level = 1), "ci_level")
  expect_error(predict(fit, ci_level = 2), "ci_level")
})

test_that("predict() works for hmm and constant", {
  skip_if_no_rstan()

  pred_hmm <- predict(get_hmm_fit())
  pred_con <- predict(get_constant_fit())

  expect_s3_class(pred_hmm, "data.frame")
  expect_s3_class(pred_con, "data.frame")
  expect_true(all(c("time", "variable", "mean", "lower", "upper") %in% names(pred_hmm)))
  expect_true(all(c("time", "variable", "mean", "lower", "upper") %in% names(pred_con)))
})

test_that("fitted() works and predict() errors for non-normal fits", {
  skip_if_no_rstan()

  fit <- get_dcvar_exponential_fit()
  fit_df <- fitted(fit)

  expect_s3_class(fit_df, "data.frame")
  expect_error(
    predict(fit),
    "only supported for normal margins"
  )
})

test_that("fitted() type='response' unstandardizes correctly", {
  skip_if_no_rstan()
  fit <- get_dcvar_fit()
  if (!isTRUE(fit$standardized)) skip("fit not standardized")

  fit_link <- fitted(fit, type = "link")
  fit_resp <- fitted(fit, type = "response")

  # Response should differ from link when standardized
  expect_false(identical(fit_link, fit_resp))

  # Response values should be on original data scale (larger magnitude)
  expect_true(any(abs(fit_resp[[fit$vars[1]]]) != abs(fit_link[[fit$vars[1]]])))
})

test_that("fitted() default type='link' is backward compatible", {
  skip_if_no_rstan()
  fit <- get_dcvar_fit()

  # Default should be same as explicit "link"
  fit_default <- fitted(fit)
  fit_link <- fitted(fit, type = "link")
  expect_equal(fit_default, fit_link)
})

test_that("predict() type='response' unstandardizes correctly", {
  skip_if_no_rstan()
  fit <- get_dcvar_fit()
  if (!isTRUE(fit$standardized)) skip("fit not standardized")

  pred_link <- predict(fit, type = "link")
  pred_resp <- predict(fit, type = "response")

  expect_false(identical(pred_link$mean, pred_resp$mean))
})

test_that("fitted() and predict() honor preserved time values", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
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

test_that("multilevel and SEM fitted()/predict() honor preserved time values", {
  skip_if_no_rstan()

  ml_fit <- get_multilevel_fit()
  attr(ml_fit$stan_data, "time_values") <- seq.Date(
    as.Date("2022-01-01"),
    by = "day",
    length.out = ml_fit$stan_data$T
  )
  ml_fitted <- fitted(ml_fit)
  ml_pred <- predict(ml_fit)
  expect_equal(ml_fitted$time, rep(attr(ml_fit$stan_data, "time_values")[-1], ml_fit$N))
  expect_equal(unique(ml_pred$time), attr(ml_fit$stan_data, "time_values")[-1])

  sem_fit <- get_sem_fit()
  attr(sem_fit$stan_data, "time_values") <- seq.Date(
    as.Date("2022-03-01"),
    by = "day",
    length.out = sem_fit$stan_data$T
  )
  sem_fitted <- fitted(sem_fit, type = "link")
  sem_pred <- predict(sem_fit, type = "response")
  expect_equal(sem_fitted$time, attr(sem_fit$stan_data, "time_values"))
  expect_equal(unique(sem_pred$time), attr(sem_fit$stan_data, "time_values"))
})

test_that("predictive methods fail clearly when required Stan outputs are missing", {
  make_stub_draws <- function(variables) {
    posterior::as_draws_array(
      array(
        seq_len(4L * length(variables)),
        dim = c(4L, 1L, length(variables)),
        dimnames = list(NULL, NULL, variables)
      )
    )
  }

  single_fit <- structure(
    list(
      fit = make_stub_draws(c("Phi[1,1]", "Phi[1,2]", "Phi[2,1]", "Phi[2,2]", "sigma_eps[1]", "sigma_eps[2]")),
      stan_data = list(T = 5, D = 2, Y = matrix(0, nrow = 5, ncol = 2)),
      model = "dcvar",
      vars = c("y1", "y2"),
      standardized = FALSE,
      margins = "normal",
      backend = "rstan",
      priors = list(),
      meta = list()
    ),
    class = c("dcvar_fit", "dcvar_model_fit")
  )
  expect_error(fitted(single_fit), "Custom Stan files must preserve")
  expect_error(predict(single_fit), "Custom Stan files must preserve")

  sem_fit <- structure(
    list(
      fit = make_stub_draws(c("rho", "Phi[1,1]", "Phi[1,2]", "Phi[2,1]", "Phi[2,2]")),
      stan_data = list(T = 4),
      model = "sem",
      vars = c("latent1", "latent2"),
      J = 2,
      lambda = c(0.8, 0.8),
      sigma_e = 0.3,
      indicators = list(latent1 = c("y1_1", "y1_2"), latent2 = c("y2_1", "y2_2")),
      backend = "rstan",
      priors = list(),
      meta = list()
    ),
    class = c("dcvar_sem_fit", "dcvar_model_fit")
  )
  expect_error(fitted(sem_fit, type = "link"), "Custom Stan files must preserve")
  expect_error(predict(sem_fit, type = "response"), "Custom Stan files must preserve")
})
