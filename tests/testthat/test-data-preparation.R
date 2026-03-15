test_that("prepare_dcvar_data works with valid input", {
  df <- data.frame(time = 1:50, y1 = rnorm(50), y2 = rnorm(50))
  result <- prepare_dcvar_data(df, vars = c("y1", "y2"))

  expect_equal(result$T, 50)
  expect_equal(result$D, 2)
  expect_equal(nrow(result$Y), 50)
  expect_equal(ncol(result$Y), 2)
  expect_true(attr(result, "standardized"))
  expect_equal(attr(result, "time_values"), 1:50)
})

test_that("prepare_dcvar_data rejects non-bivariate input", {
  df <- data.frame(time = 1:50, y1 = rnorm(50), y2 = rnorm(50), y3 = rnorm(50))
  expect_error(prepare_dcvar_data(df, vars = "y1"), "Exactly 2")
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2", "y3")), "Exactly 2")
})

test_that("prepare_dcvar_data rejects missing columns", {
  df <- data.frame(time = 1:50, y1 = rnorm(50))
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2")), "not found")
})

test_that("prepare_dcvar_data errors on interior missing values by default", {
  df <- data.frame(time = 1:50, y1 = rnorm(50), y2 = rnorm(50))
  df$y1[25] <- NA
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2")), "interior missing")
})

test_that("prepare_dcvar_data can drop interior missing values when allow_gaps = TRUE", {
  df <- data.frame(time = 1:50, y1 = rnorm(50), y2 = rnorm(50))
  df$y1[25] <- NA

  expect_warning(
    out <- prepare_dcvar_data(df, vars = c("y1", "y2"), allow_gaps = TRUE),
    "Removing 1 row"
  )
  expect_equal(out$T, 49)
  expect_equal(attr(out, "time_values"), c(1:24, 26:50))
})

test_that("prepare_dcvar_data warns on edge missing values", {
  df <- data.frame(time = 1:50, y1 = rnorm(50), y2 = rnorm(50))
  df$y1[1] <- NA

  expect_warning(
    out <- prepare_dcvar_data(df, vars = c("y1", "y2")),
    "edges only"
  )
  expect_equal(out$T, 49)
  expect_equal(attr(out, "time_values"), 2:50)
})

test_that("prepare_dcvar_data preserves standardization metadata", {
  df <- data.frame(time = 1:50, y1 = rnorm(50, 5, 2), y2 = rnorm(50, -3, 0.5))
  result <- prepare_dcvar_data(df, vars = c("y1", "y2"), standardize = TRUE)

  expect_true(!is.null(attr(result, "Y_means")))
  expect_true(!is.null(attr(result, "Y_sds")))
  expect_equal(length(attr(result, "Y_means")), 2)
})

test_that("prepare_dcvar_data skips standardization when FALSE", {
  df <- data.frame(time = 1:50, y1 = rnorm(50, 5, 2), y2 = rnorm(50))
  result <- prepare_dcvar_data(df, vars = c("y1", "y2"), standardize = FALSE)

  expect_false(attr(result, "standardized"))
  expect_null(attr(result, "Y_means"))
})

test_that("prepare_dcvar_data sets correct prior hyperparameters", {
  df <- data.frame(time = 1:50, y1 = rnorm(50), y2 = rnorm(50))
  result <- prepare_dcvar_data(df, vars = c("y1", "y2"),
                               prior_mu_sd = 3, prior_sigma_omega_rate = 0.2)

  expect_equal(result$sigma_mu_prior, 3)
  expect_equal(result$sigma_omega_prior, 0.2)
})

test_that("prepare_hmm_data adds HMM-specific fields", {
  df <- data.frame(time = 1:50, y1 = rnorm(50), y2 = rnorm(50))
  result <- prepare_hmm_data(df, vars = c("y1", "y2"), K = 3)

  expect_equal(result$K, 3L)
  expect_equal(result$kappa, 10)
  expect_equal(result$alpha_off, 1)
  expect_equal(result$z_rho_prior_sd, 1.0)
})

test_that("prepare_hmm_data validates K consistently with dcvar_hmm", {
  df <- data.frame(time = 1:50, y1 = rnorm(50), y2 = rnorm(50))
  expect_error(prepare_hmm_data(df, vars = c("y1", "y2"), K = 2.7), "integer >= 2")
})

test_that("prepare_constant_data adds copula prior", {
 df <- data.frame(time = 1:50, y1 = rnorm(50), y2 = rnorm(50))
  result <- prepare_constant_data(df, vars = c("y1", "y2"))

  expect_equal(result$z_rho_prior_sd, 1.0)
  expect_null(result$sigma_omega_prior)
})

test_that("data is sorted by time_var", {
  df <- data.frame(time = c(3, 1, 2), y1 = c(30, 10, 20), y2 = c(33, 11, 22))
  result <- prepare_dcvar_data(df, vars = c("y1", "y2"), standardize = FALSE)

  expect_equal(result$Y[1, 1], 10)
  expect_equal(result$Y[2, 1], 20)
  expect_equal(result$Y[3, 1], 30)
  expect_equal(attr(result, "time_values"), c(1, 2, 3))
})

test_that("prepare_dcvar_data rejects duplicate time values", {
  df_duplicate <- data.frame(time = c(1, 2, 2), y1 = rnorm(3), y2 = rnorm(3))
  expect_error(
    prepare_dcvar_data(df_duplicate, vars = c("y1", "y2")),
    "Duplicate"
  )
})

test_that("prepare_dcvar_data handles irregular time values consistently", {
  df_irregular <- data.frame(time = c(1, 3, 4), y1 = rnorm(3), y2 = rnorm(3))
  expect_error(
    prepare_dcvar_data(df_irregular, vars = c("y1", "y2")),
    "not evenly spaced"
  )

  expect_warning(
    out <- prepare_dcvar_data(
      df_irregular,
      vars = c("y1", "y2"),
      allow_gaps = TRUE
    ),
    "not evenly spaced"
  )
  expect_equal(attr(out, "time_values"), c(1, 3, 4))
})
