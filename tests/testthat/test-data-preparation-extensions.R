# ============================================================================
# Tests for prepare_multilevel_data() and prepare_sem_data()
# ============================================================================

# --- prepare_multilevel_data ------------------------------------------------

test_that("prepare_multilevel_data returns correct Stan data structure", {
  sim <- simulate_dcvar_multilevel(N = 5, T = 20, rho = 0.5, seed = 1)
  stan_data <- prepare_multilevel_data(sim$data, vars = c("y1", "y2"),
                                       id_var = "id")

  expect_type(stan_data, "list")
  expect_equal(stan_data$N, 5)
  expect_equal(stan_data$T, 20)
  expect_type(stan_data$y, "list")
  expect_equal(length(stan_data$y), 5)

  # Each element is a T x 2 matrix
  expect_equal(nrow(stan_data$y[[1]]), 20)
  expect_equal(ncol(stan_data$y[[1]]), 2)
})

test_that("prepare_multilevel_data includes prior hyperparameters", {
  sim <- simulate_dcvar_multilevel(N = 3, T = 15, rho = 0.5, seed = 2)
  stan_data <- prepare_multilevel_data(sim$data, vars = c("y1", "y2"),
                                       id_var = "id",
                                       prior_phi_bar_sd = 0.8,
                                       prior_tau_phi_scale = 0.3)

  expect_equal(stan_data$prior_phi_bar_sd, 0.8)
  expect_equal(stan_data$prior_tau_phi_scale, 0.3)
})

test_that("prepare_multilevel_data errors on unbalanced panels", {
  df <- data.frame(
    id = c(1, 1, 1, 2, 2),
    time = c(1, 2, 3, 1, 2),
    y1 = rnorm(5),
    y2 = rnorm(5)
  )

  expect_error(prepare_multilevel_data(df, vars = c("y1", "y2"), id_var = "id"),
               "Unbalanced")
})

test_that("prepare_multilevel_data errors on missing columns", {
  df <- data.frame(id = 1:10, time = 1:10, y1 = rnorm(10))
  expect_error(
    prepare_multilevel_data(df, vars = c("y1", "y2"), id_var = "id"),
    "not found"
  )
})

test_that("prepare_multilevel_data validates vars length", {
  sim <- simulate_dcvar_multilevel(N = 3, T = 10, rho = 0.5, seed = 5)
  expect_error(
    prepare_multilevel_data(sim$data, vars = "y1", id_var = "id"),
    "Exactly 2 variables"
  )
})

test_that("prepare_multilevel_data validates numeric columns", {
  sim <- simulate_dcvar_multilevel(N = 3, T = 10, rho = 0.5, seed = 6)
  sim$data$y1 <- as.character(sim$data$y1)
  expect_error(
    prepare_multilevel_data(sim$data, vars = c("y1", "y2"), id_var = "id"),
    "must be numeric"
  )
})

test_that("prepare_multilevel_data handles centering correctly", {
  sim <- simulate_dcvar_multilevel(N = 3, T = 20, rho = 0.5,
                                   center = FALSE, seed = 3)
  stan_centered <- prepare_multilevel_data(sim$data, vars = c("y1", "y2"),
                                           id_var = "id", center = TRUE)
  stan_uncentered <- prepare_multilevel_data(sim$data, vars = c("y1", "y2"),
                                             id_var = "id", center = FALSE)

  # Centered data person means should be near zero
  for (i in seq_len(stan_centered$N)) {
    expect_true(abs(mean(stan_centered$y[[i]][, 1])) < 1e-10)
    expect_true(abs(mean(stan_centered$y[[i]][, 2])) < 1e-10)
  }

  # Uncentered data should differ from centered
  expect_false(identical(stan_centered$y[[1]], stan_uncentered$y[[1]]))
})

test_that("prepare_multilevel_data stores attributes", {
  sim <- simulate_dcvar_multilevel(N = 4, T = 15, rho = 0.5, seed = 4)
  stan_data <- prepare_multilevel_data(sim$data, vars = c("y1", "y2"),
                                       id_var = "id")

  expect_equal(attr(stan_data, "vars"), c("y1", "y2"))
  expect_equal(length(attr(stan_data, "ids")), 4)
  expect_equal(nrow(attr(stan_data, "person_means")), 4)
  expect_equal(attr(stan_data, "time_values"), 1:15)
})

test_that("prepare_multilevel_data requires a shared regular time grid", {
  df <- data.frame(
    id = rep(c("A", "B"), each = 3),
    time = c(1, 2, 3, 2, 3, 4),
    y1 = rnorm(6),
    y2 = rnorm(6)
  )

  expect_error(
    prepare_multilevel_data(df, vars = c("y1", "y2"), id_var = "id"),
    "same.*time.*grid"
  )
})

test_that("prepare_multilevel_data rejects duplicate unit time values", {
  df <- data.frame(
    id = rep(c("A", "B"), each = 3),
    time = c(1, 2, 2, 1, 2, 3),
    y1 = rnorm(6),
    y2 = rnorm(6)
  )

  expect_error(
    prepare_multilevel_data(df, vars = c("y1", "y2"), id_var = "id"),
    "Duplicate"
  )
})


# --- prepare_sem_data -------------------------------------------------------

test_that("prepare_sem_data returns correct Stan data structure", {
  J <- 3
  sim <- simulate_dcvar_sem(T = 50, J = J,
                            lambda = rep(0.8, J),
                            rho = 0.5, seed = 10)
  indicators <- list(
    latent1 = paste0("y1_", seq_len(J)),
    latent2 = paste0("y2_", seq_len(J))
  )
  stan_data <- prepare_sem_data(sim$data, indicators = indicators,
                                J = J, lambda = rep(0.8, J),
                                sigma_e = sqrt(0.2))

  expect_type(stan_data, "list")
  expect_equal(stan_data$T, 50)
  expect_equal(stan_data$J, J)
  expect_equal(nrow(stan_data$y), 50)
  expect_equal(ncol(stan_data$y), 2 * J)
  expect_equal(stan_data$lambda, rep(0.8, J))
  expect_equal(stan_data$sigma_e, sqrt(0.2))
})

test_that("prepare_sem_data validates indicator columns exist", {
  df <- data.frame(time = 1:20, y1_1 = rnorm(20), y1_2 = rnorm(20))
  indicators <- list(
    latent1 = c("y1_1", "y1_2"),
    latent2 = c("y2_1", "y2_2")
  )

  expect_error(
    prepare_sem_data(df, indicators = indicators, J = 2,
                     lambda = c(0.8, 0.8), sigma_e = 0.5),
    "not found"
  )
})

test_that("prepare_sem_data validates numeric indicator columns", {
  df <- data.frame(time = 1:10, y1_1 = letters[1:10], y1_2 = rnorm(10),
                   y2_1 = rnorm(10), y2_2 = rnorm(10))
  indicators <- list(
    latent1 = c("y1_1", "y1_2"),
    latent2 = c("y2_1", "y2_2")
  )

  expect_error(
    prepare_sem_data(df, indicators = indicators, J = 2,
                     lambda = c(0.8, 0.8), sigma_e = 0.5),
    "must be numeric"
  )
})

test_that("prepare_sem_data rejects missing and non-finite values", {
  df <- data.frame(time = 1:10, y1_1 = rnorm(10), y1_2 = rnorm(10),
                   y2_1 = rnorm(10), y2_2 = rnorm(10))
  indicators <- list(
    latent1 = c("y1_1", "y1_2"),
    latent2 = c("y2_1", "y2_2")
  )

  df_na <- df
  df_na$y1_1[3] <- NA
  expect_error(
    prepare_sem_data(df_na, indicators = indicators, J = 2,
                     lambda = c(0.8, 0.8), sigma_e = 0.5),
    "missing values"
  )

  df_inf <- df
  df_inf$y1_1[3] <- Inf
  expect_error(
    prepare_sem_data(df_inf, indicators = indicators, J = 2,
                     lambda = c(0.8, 0.8), sigma_e = 0.5),
    "NaN or Inf"
  )
})

test_that("prepare_sem_data errors when indicators is not a list of length 2", {
  df <- data.frame(time = 1:20, y1 = rnorm(20), y2 = rnorm(20))

  expect_error(
    prepare_sem_data(df, indicators = c("y1", "y2"), J = 1,
                     lambda = 0.8, sigma_e = 0.5),
    "list of two"
  )
  expect_error(
    prepare_sem_data(df, indicators = list("y1"), J = 1,
                     lambda = 0.8, sigma_e = 0.5),
    "list of two"
  )
})

test_that("prepare_sem_data errors when lambda length mismatches J", {
  J <- 3
  sim <- simulate_dcvar_sem(T = 30, J = J, lambda = rep(0.8, J),
                            rho = 0.5, seed = 11)
  indicators <- list(
    latent1 = paste0("y1_", seq_len(J)),
    latent2 = paste0("y2_", seq_len(J))
  )

  expect_error(
    prepare_sem_data(sim$data, indicators = indicators, J = J,
                     lambda = c(0.8, 0.8), sigma_e = 0.5),
    "lambda"
  )
})

test_that("prepare_sem_data correct dimensions for indicator matrix", {
  for (J in c(2, 4)) {
    sim <- simulate_dcvar_sem(T = 40, J = J, lambda = rep(0.8, J),
                              rho = 0.5, seed = 200 + J)
    indicators <- list(
      latent1 = paste0("y1_", seq_len(J)),
      latent2 = paste0("y2_", seq_len(J))
    )
    stan_data <- prepare_sem_data(sim$data, indicators = indicators,
                                  J = J, lambda = rep(0.8, J),
                                  sigma_e = 0.4)

    expect_equal(nrow(stan_data$y), 40, info = paste("J =", J))
    expect_equal(ncol(stan_data$y), 2 * J, info = paste("J =", J))
    expect_equal(stan_data$J, J, info = paste("J =", J))
  }
})

test_that("prepare_sem_data stores latent variable name attributes", {
  J <- 2
  sim <- simulate_dcvar_sem(T = 30, J = J, lambda = rep(0.8, J),
                            rho = 0.5, seed = 12)
  indicators <- list(
    anxiety = paste0("y1_", seq_len(J)),
    depression = paste0("y2_", seq_len(J))
  )
  stan_data <- prepare_sem_data(sim$data, indicators = indicators,
                                J = J, lambda = rep(0.8, J),
                                sigma_e = 0.4)

  expect_equal(attr(stan_data, "vars"), c("anxiety", "depression"))
  expect_equal(attr(stan_data, "indicators"), indicators)
  expect_equal(attr(stan_data, "time_values"), 1:30)
})

test_that("prepare_sem_data rejects irregular or duplicate time values", {
  df <- data.frame(
    time = c(1, 3, 4),
    y1_1 = rnorm(3),
    y1_2 = rnorm(3),
    y2_1 = rnorm(3),
    y2_2 = rnorm(3)
  )
  indicators <- list(
    latent1 = c("y1_1", "y1_2"),
    latent2 = c("y2_1", "y2_2")
  )

  expect_error(
    prepare_sem_data(df, indicators = indicators, J = 2,
                     lambda = c(0.8, 0.8), sigma_e = 0.5),
    "not evenly spaced"
  )

  df$time <- c(1, 2, 2)
  expect_error(
    prepare_sem_data(df, indicators = indicators, J = 2,
                     lambda = c(0.8, 0.8), sigma_e = 0.5),
    "Duplicate"
  )
})
