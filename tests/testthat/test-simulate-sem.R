test_that("simulate_dcvar_sem returns correct structure", {
  sim <- simulate_dcvar_sem(T = 50, rho = 0.5, seed = 1)

  expect_type(sim, "list")
  expect_named(sim, c("data", "true_params", "latent_states", "innovations"))
  expect_s3_class(sim$data, "data.frame")
})

test_that("simulate_dcvar_sem data has correct dimensions", {
  J <- 3
  T_val <- 50
  sim <- simulate_dcvar_sem(T = T_val, J = J, rho = 0.5, seed = 2)

  # Data: T rows, 1 time column + 2*J indicator columns
  expect_equal(nrow(sim$data), T_val)
  expect_equal(ncol(sim$data), 1 + 2 * J)
  expect_true("time" %in% names(sim$data))
})

test_that("simulate_dcvar_sem indicator column names are correct", {
  J <- 3
  sim <- simulate_dcvar_sem(T = 40, J = J, rho = 0.3, seed = 3)

  expected_cols <- c("time",
                     paste0("y1_", seq_len(J)),
                     paste0("y2_", seq_len(J)))
  expect_equal(names(sim$data), expected_cols)
})

test_that("simulate_dcvar_sem latent_states has correct dimensions", {
  T_val <- 60
  sim <- simulate_dcvar_sem(T = T_val, rho = 0.5, seed = 4)

  expect_equal(nrow(sim$latent_states), T_val)
  expect_equal(ncol(sim$latent_states), 2)
})

test_that("simulate_dcvar_sem innovations has correct dimensions", {
  T_val <- 60
  sim <- simulate_dcvar_sem(T = T_val, rho = 0.5, seed = 5)

  expect_equal(nrow(sim$innovations), T_val)
  expect_equal(ncol(sim$innovations), 2)
})

test_that("simulate_dcvar_sem true_params contains expected elements", {
  J <- 3
  sim <- simulate_dcvar_sem(T = 50, J = J, rho = 0.4, seed = 6)

  expect_named(sim$true_params,
               c("Phi", "mu", "sigma", "rho", "lambda", "sigma_e", "J"))
  expect_equal(sim$true_params$rho, 0.4)
  expect_equal(sim$true_params$J, J)
  expect_equal(length(sim$true_params$lambda), J)
  expect_equal(dim(sim$true_params$Phi), c(2, 2))
})

test_that("simulate_dcvar_sem errors when lambda length mismatches J", {
  expect_error(
    simulate_dcvar_sem(T = 50, J = 3, lambda = c(0.8, 0.8), rho = 0.5),
    "lambda"
  )
  expect_error(
    simulate_dcvar_sem(T = 50, J = 2, lambda = c(0.8, 0.8, 0.8), rho = 0.5),
    "lambda"
  )
})

test_that("simulate_dcvar_sem validates lambda finiteness and sigma_e", {
  expect_error(
    simulate_dcvar_sem(T = 50, J = 2, lambda = c(0.8, NA_real_), rho = 0.5),
    "lambda"
  )
  expect_error(
    simulate_dcvar_sem(T = 50, J = 2, lambda = c(0.8, 0.8), sigma_e = 0, rho = 0.5),
    "sigma_e"
  )
})

test_that("simulate_dcvar_sem validates Phi and sigma shapes", {
  expect_error(
    simulate_dcvar_sem(T = 20, Phi = matrix(1:3, nrow = 1), rho = 0.5),
    "Phi"
  )
  expect_error(
    simulate_dcvar_sem(T = 20, sigma = c(1, 1, 1), rho = 0.5),
    "sigma"
  )
})

test_that("simulate_dcvar_sem works with different J values", {
  for (J in c(2, 4, 5)) {
    sim <- simulate_dcvar_sem(T = 40, J = J,
                              lambda = rep(0.8, J),
                              rho = 0.5, seed = 100 + J)
    expect_equal(ncol(sim$data), 1 + 2 * J, info = paste("J =", J))
    expect_equal(nrow(sim$data), 40, info = paste("J =", J))
    expect_equal(length(sim$true_params$lambda), J, info = paste("J =", J))
  }
})

test_that("simulate_dcvar_sem is reproducible with seed", {
  s1 <- simulate_dcvar_sem(T = 30, rho = 0.5, seed = 42)
  s2 <- simulate_dcvar_sem(T = 30, rho = 0.5, seed = 42)
  expect_identical(s1$data, s2$data)
  expect_identical(s1$latent_states, s2$latent_states)
})

test_that("simulate_dcvar_sem warns that burnin is ignored", {
  expect_warning(
    simulate_dcvar_sem(T = 20, burnin = 10, seed = 11),
    "ignored"
  )
})
