test_that("simulate_dcvar returns correct dimensions", {
  sim <- simulate_dcvar(T = 50, rho_trajectory = rho_constant(50, 0.5))

  expect_equal(nrow(sim$Y), 50)
  expect_equal(ncol(sim$Y), 2)
  expect_equal(nrow(sim$Y_df), 50)
  expect_true(all(c("time", "y1", "y2") %in% names(sim$Y_df)))
})

test_that("simulate_dcvar rejects wrong rho_trajectory length", {
  expect_error(
    simulate_dcvar(T = 50, rho_trajectory = rep(0.5, 10)),
    "length T-1"
  )
})

test_that("simulate_dcvar validates T >= 2", {
  expect_error(
    simulate_dcvar(T = 1, rho_trajectory = numeric(0)),
    "integer >= 2"
  )
})

test_that("simulate_dcvar validates sigma_eps length", {
  expect_error(
    simulate_dcvar(T = 50, rho_trajectory = rho_constant(50, 0.5),
                   sigma_eps = c(1, 1, 1)),
    "sigma_eps"
  )
})

test_that("simulate_dcvar validates gamma shape", {
  expect_error(
    simulate_dcvar(
      T = 50,
      rho_trajectory = rho_constant(50, 0.5),
      margins = "gamma",
      skew_direction = c(1, 1),
      skew_params = list(shape = 0)
    ),
    "shape"
  )
})

test_that("simulate_dcvar is reproducible with seed", {
  traj <- rho_constant(50, 0.5)
  s1 <- simulate_dcvar(T = 50, rho_trajectory = traj, seed = 42)
  s2 <- simulate_dcvar(T = 50, rho_trajectory = traj, seed = 42)
  expect_identical(s1$Y, s2$Y)
})

test_that("simulate_dcvar stores true parameters", {
  traj <- rho_decreasing(50)
  sim <- simulate_dcvar(T = 50, rho_trajectory = traj)

  expect_identical(sim$true_params$rho, traj)
  expect_equal(sim$true_params$mu, c(0, 0))
  expect_equal(dim(sim$true_params$Phi), c(2, 2))
  expect_equal(sim$true_params$sigma_eps, c(1, 1))
})

test_that("simulate_dcvar stores default non-normal margin parameters in true_params", {
  traj <- rho_constant(30, 0.2)

  gamma_sim <- simulate_dcvar(
    T = 30,
    rho_trajectory = traj,
    margins = "gamma",
    skew_direction = c(1, 1)
  )
  expect_equal(gamma_sim$true_params$skew_params$shape, 1)

  skew_sim <- simulate_dcvar(
    T = 30,
    rho_trajectory = traj,
    margins = "skew_normal"
  )
  expect_equal(skew_sim$true_params$skew_params$alpha, c(0, 0))
})

test_that("simulate_dcvar Y_df time column is 1:T", {
  sim <- simulate_dcvar(T = 30, rho_trajectory = rho_constant(30, 0.3))
  expect_equal(sim$Y_df$time, 1:30)
})
