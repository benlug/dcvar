test_that("simulate_dcVar with exponential margins (positive skew) works", {
  traj <- rho_constant(50, 0.5)
  sim <- simulate_dcVar(T = 50, rho_trajectory = traj,
                        margins = "exponential",
                        skew_direction = c(1, 1),
                        seed = 1)

  expect_type(sim, "list")
  expect_named(sim, c("Y", "Y_df", "true_params"))
  expect_equal(nrow(sim$Y), 50)
  expect_equal(ncol(sim$Y), 2)
  expect_equal(nrow(sim$Y_df), 50)
  expect_equal(sim$true_params$margins, "exponential")
  expect_equal(sim$true_params$skew_direction, c(1, 1))
})

test_that("simulate_dcVar with exponential margins (negative skew) works", {
  traj <- rho_constant(50, 0.5)
  sim <- simulate_dcVar(T = 50, rho_trajectory = traj,
                        margins = "exponential",
                        skew_direction = c(-1, -1),
                        seed = 2)

  expect_equal(nrow(sim$Y), 50)
  expect_equal(ncol(sim$Y), 2)
  expect_equal(sim$true_params$margins, "exponential")
  expect_equal(sim$true_params$skew_direction, c(-1, -1))
})

test_that("simulate_dcVar with gamma margins works", {
  traj <- rho_constant(50, 0.5)
  sim <- simulate_dcVar(T = 50, rho_trajectory = traj,
                        margins = "gamma",
                        skew_direction = c(1, 1),
                        skew_params = list(shape = 2),
                        seed = 3)

  expect_equal(nrow(sim$Y), 50)
  expect_equal(ncol(sim$Y), 2)
  expect_equal(sim$true_params$margins, "gamma")
  expect_equal(sim$true_params$skew_direction, c(1, 1))
  expect_equal(sim$true_params$skew_params$shape, 2)
})

test_that("simulate_dcVar with skew_normal margins works", {
  skip_if_not_installed("sn")

  traj <- rho_constant(50, 0.5)
  sim <- simulate_dcVar(T = 50, rho_trajectory = traj,
                        margins = "skew_normal",
                        skew_params = list(alpha = c(3, -3)),
                        seed = 4)

  expect_equal(nrow(sim$Y), 50)
  expect_equal(ncol(sim$Y), 2)
  expect_equal(sim$true_params$margins, "skew_normal")
  expect_equal(sim$true_params$skew_params$alpha, c(3, -3))
})

test_that("simulate_dcVar with non-normal margins has correct Y_df columns", {
  traj <- rho_constant(30, 0.3)
  sim <- simulate_dcVar(T = 30, rho_trajectory = traj,
                        margins = "exponential",
                        skew_direction = c(1, -1),
                        seed = 5)

  expect_true(all(c("time", "y1", "y2") %in% names(sim$Y_df)))
  expect_equal(sim$Y_df$time, 1:30)
})

test_that("simulate_dcVar with non-normal margins is reproducible with seed", {
  traj <- rho_constant(40, 0.5)
  s1 <- simulate_dcVar(T = 40, rho_trajectory = traj,
                       margins = "exponential",
                       skew_direction = c(1, 1), seed = 42)
  s2 <- simulate_dcVar(T = 40, rho_trajectory = traj,
                       margins = "exponential",
                       skew_direction = c(1, 1), seed = 42)
  expect_identical(s1$Y, s2$Y)
})

test_that("simulate_dcVar normal margins stores sigma_eps in true_params", {
  traj <- rho_constant(30, 0.5)
  sim <- simulate_dcVar(T = 30, rho_trajectory = traj,
                        margins = "normal", seed = 6)

  expect_equal(sim$true_params$margins, "normal")
  expect_equal(sim$true_params$sigma_eps, c(1, 1))
})
