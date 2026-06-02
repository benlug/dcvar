test_that("simulate_dcvar with exponential margins (positive skew) works", {
  traj <- rho_constant(50, 0.5)
  sim <- simulate_dcvar(n_time = 50, rho_trajectory = traj,
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

test_that("simulate_dcvar with exponential margins (negative skew) works", {
  traj <- rho_constant(50, 0.5)
  sim <- simulate_dcvar(n_time = 50, rho_trajectory = traj,
                        margins = "exponential",
                        skew_direction = c(-1, -1),
                        seed = 2)

  expect_equal(nrow(sim$Y), 50)
  expect_equal(ncol(sim$Y), 2)
  expect_equal(sim$true_params$margins, "exponential")
  expect_equal(sim$true_params$skew_direction, c(-1, -1))
})

test_that("simulate_dcvar with gamma margins works", {
  traj <- rho_constant(50, 0.5)
  sim <- simulate_dcvar(n_time = 50, rho_trajectory = traj,
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

test_that("simulate_dcvar with skew_normal margins works", {
  skip_if_not_installed("sn")

  traj <- rho_constant(50, 0.5)
  sim <- simulate_dcvar(n_time = 50, rho_trajectory = traj,
                        margins = "skew_normal",
                        skew_params = list(alpha = c(3, -3)),
                        seed = 4)

  expect_equal(nrow(sim$Y), 50)
  expect_equal(ncol(sim$Y), 2)
  expect_equal(sim$true_params$margins, "skew_normal")
  expect_equal(sim$true_params$skew_params$alpha, c(3, -3))
})

test_that("simulate_dcvar with non-normal margins has correct Y_df columns", {
  traj <- rho_constant(30, 0.3)
  sim <- simulate_dcvar(n_time = 30, rho_trajectory = traj,
                        margins = "exponential",
                        skew_direction = c(1, -1),
                        seed = 5)

  expect_true(all(c("time", "y1", "y2") %in% names(sim$Y_df)))
  expect_equal(sim$Y_df$time, 1:30)
})

test_that("simulate_dcvar with non-normal margins is reproducible with seed", {
  traj <- rho_constant(40, 0.5)
  s1 <- simulate_dcvar(n_time = 40, rho_trajectory = traj,
                       margins = "exponential",
                       skew_direction = c(1, 1), seed = 42)
  s2 <- simulate_dcvar(n_time = 40, rho_trajectory = traj,
                       margins = "exponential",
                       skew_direction = c(1, 1), seed = 42)
  expect_identical(s1$Y, s2$Y)
})

test_that("simulate_dcvar normal margins stores sigma_eps in true_params", {
  traj <- rho_constant(30, 0.5)
  sim <- simulate_dcvar(n_time = 30, rho_trajectory = traj,
                        margins = "normal", seed = 6)

  expect_equal(sim$true_params$margins, "normal")
  expect_equal(sim$true_params$sigma_eps, c(1, 1))
})

test_that("simulate_dcvar supports mixed (per-variable) margins", {
  traj <- rho_constant(40, 0.5)
  sim <- simulate_dcvar(n_time = 40, rho_trajectory = traj,
                        margins = c("normal", "exponential"),
                        skew_direction = c(1, 1), seed = 7)

  expect_equal(ncol(sim$Y), 2)
  expect_equal(sim$true_params$margins, c("normal", "exponential"))
  expect_equal(sim$true_params$sigma_eps, c(1, 1))
  expect_equal(sim$true_params$skew_direction, c(1, 1))
})

test_that("simulate_dcvar mixed margins are reproducible with a seed", {
  traj <- rho_constant(40, 0.4)
  s1 <- simulate_dcvar(n_time = 40, rho_trajectory = traj,
                       margins = c("normal", "exponential"),
                       skew_direction = c(1, -1), seed = 99)
  s2 <- simulate_dcvar(n_time = 40, rho_trajectory = traj,
                       margins = c("normal", "exponential"),
                       skew_direction = c(1, -1), seed = 99)
  expect_identical(s1$Y, s2$Y)
})

test_that("simulate_dcvar collapses identical margin vectors to scalar", {
  traj <- rho_constant(30, 0.5)
  sim <- simulate_dcvar(n_time = 30, rho_trajectory = traj,
                        margins = c("normal", "normal"), seed = 8)
  expect_identical(sim$true_params$margins, "normal")

  # An all-identical vector must reproduce the scalar simulation exactly.
  sim_scalar <- simulate_dcvar(n_time = 30, rho_trajectory = traj,
                               margins = "normal", seed = 8)
  expect_identical(sim$Y, sim_scalar$Y)
})

test_that("mixed-margin simulation feeds the single-level prep functions", {
  # Mixed margins are supported by the constant, dynamic, and HMM prep paths;
  # each builds a per-dimension family array.
  traj <- rho_constant(20, 0.5)
  sim <- simulate_dcvar(n_time = 20, rho_trajectory = traj,
                        margins = c("normal", "exponential"),
                        skew_direction = c(1, 1), seed = 9)
  for (prep in list(
    prepare_constant_data(sim$Y_df, vars = c("y1", "y2"),
                          margins = c("normal", "exponential"),
                          skew_direction = c(1, 1)),
    prepare_dcvar_data(sim$Y_df, vars = c("y1", "y2"),
                       margins = c("normal", "exponential"),
                       skew_direction = c(1, 1)),
    prepare_hmm_data(sim$Y_df, vars = c("y1", "y2"), K = 2,
                     margins = c("normal", "exponential"),
                     skew_direction = c(1, 1))
  )) {
    expect_equal(prep$family, c(1L, 2L))
  }

  # Multilevel still requires a single margin family.
  dfm <- data.frame(id = rep(1:2, each = 10), time = rep(1:10, 2),
                    y1 = rnorm(20), y2 = rnorm(20))
  expect_error(
    prepare_multilevel_data(dfm, vars = c("y1", "y2"), id_var = "id",
                            margins = c("normal", "exponential"),
                            skew_direction = c(1, 1)),
    "not supported by"
  )
})
