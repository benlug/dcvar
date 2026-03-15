test_that("simulate_breakpoint_data() works for single breakpoint", {
  sim <- simulate_breakpoint_data(T = 50, type = "single", seed = 42)

  expect_equal(nrow(sim$Y), 50)
  expect_equal(ncol(sim$Y), 2)
  expect_equal(length(sim$true_params$rho), 49)
  expect_true(all(c("time", "y1", "y2") %in% names(sim$Y_df)))
})

test_that("simulate_breakpoint_data() works for double breakpoint", {
  sim <- simulate_breakpoint_data(T = 50, type = "double", seed = 42)

  expect_equal(nrow(sim$Y), 50)
  expect_equal(length(sim$true_params$rho), 49)
})

test_that("simulate_breakpoint_data() is reproducible with seed", {
  s1 <- simulate_breakpoint_data(T = 50, type = "single", seed = 123)
  s2 <- simulate_breakpoint_data(T = 50, type = "single", seed = 123)
  expect_identical(s1$Y, s2$Y)
})

test_that("simulate_breakpoint_data() respects custom parameters", {
  sim <- simulate_breakpoint_data(
    T = 50, type = "single",
    rho_before = 0.9, rho_after = 0.1,
    mu = c(1, -1),
    sigma_eps = c(0.5, 0.5),
    seed = 42
  )

  expect_equal(sim$true_params$mu, c(1, -1))
  expect_equal(sim$true_params$sigma_eps, c(0.5, 0.5))
})
