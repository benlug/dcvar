test_that("rho_constant returns correct length and value", {
  result <- rho_constant(100, rho = 0.5)
  expect_equal(length(result), 99)
  expect_true(all(result == 0.5))
})

test_that("rho_constant rejects out-of-range rho", {
  expect_error(rho_constant(100, rho = 1.5))
  expect_error(rho_constant(100, rho = -1.5))
})

test_that("rho_decreasing returns correct length and trend", {
  result <- rho_decreasing(100)
  expect_equal(length(result), 99)
  expect_true(result[1] > result[99])
})

test_that("rho_increasing returns correct length and trend", {
  result <- rho_increasing(100)
  expect_equal(length(result), 99)
  expect_true(result[1] < result[99])
})

test_that("rho_random_walk returns correct length", {
  result <- rho_random_walk(100, seed = 42)
  expect_equal(length(result), 99)
  expect_true(all(result > -1 & result < 1))
})

test_that("rho_random_walk is reproducible with seed", {
  r1 <- rho_random_walk(50, seed = 123)
  r2 <- rho_random_walk(50, seed = 123)
  expect_identical(r1, r2)
})

test_that("rho_step returns correct length", {
  result <- rho_step(100)
  expect_equal(length(result), 99)
})

test_that("rho_step has correct levels before and after breakpoint", {
  result <- rho_step(101, rho_before = 0.8, rho_after = 0.2, breakpoint = 0.5)
  expect_equal(result[1], 0.8)
  expect_equal(result[99], 0.2)
})

test_that("rho_step validates rho bounds and breakpoint arguments", {
  expect_error(rho_step(100, rho_before = 1.2), "\\[-1, 1\\]")
  expect_error(rho_step(100, rho_after = -1.2), "\\[-1, 1\\]")
  expect_error(rho_step(100, breakpoint = -0.1), "breakpoint")
  expect_error(rho_step(100, breakpoint = 1.4), "breakpoint")
  expect_error(rho_step(100, transition_width = -1), "transition_width")
})

test_that("rho_double_step returns correct length", {
  result <- rho_double_step(100)
  expect_equal(length(result), 99)
})

test_that("rho_double_step has three phases", {
  result <- rho_double_step(100, rho_levels = c(0.8, 0.2, 0.8))
  expect_equal(result[1], 0.8)
  expect_equal(result[50], 0.2)
  expect_equal(result[99], 0.8)
})

test_that("rho_double_step validates levels and breakpoints", {
  expect_error(rho_double_step(100, rho_levels = c(0.8, 1.2, 0.7)), "\\[-1, 1\\]")
  expect_error(rho_double_step(100, breakpoints = c(0.7, 0.2)), "strictly increasing")
  expect_error(rho_double_step(100, breakpoints = c(-0.1, 0.2)), "breakpoints")
  expect_error(rho_double_step(100, transition_width = -0.5), "transition_width")
})

test_that("rho_scenario works for all named scenarios", {
  scenarios <- c("constant", "decreasing", "increasing", "random_walk",
                 "single_middle", "large_change", "small_change",
                 "increase", "double_relapse")

  for (s in scenarios) {
    result <- rho_scenario(s, T = 50)
    expect_equal(length(result), 49, info = paste("scenario:", s))
    expect_true(all(result >= -1 & result <= 1), info = paste("scenario:", s))
  }
})

test_that("rho_scenario forwards overrides to scenario generators", {
  early <- rho_scenario("single_middle", T = 10, breakpoint = 0.1)
  late <- rho_scenario("single_middle", T = 10, breakpoint = 0.9)

  expect_false(identical(early, late))
})

test_that("rho_scenario errors for unknown scenario", {
  expect_error(rho_scenario("nonexistent", T = 50), "Unknown scenario")
})
