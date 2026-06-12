test_that("compute_rho_metrics returns all expected fields", {
  n <- 50
  rho_true <- seq(0.7, 0.3, length.out = n)
  rho_est <- rho_true + rnorm(n, 0, 0.01)
  rho_lower <- rho_est - 0.1
  rho_upper <- rho_est + 0.1

  result <- compute_rho_metrics(rho_true, rho_est, rho_lower, rho_upper)

  expect_true(is.list(result))
  expected_names <- c("bias", "relative_bias", "coverage", "interval_width",
                      "correlation", "bias_start", "bias_end",
                      "coverage_start", "coverage_end")
  expect_true(all(expected_names %in% names(result)))
})

test_that("compute_rho_metrics rejects malformed inputs", {
  expect_error(
    compute_rho_metrics(1:3, 1:2, 1:3, 1:3),
    "same length"
  )
  expect_error(
    compute_rho_metrics(c(1, NA_real_), c(1, 2), c(0, 0), c(2, 2)),
    "finite"
  )
  expect_error(
    compute_rho_metrics(c(1, 2), c(1, 2), c(0.8, 0.9), c(0.7, 1.0)),
    "must not exceed"
  )
})

test_that("compute_rho_metrics has 100% coverage when CI is wide", {
  n <- 20
  rho_true <- rep(0.5, n)
  rho_est <- rep(0.5, n)
  rho_lower <- rep(-1, n)
  rho_upper <- rep(1, n)

  result <- compute_rho_metrics(rho_true, rho_est, rho_lower, rho_upper)
  expect_equal(result$coverage, 1)
})

test_that("compute_rho_metrics has 0% coverage when CI misses", {
  n <- 20
  rho_true <- rep(0.5, n)
  rho_est <- rep(-0.5, n)
  rho_lower <- rep(-0.9, n)
  rho_upper <- rep(-0.1, n)

  result <- compute_rho_metrics(rho_true, rho_est, rho_lower, rho_upper)
  expect_equal(result$coverage, 0)
})

test_that("compute_param_metrics works for scalars", {
  result <- compute_param_metrics(
    true_value = 0.3,
    est_mean = 0.32,
    est_lower = 0.25,
    est_upper = 0.39
  )

  expect_equal(result$bias, 0.02)
  expect_equal(result$covered, 1)
  expect_equal(result$interval_width, 0.14)
})

test_that("compute_param_metrics detects non-coverage", {
  result <- compute_param_metrics(
    true_value = 0.5,
    est_mean = 0.1,
    est_lower = 0.05,
    est_upper = 0.15
  )
  expect_equal(result$covered, 0)
})

test_that("compute_param_metrics rejects malformed inputs", {
  expect_error(
    compute_param_metrics(true_value = 1:2, est_mean = 0.1, est_lower = 0, est_upper = 1),
    "single finite numeric"
  )
  expect_error(
    compute_param_metrics(true_value = 1, est_mean = 0.1, est_lower = 0.5, est_upper = 0.1),
    "must not exceed"
  )
})

test_that("aggregate_metrics returns correct structure", {
  rho_true <- seq(0.7, 0.3, length.out = 10)
  metrics_list <- lapply(1:5, function(i) {
    rho_est <- rho_true + rnorm(length(rho_true), 0, 0.01)
    rho_lower <- rho_est - 0.1
    rho_upper <- rho_est + 0.1
    list(rho = compute_rho_metrics(rho_true, rho_est, rho_lower, rho_upper))
  })

  result <- aggregate_metrics(metrics_list)

  expect_equal(result$n_reps, 5)
  expect_true(is.data.frame(result$rho))
  expect_true("metric" %in% names(result$rho))
  expect_true("mean" %in% names(result$rho))
  expect_equal(nrow(result$rho), 9)
  expect_equal(
    result$rho$metric,
    c(
      "bias", "relative_bias", "coverage", "interval_width",
      "correlation", "bias_start", "bias_end",
      "coverage_start", "coverage_end"
    )
  )
})

test_that("aggregate_metrics rejects empty input", {
  expect_error(aggregate_metrics(list()), "at least one replication")
})

test_that("aggregate_metrics rejects malformed rho elements", {
  bad_metrics <- list(
    list(rho = list(
      bias = 0,
      relative_bias = 0,
      coverage = 1,
      interval_width = 0.1,
      correlation = 1,
      bias_start = 0,
      bias_end = 0,
      coverage_start = 1
    ))
  )

  expect_error(
    aggregate_metrics(bad_metrics),
    "full set"
  )
})

test_that(".relative_bias preserves the sign convention and handles edge cases", {
  # all-positive truth: ratio-of-means
  expect_equal(.relative_bias(rep(0.6, 5), rep(0.5, 5)), 20)
  # all-negative truth: attenuation toward zero stays negative (conventional)
  expect_equal(.relative_bias(-0.4, -0.5), -20)
  # mixed-sign truth: magnitude-normalized, no explosion at zero crossings
  rb <- .relative_bias(c(-0.4, 0, 0.6), c(-0.5, 0, 0.5))
  expect_true(is.finite(rb))
  # all-zero truth: NA with a warning
  expect_warning(rb0 <- .relative_bias(c(0.1, -0.1), c(0, 0)), "undefined")
  expect_true(is.na(rb0))
})

test_that("aggregate_metrics tolerates NA relative_bias from null conditions", {
  m <- suppressWarnings(compute_rho_metrics(
    rho_true = rep(0, 10),
    rho_est = rnorm(10, 0, 0.01),
    rho_lower = rep(-0.2, 10),
    rho_upper = rep(0.2, 10)
  ))
  expect_true(is.na(m$relative_bias))
  agg <- aggregate_metrics(list(list(rho = m), list(rho = m)))
  expect_s3_class(agg$rho, "data.frame")
  expect_true(is.nan(agg$rho$mean[agg$rho$metric == "relative_bias"]) ||
                is.na(agg$rho$mean[agg$rho$metric == "relative_bias"]))
})
