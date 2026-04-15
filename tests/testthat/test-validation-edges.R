# ============================================================================
# Tests for edge-case validation behaviors
# ============================================================================

# --- NA handling in prepare_dcvar_data ---------------------------------------

test_that("interior NA errors by default", {
  df <- data.frame(time = 1:10, y1 = c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10), y2 = rnorm(10))
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2")), "interior|adjacency|missing")
})

test_that("interior NA warns when allow_gaps = TRUE", {
  df <- data.frame(time = 1:10, y1 = c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10), y2 = rnorm(10))
  expect_warning(prepare_dcvar_data(df, vars = c("y1", "y2"), allow_gaps = TRUE), "missing")
})

test_that("edge-only NA still just warns", {
  df <- data.frame(time = 1:10, y1 = c(NA, 2, 3, 4, 5, 6, 7, 8, 9, 10), y2 = rnorm(10))
  expect_warning(prepare_dcvar_data(df, vars = c("y1", "y2")), "missing")
})

test_that("all-NA data errors", {
  df <- data.frame(time = 1:5, y1 = rep(NA_real_, 5), y2 = rep(NA_real_, 5))
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2")), "missing|observations|NA")
})

# --- Minimum observation requirements ---------------------------------------

test_that("T=2 errors (too few observations)", {
  df <- data.frame(time = 1:2, y1 = c(1, 2), y2 = c(3, 4))
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2")), "3 observations")
})

test_that("T=3 succeeds (minimum viable)", {
  df <- data.frame(time = 1:3, y1 = c(1, 2, 3), y2 = c(4, 5, 6))
  result <- prepare_dcvar_data(df, vars = c("y1", "y2"))
  expect_equal(result$T, 3)
})

# --- Data type validation ---------------------------------------------------

test_that("NaN or Inf values error", {
  df_nan <- data.frame(time = 1:5, y1 = c(1, NaN, 3, 4, 5), y2 = rnorm(5))
  expect_error(prepare_dcvar_data(df_nan, vars = c("y1", "y2")), "NaN|Inf")

  df_inf <- data.frame(time = 1:5, y1 = c(1, Inf, 3, 4, 5), y2 = rnorm(5))
  expect_error(prepare_dcvar_data(df_inf, vars = c("y1", "y2")), "NaN|Inf")
})

test_that("non-numeric column errors", {
  df <- data.frame(time = 1:5, y1 = letters[1:5], y2 = rnorm(5))
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2")), "numeric")
})

test_that("zero-variance column errors on standardize", {
  df <- data.frame(time = 1:10, y1 = rep(5, 10), y2 = rnorm(10))
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2")), "zero|variance")
})

test_that("zero-variance column after NA removal errors", {
  df <- data.frame(time = 1:10, y1 = c(5, 5, NA, 5, 5, 5, 5, 5, 5, 5), y2 = rnorm(10))
  expect_warning(
    expect_error(prepare_dcvar_data(df, vars = c("y1", "y2"),
                                    allow_gaps = TRUE), "zero|variance"),
    "Removing 1 row"
  )
})

# --- Margin validation ------------------------------------------------------

test_that("exponential margins require skew_direction", {
  df <- data.frame(time = 1:20, y1 = abs(rnorm(20)), y2 = abs(rnorm(20)))
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2"),
                                  margins = "exponential"), "skew_direction")
})

test_that("gamma margins require skew_direction", {
  df <- data.frame(time = 1:20, y1 = abs(rnorm(20)), y2 = abs(rnorm(20)))
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2"),
                                  margins = "gamma"), "skew_direction")
})

test_that("invalid margin type errors", {
  df <- data.frame(time = 1:20, y1 = rnorm(20), y2 = rnorm(20))
  expect_error(prepare_dcvar_data(df, vars = c("y1", "y2"),
                                  margins = "invalid_margin"), "margin")
})

test_that("normal margins work without skew_direction", {
  df <- data.frame(time = 1:20, y1 = rnorm(20), y2 = rnorm(20))
  result <- prepare_dcvar_data(df, vars = c("y1", "y2"), margins = "normal")
  expect_equal(result$T, 20)
})

test_that("prepare_dcvar_data validates scalar logical flags", {
  df <- data.frame(time = 1:20, y1 = rnorm(20), y2 = rnorm(20))
  expect_error(
    prepare_dcvar_data(df, vars = c("y1", "y2"), standardize = 1),
    "single logical"
  )
  expect_error(
    prepare_dcvar_data(df, vars = c("y1", "y2"), allow_gaps = c(TRUE, FALSE)),
    "single logical"
  )
  expect_error(
    prepare_dcvar_data(df, vars = c("y1", "y2"), standardize = NA),
    "single logical"
  )
})

# --- HMM K validation -------------------------------------------------------

test_that("K must be integer-like", {
  skip_if_no_rstan()
  sim <- simulate_dcvar(T = 30, rho_trajectory = rho_step(30), seed = 1)
  expect_error(dcvar_hmm(sim$Y_df, vars = c("y1", "y2"), K = 2.5), "integer")
})

test_that("K < 2 errors", {
  skip_if_no_rstan()
  sim <- simulate_dcvar(T = 30, rho_trajectory = rho_step(30), seed = 1)
  expect_error(dcvar_hmm(sim$Y_df, vars = c("y1", "y2"), K = 1), "integer.*>= 2")
})

test_that("non-numeric K errors", {
  skip_if_no_rstan()
  sim <- simulate_dcvar(T = 30, rho_trajectory = rho_step(30), seed = 1)
  expect_error(dcvar_hmm(sim$Y_df, vars = c("y1", "y2"), K = "two"), "integer")
})

# --- Multilevel validation ---------------------------------------------------

test_that("multilevel errors on missing columns", {
  df <- data.frame(time = 1:10, y1 = rnorm(10), y2 = rnorm(10), id = rep(1, 10))
  expect_error(prepare_multilevel_data(df, vars = c("y1", "missing_var"), id_var = "id"),
               "not found")
})

test_that("multilevel errors on unbalanced panels", {
  df <- data.frame(
    time = c(1:5, 1:3),
    y1 = rnorm(8), y2 = rnorm(8),
    id = c(rep(1, 5), rep(2, 3))
  )
  expect_error(prepare_multilevel_data(df, vars = c("y1", "y2"), id_var = "id"),
               "Unbalanced|balanced")
})

test_that("multilevel errors on non-dataframe input", {
  expect_error(prepare_multilevel_data(list(a = 1), vars = c("y1", "y2"), id_var = "id"),
               "data frame")
})

test_that("multilevel rejects missing unit ids", {
  df <- data.frame(
    time = c(1, 2, 1, 2),
    y1 = rnorm(4),
    y2 = rnorm(4),
    id = c("A", NA, "B", "B")
  )
  expect_error(
    prepare_multilevel_data(df, vars = c("y1", "y2"), id_var = "id"),
    "missing values"
  )
})

test_that("multilevel drops unused factor levels in ids", {
  df <- data.frame(
    time = c(1, 2, 1, 2),
    y1 = rnorm(4),
    y2 = rnorm(4),
    id = factor(c("A", "A", "B", "B"), levels = c("A", "B", "C"))
  )
  out <- prepare_multilevel_data(df, vars = c("y1", "y2"), id_var = "id")
  expect_equal(out$N, 2)
  expect_equal(as.character(attr(out, "ids")), c("A", "B"))
  expect_equal(levels(attr(out, "ids")), c("A", "B"))
})

test_that("multilevel errors on NAs in unit data", {
  df <- data.frame(
    time = rep(1:5, 2),
    y1 = c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10),
    y2 = rnorm(10),
    id = rep(1:2, each = 5)
  )
  expect_error(prepare_multilevel_data(df, vars = c("y1", "y2"), id_var = "id"),
               "missing|NA")
})

test_that("multilevel errors on Inf values", {
  df <- data.frame(
    time = rep(1:5, 2),
    y1 = c(1, 2, Inf, 4, 5, 6, 7, 8, 9, 10),
    y2 = rnorm(10),
    id = rep(1:2, each = 5)
  )
  expect_error(prepare_multilevel_data(df, vars = c("y1", "y2"), id_var = "id"),
               "NaN|Inf")
})

# --- SEM validation ----------------------------------------------------------

test_that("SEM errors on wrong indicator count", {
  skip_if_no_rstan()
  expect_error(
    dcvar_sem(data.frame(), indicators = list(a = "x1", b = c("x2", "x3")),
              J = 2, lambda = c(1, 1), sigma_e = 1),
    "indicator"
  )
})

test_that("SEM errors on mismatched lambda length", {
  df <- data.frame(time = 1:10, y1_1 = rnorm(10), y1_2 = rnorm(10),
                   y2_1 = rnorm(10), y2_2 = rnorm(10))
  expect_error(
    prepare_sem_data(df, indicators = list(a = c("y1_1", "y1_2"), b = c("y2_1", "y2_2")),
                     J = 2, lambda = c(1, 1, 1), sigma_e = 1),
    "lambda.*length|length.*2"
  )
})

test_that("J must be a positive integer", {
  skip_if_no_rstan()
  expect_error(dcvar_sem(data.frame(), indicators = list(), J = 0,
                         lambda = numeric(), sigma_e = 1), "positive integer|J")
  expect_error(dcvar_sem(data.frame(), indicators = list(), J = -1,
                         lambda = numeric(), sigma_e = 1), "positive integer|J")
})

# --- Probs validation -------------------------------------------------------

test_that("probs must be in [0, 1]", {
  skip_if_no_rstan()
  fit <- get_dcvar_fit()
  expect_error(rho_trajectory(fit, probs = c(-0.5, 0.5)), "probs")
  expect_error(rho_trajectory(fit, probs = c(0.5, 1.5)), "probs")
})
