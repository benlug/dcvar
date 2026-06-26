# Logic tests for applicability_check(). The two posterior accessors
# (var_params, dcvar_diagnostics) are mocked so the verdict logic can be driven
# deterministically without a Stan fit; rstan-gated smoke tests at the bottom
# exercise the real accessors end to end.

.mock_var_params <- function(object, ...) {
  D <- ncol(object$stan_data$Y)
  out <- list(
    Phi = data.frame(
      variable = sprintf("Phi[%d,%d]", seq_len(D), seq_len(D)),
      mean = object$.mock_phi,
      stringsAsFactors = FALSE
    )
  )
  if (!is.null(object$.mock_delta)) {
    k <- length(object$.mock_delta)
    out$delta <- data.frame(
      variable = sprintf("delta[%d]", seq_len(k)),
      mean = object$.mock_delta,
      stringsAsFactors = FALSE
    )
  }
  out
}

.mock_diagnostics <- function(object, ...) object$.mock_diag

mock_constant_fit <- function(Y, margins = "normal", phi_self,
                              delta = NULL, n_divergent = 0L, max_rhat = 1.0,
                              vars = NULL) {
  D <- ncol(Y)
  structure(
    list(
      stan_data = list(Y = Y),
      vars = if (is.null(vars)) paste0("v", seq_len(D)) else vars,
      margins = margins,
      backend = "rstan",
      .mock_phi = phi_self,
      .mock_delta = delta,
      .mock_diag = list(n_divergent = n_divergent, max_rhat = max_rhat)
    ),
    class = c("dcvar_constant_fit", "dcvar_model_fit")
  )
}

# A clean AR(1) series with no ties at the bounds.
ar1_series <- function(n, phi, sd = 0.3, seed = 1) {
  set.seed(seed)
  y <- numeric(n)
  for (t in 2:n) y[t] <- phi * y[t - 1] + stats::rnorm(1, sd = sd)
  y
}

# An AR(1) series with a floor pile-up (a genuine boundary atom).
ar1_with_floor <- function(n, phi, frac = 0.15, seed = 2) {
  y <- ar1_series(n, phi, seed = seed)
  thr <- stats::quantile(y, frac)
  y[y < thr] <- min(y)
  y
}

test_that("clean continuous series with converged sampler is suitable", {
  Y <- cbind(ar1_series(200, 0.6, seed = 1), ar1_series(200, 0.6, seed = 2))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(mock_constant_fit(Y, "normal", phi_self = c(0.58, 0.6)))
  expect_s3_class(res, "dcvar_applicability")
  expect_identical(res$verdict, "suitable")
  expect_length(res$reasons, 0L)
})

test_that("flexible margin with atom and collapse is unsuitable", {
  Y <- cbind(ar1_with_floor(200, 0.6), ar1_series(200, 0.6, seed = 3))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(
    mock_constant_fit(Y, "skew_normal", phi_self = c(0.02, 0.03), delta = c(0.99, 0.1))
  )
  expect_identical(res$verdict, "unsuitable")
  expect_true(any(grepl("boundary atom", res$reasons)))
  expect_true(any(grepl("collapse", res$reasons)))
  expect_true(any(grepl("slant at its boundary", res$reasons)))
})

test_that("convergence problems block a 'suitable' verdict (no self-contradiction)", {
  Y <- cbind(ar1_series(200, 0.6, seed = 1), ar1_series(200, 0.6, seed = 2))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(
    mock_constant_fit(Y, "normal", phi_self = c(0.58, 0.6),
                      n_divergent = 12L, max_rhat = 1.2)
  )
  expect_false(res$verdict == "suitable")
  expect_true(any(grepl("divergent", res$reasons)))
  expect_true(any(grepl("Rhat", res$reasons)))
})

test_that("collapse is detected for negative autoregression (two-sided anchor)", {
  Y <- cbind(ar1_series(200, -0.55, seed = 4), ar1_series(200, -0.55, seed = 5))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(
    mock_constant_fit(Y, "exponential", phi_self = c(0.01, 0.02))
  )
  expect_true(all(res$ols_self_lag < -0.2))
  expect_true(any(grepl("collapse", res$reasons)))
})

test_that("a lone extreme value on a short series is not counted as an atom", {
  set.seed(99)
  Y <- cbind(stats::rnorm(15), stats::rnorm(15))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(mock_constant_fit(Y, "normal", phi_self = c(0.1, 0.1)))
  expect_equal(unname(res$atom), c(0, 0))
  expect_identical(res$verdict, "suitable")
})

test_that("a genuine tie at the bound is counted as an atom", {
  set.seed(7)
  col <- stats::rnorm(40)
  col <- col - min(col) + 1     # shift so the series has a unique minimum at 1
  col[c(1, 2, 3)] <- 0          # introduce a clean three-way tie at a new floor
  Y <- cbind(col, stats::rnorm(40))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(mock_constant_fit(Y, "normal", phi_self = c(0.1, 0.1)))
  expect_equal(unname(res$atom[1]), 3 / 40)
})

test_that("undefined OLS anchor blocks a suitable verdict", {
  Y <- cbind(c(0.1, 0.3, -0.2), c(0.4, -0.1, 0.2))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(mock_constant_fit(Y, "normal", phi_self = c(0.5, 0.5)))
  expect_identical(res$verdict, "caution")
  expect_true(any(is.na(res$ols_self_lag)))
  expect_true(any(grepl("OLS VAR\\(1\\) anchor self-lag undefined", res$reasons)))
})

test_that("delta rail is checked for mixed margins and named by the skew dimension", {
  Y <- cbind(ar1_with_floor(200, 0.6), ar1_series(200, 0.6, seed = 8))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(
    mock_constant_fit(Y, c("skew_normal", "gamma"),
                      phi_self = c(0.02, 0.55), delta = c(0.99))
  )
  expect_named(res$delta, "v1")
  expect_true(any(grepl("slant at its boundary", res$reasons)))
})

test_that("non-finite delta summary produces caution instead of an error", {
  Y <- cbind(ar1_series(200, 0.6, seed = 1), ar1_series(200, 0.6, seed = 2))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(
    mock_constant_fit(Y, "skew_normal", phi_self = c(0.58, 0.6),
                      delta = c(NA_real_, 0.1))
  )
  expect_identical(res$verdict, "caution")
  expect_true(any(grepl("slant summary is non-finite", res$reasons)))
})

test_that("reference fit participates in collapse detection and is reported", {
  Y <- cbind(ar1_series(200, 0.6, seed = 1), ar1_series(200, 0.6, seed = 2))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  reference <- mock_constant_fit(Y, "normal", phi_self = c(0.6, 0.6))
  fit <- mock_constant_fit(Y, "skew_normal", phi_self = c(0.01, 0.02), delta = c(0.2, 0.2))
  res <- applicability_check(fit, reference = reference)
  expect_false(is.null(res$reference_self_lag))
  expect_equal(unname(res$reference_self_lag), c(0.6, 0.6))
  expect_true(any(grepl("collapse", res$reasons)))
})

test_that("reference is validated for class and dimensionality", {
  Y <- cbind(ar1_series(60, 0.5, seed = 1), ar1_series(60, 0.5, seed = 2))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  fit <- mock_constant_fit(Y, "normal", phi_self = c(0.5, 0.5))
  expect_error(applicability_check(fit, reference = list(1)), "reference")
  bad_dim <- mock_constant_fit(cbind(ar1_series(60, 0.5, seed = 3)), "normal", phi_self = 0.5)
  expect_error(applicability_check(fit, reference = bad_dim), "same number of variables")
  bad_n <- mock_constant_fit(
    cbind(ar1_series(50, 0.5, seed = 3), ar1_series(50, 0.5, seed = 4)),
    "normal",
    phi_self = c(0.5, 0.5)
  )
  expect_error(applicability_check(fit, reference = bad_n), "same number of observations")
  bad_vars <- mock_constant_fit(Y, "normal", phi_self = c(0.5, 0.5), vars = c("x", "y"))
  expect_error(applicability_check(fit, reference = bad_vars), "same variables")
  bad_margin <- mock_constant_fit(Y, "gamma", phi_self = c(0.5, 0.5))
  expect_error(applicability_check(fit, reference = bad_margin), "normal-margin")
})

test_that("ols_clear_tol is a tunable threshold and is reported", {
  Y <- cbind(ar1_series(200, 0.6, seed = 1), ar1_series(200, 0.6, seed = 2))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(mock_constant_fit(Y, "normal", phi_self = c(0.58, 0.6)))
  expect_true("ols_clear_tol" %in% names(res$thresholds))
  expect_equal(res$thresholds$ols_clear_tol, 0.20)
})

test_that("missing values in the data matrix abort with an informative error", {
  Y <- cbind(ar1_series(50, 0.5, seed = 1), ar1_series(50, 0.5, seed = 2))
  Y[5, 1] <- NA
  expect_error(
    applicability_check(mock_constant_fit(Y, "normal", phi_self = c(0.5, 0.5))),
    "missing values"
  )
})

test_that("non-constant-fit input is rejected", {
  expect_error(
    applicability_check(structure(list(), class = "dcvar_hmm_fit")),
    "dcvar_constant_fit"
  )
})

test_that("print method returns its input invisibly", {
  Y <- cbind(ar1_series(120, 0.6, seed = 1), ar1_series(120, 0.6, seed = 2))
  local_mocked_bindings(
    var_params = .mock_var_params, dcvar_diagnostics = .mock_diagnostics,
    .package = "dcvar"
  )
  res <- applicability_check(mock_constant_fit(Y, "normal", phi_self = c(0.58, 0.6)))
  expect_output(out <- withVisible(print(res)), "SUITABLE")
  expect_false(out$visible)
  expect_identical(out$value, res)
})

test_that("applicability_check runs on a real constant skew-normal fit", {
  skip_if_no_rstan()

  fit <- get_constant_skew_normal_fit()
  res <- applicability_check(fit)
  expect_s3_class(res, "dcvar_applicability")
  expect_true(res$verdict %in% c("suitable", "caution", "unsuitable"))
  expect_length(res$fit_self_lag, ncol(as.matrix(fit$stan_data$Y)))
  expect_output(print(res), "Applicability check")
})
