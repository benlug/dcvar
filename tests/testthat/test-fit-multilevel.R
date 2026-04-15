# ============================================================================
# Tests for dcvar_multilevel() fit and methods
# ============================================================================

# --- Fit object structure ----------------------------------------------------

test_that("dcvar_multilevel returns correct class", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  expect_s3_class(fit, "dcvar_multilevel_fit")
  expect_s3_class(fit, "dcvar_model_fit")
})

test_that("dcvar_multilevel rejects center = FALSE with bundled Stan model", {
  df <- data.frame(
    id = rep(1:2, each = 4),
    time = rep(1:4, times = 2),
    y1 = rnorm(8),
    y2 = rnorm(8)
  )

  expect_error(
    dcvar_multilevel(df, vars = c("y1", "y2"), center = FALSE),
    "not supported by the bundled multilevel model"
  )
})

test_that("dcvar_multilevel rejects explicit bundled stan_file with center = FALSE", {
  df <- data.frame(
    id = rep(1:2, each = 4),
    time = rep(1:4, times = 2),
    y1 = rnorm(8),
    y2 = rnorm(8)
  )

  expect_error(
    dcvar_multilevel(
      df,
      vars = c("y1", "y2"),
      center = FALSE,
      stan_file = dcvar_stan_path("multilevel")
    ),
    "not supported by the bundled multilevel model"
  )
})

test_that("multilevel coef() returns expected structure", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  co <- coef(fit)

  expect_type(co, "list")
  expect_named(co, c("phi_bar", "tau_phi", "sigma", "rho"))
  expect_true(length(co$phi_bar) == 4)
  expect_true(length(co$tau_phi) == 4)
  expect_true(length(co$sigma) == 2)
  expect_true(length(co$rho) == 1)
})

test_that("multilevel summary() returns correct class", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  s <- summary(fit)

  expect_s3_class(s, "dcvar_multilevel_summary")
  expect_true(all(c("var_params", "random_effects", "diagnostics") %in% names(s)))
})

test_that("multilevel print() works", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  expect_output(print(fit), "multilevel|Multilevel")
})

# --- Extraction functions ----------------------------------------------------

test_that("multilevel rho_trajectory() returns expected structure", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  rho_df <- rho_trajectory(fit)

  expect_s3_class(rho_df, "data.frame")
  expect_true(all(c("time", "mean", "sd") %in% names(rho_df)))
  # Constant rho: all means should be identical
  expect_equal(length(unique(rho_df$mean)), 1)
})

test_that("multilevel var_params() returns expected structure", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  vp <- var_params(fit)

  expect_type(vp, "list")
  expect_named(vp, c("phi_bar", "tau_phi", "sigma", "rho"))
  expect_s3_class(vp$phi_bar, "data.frame")
  expect_true(all(c("variable", "mean", "sd", "q2.5", "q97.5") %in% names(vp$phi_bar)))
})

test_that("multilevel random_effects() returns expected structure", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  re <- random_effects(fit)

  expect_s3_class(re, "data.frame")
  expect_true(all(c("unit", "parameter", "mean", "sd", "q2.5", "q97.5") %in% names(re)))
  expect_true(all(re$parameter %in% c("phi11", "phi12", "phi21", "phi22")))
  expect_equal(nrow(re), fit$N * 4)
})

# --- Plotting ----------------------------------------------------------------

test_that("multilevel plot random_effects returns ggplot", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  p <- plot(fit, type = "random_effects")
  expect_s3_class(p, "ggplot")
})

test_that("multilevel plot diagnostics returns ggplot", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  p <- plot(fit, type = "diagnostics")
  expect_s3_class(p, "ggplot")
})

# --- Fitted values and prediction --------------------------------------------

test_that("multilevel fitted() returns unit-specific trajectories", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()

  fit_link <- fitted(fit, type = "link")
  fit_resp <- fitted(fit, type = "response")

  expect_s3_class(fit_link, "data.frame")
  expect_named(fit_link, c("unit", "time", "y1", "y2"))
  expect_equal(nrow(fit_link), fit$N * (fit$stan_data$T - 1))
  expect_equal(length(unique(fit_link$unit)), fit$N)
  expect_false(identical(fit_link[, fit$vars], fit_resp[, fit$vars]))
})

test_that("multilevel predict() returns unit-specific intervals", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()

  pred <- predict(fit)

  expect_s3_class(pred, "data.frame")
  expect_named(pred, c("unit", "time", "variable", "mean", "lower", "upper"))
  expect_equal(nrow(pred), fit$N * (fit$stan_data$T - 1) * length(fit$vars))
  expect_equal(sort(unique(pred$variable)), sort(fit$vars))
  expect_true(all(pred$lower <= pred$mean))
  expect_true(all(pred$upper >= pred$mean))
})

test_that("multilevel pit_values() errors informatively", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  expect_error(pit_values(fit), "not.*supported")
})

# --- Diagnostics -------------------------------------------------------------

test_that("multilevel diagnostics are finite", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  diag <- dcvar_diagnostics(fit)

  expect_lte(diag$n_divergent, 1)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.25)
  expect_true(diag$min_ess_bulk > 8)
  expect_true(diag$min_ess_tail > 8)
})

test_that("multilevel fit cache emits only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_multilevel_fit_warnings(), "multilevel")
})
