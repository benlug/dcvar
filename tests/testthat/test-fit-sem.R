# ============================================================================
# Tests for dcVar_sem() fit and methods
# ============================================================================

# --- Fit object structure ----------------------------------------------------

test_that("dcVar_sem returns correct class", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  expect_s3_class(fit, "dcVar_sem_fit")
  expect_s3_class(fit, "dcVar_model_fit")
})

test_that("SEM coef() returns expected structure", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  co <- coef(fit)

  expect_type(co, "list")
  expect_named(co, c("mu", "Phi", "sigma", "rho"))
  expect_true(length(co$mu) == 2)
  expect_true(length(co$Phi) == 4)
  expect_true(length(co$sigma) == 2)
  expect_true(length(co$rho) == 1)
})

test_that("SEM summary() returns correct class", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  s <- summary(fit)

  expect_s3_class(s, "dcVar_sem_summary")
  expect_true(all(c("var_params", "diagnostics") %in% names(s)))
})

test_that("SEM print() works", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  expect_output(print(fit), "SEM|sem")
})

# --- Extraction functions ----------------------------------------------------

test_that("SEM rho_trajectory() returns expected structure", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  rho_df <- rho_trajectory(fit)

  expect_s3_class(rho_df, "data.frame")
  expect_true(all(c("time", "mean", "sd") %in% names(rho_df)))
  # SEM has constant rho: all means should be identical
  expect_equal(length(unique(rho_df$mean)), 1)
  # rho should be bounded within (-0.97, 0.97)
  expect_true(all(rho_df$mean > -1 & rho_df$mean < 1))
})

test_that("SEM latent_states() returns expected structure", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  ls <- latent_states(fit)

  expect_s3_class(ls, "data.frame")
  expect_true(all(c("time", "variable", "mean", "sd") %in% names(ls)))
  expect_true(all(fit$vars %in% ls$variable))
  # Should have T rows per variable (2 variables)
  expect_equal(nrow(ls), fit$stan_data$T * 2)
})

test_that("SEM var_params() returns expected structure", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  vp <- var_params(fit)

  expect_type(vp, "list")
  expect_named(vp, c("mu", "Phi", "sigma", "rho"))
  expect_s3_class(vp$mu, "data.frame")
  expect_true(all(c("variable", "mean", "sd", "q2.5", "q97.5") %in% names(vp$mu)))
})

# --- Plotting ----------------------------------------------------------------

test_that("SEM plot latent_states returns ggplot", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  p <- plot(fit, type = "latent_states")
  expect_s3_class(p, "ggplot")
})

test_that("SEM plot diagnostics returns ggplot", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  p <- plot(fit, type = "diagnostics")
  expect_s3_class(p, "ggplot")
})

# --- Unsupported methods abort gracefully ------------------------------------

test_that("SEM fitted() errors informatively", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  expect_error(fitted(fit), "not yet implemented")
})

test_that("SEM predict() errors informatively", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  expect_error(predict(fit), "not yet implemented")
})

# --- Diagnostics -------------------------------------------------------------

test_that("SEM diagnostics are finite", {
  skip_if_no_cmdstanr()
  fit <- get_sem_fit()
  diag <- dcVar_diagnostics(fit)

  expect_equal(diag$n_divergent, 0)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.20)
  expect_true(diag$min_ess_bulk > 10)
  expect_true(diag$min_ess_tail > 10)
})
