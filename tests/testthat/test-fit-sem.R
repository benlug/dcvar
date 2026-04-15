# ============================================================================
# Tests for dcvar_sem() fit and methods
# ============================================================================

# --- Fit object structure ----------------------------------------------------

test_that("dcvar_sem returns correct class", {
  skip_if_no_rstan()
  fit <- get_sem_fit()
  expect_s3_class(fit, "dcvar_sem_fit")
  expect_s3_class(fit, "dcvar_model_fit")
  expect_equal(fit$margins, "normal")
})

test_that("SEM coef() returns expected structure", {
  skip_if_no_rstan()
  fit <- get_sem_fit()
  co <- coef(fit)

  expect_type(co, "list")
  expect_named(co, c("mu", "Phi", "sigma", "rho"))
  expect_true(length(co$mu) == 2)
  expect_true(length(co$Phi) == 4)
  expect_true(length(co$sigma) == 2)
  expect_true(length(co$rho) == 1)
})

test_that("SEM exponential margins return expected structure", {
  skip_if_no_rstan()
  fit <- get_sem_exponential_fit()
  co <- coef(fit)
  vp <- var_params(fit)

  expect_equal(fit$margins, "exponential")
  expect_equal(fit$skew_direction, c(1, -1))
  expect_named(co, c("mu", "Phi", "sigma_exp", "rho"))
  expect_true(length(co$sigma_exp) == 2)
  expect_named(vp, c("mu", "Phi", "sigma_exp", "rho"))
  expect_s3_class(vp$sigma_exp, "data.frame")
})

test_that("SEM summary() returns correct class", {
  skip_if_no_rstan()
  fit <- get_sem_fit()
  s <- summary(fit)

  expect_s3_class(s, "dcvar_sem_summary")
  expect_true(all(c("var_params", "diagnostics") %in% names(s)))
})

test_that("SEM print() works", {
  skip_if_no_rstan()
  fit <- get_sem_fit()
  expect_output(print(fit), "SEM|sem")
})

# --- Extraction functions ----------------------------------------------------

test_that("SEM rho_trajectory() returns expected structure", {
  skip_if_no_rstan()
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
  skip_if_no_rstan()
  fit <- get_sem_fit()
  ls <- latent_states(fit)

  expect_s3_class(ls, "data.frame")
  expect_true(all(c("time", "variable", "mean", "sd") %in% names(ls)))
  expect_true(all(fit$vars %in% ls$variable))
  # Should have T rows per variable (2 variables)
  expect_equal(nrow(ls), fit$stan_data$T * 2)
})

test_that("SEM var_params() returns expected structure", {
  skip_if_no_rstan()
  fit <- get_sem_fit()
  vp <- var_params(fit)

  expect_type(vp, "list")
  expect_named(vp, c("mu", "Phi", "sigma", "rho"))
  expect_s3_class(vp$mu, "data.frame")
  expect_true(all(c("variable", "mean", "sd", "q2.5", "q97.5") %in% names(vp$mu)))
})

test_that("SEM exponential fitted() and predict() return expected structures", {
  skip_if_no_rstan()
  fit <- get_sem_exponential_fit()

  fit_link <- fitted(fit, type = "link")
  fit_resp <- fitted(fit, type = "response")
  pred_resp <- predict(fit, type = "response")
  indicator_names <- unlist(fit$indicators, use.names = FALSE)

  expect_s3_class(fit_link, "data.frame")
  expect_named(fit_link, c("time", fit$vars))
  expect_s3_class(fit_resp, "data.frame")
  expect_named(fit_resp, c("time", indicator_names))
  expect_s3_class(pred_resp, "data.frame")
  expect_named(pred_resp, c("time", "variable", "mean", "lower", "upper"))
})

# --- Plotting ----------------------------------------------------------------

test_that("SEM plot latent_states returns ggplot", {
  skip_if_no_rstan()
  fit <- get_sem_fit()
  p <- plot(fit, type = "latent_states")
  expect_s3_class(p, "ggplot")
})

test_that("SEM plot diagnostics returns ggplot", {
  skip_if_no_rstan()
  fit <- get_sem_fit()
  p <- plot(fit, type = "diagnostics")
  expect_s3_class(p, "ggplot")
})

# --- Fitted values and prediction --------------------------------------------

test_that("SEM fitted() returns latent and indicator trajectories", {
  skip_if_no_rstan()
  fit <- get_sem_fit()

  fit_link <- fitted(fit, type = "link")
  fit_resp <- fitted(fit, type = "response")
  indicator_names <- unlist(fit$indicators, use.names = FALSE)

  expect_s3_class(fit_link, "data.frame")
  expect_named(fit_link, c("time", fit$vars))
  expect_equal(nrow(fit_link), fit$stan_data$T)

  expect_s3_class(fit_resp, "data.frame")
  expect_named(fit_resp, c("time", indicator_names))
  expect_equal(nrow(fit_resp), fit$stan_data$T)
})

test_that("SEM predict() returns latent and indicator intervals", {
  skip_if_no_rstan()
  fit <- get_sem_fit()

  pred_link <- predict(fit, type = "link")
  pred_resp <- predict(fit, type = "response")
  indicator_names <- unlist(fit$indicators, use.names = FALSE)

  expect_s3_class(pred_link, "data.frame")
  expect_named(pred_link, c("time", "variable", "mean", "lower", "upper"))
  expect_equal(nrow(pred_link), fit$stan_data$T * length(fit$vars))
  expect_equal(sort(unique(pred_link$variable)), sort(fit$vars))

  expect_s3_class(pred_resp, "data.frame")
  expect_named(pred_resp, c("time", "variable", "mean", "lower", "upper"))
  expect_equal(nrow(pred_resp), fit$stan_data$T * length(indicator_names))
  expect_equal(sort(unique(pred_resp$variable)), sort(indicator_names))
  expect_true(all(pred_resp$lower <= pred_resp$mean))
  expect_true(all(pred_resp$upper >= pred_resp$mean))
})

# --- Diagnostics -------------------------------------------------------------

test_that("SEM diagnostics are finite", {
  skip_if_no_rstan()
  fit <- get_sem_fit()
  diag <- dcvar_diagnostics(fit)

  expect_lte(diag$n_divergent, 1)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.30)
  expect_true(diag$min_ess_bulk > 8)
  expect_true(diag$min_ess_tail > 8)
})

test_that("SEM fit cache emits only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_sem_fit_warnings(), "SEM")
})

test_that("SEM exponential fit cache emits only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_sem_exponential_fit_warnings(), "SEM exponential")
})
