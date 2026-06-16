# ============================================================================
# Tests for dcvar_sem() fit and methods
# ============================================================================

# --- Fit object structure ----------------------------------------------------

make_sem_mixed_engine_stub_fit <- function(margins = c("skew_normal", "gamma")) {
  margins <- match.arg(margins)
  n_time <- 4
  variables <- c(
    "mu[1]", "mu[2]",
    "phi11", "phi12", "phi21", "phi22",
    "Phi[1,1]", "Phi[1,2]", "Phi[2,1]", "Phi[2,2]",
    "rho_raw", "rho",
    paste0("state[", rep(seq_len(n_time), each = 2), ",", rep(1:2, times = n_time), "]"),
    paste0("log_lik[", seq_len(n_time), "]")
  )
  if (identical(margins, "skew_normal")) {
    variables <- c(variables, "omega[1]", "omega[2]", "delta[1]", "delta[2]")
  } else {
    variables <- c(
      variables,
      "eta[1]", "eta[2]",
      "sigma_gam[1]", "sigma_gam[2]",
      "shape_gam[1]", "shape_gam[2]"
    )
  }

  draws <- array(
    rnorm(40 * length(variables), 0, 0.1),
    dim = c(40, 1, length(variables)),
    dimnames = list(NULL, NULL, variables)
  )
  draws[, , "rho"] <- runif(40, -0.5, 0.5)
  draws[, , paste0("log_lik[", seq_len(n_time), "]")] <-
    -abs(draws[, , paste0("log_lik[", seq_len(n_time), "]")])
  if (identical(margins, "skew_normal")) {
    draws[, , c("omega[1]", "omega[2]")] <- runif(40 * 2, 0.5, 1.5)
    draws[, , c("delta[1]", "delta[2]")] <- runif(40 * 2, -0.5, 0.5)
  } else {
    draws[, , c("sigma_gam[1]", "sigma_gam[2]")] <- runif(40 * 2, 0.5, 1.5)
    draws[, , c("shape_gam[1]", "shape_gam[2]")] <- runif(40 * 2, 0.5, 2.0)
  }

  indicators <- list(
    latent1 = c("y1_1", "y1_2"),
    latent2 = c("y2_1", "y2_2")
  )
  stan_data <- structure(
    list(
      n_time = n_time,
      J = 2,
      y = matrix(rnorm(n_time * 4), n_time, 4),
      lambda = c(0.8, 0.8),
      sigma_e = 0.5
    ),
    time_values = 21:24,
    indicators = indicators,
    vars = names(indicators),
    margins = margins,
    method = "indicator",
    J = 2,
    skew_direction = if (identical(margins, "gamma")) c(1, -1) else NULL
  )

  new_dcvar_sem_fit(
    fit = posterior::as_draws_array(draws),
    stan_data = stan_data,
    vars = names(indicators),
    J = 2,
    lambda = c(0.8, 0.8),
    sigma_e = 0.5,
    indicators = indicators,
    margins = margins,
    method = "indicator",
    skew_direction = attr(stan_data, "skew_direction"),
    backend = "rstan",
    priors = list(),
    meta = list(chains = 1, iter_warmup = 10, iter_sampling = 40)
  )
}

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

test_that("SEM homogeneous skew_normal and gamma stubs report mixed-engine params", {
  skew_fit <- make_sem_mixed_engine_stub_fit("skew_normal")
  skew_co <- coef(skew_fit)
  skew_vp <- var_params(skew_fit)
  expect_named(skew_co, c("mu", "Phi", "omega", "delta", "rho"))
  expect_equal(names(skew_co$omega), c("omega[1]", "omega[2]"))
  expect_true(all(c("omega", "delta") %in% names(skew_vp)))

  gamma_fit <- make_sem_mixed_engine_stub_fit("gamma")
  gamma_co <- coef(gamma_fit)
  gamma_vp <- var_params(gamma_fit)
  expect_named(gamma_co, c("mu", "Phi", "sigma_gam", "shape_gam", "rho"))
  expect_equal(names(gamma_co$shape_gam), c("shape_gam[1]", "shape_gam[2]"))
  expect_true(all(c("sigma_gam", "shape_gam") %in% names(gamma_vp)))
})

test_that("SEM homogeneous skew_normal and gamma slow fits run", {
  skip_if_no_rstan()
  skip_if_not_slow()

  for (margin in c("skew_normal", "gamma")) {
    J <- 2
    skew_direction <- if (identical(margin, "gamma")) c(1, 1) else NULL
    sim <- simulate_dcvar_sem(
      n_time = 80, J = J, lambda = rep(0.8, J), sigma_e = sqrt(0.2),
      margins = margin,
      skew_direction = skew_direction,
      skew_params = if (identical(margin, "skew_normal")) list(alpha = c(2, 2)) else list(shape = 2),
      rho = 0.5,
      seed = 500 + match(margin, c("skew_normal", "gamma"))
    )
    fit <- dcvar_sem(
      sim$data,
      indicators = list(latent1 = paste0("y1_", seq_len(J)), latent2 = paste0("y2_", seq_len(J))),
      J = J,
      lambda = rep(0.8, J),
      sigma_e = sqrt(0.2),
      margins = margin,
      skew_direction = skew_direction,
      chains = 2, iter_warmup = 250, iter_sampling = 250,
      adapt_delta = 0.999, max_treedepth = 13, refresh = 0, seed = 123
    )

    expect_equal(fit$margins, margin)
    expect_equal(fit$stan_data$family, rep(unname(.family_codes[[margin]]), 2))
    expect_true(coef(fit)$rho > -0.2)
  }
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
  expect_equal(nrow(ls), fit$stan_data$n_time * 2)
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
  expect_equal(nrow(fit_link), fit$stan_data$n_time)

  expect_s3_class(fit_resp, "data.frame")
  expect_named(fit_resp, c("time", indicator_names))
  expect_equal(nrow(fit_resp), fit$stan_data$n_time)
})

test_that("SEM predict() returns latent and indicator intervals", {
  skip_if_no_rstan()
  fit <- get_sem_fit()

  pred_link <- predict(fit, type = "link")
  pred_resp <- predict(fit, type = "response")
  indicator_names <- unlist(fit$indicators, use.names = FALSE)

  expect_s3_class(pred_link, "data.frame")
  expect_named(pred_link, c("time", "variable", "mean", "lower", "upper"))
  expect_equal(nrow(pred_link), fit$stan_data$n_time * length(fit$vars))
  expect_equal(sort(unique(pred_link$variable)), sort(fit$vars))

  expect_s3_class(pred_resp, "data.frame")
  expect_named(pred_resp, c("time", "variable", "mean", "lower", "upper"))
  expect_equal(nrow(pred_resp), fit$stan_data$n_time * length(indicator_names))
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
