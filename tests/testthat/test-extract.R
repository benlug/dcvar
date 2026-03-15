test_that("rho_trajectory() returns correct structure for dcVar", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  rho_df <- rho_trajectory(fit)

  expect_s3_class(rho_df, "data.frame")
  expect_true(all(c("time", "mean", "sd") %in% names(rho_df)))
  expect_true("q2.5" %in% names(rho_df))
  expect_true("q97.5" %in% names(rho_df))
  expect_true(all(rho_df$mean >= -1 & rho_df$mean <= 1))
  expect_equal(rho_df$time, attr(fit$stan_data, "time_values")[-1])
})

test_that("rho_trajectory() returns correct structure for hmm", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  rho_df <- rho_trajectory(fit)

  expect_s3_class(rho_df, "data.frame")
  expect_true(all(c("time", "mean", "sd") %in% names(rho_df)))
  expect_equal(rho_df$time, attr(fit$stan_data, "time_values")[-1])
})

test_that("rho_trajectory() returns full trajectory for constant", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  rho_df <- rho_trajectory(fit)

  expect_s3_class(rho_df, "data.frame")
  # Constant model now expands rho to all T-1 time points
  expect_equal(nrow(rho_df), fit$stan_data$T - 1)
  expect_true(all(c("time", "mean", "sd") %in% names(rho_df)))
  # All means should be identical (constant)
  expect_equal(length(unique(rho_df$mean)), 1)
  expect_true(all(rho_df$mean >= -1 & rho_df$mean <= 1))
  expect_equal(rho_df$time, attr(fit$stan_data, "time_values")[-1])
})

test_that("var_params() returns correct structure", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  vp <- var_params(fit)

  expect_type(vp, "list")
  expect_true(all(c("mu", "Phi", "sigma_eps") %in% names(vp)))
  expect_true("sigma_omega" %in% names(vp))
  expect_true(all(c("variable", "mean", "sd", "q2.5", "q97.5") %in% names(vp$mu)))
})

test_that("extractors work for non-normal fitted objects", {
  skip_if_no_cmdstanr()

  rho_df <- rho_trajectory(get_dcVar_exponential_fit())
  expect_s3_class(rho_df, "data.frame")
  expect_true(all(rho_df$mean >= -1 & rho_df$mean <= 1))

  gamma_params <- var_params(get_constant_gamma_fit())
  expect_true("sigma_gam" %in% names(gamma_params))
  expect_true("shape_gam" %in% names(gamma_params))

  skip_if_not_installed("sn")
  skew_params <- var_params(get_constant_skew_normal_fit())
  expect_true("omega" %in% names(skew_params))
  expect_true("delta" %in% names(skew_params))
})

test_that("hmm_states() returns correct structure", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  states <- hmm_states(fit)

  expect_type(states, "list")
  expect_named(states, c("gamma", "viterbi", "rho_state", "A", "rho_hmm"))
  expect_equal(ncol(states$gamma), fit$K)
  expect_equal(dim(states$A), c(fit$K, fit$K))
})

test_that("hmm_states() returns a Viterbi path observed in the posterior draws", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  states <- hmm_states(fit)
  viterbi_draws <- posterior::as_draws_matrix(fit$fit$draws("viterbi_state"))
  draw_paths <- apply(viterbi_draws, 1, paste, collapse = ",")

  expect_true(paste(states$viterbi, collapse = ",") %in% draw_paths)
})

test_that("draws() returns posterior draws", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  d <- draws(fit, variable = "mu")
  expect_true(inherits(d, "draws_array"))
})

test_that("dcVar_diagnostics() returns correct fields", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  diag <- dcVar_diagnostics(fit)

  expect_type(diag, "list")
  expect_named(diag, c("n_divergent", "n_max_treedepth", "max_rhat",
                        "min_ess_bulk", "min_ess_tail", "mean_accept_prob"))
  expect_true(diag$max_rhat >= 1)
})

test_that("rho_trajectory() supports a single quantile request", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  rho_df <- rho_trajectory(fit, probs = 0.5)

  expect_true("q50" %in% names(rho_df))
})

test_that("latent_states() supports a single quantile request", {
  skip_if_no_cmdstanr()

  fit <- get_sem_fit()
  states <- latent_states(fit, probs = 0.5)

  expect_true("q50" %in% names(states))
})

test_that("rho_trajectory() and latent_states() honor preserved time values", {
  skip_if_no_cmdstanr()

  fit <- get_dcVar_fit()
  attr(fit$stan_data, "time_values") <- 101:(100 + fit$stan_data$T)
  expect_equal(rho_trajectory(fit)$time, 102:(100 + fit$stan_data$T))

  sem_fit <- get_sem_fit()
  attr(sem_fit$stan_data, "time_values") <- seq.Date(
    as.Date("2020-01-01"),
    by = "day",
    length.out = sem_fit$stan_data$T
  )
  sem_rho <- rho_trajectory(sem_fit)
  sem_states <- latent_states(sem_fit)
  expect_equal(sem_rho$time, attr(sem_fit$stan_data, "time_values")[-1])
  expect_equal(
    unique(sem_states$time),
    attr(sem_fit$stan_data, "time_values")
  )
})

test_that("random_effects() preserves original unit ids", {
  skip_if_no_cmdstanr()

  fit <- get_multilevel_fit()
  attr(fit$stan_data, "ids") <- c("A", "B", "C")
  re <- random_effects(fit)

  expect_equal(unique(re$unit), c("A", "B", "C"))
})
