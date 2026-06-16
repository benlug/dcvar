# ============================================================================
# Real-MCMC recovery tests for the time-varying SEM and multilevel engines
# (v0.9.0). The other TV test files exercise routing / data-prep / extractor
# contracts on stub fits; these tests actually sample sem_tv_mixed.stan and
# multilevel_tv_mixed.stan and check that the time-varying trajectories are
# recovered. Gated as slow MCMC fits (set DCVAR_SLOW_TESTS=true to run).
# ============================================================================

# Correlation of a posterior-mean trajectory (ordered by time) with the truth.
.traj_corr <- function(df, key_col, key, truth) {
  sub <- df[df[[key_col]] == key, , drop = FALSE]
  sub <- sub[order(sub$time), , drop = FALSE]
  stats::cor(sub$mean, truth)
}

test_that("time-varying SEM recovers AR / scale / coupling trajectories (MCMC)", {
  skip_if_not_slow()
  skip_if_no_rstan()

  J <- 2L
  n_time <- 200L
  n_eff <- n_time - 1L
  # phi11 (AR) ramps 0.0 -> 0.5; phi22 held constant (nominally TV under "ar");
  # sigma_1 ramps up; the Gaussian-copula correlation declines.
  phi_traj <- list(
    seq(0.0, 0.5, length.out = n_eff), # phi11 (varying)
    rep(0.10, n_eff),                  # phi12
    rep(0.05, n_eff),                  # phi21
    rep(0.30, n_eff)                   # phi22 (constant)
  )
  sigma_traj <- cbind(seq(0.8, 1.4, length.out = n_eff), rep(1.0, n_eff))
  sim <- simulate_dcvar_sem(
    n_time = n_time, J = J, lambda = rep(sqrt(0.8), J), sigma_e = sqrt(0.2),
    rho_trajectory = rho_decreasing(n_time, 0.7, 0.1),
    phi_trajectory = phi_traj, sigma_trajectory = sigma_traj, seed = 42
  )
  indicators <- list(
    latent1 = paste0("y1_", seq_len(J)),
    latent2 = paste0("y2_", seq_len(J))
  )

  fit <- dcvar_sem(
    sim$data, indicators = indicators, J = J,
    lambda = rep(sqrt(0.8), J), sigma_e = sqrt(0.2),
    tv_phi = "ar", tv_sigma = TRUE,
    chains = 2, iter_warmup = 500, iter_sampling = 500,
    adapt_delta = 0.995, max_treedepth = 12, refresh = 0, seed = 123
  )

  expect_s3_class(fit, "dcvar_sem_tv_fit")

  tp <- sim$true_params
  phi_df <- phi_trajectory(fit)
  sigma_df <- sigma_trajectory(fit)
  rho_df <- rho_trajectory(fit)

  # Varying AR coefficient is recovered; the nominally-TV constant one stays flat.
  expect_gt(.traj_corr(phi_df, "coefficient", "phi11", tp$Phi[, 1]), 0.5)
  expect_lt(stats::sd(phi_df$mean[phi_df$coefficient == "phi22"]), 0.15)
  # Varying innovation scale and time-varying coupling are recovered.
  svar <- sort(unique(sigma_df$variable))
  expect_gt(.traj_corr(sigma_df, "variable", svar[1], tp$sigma[, 1]), 0.6)
  expect_gt(stats::cor(rho_df$mean[order(rho_df$time)], tp$rho), 0.7)

  # No catastrophic non-convergence.
  expect_lt(dcvar_diagnostics(fit)$max_rhat, 1.15)
})

test_that("time-varying multilevel recovers cross-lag / scale trajectories (MCMC)", {
  skip_if_not_slow()
  skip_if_no_rstan()

  N <- 8L
  n_time <- 80L
  n_eff <- n_time - 1L
  # phi12 (cross-lag) population drift ramps 0.0 -> 0.3; phi21 held at 0;
  # sigma_1 ramps up; rho is constant (multilevel copula correlation is global).
  phi_traj <- list(
    rep(0.30, n_eff),                  # phi11 (constant baseline)
    seq(0.0, 0.30, length.out = n_eff),# phi12 (varying, "cross")
    rep(0.0, n_eff),                   # phi21
    rep(0.30, n_eff)                   # phi22
  )
  sigma_traj <- cbind(seq(0.8, 1.3, length.out = n_eff), rep(1.0, n_eff))
  sim <- simulate_dcvar_multilevel(
    N = N, n_time = n_time, rho = 0.4,
    phi_trajectory = phi_traj, sigma_trajectory = sigma_traj, seed = 42
  )

  fit <- dcvar_multilevel(
    sim$data, vars = c("y1", "y2"), id_var = "id",
    tv_phi = "cross", tv_sigma = TRUE,
    chains = 2, iter_warmup = 500, iter_sampling = 500,
    adapt_delta = 0.95, max_treedepth = 13, refresh = 0, seed = 123
  )

  expect_s3_class(fit, "dcvar_multilevel_tv_fit")

  tp <- sim$true_params
  phi_pop <- phi_trajectory(fit) # population drift (unit = NULL)
  sigma_df <- sigma_trajectory(fit)

  # Varying cross-lag (population path) recovered; constant baseline stays flat.
  expect_gt(.traj_corr(phi_pop, "coefficient", "phi12", tp$Phi_population[, 2]), 0.5)
  expect_lt(stats::sd(phi_pop$mean[phi_pop$coefficient == "phi11"]), 0.1)
  # Varying innovation scale recovered.
  svar <- sort(unique(sigma_df$variable))
  expect_gt(.traj_corr(sigma_df, "variable", svar[1], tp$sigma[, 1]), 0.6)

  expect_lt(dcvar_diagnostics(fit)$max_rhat, 1.15)
})
