# ============================================================================
# Convergence Diagnostics and Parameter Recovery Tests
# ============================================================================
# These tests validate that fitted models produce reasonable MCMC diagnostics
# and that parameter estimates are in the right ballpark. They use the cached
# minimal fits from helper-fits.R (small iter counts), so thresholds are
# deliberately lenient.

# --- Convergence diagnostics ------------------------------------------------

test_that("cached baseline fits emit only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_dcvar_fit_warnings(), "dcvar")
  expect_known_fit_warnings(get_hmm_fit_warnings(), "hmm")
  expect_known_fit_warnings(get_constant_fit_warnings(), "constant")
})

test_that("gamma fits emit only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_dcvar_gamma_fit_warnings(), "dcvar gamma")
  expect_known_fit_warnings(get_hmm_gamma_fit_warnings(), "hmm gamma")
  expect_known_fit_warnings(get_constant_gamma_fit_warnings(), "constant gamma")
})

.expect_gamma_support_consistent <- function(fit, co, fit_name) {
  required_mean <- pmax(.gamma_support_bound(fit), 0)
  implied_mean <- sqrt(unname(co$shape_gam[[1]])) * unname(co$sigma_gam)

  expect_true(
    all(implied_mean + 1e-8 >= required_mean),
    info = paste(
      fit_name,
      "gamma support violated:",
      paste(
        sprintf(
          "dim %d mean=%.4f lb=%.4f",
          seq_along(implied_mean),
          implied_mean,
          required_mean
        ),
        collapse = ", "
      )
    )
  )
}

test_that("dcvar fit has no divergences and finite Rhat", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
  diag <- dcvar_diagnostics(fit)

  expect_equal(diag$n_divergent, 0)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.25)
  expect_true(diag$min_ess_bulk > 10)
  expect_true(diag$min_ess_tail > 10)
  expect_true(diag$mean_accept_prob > 0 && diag$mean_accept_prob < 1)
})

test_that("hmm fit has finite Rhat", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  diag <- dcvar_diagnostics(fit)

  expect_lte(diag$n_divergent, 1)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.30)
  expect_true(diag$min_ess_bulk > 8)
  expect_true(diag$min_ess_tail > 8)
})

test_that("constant fit has finite Rhat", {
  skip_if_no_rstan()

  fit <- get_constant_fit()
  diag <- dcvar_diagnostics(fit)

  expect_equal(diag$n_divergent, 0)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.25)
  expect_true(diag$min_ess_bulk > 10)
  expect_true(diag$min_ess_tail > 10)
})

test_that("dcvar exponential fit has usable diagnostics", {
  skip_if_no_rstan()

  fit <- get_dcvar_exponential_fit()
  diag <- dcvar_diagnostics(fit)
  co <- coef(fit)

  expect_equal(diag$n_divergent, 0)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.30)
  expect_true(all(co$sigma_exp > 0))
})

test_that("dcvar gamma fit has usable diagnostics", {
  skip_if_no_rstan()

  fit <- get_dcvar_gamma_fit()
  diag <- dcvar_diagnostics(fit)
  co <- coef(fit)

  expect_equal(diag$n_divergent, 0)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.35, info = paste("max_rhat =", signif(diag$max_rhat, 4)))
  expect_true(all(co$sigma_gam > 0))
  expect_true(co$shape_gam > 0)
  .expect_gamma_support_consistent(fit, co, "dcvar gamma")
})

test_that("hmm exponential fit has usable diagnostics", {
  skip_if_no_rstan()

  fit <- get_hmm_exponential_fit()
  diag <- dcvar_diagnostics(fit)
  co <- coef(fit)

  expect_equal(diag$n_divergent, 0)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.30)
  expect_true(all(co$sigma_exp > 0))
})

test_that("hmm gamma fit has usable diagnostics", {
  skip_if_no_rstan()

  fit <- get_hmm_gamma_fit()
  diag <- dcvar_diagnostics(fit)
  co <- coef(fit)

  expect_equal(diag$n_divergent, 0)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.35, info = paste("max_rhat =", signif(diag$max_rhat, 4)))
  expect_true(all(co$sigma_gam > 0))
  expect_true(co$shape_gam > 0)
  .expect_gamma_support_consistent(fit, co, "hmm gamma")
})

test_that("constant gamma fit has usable diagnostics", {
  skip_if_no_rstan()

  fit <- get_constant_gamma_fit()
  diag <- dcvar_diagnostics(fit)
  co <- coef(fit)

  expect_equal(diag$n_divergent, 0)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.35, info = paste("max_rhat =", signif(diag$max_rhat, 4)))
  expect_true(all(co$sigma_gam > 0))
  expect_true(co$shape_gam > 0)
  .expect_gamma_support_consistent(fit, co, "constant gamma")
})

# --- Parameter recovery: rho estimates are bounded ----------------------------

test_that("dcvar rho trajectory is within valid bounds", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
  rho_df <- rho_trajectory(fit)

  # All posterior means must be in (-1, 1)
  expect_true(all(rho_df$mean > -1 & rho_df$mean < 1))
  # Credible intervals must not extend beyond bounds
  expect_true(all(rho_df$q2.5 >= -1))
  expect_true(all(rho_df$q97.5 <= 1))
  # SDs must be positive

  expect_true(all(rho_df$sd > 0))
})

test_that("hmm rho_state values are ordered and bounded", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  states <- hmm_states(fit)

  rho_means <- states$rho_state$mean
  # Ordered constraint: state 1 rho < state 2 rho (for K=2)
  expect_true(all(diff(rho_means) >= 0))
  # Bounded in (-1, 1)
  expect_true(all(rho_means > -1 & rho_means < 1))
})

test_that("constant rho is bounded", {
  skip_if_no_rstan()

  fit <- get_constant_fit()
  rho_df <- rho_trajectory(fit)

  expect_true(all(rho_df$mean > -1 & rho_df$mean < 1))
  # All values should be identical (constant model)
  expect_equal(length(unique(rho_df$mean)), 1)
})

# --- Parameter recovery: VAR coefficients -----------------------------------

test_that("dcvar sigma_eps estimates are positive", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
  co <- coef(fit)

  expect_true(all(co$sigma_eps > 0))
  expect_true(co$sigma_omega > 0)
})

test_that("hmm sigma_eps estimates are positive", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  co <- coef(fit)

  expect_true(all(co$sigma_eps > 0))
})

test_that("constant sigma_eps estimates are positive", {
  skip_if_no_rstan()

  fit <- get_constant_fit()
  co <- coef(fit)

  expect_true(all(co$sigma_eps > 0))
})

# --- Parameter recovery: HMM transition matrix is valid ----------------------

test_that("hmm transition matrix rows sum to 1", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  states <- hmm_states(fit)

  row_sums <- rowSums(states$A)
  expect_equal(row_sums, rep(1, fit$K), tolerance = 1e-6)
  # All entries non-negative
  expect_true(all(states$A >= 0))
})

test_that("hmm gamma posteriors sum to 1 at each time point", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  states <- hmm_states(fit)

  row_sums <- rowSums(states$gamma)
  expect_equal(row_sums, rep(1, nrow(states$gamma)), tolerance = 1e-4)
})

test_that("hmm viterbi states are valid integers", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  states <- hmm_states(fit)

  expect_true(all(states$viterbi %in% 1:fit$K))
  expect_equal(length(states$viterbi), fit$stan_data$T - 1)
})

# --- simulate_dcvar input validation ----------------------------------------

test_that("simulate_dcvar rejects D != 2", {
  expect_error(
    simulate_dcvar(T = 50, rho_trajectory = rho_constant(50, 0.5),
                   mu = c(0, 0, 0)),
    "bivariate"
  )
})

test_that("simulate_dcvar rejects mismatched Phi dimensions", {
  expect_error(
    simulate_dcvar(T = 50, rho_trajectory = rho_constant(50, 0.5),
                   Phi = matrix(0, 3, 3)),
    "matrix"
  )
})

test_that("simulate_dcvar rejects mismatched sigma_eps length", {
  expect_error(
    simulate_dcvar(T = 50, rho_trajectory = rho_constant(50, 0.5),
                   sigma_eps = c(1, 1, 1)),
    "length"
  )
})
