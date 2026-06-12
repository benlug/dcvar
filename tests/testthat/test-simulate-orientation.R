# Regression tests for the copula orientation of simulated non-normal margins.
#
# The Stan likelihoods use u = 1 - F(x_shifted) on left-skewed dimensions
# (skew_direction = -1), so simulated data must imply the recorded rho -- not
# its negative -- when residuals are pushed through that convention. Before
# the fix, any configuration with exactly one negative skew_direction silently
# recovered -rho.

# Map simulated residuals to copula-scale z-scores using the likelihood's
# PIT convention (standardized margins: sigma_exp = 1, sd-1 gamma).
likelihood_copula_z <- function(eps, margin, skew, shape = 1) {
  u <- switch(margin,
    exponential = stats::pexp(1 + skew * eps, rate = 1),
    gamma = stats::pgamma(sqrt(shape) + skew * eps, shape = shape, rate = sqrt(shape))
  )
  if (skew < 0) u <- 1 - u
  stats::qnorm(u)
}

test_that(".sim_marginal_quantile preserves copula orientation for all skew_direction combinations", {
  set.seed(42)
  n <- 20000
  rho <- 0.6
  w1 <- rnorm(n)
  w2 <- rho * w1 + sqrt(1 - rho^2) * rnorm(n)

  for (margin in c("exponential", "gamma")) {
    shape <- if (margin == "gamma") 2 else 1
    skew_params <- if (margin == "gamma") list(shape = shape) else list()
    for (skew in list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))) {
      eps <- vapply(
        seq_len(n),
        function(i) {
          .sim_marginal_quantile(c(w1[i], w2[i]), margin, c(1, 1), skew, skew_params)
        },
        numeric(2)
      )
      z1 <- likelihood_copula_z(eps[1, ], margin, skew[1], shape)
      z2 <- likelihood_copula_z(eps[2, ], margin, skew[2], shape)
      expect_equal(
        cor(z1, z2), rho,
        tolerance = 0.05,
        info = sprintf("margin=%s skew=(%d,%d)", margin, skew[1], skew[2])
      )
    }
  }
})

test_that("simulate_dcvar data implies the recorded rho under the likelihood convention", {
  n_time <- 4000
  sim <- simulate_dcvar(
    n_time = n_time,
    rho_trajectory = rho_constant(n_time, 0.6),
    margins = "exponential",
    skew_direction = c(1, -1),
    seed = 7
  )
  tp <- sim$true_params
  Y <- sim$Y
  n_eff <- nrow(Y) - 1L
  eps <- matrix(NA_real_, n_eff, 2L)
  for (t in seq_len(n_eff)) {
    eps[t, ] <- Y[t + 1L, ] - (tp$mu + tp$Phi %*% (Y[t, ] - tp$mu))
  }

  z1 <- likelihood_copula_z(eps[, 1], "exponential", 1)
  z2 <- likelihood_copula_z(eps[, 2], "exponential", -1)
  expect_equal(cor(z1, z2), 0.6, tolerance = 0.05)
})

test_that("simulate_dcvar_sem exponential innovations imply the recorded rho", {
  sigma_exp <- c(1.5, 2)
  sim <- simulate_dcvar_sem(
    n_time = 4000,
    margins = "exponential",
    sigma_exp = sigma_exp,
    skew_direction = c(1, -1),
    rho = 0.6,
    seed = 11
  )
  zeta <- sim$innovations

  pit_sem <- function(zeta_i, sigma_exp_i, skew) {
    u <- stats::pexp(sigma_exp_i + skew * zeta_i, rate = 1 / sigma_exp_i)
    if (skew < 0) u <- 1 - u
    u
  }
  z1 <- stats::qnorm(pit_sem(zeta[, 1], sigma_exp[1], 1))
  z2 <- stats::qnorm(pit_sem(zeta[, 2], sigma_exp[2], -1))
  expect_equal(cor(z1, z2), 0.6, tolerance = 0.05)
})
