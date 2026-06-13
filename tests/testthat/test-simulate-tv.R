# ============================================================================
# Time-varying simulator: trajectory resolution, Phi reshaping, and the
# copula orientation under per-t scales (the bug class fixed in 0.4.0).
# ============================================================================

test_that(".tv_resolve_trajectory accepts equivalent input forms", {
  n <- 21
  m <- cbind(rep(0.3, 20), rep(0.1, 20), rep(-0.2, 20), rep(0.4, 20))

  from_matrix <- .tv_resolve_trajectory(m, n, 4L, "phi_trajectory")
  from_list <- .tv_resolve_trajectory(
    list(rep(0.3, 20), rep(0.1, 20), rep(-0.2, 20), rep(0.4, 20)), n, 4L, "phi_trajectory"
  )
  from_const <- .tv_resolve_trajectory(
    matrix(c(0.3, -0.2, 0.1, 0.4), 2, 2), n, 4L, "phi_trajectory"
  )
  expect_equal(from_matrix, from_list)
  # Constant 2x2 expands row-major: (phi11, phi12, phi21, phi22)
  expect_equal(from_matrix, from_const)

  expect_error(.tv_resolve_trajectory(m[, 1:3], n, 4L, "phi_trajectory"), "matrix")
  expect_error(.tv_resolve_trajectory(list(rep(0.3, 19)), n, 4L, "phi_trajectory"), "4")
  expect_error(.tv_resolve_trajectory(cbind(rep(0.3, 20), rep(NA, 20)), n, 2L, "sigma_trajectory"), "finite")
})

test_that("simulate_dcvar with a constant phi_trajectory reproduces the constant-Phi path", {
  Phi <- matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2)
  s_const <- simulate_dcvar(30, rho_constant(30, 0.4), seed = 11)
  s_traj <- simulate_dcvar(30, rho_constant(30, 0.4), seed = 11, phi_trajectory = Phi)
  expect_equal(s_const$Y, s_traj$Y)
})

test_that("asymmetric Phi paths are applied in row-major order", {
  # phi12 = 0.4, phi21 = -0.2: a transposed reshape would swap them.
  n <- 4000
  phi_const <- matrix(c(0.2, 0.4, -0.2, 0.3), 2, 2, byrow = TRUE)
  sim <- simulate_dcvar(n, rho_constant(n, 0), seed = 21, phi_trajectory = phi_const)
  Y <- sim$Y

  # OLS VAR(1) recovery of the coefficient matrix
  X <- Y[-n, ]
  resp <- Y[-1, ]
  bhat <- solve(crossprod(X), crossprod(X, resp))  # 2x2: column d = coefficients of y_d
  expect_equal(unname(bhat[2, 1]), 0.4, tolerance = 0.1)   # phi12: y2's effect on y1
  expect_equal(unname(bhat[1, 2]), -0.2, tolerance = 0.1)  # phi21: y1's effect on y2
})

test_that("exp/gamma dimensions reject non-constant scale paths", {
  expect_error(
    simulate_dcvar(30, rho_constant(30, 0.4),
                   margins = c("normal", "exponential"),
                   skew_direction = c(1, 1),
                   sigma_trajectory = cbind(rep(1, 29), seq(0.5, 2, length.out = 29))),
    "constant on dimension 2"
  )
  # Constant non-unit scales on exp dims are allowed
  sim <- simulate_dcvar(30, rho_constant(30, 0.4),
                        margins = c("normal", "exponential"),
                        skew_direction = c(1, 1),
                        sigma_trajectory = cbind(seq(0.5, 2, length.out = 29), rep(1.5, 29)))
  expect_equal(unname(sim$true_params$sigma[1, ]), c(0.5, 1.5))
})

# Map simulated residuals to copula-scale z-scores using the likelihood's
# per-t PIT convention (u = F(scale + skew * eps), flipped on left-skew dims).
.tv_likelihood_copula_z <- function(eps, margin, skew, scale, shape = 1) {
  u <- switch(margin,
    normal = stats::pnorm(eps / scale),
    exponential = stats::pexp(scale + skew * eps, rate = 1 / scale),
    gamma = stats::pgamma(sqrt(shape) * scale + skew * eps,
                          shape = shape, rate = sqrt(shape) / scale)
  )
  if (margin != "normal" && skew < 0) u <- 1 - u
  stats::qnorm(u)
}

test_that("TV simulation preserves the copula orientation for every skew combination", {
  n <- 20000
  rho <- 0.6
  scale_path <- seq(0.5, 2, length.out = n - 1)  # non-constant scale on the normal dim

  for (margin2 in c("exponential", "gamma")) {
    shape <- if (margin2 == "gamma") 2 else 1
    skew_params <- if (margin2 == "gamma") list(shape = shape) else NULL
    for (skew2 in c(1, -1)) {
      sim <- simulate_dcvar(
        n, rho_constant(n, rho),
        mu = c(0, 0),
        margins = c("normal", margin2),
        skew_direction = c(1, skew2),
        skew_params = skew_params,
        phi_trajectory = matrix(c(0.2, 0.4, -0.2, 0.3), 2, 2, byrow = TRUE),
        sigma_trajectory = cbind(scale_path, rep(1.5, n - 1)),
        seed = 100 + skew2 + (margin2 == "gamma") * 10
      )

      # Recompute the residuals with the TRUE per-t coefficients
      tp <- sim$true_params
      Y <- sim$Y
      eps <- matrix(NA_real_, n - 1, 2)
      for (t in seq_len(n - 1)) {
        Phi_t <- matrix(tp$Phi[t, ], 2, 2, byrow = TRUE)
        eps[t, ] <- Y[t + 1, ] - (tp$mu + Phi_t %*% (Y[t, ] - tp$mu))
      }

      # Push through the likelihood PIT with the TRUE per-t scales
      z1 <- .tv_likelihood_copula_z(eps[, 1], "normal", 1, tp$sigma[, 1])
      z2 <- .tv_likelihood_copula_z(eps[, 2], margin2, skew2, tp$sigma[, 2], shape)
      expect_equal(
        cor(z1, z2), rho,
        tolerance = 0.05,
        info = sprintf("margin2=%s skew2=%d", margin2, skew2)
      )
    }
  }
})

test_that("TV simulation with time-varying rho recovers the rho path shape", {
  n <- 30000
  rho_path <- rho_decreasing(n, 0.7, 0.1)
  sim <- simulate_dcvar(
    n, rho_path,
    margins = "normal",
    phi_trajectory = list(rho_constant(n, 0.3), rho_constant(n, 0.1),
                          rho_constant(n, 0.1), rho_constant(n, 0.3)),
    sigma_trajectory = cbind(seq(1, 2, length.out = n - 1), seq(2, 1, length.out = n - 1)),
    seed = 77
  )
  tp <- sim$true_params
  Y <- sim$Y
  eps <- matrix(NA_real_, n - 1, 2)
  for (t in seq_len(n - 1)) {
    Phi_t <- matrix(tp$Phi[t, ], 2, 2, byrow = TRUE)
    eps[t, ] <- Y[t + 1, ] - (tp$mu + Phi_t %*% (Y[t, ] - tp$mu))
  }
  z <- eps / tp$sigma

  # Windowed empirical correlations should track the declining rho path
  early <- cor(z[1:5000, 1], z[1:5000, 2])
  late <- cor(z[(n - 5001):(n - 1), 1], z[(n - 5001):(n - 1), 2])
  expect_equal(early, mean(rho_path[1:5000]), tolerance = 0.07)
  expect_equal(late, mean(rho_path[(n - 5001):(n - 1)]), tolerance = 0.07)
  expect_gt(early, late)
})

# Soft-barrier generative model for time-varying exponential/gamma scales.

test_that("non-constant exp/gamma scale needs tv_sigma_k, then simulates", {
  n <- 60
  # Without tv_sigma_k a non-constant exp scale is rejected
  expect_error(
    simulate_dcvar(n, rho_constant(n, 0.4),
                   margins = c("normal", "exponential"),
                   skew_direction = c(1, 1),
                   sigma_trajectory = cbind(rep(1, n - 1), seq(0.6, 1.8, length.out = n - 1))),
    "tv_sigma_k"
  )
  # With tv_sigma_k it is allowed
  sim <- simulate_dcvar(n, rho_constant(n, 0.4),
                        margins = c("normal", "exponential"),
                        skew_direction = c(1, 1),
                        sigma_trajectory = cbind(rep(1, n - 1), seq(0.6, 1.8, length.out = n - 1)),
                        tv_sigma_k = 8, seed = 5)
  expect_equal(nrow(sim$Y), n)
  expect_equal(unname(sim$true_params$sigma[, 2]), seq(0.6, 1.8, length.out = n - 1))
})

test_that("soft-barrier simulation preserves copula orientation under a TV exp scale", {
  n <- 20000
  rho <- 0.6
  k <- 8
  softplus <- function(a) log1p(exp(k * a)) / k
  scale_path <- seq(0.6, 1.8, length.out = n - 1)

  for (skew2 in c(1, -1)) {
    sim <- simulate_dcvar(
      n, rho_constant(n, rho),
      margins = c("normal", "exponential"),
      skew_direction = c(1, skew2),
      phi_trajectory = matrix(c(0.2, 0.4, -0.2, 0.3), 2, 2, byrow = TRUE),
      sigma_trajectory = cbind(rep(1, n - 1), scale_path),
      tv_sigma_k = k,
      seed = 200 + skew2
    )
    tp <- sim$true_params
    Y <- sim$Y
    eps <- matrix(NA_real_, n - 1, 2)
    for (t in seq_len(n - 1)) {
      P <- matrix(tp$Phi[t, ], 2, 2, byrow = TRUE)
      eps[t, ] <- Y[t + 1, ] - (tp$mu + P %*% (Y[t, ] - tp$mu))
    }
    z1 <- stats::qnorm(stats::pnorm(eps[, 1]))            # normal dim, scale 1
    m <- scale_path                                       # exp scale path
    u2 <- stats::pexp(softplus(m + skew2 * eps[, 2]), rate = 1 / m)
    if (skew2 < 0) u2 <- 1 - u2
    z2 <- stats::qnorm(u2)
    expect_equal(cor(z1, z2), rho, tolerance = 0.05,
                 info = sprintf("skew2=%d", skew2))
  }
})

test_that("soft-barrier exp/gamma change-of-variables density is normalized (Jacobian present and correct)", {
  # The eps-density of the soft-barrier margin is
  #   exp/gamma_lpdf(softplus_k(m + skew*eps) | rate) + log_inv_logit(k*(m + skew*eps)).
  # A missing, doubled, or sign-flipped Jacobian would make this integrate to
  # something other than 1. This is the falsifiable guard the orientation
  # round-trips (which only check the copula rank correlation) cannot provide.
  k <- 8
  softplus <- function(a) log1p(exp(k * a)) / k
  log_jac <- function(a) stats::plogis(k * a, log.p = TRUE)   # log_inv_logit(k*a)

  exp_eps_density <- function(eps, m, skew) {
    arg <- m + skew * eps
    exp(stats::dexp(softplus(arg), rate = 1 / m, log = TRUE) + log_jac(arg))
  }
  gam_eps_density <- function(eps, m, skew, shape) {
    arg <- m + skew * eps
    exp(stats::dgamma(softplus(arg), shape = shape, rate = shape / m, log = TRUE) + log_jac(arg))
  }

  for (skew in c(1, -1)) {
    for (m in c(0.7, 1.3)) {
      area_exp <- stats::integrate(exp_eps_density, -Inf, Inf, m = m, skew = skew,
                                   rel.tol = 1e-8)$value
      expect_equal(area_exp, 1, tolerance = 1e-4,
                   info = sprintf("exp m=%.1f skew=%d", m, skew))
      area_gam <- stats::integrate(gam_eps_density, -Inf, Inf, m = m, skew = skew,
                                   shape = 2, rel.tol = 1e-8)$value
      expect_equal(area_gam, 1, tolerance = 1e-4,
                   info = sprintf("gamma m=%.1f skew=%d", m, skew))
    }
  }
})

test_that("the soft-barrier log-Jacobian is the exact transform Jacobian", {
  # log_inv_logit(k*a) must equal log|d softplus_k(a)/da|; a sign flip or a
  # constant offset (the easy regressions) would break this finite-difference
  # identity that the Stan model relies on.
  k <- 8
  softplus <- function(a) log1p(exp(k * a)) / k
  a <- c(-2, -0.5, 0, 0.5, 2, 5)
  num <- (softplus(a + 1e-6) - softplus(a - 1e-6)) / 2e-6
  ana <- exp(stats::plogis(k * a, log.p = TRUE))
  expect_equal(ana, num, tolerance = 1e-5)
})
