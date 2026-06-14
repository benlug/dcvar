# ============================================================================
# Unified dynamic engine (dcvar_dynamic.stan): routing + density equivalence
#
# The DC-VAR time-varying path and the drift covariate model are both served by
# the single dcvar_dynamic.stan engine. These tests pin the routing and prove
# the engine reproduces the legacy densities: an independent R recomputation of
# the covariate log-likelihood (the z-score copula fast path), the beta_0 alias
# identity, and (gated, slow) posterior equivalence with the legacy covariate
# Stan file fitted side by side.
# ============================================================================

# --- routing (fast, no fit) -------------------------------------------------

test_that("the dynamic engine is reachable and public paths remain prepare-compatible", {
  dyn <- dcvar_stan_path("dcvar_dynamic")
  expect_true(file.exists(dyn))
  expect_match(dyn, "dcvar_dynamic\\.stan$")

  # Public model-specific paths match the exported prepare_*() data contracts;
  # wrappers route bundled TV / drift-covariate fits to dcvar_dynamic internally.
  expect_match(dcvar_stan_path("dcvar_tv"), "dcvar_tv_mixed\\.stan$")
  expect_match(dcvar_stan_path("dcvar_covariate"), "dcvar_covariate_ncp\\.stan$")
  expect_match(dcvar_stan_path("dcvar_covariate_nodrift"), "dcvar_covariate_nodrift\\.stan$")

  expect_true(.uses_dynamic_stan_file(NULL))
  expect_true(.uses_dynamic_stan_file(dyn))
  expect_false(.uses_dynamic_stan_file(dcvar_stan_path("dcvar_tv")))
  expect_false(.uses_dynamic_stan_file(dcvar_stan_path("dcvar_covariate")))
})

test_that(".as_dynamic_stan_data fills the engine block for both entry shapes", {
  df <- data.frame(time = 1:30, y1 = rnorm(30), y2 = rnorm(30))

  # TV entry: prepare_dcvar_data output (no covariate fields) -> P = 0 predictor.
  tv <- prepare_dcvar_data(df, c("y1", "y2"), tv_phi = "ar")
  tv_aug <- .as_dynamic_stan_data(tv)
  expect_identical(tv_aug$P, 0L)
  expect_identical(dim(tv_aug$X), c(30L, 0L))
  expect_identical(tv_aug$zero_init_eta, 0L)
  expect_identical(tv_aug$copula_normal_fastpath, 0L)
  # TV settings are preserved, not overwritten.
  expect_identical(tv_aug$tv_phi, 1L)
  expect_identical(unname(tv_aug$phi_tv_mask), c(1L, 0L, 0L, 1L))

  # Covariate entry: prepare_dcvar_covariate_data output -> margin/TV defaults.
  dfc <- df
  dfc$phase <- as.numeric(dfc$time > 15)
  cov <- prepare_dcvar_covariate_data(dfc, c("y1", "y2"), covariates = "phase")
  cov$copula_normal_fastpath <- 1L
  cov_aug <- .as_dynamic_stan_data(cov)
  expect_identical(cov_aug$family, c(1L, 1L))
  expect_identical(cov_aug$tv_phi, 0L)
  expect_identical(cov_aug$tv_sigma, 0L)
  expect_identical(cov_aug$P, 1L)
  expect_identical(cov_aug$copula_normal_fastpath, 1L)
  # The full superset block must be present (every engine-declared field).
  needed <- c("family", "skew_direction", "tv_phi", "phi_tv_mask", "tv_sigma",
              "tau_phi_prior", "tau_sigma_prior", "barrier_k", "sigma_eps_prior",
              "P", "X", "sigma_beta_prior", "zero_init_eta", "copula_normal_fastpath")
  expect_true(all(needed %in% names(cov_aug)))
})

test_that(".init_dcvar_dynamic_params omits inactive components and shapes length-1 vectors", {
  # P = 0, drift first state non-zero, no TV: beta and TV containers absent.
  i0 <- .init_dcvar_dynamic_params(2, 20, "normal", FALSE, FALSE, P = 0L, zero_init_eta = FALSE)
  expect_null(i0$beta)
  expect_null(i0$tau_phi)
  expect_null(i0$sigma_raw)
  expect_false("beta_0" %in% names(i0))      # the engine samples z_rho_init, not beta_0
  expect_true("z_rho_init" %in% names(i0))
  expect_length(i0$omega_raw, 19L)           # n_eff - 0

  # P = 1, zero_init_eta = TRUE: beta is length-1 array, omega_raw is n_eff - 1.
  i1 <- .init_dcvar_dynamic_params(2, 20, "normal", FALSE, FALSE, P = 1L, zero_init_eta = TRUE)
  expect_identical(dim(i1$beta), 1L)
  expect_length(i1$omega_raw, 18L)
})

# --- density equivalence (fast: one compiled fit, deterministic recomputation) --

test_that("covariate engine exposes beta_0 as a draw-for-draw alias of z_rho_init", {
  skip_if_no_rstan()

  fit <- get_covariate_fit(drift = TRUE)
  dm <- posterior::as_draws_matrix(draws(fit))
  expect_true(all(c("beta_0", "z_rho_init") %in% posterior::variables(dm)))
  expect_equal(as.numeric(dm[, "beta_0"]), as.numeric(dm[, "z_rho_init"]))
  expect_true(isTRUE(fit$dynamic_engine))
  monitored <- .diagnostic_parameter_variables(fit)
  expect_true("z_rho_init" %in% monitored)
  expect_false("beta_0" %in% monitored)
})

test_that("covariate engine log_lik matches an independent R recomputation", {
  skip_if_no_rstan()

  # Recompute log_lik[t] for one posterior draw using the documented covariate
  # density: independent normal margins + the z-score Gaussian copula on the
  # standardized residuals (the engine's fast path). A wrong fast path, scale, or
  # rho assembly would make these disagree.
  fit <- get_covariate_fit(drift = TRUE)
  dm <- posterior::as_draws_matrix(draws(fit))
  s <- 1L
  g <- function(nm) as.numeric(dm[s, nm])

  n_eff <- fit$stan_data$n_time - 1L
  sigma_eps <- c(g("sigma_eps[1]"), g("sigma_eps[2]"))

  eps <- matrix(NA_real_, n_eff, 2L)
  for (d in 1:2) {
    eps[, d] <- vapply(seq_len(n_eff), function(t) g(sprintf("eps[%d,%d]", t, d)), numeric(1))
  }
  rho <- vapply(seq_len(n_eff), function(t) g(sprintf("rho[%d]", t)), numeric(1))

  gauss_copula_z <- function(z1, z2, r) {
    rs <- r^2
    -0.5 * log1p(-rs) - (rs * (z1^2 + z2^2) - 2 * r * z1 * z2) / (2 * (1 - rs))
  }
  ll_R <- vapply(seq_len(n_eff), function(t) {
    z1 <- eps[t, 1] / sigma_eps[1]
    z2 <- eps[t, 2] / sigma_eps[2]
    stats::dnorm(eps[t, 1], 0, sigma_eps[1], log = TRUE) +
      stats::dnorm(eps[t, 2], 0, sigma_eps[2], log = TRUE) +
      gauss_copula_z(z1, z2, rho[t])
  }, numeric(1))

  ll_stan <- vapply(seq_len(n_eff), function(t) g(sprintf("log_lik[%d]", t)), numeric(1))
  expect_equal(ll_R, ll_stan, tolerance = 1e-8)
})

# --- engine coverage: TV path, multiple covariates, alt drift, methods ------

test_that("the TV path runs on the unified engine and emits its outputs (incl. beta_0 alias)", {
  skip_if_no_rstan()

  fit <- get_dcvar_tv_fit()
  expect_s3_class(fit, "dcvar_tv_fit")
  v <- posterior::variables(draws(fit))
  # The engine exposes beta_0 (= z_rho_init) even on the P = 0 TV path.
  expect_true(all(c("z_rho_init", "beta_0", "rho[1]", "log_lik[1]") %in% v))
  expect_identical(setdiff(.diagnostic_parameter_variables(fit), v), character(0))
})

test_that("explicit dcvar_dynamic stan_file uses the engine data block", {
  skip_if_no_rstan()

  dyn <- dcvar_stan_path("dcvar_dynamic")

  sim_tv <- simulate_dcvar(n_time = 24, rho_trajectory = rho_decreasing(24), seed = 21)
  fit_tv <- dcvar(
    sim_tv$Y_df, vars = c("y1", "y2"), tv_phi = "ar",
    stan_file = dyn,
    chains = 1, iter_warmup = 40, iter_sampling = 40, refresh = 0, seed = 22
  )
  expect_s3_class(fit_tv, "dcvar_tv_fit")
  expect_true("beta_0" %in% posterior::variables(draws(fit_tv)))

  sim_cov <- simulate_dcvar(n_time = 24, rho_trajectory = rho_step(24), seed = 23)
  df <- sim_cov$Y_df
  df$phase <- as.numeric(df$time > 12)
  fit_cov <- dcvar_covariate(
    df, vars = c("y1", "y2"), covariates = "phase",
    stan_file = dyn,
    chains = 1, iter_warmup = 40, iter_sampling = 40, refresh = 0, seed = 24
  )
  expect_s3_class(fit_cov, "dcvar_covariate_fit")
  expect_true(isTRUE(fit_cov$dynamic_engine))
  expect_true("z_rho_init" %in% .diagnostic_parameter_variables(fit_cov))
})

test_that("covariate engine fits with multiple covariates and zero_init_eta = FALSE", {
  skip_if_no_rstan()

  sim <- simulate_dcvar(
    n_time = 35,
    rho_trajectory = rho_step(35, rho_before = 0.6, rho_after = 0.2, breakpoint = 0.5),
    seed = 7
  )
  df <- sim$Y_df
  df$phase <- as.numeric(df$time > 17)
  df$trend <- as.numeric(scale(df$time))

  # P = 2 exercises the dot_product(X, beta) predictor; zero_init_eta = FALSE
  # exercises the compute_rw_ncp drift branch of the engine.
  fit <- dcvar_covariate(
    df, vars = c("y1", "y2"), covariates = c("phase", "trend"),
    drift = TRUE, zero_init_eta = FALSE,
    chains = 1, iter_warmup = 60, iter_sampling = 60, refresh = 0, seed = 11
  )
  expect_s3_class(fit, "dcvar_covariate_fit")
  expect_identical(fit$stan_data$P, 2L)
  expect_identical(fit$stan_data$zero_init_eta, 0L)

  ce <- covariate_effects(fit)
  expect_true(all(c("phase", "trend") %in% ce$term))
  expect_length(coef(fit)$beta, 2L)
  expect_equal(nrow(rho_trajectory(fit)), fit$stan_data$n_time - 1L)
  # Name-pin holds for P = 2, drift, zero_init_eta = FALSE (n_omega = n_eff).
  expect_identical(
    setdiff(.diagnostic_parameter_variables(fit), posterior::variables(draws(fit))),
    character(0)
  )
})

test_that("covariate engine fit supports predict / fitted / pit_values", {
  skip_if_no_rstan()

  fit <- get_covariate_fit(drift = TRUE)
  pr <- predict(fit)
  expect_true(all(pr$lower <= pr$upper))
  ft <- fitted(fit)
  expect_false(is.null(ft))
  pv <- pit_values(fit)
  expect_true(all(pv$pit >= 0 & pv$pit <= 1, na.rm = TRUE))
})

test_that("covariate engine fits under the cmdstanr backend", {
  skip_if_no_cmdstanr_backend()

  # Validates the engine's zero-extent init wrapping on the cmdstanr serializer
  # for the covariate (drift) path specifically.
  sim <- simulate_dcvar(
    n_time = 24,
    rho_trajectory = rho_step(24, rho_before = 0.6, rho_after = 0.2, breakpoint = 0.5),
    seed = 9
  )
  df <- sim$Y_df
  df$phase <- as.numeric(df$time > 12)
  fit <- dcvar_covariate(
    df, vars = c("y1", "y2"), covariates = "phase", drift = TRUE,
    chains = 1, iter_warmup = 50, iter_sampling = 50, refresh = 0,
    seed = 13, backend = "cmdstanr"
  )
  expect_s3_class(fit, "dcvar_covariate_fit")
  expect_identical(fit$backend, "cmdstanr")
  expect_true("beta_0" %in% posterior::variables(draws(fit)))
  expect_equal(nrow(rho_trajectory(fit)), fit$stan_data$n_time - 1L)
})

test_that("a custom stan_file bypasses the dynamic engine on the TV path", {
  skip_if_no_rstan()
  skip_if_not_slow()

  # With a user-supplied legacy TV Stan file, the data is NOT augmented (no
  # covariate predictor / drift fields) and the legacy parameter block is used.
  sim <- simulate_dcvar(
    n_time = 30, rho_trajectory = rho_decreasing(30, 0.6, 0.2), seed = 5
  )
  legacy_tv <- dcvar_stan_path("dcvar_tv")
  expect_true(nzchar(legacy_tv))
  fit <- dcvar(
    sim$Y_df, vars = c("y1", "y2"), tv_phi = "ar",
    stan_file = legacy_tv,
    chains = 1, iter_warmup = 60, iter_sampling = 60, refresh = 0, seed = 6
  )
  expect_s3_class(fit, "dcvar_tv_fit")
  v <- posterior::variables(draws(fit))
  # The legacy TV file has no covariate predictor, so beta_0 is absent.
  expect_false("beta_0" %in% v)
  expect_true(all(c("z_rho_init", "rho[1]") %in% v))
})

# --- posterior equivalence vs the legacy covariate Stan file (slow, gated) ---

test_that("covariate engine reproduces the legacy covariate posterior", {
  skip_if_no_rstan()
  skip_if_not_slow()

  sim <- simulate_dcvar(
    n_time = 50,
    rho_trajectory = rho_step(50, rho_before = 0.6, rho_after = 0.2, breakpoint = 0.5),
    seed = 42
  )
  df <- sim$Y_df
  df$phase <- as.numeric(df$time > 25)

  common <- list(
    data = df, vars = c("y1", "y2"), covariates = "phase",
    chains = 2, iter_warmup = 500, iter_sampling = 500,
    adapt_delta = 0.99, refresh = 0, seed = 999
  )

  fit_engine <- do.call(dcvar_covariate, common)
  legacy_path <- dcvar_stan_path("dcvar_covariate")
  expect_true(nzchar(legacy_path))
  fit_legacy <- do.call(dcvar_covariate, c(common, list(stan_file = legacy_path)))

  ce <- coef(fit_engine)
  cl <- coef(fit_legacy)
  expect_equal(unname(ce$mu), unname(cl$mu), tolerance = 0.05)
  expect_equal(unname(ce$Phi), unname(cl$Phi), tolerance = 0.06)
  expect_equal(unname(ce$sigma_eps), unname(cl$sigma_eps), tolerance = 0.05)
  expect_equal(unname(ce$beta_0), unname(cl$beta_0), tolerance = 0.08)
  expect_equal(unname(ce$beta), unname(cl$beta), tolerance = 0.08)
  expect_equal(unname(ce$sigma_omega), unname(cl$sigma_omega), tolerance = 0.06)

  # The dependence trajectories should agree across the whole series.
  re <- rho_trajectory(fit_engine)$mean
  rl <- rho_trajectory(fit_legacy)$mean
  expect_equal(re, rl, tolerance = 0.05)
})
