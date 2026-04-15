# ============================================================================
# Shared test helper: fit minimal models once and cache for reuse
# ============================================================================

# Skip helpers
skip_if_no_cmdstanr <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) skip("CmdStan not found")
  )
}

# Cache fitted objects in the test environment
cache_env <- new.env(parent = emptyenv())

core_iter_warmup <- 150
core_iter_sampling <- 150
hier_iter_warmup <- 200
hier_iter_sampling <- 200
smoke_iter_warmup <- 75
smoke_iter_sampling <- 75
margin_iter_warmup <- 150
margin_iter_sampling <- 150

get_dcvar_fit <- function() {
  if (is.null(cache_env$dcvar_fit)) {
    sim <- simulate_dcvar(T = 50, rho_trajectory = rho_decreasing(50), seed = 42)
    cache_env$dcvar_fit <- dcvar(
      sim$Y_df, vars = c("y1", "y2"),
      chains = 2,
      iter_warmup = core_iter_warmup,
      iter_sampling = core_iter_sampling,
      refresh = 0, seed = 123
    )
  }
  cache_env$dcvar_fit
}

get_hmm_fit <- function() {
  if (is.null(cache_env$hmm_fit)) {
    sim <- simulate_dcvar(T = 50, rho_trajectory = rho_step(50), seed = 42)
    cache_env$hmm_fit <- dcvar_hmm(
      sim$Y_df, vars = c("y1", "y2"), K = 2,
      chains = 2,
      iter_warmup = core_iter_warmup,
      iter_sampling = core_iter_sampling,
      refresh = 0, seed = 123
    )
  }
  cache_env$hmm_fit
}

get_constant_fit <- function() {
  if (is.null(cache_env$constant_fit)) {
    sim <- simulate_dcvar(T = 50, rho_trajectory = rho_constant(50, 0.5), seed = 42)
    cache_env$constant_fit <- dcvar_constant(
      sim$Y_df, vars = c("y1", "y2"),
      chains = 2,
      iter_warmup = core_iter_warmup,
      iter_sampling = core_iter_sampling,
      refresh = 0, seed = 123
    )
  }
  cache_env$constant_fit
}

get_multilevel_fit <- function() {
  if (is.null(cache_env$multilevel_fit)) {
    sim <- simulate_dcvar_multilevel(N = 3, T = 20, rho = 0.5, seed = 42)
    cache_env$multilevel_fit <- dcvar_multilevel(
      sim$data, vars = c("y1", "y2"), id_var = "id",
      chains = 2,
      iter_warmup = hier_iter_warmup,
      iter_sampling = hier_iter_sampling,
      adapt_delta = 0.95,
      refresh = 0, seed = 123
    )
  }
  cache_env$multilevel_fit
}

get_sem_fit <- function() {
  if (is.null(cache_env$sem_fit)) {
    J <- 2
    sim <- simulate_dcvar_sem(T = 50, J = J, lambda = rep(0.8, J),
                              rho = 0.5, seed = 42)
    indicators <- list(
      latent1 = paste0("y1_", seq_len(J)),
      latent2 = paste0("y2_", seq_len(J))
    )
    cache_env$sem_fit <- dcvar_sem(
      sim$data, indicators = indicators, J = J,
      lambda = rep(0.8, J), sigma_e = sqrt(0.2),
      chains = 2,
      iter_warmup = hier_iter_warmup,
      iter_sampling = hier_iter_sampling,
      adapt_delta = 0.99,
      refresh = 0, seed = 123
    )
  }
  cache_env$sem_fit
}

get_dcvar_exponential_fit <- function() {
  if (is.null(cache_env$dcvar_fit_exponential)) {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_decreasing(30),
      margins = "exponential",
      skew_direction = c(1, -1),
      seed = 42
    )
    cache_env$dcvar_fit_exponential <- dcvar(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "exponential",
      skew_direction = c(1, -1),
      chains = 2,
      iter_warmup = margin_iter_warmup,
      iter_sampling = margin_iter_sampling,
      adapt_delta = 0.995,
      refresh = 0,
      seed = 123
    )
  }
  cache_env$dcvar_fit_exponential
}

get_dcvar_gamma_fit <- function() {
  if (is.null(cache_env$dcvar_fit_gamma)) {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_decreasing(30),
      margins = "gamma",
      skew_direction = c(1, 1),
      skew_params = list(shape = 2),
      seed = 42
    )
    cache_env$dcvar_fit_gamma <- dcvar(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "gamma",
      skew_direction = c(1, 1),
      chains = 2,
      iter_warmup = margin_iter_warmup,
      iter_sampling = margin_iter_sampling,
      refresh = 0,
      seed = 123
    )
  }
  cache_env$dcvar_fit_gamma
}

get_dcvar_skew_normal_fit <- function() {
  skip_if_not_installed("sn")

  if (is.null(cache_env$dcvar_fit_skew_normal)) {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_decreasing(30),
      margins = "skew_normal",
      skew_params = list(alpha = c(3, -3)),
      seed = 42
    )
    cache_env$dcvar_fit_skew_normal <- dcvar(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "skew_normal",
      chains = 1,
      iter_warmup = smoke_iter_warmup,
      iter_sampling = smoke_iter_sampling,
      refresh = 0,
      seed = 123
    )
  }
  cache_env$dcvar_fit_skew_normal
}

get_hmm_exponential_fit <- function() {
  if (is.null(cache_env$hmm_fit_exponential)) {
    # Use a strongly separated step trajectory so the tiny test fit remains
    # identifiable across platforms and avoids pathological treedepth behavior.
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_step(30, rho_before = 0.9, rho_after = 0.1),
      margins = "exponential",
      skew_direction = c(1, -1),
      seed = 42
    )
    cache_env$hmm_fit_exponential <- dcvar_hmm(
      sim$Y_df,
      vars = c("y1", "y2"),
      K = 2,
      margins = "exponential",
      skew_direction = c(1, -1),
      chains = 2,
      iter_warmup = margin_iter_warmup,
      iter_sampling = margin_iter_sampling,
      adapt_delta = 0.999,
      max_treedepth = 14,
      refresh = 0,
      seed = 123
    )
  }
  cache_env$hmm_fit_exponential
}

get_hmm_gamma_fit <- function() {
  if (is.null(cache_env$hmm_fit_gamma)) {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_step(30, rho_before = 0.8, rho_after = 0.2),
      margins = "gamma",
      skew_direction = c(1, 1),
      skew_params = list(shape = 2),
      seed = 42
    )
    cache_env$hmm_fit_gamma <- dcvar_hmm(
      sim$Y_df,
      vars = c("y1", "y2"),
      K = 2,
      margins = "gamma",
      skew_direction = c(1, 1),
      chains = 1,
      iter_warmup = smoke_iter_warmup,
      iter_sampling = smoke_iter_sampling,
      adapt_delta = 0.99,
      max_treedepth = 14,
      refresh = 0,
      seed = 123
    )
  }
  cache_env$hmm_fit_gamma
}

get_hmm_skew_normal_fit <- function() {
  skip_if_not_installed("sn")

  if (is.null(cache_env$hmm_fit_skew_normal)) {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_step(30, rho_before = 0.8, rho_after = 0.2),
      margins = "skew_normal",
      skew_params = list(alpha = c(3, -3)),
      seed = 42
    )
    cache_env$hmm_fit_skew_normal <- dcvar_hmm(
      sim$Y_df,
      vars = c("y1", "y2"),
      K = 2,
      margins = "skew_normal",
      chains = 1,
      iter_warmup = smoke_iter_warmup,
      iter_sampling = smoke_iter_sampling,
      adapt_delta = 0.99,
      max_treedepth = 14,
      refresh = 0,
      seed = 123
    )
  }
  cache_env$hmm_fit_skew_normal
}

get_constant_gamma_fit <- function() {
  if (is.null(cache_env$constant_fit_gamma)) {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_constant(30, 0.4),
      margins = "gamma",
      skew_direction = c(1, 1),
      skew_params = list(shape = 2),
      seed = 42
    )
    cache_env$constant_fit_gamma <- dcvar_constant(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "gamma",
      skew_direction = c(1, 1),
      chains = 2,
      iter_warmup = margin_iter_warmup,
      iter_sampling = margin_iter_sampling,
      refresh = 0,
      seed = 123
    )
  }
  cache_env$constant_fit_gamma
}

get_constant_skew_normal_fit <- function() {
  skip_if_not_installed("sn")

  if (is.null(cache_env$constant_fit_skew_normal)) {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_constant(30, 0.4),
      margins = "skew_normal",
      skew_params = list(alpha = c(3, -3)),
      seed = 42
    )
    cache_env$constant_fit_skew_normal <- dcvar_constant(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "skew_normal",
      chains = 1,
      iter_warmup = smoke_iter_warmup,
      iter_sampling = smoke_iter_sampling,
      refresh = 0,
      seed = 123
    )
  }
  cache_env$constant_fit_skew_normal
}
