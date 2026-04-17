# ============================================================================
# Shared test helper: fit minimal models once and cache for reuse
# ============================================================================

# Skip helpers
skip_if_no_rstan <- function() {
  skip_if_no_backend("rstan")
}

skip_if_no_cmdstanr_toolchain <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) skip("CmdStan not found")
  )
}

skip_if_no_backend <- function(backend = "rstan") {
  skip_on_cran()
  if (backend == "cmdstanr") {
    skip_if_not_installed("cmdstanr")
    skip_if(
      is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
      "CmdStan not installed"
    )
  }
  # rstan is always available (in Imports), no skip needed
}

skip_if_no_cmdstanr_backend <- function() {
  skip_if_no_backend("cmdstanr")
}

# Cache fitted objects in the test environment
cache_env <- new.env(parent = emptyenv())

.capture_fit_warnings <- function(expr) {
  warnings <- character()
  fit <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = warnings)
}

.cache_fit_result <- function(key, builder) {
  if (is.null(cache_env[[key]])) {
    cache_env[[key]] <- .capture_fit_warnings(builder())
  }
  cache_env[[key]]
}

expect_known_fit_warnings <- function(warnings, fit_name) {
  known_warning_patterns <- c(
    "divergent transitions after warmup",
    "Examine the pairs\\(\\) plot",
    "Bulk Effective Samples Size \\(ESS\\) is too low",
    "Tail Effective Samples Size \\(ESS\\) is too low",
    "largest R-hat is",
    "R-hat statistic is larger than",
    "sampling finished with diagnostic issues",
    "Inspect `dcvar_diagnostics\\(\\)`",
    "maximum treedepth",
    "There were [0-9]+ divergent transitions after warmup",
    "Bayesian Fraction of Missing Information",
    "E-BFMI",
    "BFMI"
  )

  unexpected_warnings <- warnings[!vapply(
    warnings,
    function(warning_text) any(vapply(
      known_warning_patterns,
      function(pattern) grepl(pattern, warning_text, perl = TRUE),
      logical(1)
    )),
    logical(1)
  )]

  expect_true(
    length(unexpected_warnings) == 0,
    info = paste0(
      fit_name,
      " fit emitted unexpected warnings: ",
      paste(unexpected_warnings, collapse = " | ")
    )
  )
}

core_iter_warmup <- 250
core_iter_sampling <- 250
hier_iter_warmup <- 300
hier_iter_sampling <- 300
smoke_iter_warmup <- 75
smoke_iter_sampling <- 75
margin_iter_warmup <- 300
margin_iter_sampling <- 300

.gamma_support_bound <- function(fit) {
  eps_draws <- posterior::as_draws_matrix(.fit_draws(
    fit$fit, "eps", backend = fit$backend,
    required = .stan_output_group_pattern("eps"),
    required_type = "pattern",
    context = ".gamma_support_bound()",
    output_type = "transformed parameter group"
  ))

  T_eff <- fit$stan_data$T - 1
  D <- fit$stan_data$D
  eps_mean <- matrix(NA_real_, nrow = T_eff, ncol = D)

  for (d in seq_len(D)) {
    cols <- paste0("eps[", seq_len(T_eff), ",", d, "]")
    eps_mean[, d] <- colMeans(eps_draws[, cols, drop = FALSE])
  }

  vapply(seq_len(D), function(d) {
    max(-fit$skew_direction[d] * eps_mean[, d])
  }, numeric(1))
}

get_dcvar_fit <- function() {
  .cache_fit_result("dcvar_fit", function() {
    sim <- simulate_dcvar(T = 50, rho_trajectory = rho_decreasing(50), seed = 42)
    dcvar(
      sim$Y_df, vars = c("y1", "y2"),
      chains = 2,
      iter_warmup = core_iter_warmup,
      iter_sampling = core_iter_sampling,
      adapt_delta = 0.995,
      refresh = 0, seed = 123
    )
  })$fit
}

get_dcvar_fit_warnings <- function() {
  .cache_fit_result("dcvar_fit", function() {
    sim <- simulate_dcvar(T = 50, rho_trajectory = rho_decreasing(50), seed = 42)
    dcvar(
      sim$Y_df, vars = c("y1", "y2"),
      chains = 2,
      iter_warmup = core_iter_warmup,
      iter_sampling = core_iter_sampling,
      adapt_delta = 0.995,
      refresh = 0, seed = 123
    )
  })$warnings
}

get_hmm_fit <- function() {
  .cache_fit_result("hmm_fit", function() {
    sim <- simulate_dcvar(T = 50, rho_trajectory = rho_step(50), seed = 42)
    dcvar_hmm(
      sim$Y_df, vars = c("y1", "y2"), K = 2,
      chains = 2,
      iter_warmup = max(core_iter_warmup, 250),
      iter_sampling = max(core_iter_sampling, 250),
      adapt_delta = 0.995,
      max_treedepth = 14,
      refresh = 0, seed = 123
    )
  })$fit
}

get_hmm_fit_warnings <- function() {
  .cache_fit_result("hmm_fit", function() {
    sim <- simulate_dcvar(T = 50, rho_trajectory = rho_step(50), seed = 42)
    dcvar_hmm(
      sim$Y_df, vars = c("y1", "y2"), K = 2,
      chains = 2,
      iter_warmup = max(core_iter_warmup, 250),
      iter_sampling = max(core_iter_sampling, 250),
      adapt_delta = 0.995,
      max_treedepth = 14,
      refresh = 0, seed = 123
    )
  })$warnings
}

get_constant_fit <- function() {
  .cache_fit_result("constant_fit", function() {
    sim <- simulate_dcvar(T = 50, rho_trajectory = rho_constant(50, 0.5), seed = 42)
    dcvar_constant(
      sim$Y_df, vars = c("y1", "y2"),
      chains = 2,
      iter_warmup = core_iter_warmup,
      iter_sampling = core_iter_sampling,
      adapt_delta = 0.995,
      refresh = 0, seed = 123
    )
  })$fit
}

get_constant_fit_warnings <- function() {
  .cache_fit_result("constant_fit", function() {
    sim <- simulate_dcvar(T = 50, rho_trajectory = rho_constant(50, 0.5), seed = 42)
    dcvar_constant(
      sim$Y_df, vars = c("y1", "y2"),
      chains = 2,
      iter_warmup = core_iter_warmup,
      iter_sampling = core_iter_sampling,
      adapt_delta = 0.995,
      refresh = 0, seed = 123
    )
  })$warnings
}

get_multilevel_fit <- function() {
  .cache_fit_result("multilevel_fit", function() {
    sim <- simulate_dcvar_multilevel(N = 3, T = 20, rho = 0.5, seed = 42)
    dcvar_multilevel(
      sim$data, vars = c("y1", "y2"), id_var = "id",
      chains = 2,
      iter_warmup = hier_iter_warmup,
      iter_sampling = hier_iter_sampling,
      adapt_delta = 0.995,
      max_treedepth = 14,
      refresh = 0, seed = 123
    )
  })$fit
}

get_multilevel_fit_warnings <- function() {
  .cache_fit_result("multilevel_fit", function() {
    sim <- simulate_dcvar_multilevel(N = 3, T = 20, rho = 0.5, seed = 42)
    dcvar_multilevel(
      sim$data, vars = c("y1", "y2"), id_var = "id",
      chains = 2,
      iter_warmup = hier_iter_warmup,
      iter_sampling = hier_iter_sampling,
      adapt_delta = 0.995,
      max_treedepth = 14,
      refresh = 0, seed = 123
    )
  })$warnings
}

get_sem_fit <- function() {
  .cache_fit_result("sem_fit", function() {
    J <- 2
    sim <- simulate_dcvar_sem(T = 50, J = J, lambda = rep(0.8, J),
                              rho = 0.5, seed = 42)
    indicators <- list(
      latent1 = paste0("y1_", seq_len(J)),
      latent2 = paste0("y2_", seq_len(J))
    )
    dcvar_sem(
      sim$data, indicators = indicators, J = J,
      lambda = rep(0.8, J), sigma_e = sqrt(0.2),
      chains = 2,
      iter_warmup = max(hier_iter_warmup, 350),
      iter_sampling = max(hier_iter_sampling, 350),
      adapt_delta = 0.999,
      max_treedepth = 14,
      refresh = 0, seed = 123
    )
  })$fit
}

get_sem_fit_warnings <- function() {
  .cache_fit_result("sem_fit", function() {
    J <- 2
    sim <- simulate_dcvar_sem(T = 50, J = J, lambda = rep(0.8, J),
                              rho = 0.5, seed = 42)
    indicators <- list(
      latent1 = paste0("y1_", seq_len(J)),
      latent2 = paste0("y2_", seq_len(J))
    )
    dcvar_sem(
      sim$data, indicators = indicators, J = J,
      lambda = rep(0.8, J), sigma_e = sqrt(0.2),
      chains = 2,
      iter_warmup = max(hier_iter_warmup, 350),
      iter_sampling = max(hier_iter_sampling, 350),
      adapt_delta = 0.999,
      max_treedepth = 14,
      refresh = 0, seed = 123
    )
  })$warnings
}

get_sem_exponential_fit <- function() {
  .cache_fit_result("sem_fit_exponential", function() {
    J <- 2
    sim <- simulate_dcvar_sem(
      T = 50,
      J = J,
      lambda = rep(0.8, J),
      margins = "exponential",
      sigma_exp = c(0.8, 1.1),
      skew_direction = c(1, -1),
      rho = 0.5,
      seed = 42
    )
    indicators <- list(
      latent1 = paste0("y1_", seq_len(J)),
      latent2 = paste0("y2_", seq_len(J))
    )
    dcvar_sem(
      sim$data,
      indicators = indicators,
      J = J,
      lambda = rep(0.8, J),
      sigma_e = sqrt(0.2),
      margins = "exponential",
      skew_direction = c(1, -1),
      chains = 2,
      iter_warmup = max(hier_iter_warmup, 350),
      iter_sampling = max(hier_iter_sampling, 350),
      adapt_delta = 0.999,
      max_treedepth = 14,
      refresh = 0,
      seed = 123
    )
  })$fit
}

get_sem_exponential_fit_warnings <- function() {
  .cache_fit_result("sem_fit_exponential", function() {
    J <- 2
    sim <- simulate_dcvar_sem(
      T = 50,
      J = J,
      lambda = rep(0.8, J),
      margins = "exponential",
      sigma_exp = c(0.8, 1.1),
      skew_direction = c(1, -1),
      rho = 0.5,
      seed = 42
    )
    indicators <- list(
      latent1 = paste0("y1_", seq_len(J)),
      latent2 = paste0("y2_", seq_len(J))
    )
    dcvar_sem(
      sim$data,
      indicators = indicators,
      J = J,
      lambda = rep(0.8, J),
      sigma_e = sqrt(0.2),
      margins = "exponential",
      skew_direction = c(1, -1),
      chains = 2,
      iter_warmup = max(hier_iter_warmup, 350),
      iter_sampling = max(hier_iter_sampling, 350),
      adapt_delta = 0.999,
      max_treedepth = 14,
      refresh = 0,
      seed = 123
    )
  })$warnings
}

get_dcvar_exponential_fit <- function() {
  .cache_fit_result("dcvar_fit_exponential", function() {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_decreasing(30),
      margins = "exponential",
      skew_direction = c(1, -1),
      seed = 42
    )
    dcvar(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "exponential",
      skew_direction = c(1, -1),
      chains = 2,
      iter_warmup = margin_iter_warmup,
      iter_sampling = margin_iter_sampling,
      adapt_delta = 0.999,
      refresh = 0,
      seed = 123
    )
  })$fit
}

get_dcvar_exponential_fit_warnings <- function() {
  .cache_fit_result("dcvar_fit_exponential", function() {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_decreasing(30),
      margins = "exponential",
      skew_direction = c(1, -1),
      seed = 42
    )
    dcvar(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "exponential",
      skew_direction = c(1, -1),
      chains = 2,
      iter_warmup = margin_iter_warmup,
      iter_sampling = margin_iter_sampling,
      adapt_delta = 0.999,
      refresh = 0,
      seed = 123
    )
  })$warnings
}

get_dcvar_gamma_fit <- function() {
  .cache_fit_result("dcvar_fit_gamma", function() {
    sim <- simulate_dcvar(
      T = 100,
      rho_trajectory = rho_decreasing(100),
      margins = "gamma",
      skew_direction = c(1, 1),
      skew_params = list(shape = 2),
      seed = 42
    )
    dcvar(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "gamma",
      skew_direction = c(1, 1),
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      adapt_delta = 0.999,
      max_treedepth = 15,
      refresh = 0,
      seed = 123
    )
  })$fit
}

get_dcvar_gamma_fit_warnings <- function() {
  .cache_fit_result("dcvar_fit_gamma", function() {
    sim <- simulate_dcvar(
      T = 100,
      rho_trajectory = rho_decreasing(100),
      margins = "gamma",
      skew_direction = c(1, 1),
      skew_params = list(shape = 2),
      seed = 42
    )
    dcvar(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "gamma",
      skew_direction = c(1, 1),
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      adapt_delta = 0.999,
      max_treedepth = 15,
      refresh = 0,
      seed = 123
    )
  })$warnings
}

get_dcvar_skew_normal_fit <- function() {
  skip_if_not_installed("sn")

  .cache_fit_result("dcvar_fit_skew_normal", function() {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_decreasing(30),
      margins = "skew_normal",
      skew_params = list(alpha = c(3, -3)),
      seed = 42
    )
    dcvar(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "skew_normal",
      chains = 1,
      iter_warmup = smoke_iter_warmup,
      iter_sampling = smoke_iter_sampling,
      refresh = 0,
      seed = 123
    )
  })$fit
}

get_dcvar_skew_normal_fit_warnings <- function() {
  skip_if_not_installed("sn")

  .cache_fit_result("dcvar_fit_skew_normal", function() {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_decreasing(30),
      margins = "skew_normal",
      skew_params = list(alpha = c(3, -3)),
      seed = 42
    )
    dcvar(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "skew_normal",
      chains = 1,
      iter_warmup = smoke_iter_warmup,
      iter_sampling = smoke_iter_sampling,
      refresh = 0,
      seed = 123
    )
  })$warnings
}

get_hmm_exponential_fit <- function() {
  .cache_fit_result("hmm_fit_exponential", function() {
    # Use a strongly separated step trajectory so the tiny test fit remains
    # identifiable across platforms and avoids pathological treedepth behavior.
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_step(30, rho_before = 0.9, rho_after = 0.1),
      margins = "exponential",
      skew_direction = c(1, -1),
      seed = 42
    )
    dcvar_hmm(
      sim$Y_df,
      vars = c("y1", "y2"),
      K = 2,
      margins = "exponential",
      skew_direction = c(1, -1),
      chains = 2,
      iter_warmup = max(margin_iter_warmup, 350),
      iter_sampling = max(margin_iter_sampling, 350),
      adapt_delta = 0.999,
      max_treedepth = 15,
      refresh = 0,
      seed = 123
    )
  })$fit
}

get_hmm_exponential_fit_warnings <- function() {
  .cache_fit_result("hmm_fit_exponential", function() {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_step(30, rho_before = 0.9, rho_after = 0.1),
      margins = "exponential",
      skew_direction = c(1, -1),
      seed = 42
    )
    dcvar_hmm(
      sim$Y_df,
      vars = c("y1", "y2"),
      K = 2,
      margins = "exponential",
      skew_direction = c(1, -1),
      chains = 2,
      iter_warmup = max(margin_iter_warmup, 350),
      iter_sampling = max(margin_iter_sampling, 350),
      adapt_delta = 0.999,
      max_treedepth = 15,
      refresh = 0,
      seed = 123
    )
  })$warnings
}

get_hmm_gamma_fit <- function() {
  .cache_fit_result("hmm_fit_gamma", function() {
    sim <- simulate_dcvar(
      T = 100,
      rho_trajectory = rho_step(100, rho_before = 0.9, rho_after = 0.1),
      margins = "gamma",
      skew_direction = c(1, 1),
      skew_params = list(shape = 2),
      seed = 42
    )
    dcvar_hmm(
      sim$Y_df,
      vars = c("y1", "y2"),
      K = 2,
      margins = "gamma",
      skew_direction = c(1, 1),
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      adapt_delta = 0.999,
      max_treedepth = 15,
      refresh = 0,
      seed = 123
    )
  })$fit
}

get_hmm_gamma_fit_warnings <- function() {
  .cache_fit_result("hmm_fit_gamma", function() {
    sim <- simulate_dcvar(
      T = 100,
      rho_trajectory = rho_step(100, rho_before = 0.9, rho_after = 0.1),
      margins = "gamma",
      skew_direction = c(1, 1),
      skew_params = list(shape = 2),
      seed = 42
    )
    dcvar_hmm(
      sim$Y_df,
      vars = c("y1", "y2"),
      K = 2,
      margins = "gamma",
      skew_direction = c(1, 1),
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      adapt_delta = 0.999,
      max_treedepth = 15,
      refresh = 0,
      seed = 123
    )
  })$warnings
}

get_hmm_skew_normal_fit <- function() {
  skip_if_not_installed("sn")

  .cache_fit_result("hmm_fit_skew_normal", function() {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_step(30, rho_before = 0.8, rho_after = 0.2),
      margins = "skew_normal",
      skew_params = list(alpha = c(3, -3)),
      seed = 42
    )
    dcvar_hmm(
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
  })$fit
}

get_hmm_skew_normal_fit_warnings <- function() {
  skip_if_not_installed("sn")

  .cache_fit_result("hmm_fit_skew_normal", function() {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_step(30, rho_before = 0.8, rho_after = 0.2),
      margins = "skew_normal",
      skew_params = list(alpha = c(3, -3)),
      seed = 42
    )
    dcvar_hmm(
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
  })$warnings
}

get_constant_gamma_fit <- function() {
  .cache_fit_result("constant_fit_gamma", function() {
    # Gamma-margin smoke fits need a bit more adaptation headroom than the
    # normal/exponential cases to stay stable across older R/rstan toolchains.
    sim <- simulate_dcvar(
      T = 100,
      rho_trajectory = rho_constant(100, 0.4),
      margins = "gamma",
      skew_direction = c(1, 1),
      skew_params = list(shape = 2),
      seed = 42
    )
    dcvar_constant(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "gamma",
      skew_direction = c(1, 1),
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      adapt_delta = 0.999,
      max_treedepth = 15,
      refresh = 0,
      seed = 123
    )
  })$fit
}

get_constant_gamma_fit_warnings <- function() {
  .cache_fit_result("constant_fit_gamma", function() {
    sim <- simulate_dcvar(
      T = 100,
      rho_trajectory = rho_constant(100, 0.4),
      margins = "gamma",
      skew_direction = c(1, 1),
      skew_params = list(shape = 2),
      seed = 42
    )
    dcvar_constant(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "gamma",
      skew_direction = c(1, 1),
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      adapt_delta = 0.999,
      max_treedepth = 15,
      refresh = 0,
      seed = 123
    )
  })$warnings
}

get_constant_skew_normal_fit <- function() {
  skip_if_not_installed("sn")

  .cache_fit_result("constant_fit_skew_normal", function() {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_constant(30, 0.4),
      margins = "skew_normal",
      skew_params = list(alpha = c(3, -3)),
      seed = 42
    )
    dcvar_constant(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "skew_normal",
      chains = 1,
      iter_warmup = smoke_iter_warmup,
      iter_sampling = smoke_iter_sampling,
      refresh = 0,
      seed = 123
    )
  })$fit
}

get_constant_skew_normal_fit_warnings <- function() {
  skip_if_not_installed("sn")

  .cache_fit_result("constant_fit_skew_normal", function() {
    sim <- simulate_dcvar(
      T = 30,
      rho_trajectory = rho_constant(30, 0.4),
      margins = "skew_normal",
      skew_params = list(alpha = c(3, -3)),
      seed = 42
    )
    dcvar_constant(
      sim$Y_df,
      vars = c("y1", "y2"),
      margins = "skew_normal",
      chains = 1,
      iter_warmup = smoke_iter_warmup,
      iter_sampling = smoke_iter_sampling,
      refresh = 0,
      seed = 123
    )
  })$warnings
}

get_dcvar_cmdstanr_fit <- function() {
  .cache_fit_result("dcvar_cmdstanr_fit", function() {
    skip_if_no_cmdstanr_backend()

    sim <- simulate_dcvar(T = 24, rho_trajectory = rho_decreasing(24), seed = 42)
    dcvar(
      sim$Y_df, vars = c("y1", "y2"),
      chains = 1,
      iter_warmup = 50,
      iter_sampling = 50,
      adapt_delta = 0.95,
      max_treedepth = 10,
      refresh = 0, seed = 123,
      backend = "cmdstanr"
    )
  })$fit
}

get_hmm_cmdstanr_fit <- function() {
  .cache_fit_result("hmm_cmdstanr_fit", function() {
    skip_if_no_cmdstanr_backend()

    sim <- simulate_dcvar(T = 24, rho_trajectory = rho_step(24), seed = 42)
    dcvar_hmm(
      sim$Y_df, vars = c("y1", "y2"), K = 2,
      chains = 1,
      iter_warmup = 50,
      iter_sampling = 50,
      adapt_delta = 0.95,
      max_treedepth = 10,
      refresh = 0, seed = 123,
      backend = "cmdstanr"
    )
  })$fit
}

get_constant_cmdstanr_fit <- function() {
  .cache_fit_result("constant_cmdstanr_fit", function() {
    skip_if_no_cmdstanr_backend()

    sim <- simulate_dcvar(T = 24, rho_trajectory = rho_constant(24, 0.5), seed = 42)
    dcvar_constant(
      sim$Y_df, vars = c("y1", "y2"),
      chains = 1,
      iter_warmup = 50,
      iter_sampling = 50,
      adapt_delta = 0.95,
      max_treedepth = 10,
      refresh = 0, seed = 123,
      backend = "cmdstanr"
    )
  })$fit
}
