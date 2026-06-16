# ============================================================================
# Tests for time-varying SEM support
# ============================================================================

make_sem_tv_stub_fit <- function(tv_phi = TRUE, tv_sigma = TRUE) {
  n_draw <- 30
  n_time <- 4L
  n_eff <- n_time - 1L
  phi_mask <- .resolve_phi_tv_mask(if (isTRUE(tv_phi)) "ar" else FALSE)
  n_phi_tv <- sum(phi_mask)
  indicators <- list(latent1 = c("y1_1", "y1_2"), latent2 = c("y2_1", "y2_2"))

  variables <- c(
    "mu[1]", "mu[2]",
    "Phi[1,1]", "Phi[1,2]", "Phi[2,1]", "Phi[2,2]",
    "sigma_eps[1]", "sigma_eps[2]",
    "z_rho_init", "sigma_omega",
    paste0("omega_raw[", seq_len(n_eff), "]"),
    paste0("rho[", seq_len(n_eff), "]"),
    paste0(
      "state[",
      rep(seq_len(n_time), each = 2),
      ",",
      rep(seq_len(2), times = n_time),
      "]"
    ),
    paste0(
      "zeta[",
      rep(seq_len(n_time), each = 2),
      ",",
      rep(seq_len(2), times = n_time),
      "]"
    )
  )
  if (isTRUE(tv_phi)) {
    variables <- c(
      variables,
      paste0("tau_phi[", seq_len(n_phi_tv), "]"),
      paste0(
        "phi_raw[",
        rep(seq_len(n_eff), each = n_phi_tv),
        ",",
        rep(seq_len(n_phi_tv), times = n_eff),
        "]"
      ),
      paste0(
        "phi_t[",
        rep(seq_len(n_eff), each = 4),
        ",",
        rep(seq_len(4), times = n_eff),
        "]"
      )
    )
  }
  if (isTRUE(tv_sigma)) {
    variables <- c(
      variables,
      "tau_sigma[1]", "tau_sigma[2]",
      paste0(
        "sigma_raw[",
        rep(seq_len(n_eff), each = 2),
        ",",
        rep(seq_len(2), times = n_eff),
        "]"
      ),
      paste0(
        "sigma_t[",
        rep(seq_len(n_eff), each = 2),
        ",",
        rep(seq_len(2), times = n_eff),
        "]"
      )
    )
  }

  draws <- array(
    0,
    dim = c(n_draw, 1, length(variables)),
    dimnames = list(NULL, NULL, variables)
  )
  fill <- function(name, value) {
    draws[, , name] <<- rep(value, length.out = n_draw)
  }

  fill("mu[1]", 0.1)
  fill("mu[2]", -0.2)
  fill("Phi[1,1]", 0.30)
  fill("Phi[1,2]", 0.10)
  fill("Phi[2,1]", -0.05)
  fill("Phi[2,2]", 0.40)
  fill("sigma_eps[1]", 0.8)
  fill("sigma_eps[2]", 1.2)
  fill("z_rho_init", 0.0)
  fill("sigma_omega", 0.08)
  for (t in seq_len(n_eff)) {
    fill(paste0("omega_raw[", t, "]"), 0)
    fill(paste0("rho[", t, "]"), -0.2 + 0.2 * t)
  }
  for (t in seq_len(n_time)) {
    fill(paste0("state[", t, ",1]"), 0.1 * t)
    fill(paste0("state[", t, ",2]"), -0.1 * t)
    fill(paste0("zeta[", t, ",1]"), 0.05 * t)
    fill(paste0("zeta[", t, ",2]"), -0.05 * t)
  }
  if (isTRUE(tv_phi)) {
    for (k in seq_len(n_phi_tv)) fill(paste0("tau_phi[", k, "]"), 0.03 + 0.01 * k)
    for (t in seq_len(n_eff)) {
      for (k in seq_len(n_phi_tv)) {
        fill(paste0("phi_raw[", t, ",", k, "]"), 0.01 * t)
      }
    }
    base_phi <- c(0.30, 0.10, -0.05, 0.40)
    for (t in seq_len(n_eff)) {
      for (k in seq_len(4)) {
        fill(paste0("phi_t[", t, ",", k, "]"), base_phi[k] + 0.02 * t)
      }
    }
  }
  if (isTRUE(tv_sigma)) {
    fill("tau_sigma[1]", 0.04)
    fill("tau_sigma[2]", 0.05)
    for (t in seq_len(n_eff)) {
      for (d in seq_len(2)) {
        fill(paste0("sigma_raw[", t, ",", d, "]"), 0.01 * d)
      }
    }
    for (t in seq_len(n_eff)) {
      fill(paste0("sigma_t[", t, ",1]"), 0.8 + 0.05 * t)
      fill(paste0("sigma_t[", t, ",2]"), 1.2 - 0.03 * t)
    }
  }

  stan_data <- structure(
    list(
      n_time = n_time,
      J = 2L,
      y = matrix(rnorm(n_time * 4), n_time, 4),
      lambda = c(0.8, 0.8),
      sigma_e = 0.5,
      family = c(1L, 1L),
      tv_phi = as.integer(isTRUE(tv_phi)),
      phi_tv_mask = unname(phi_mask),
      tv_sigma = as.integer(isTRUE(tv_sigma)),
      tau_phi_prior = 0.025,
      tau_sigma_prior = 0.05,
      sigma_omega_prior = 0.1,
      barrier_k = 8
    ),
    time_values = 11:14,
    indicators = indicators,
    vars = names(indicators),
    margins = "normal",
    method = "indicator",
    J = 2L,
    tv_phi = isTRUE(tv_phi),
    phi_tv_mask = phi_mask,
    tv_sigma = isTRUE(tv_sigma)
  )

  new_dcvar_sem_tv_fit(
    fit = posterior::as_draws_array(draws),
    stan_data = stan_data,
    vars = names(indicators),
    J = 2L,
    lambda = c(0.8, 0.8),
    sigma_e = 0.5,
    indicators = indicators,
    margins = "normal",
    method = "indicator",
    tv_phi = isTRUE(tv_phi),
    phi_tv_mask = phi_mask,
    tv_sigma = isTRUE(tv_sigma),
    backend = "rstan",
    priors = list(
      mu_sd = 0.25,
      phi_sd = 0.5,
      sigma_sd = 0.5,
      sigma_omega_rate = 0.1,
      rho_sd = 0.75,
      tau_phi_rate = 0.025,
      tau_sigma_rate = 0.05
    ),
    meta = list(chains = 1, iter_warmup = 10, iter_sampling = n_draw)
  )
}

test_that("SEM TV routing and prepared data expose the mixed engine contract", {
  expect_identical(.margin_stan_file("sem_tv", "normal"), "sem_tv_mixed.stan")
  expect_identical(.margin_stan_file("sem_tv", c("normal", "gamma")), "sem_tv_mixed.stan")
  expect_match(dcvar_stan_path("sem_tv"), "sem_tv_mixed\\.stan$")
  expect_true(file.exists(dcvar_stan_path("sem_tv")))

  sim <- simulate_dcvar_sem(
    n_time = 6,
    J = 2,
    lambda = rep(sqrt(0.8), 2),
    sigma_e = sqrt(0.2),
    rho_trajectory = rho_decreasing(6),
    phi_trajectory = list(rep(0.3, 5), rep(0.1, 5), rep(-0.05, 5), seq(0.2, 0.4, length.out = 5)),
    sigma_trajectory = cbind(seq(0.8, 1.1, length.out = 5), rep(1, 5)),
    seed = 101
  )
  expect_equal(dim(sim$true_params$Phi), c(5L, 4L))
  expect_equal(dim(sim$true_params$sigma), c(5L, 2L))
  expect_length(sim$true_params$rho, 5L)

  sd <- prepare_sem_data(
    sim$data,
    indicators = list(latent1 = c("y1_1", "y1_2"), latent2 = c("y2_1", "y2_2")),
    J = 2,
    lambda = rep(sqrt(0.8), 2),
    sigma_e = sqrt(0.2),
    tv_phi = "ar",
    tv_sigma = TRUE
  )
  expect_identical(sd$tv_phi, 1L)
  expect_identical(sd$phi_tv_mask, c(1L, 0L, 0L, 1L))
  expect_identical(sd$tv_sigma, 1L)
  expect_identical(sd$family, c(1L, 1L))
  expect_identical(sd$skew_direction, c(1, 1))
  expect_true(attr(sd, "tv_phi"))
  expect_true(attr(sd, "tv_sigma"))

  expect_error(
    prepare_sem_data(
      sim$data,
      indicators = list(latent1 = c("y1_1", "y1_2"), latent2 = c("y2_1", "y2_2")),
      J = 2,
      method = "naive",
      tv_phi = TRUE
    ),
    "not implemented"
  )
})

test_that("SEM TV initial values size active random-walk containers only", {
  init <- .init_sem_tv_params(6, "normal", tv_phi = "cross", tv_sigma = TRUE)
  expect_length(init$tau_phi, 2)
  expect_identical(dim(init$phi_raw), c(5L, 2L))
  expect_length(init$tau_sigma, 2)
  expect_identical(dim(init$sigma_raw), c(5L, 2L))
  expect_length(init$omega_raw, 5)

  constant <- .init_sem_tv_params(6, "normal", tv_phi = FALSE, tv_sigma = FALSE)
  expect_null(constant$tau_phi)
  expect_null(constant$phi_raw)
  expect_null(constant$tau_sigma)
  expect_null(constant$sigma_raw)
})

test_that("SEM TV fit object methods report trajectories and walk parameters", {
  fit <- make_sem_tv_stub_fit(tv_phi = TRUE, tv_sigma = TRUE)
  expect_s3_class(fit, "dcvar_sem_tv_fit")
  expect_identical(class(fit), c("dcvar_sem_tv_fit", "dcvar_sem_fit", "dcvar_model_fit"))
  expect_identical(fit$model, "sem_tv")
  expect_true(fit$tv_phi)
  expect_true(fit$tv_sigma)
  expect_identical(unname(fit$phi_tv_mask), c(1L, 0L, 0L, 1L))

  co <- coef(fit)
  expect_true(all(c("mu", "Phi", "sigma_eps", "sigma_omega", "rho", "tau_phi", "tau_sigma") %in% names(co)))
  expect_identical(names(co$tau_phi), c("tau_phi[phi11]", "tau_phi[phi22]"))
  expect_length(co$rho, fit$stan_data$n_time - 1L)

  vp <- var_params(fit)
  expect_true(all(c("mu", "Phi", "sigma_eps", "sigma_omega", "tau_phi", "tau_sigma", "rho") %in% names(vp)))
  expect_identical(vp$tau_phi$variable, c("tau_phi[phi11]", "tau_phi[phi22]"))

  rho_df <- rho_trajectory(fit)
  expect_equal(nrow(rho_df), 3L)
  expect_identical(rho_df$time, 12:14)

  phi_df <- phi_trajectory(fit)
  expect_equal(nrow(phi_df), 12L)
  expect_identical(sort(unique(phi_df$coefficient)), c("phi11", "phi12", "phi21", "phi22"))
  expect_gt(length(unique(phi_df$mean[phi_df$coefficient == "phi11"])), 1L)

  sigma_df <- sigma_trajectory(fit)
  expect_equal(nrow(sigma_df), 6L)
  expect_identical(sort(unique(sigma_df$variable)), sort(fit$vars))
  expect_true(all(sigma_df$mean > 0))

  dep_df <- dependence_summary(fit)
  expect_equal(nrow(dep_df), 3L)
  expect_true(all(abs(dep_df$mean) <= 1))

  expect_output(print(fit), "TV SEM")
  s <- summary(fit)
  expect_s3_class(s, "dcvar_sem_tv_summary")
  expect_false(is.null(s$phi_summary))
  expect_false(is.null(s$sigma_summary))
  expect_silent(invisible(capture.output(print(s))))
})

test_that("SEM TV phi_trajectory tiles the baseline Phi when tv_phi is off", {
  fit <- make_sem_tv_stub_fit(tv_phi = FALSE, tv_sigma = TRUE)
  expect_false(fit$tv_phi)

  phi_df <- phi_trajectory(fit)
  expect_equal(nrow(phi_df), 12L)
  per_coef <- tapply(phi_df$mean, phi_df$coefficient, function(x) length(unique(x)))
  expect_true(all(per_coef == 1L))

  sigma_df <- sigma_trajectory(fit)
  expect_equal(nrow(sigma_df), 6L)
  expect_true(all(sigma_df$mean > 0))
})
