# ============================================================================
# Routing, data preparation, and fit-object tests for the time-varying DC-VAR
# ============================================================================

test_that("dcvar_tv routes to the single generic Stan file for all margins", {
  expect_identical(.margin_stan_file("dcvar_tv", "normal"), "dcvar_tv_mixed.stan")
  expect_identical(.margin_stan_file("dcvar_tv", c("normal", "exponential")), "dcvar_tv_mixed.stan")
  # Homogeneous non-normal also routes to the generic file (no specialised zoo)
  expect_identical(.margin_stan_file("dcvar_tv", c("gamma", "gamma")), "dcvar_tv_mixed.stan")
  expect_error(.margin_stan_file("dcvar_tv", "normal", copula = "clayton"), "Gaussian")

  expect_identical(
    .margin_cache_key("dcvar_tv", c("normal", "exponential")),
    "dcvar_tv_mixed12_model"
  )
  expect_identical(.margin_cache_key("dcvar_tv", "normal"), "dcvar_tv_mixed11_model")

  path <- dcvar_stan_path("dcvar_tv")
  expect_true(file.exists(path))
  expect_match(path, "dcvar_tv_mixed\\.stan$")
})

test_that(".resolve_phi_tv_mask maps logicals and selectors to a row-major mask", {
  expect_identical(unname(.resolve_phi_tv_mask(TRUE)), c(1L, 1L, 1L, 1L))
  expect_identical(unname(.resolve_phi_tv_mask(FALSE)), c(0L, 0L, 0L, 0L))
  # AR = phi11, phi22 (positions 1, 4); cross = phi12, phi21 (positions 2, 3)
  expect_identical(unname(.resolve_phi_tv_mask("ar")), c(1L, 0L, 0L, 1L))
  expect_identical(unname(.resolve_phi_tv_mask("cross")), c(0L, 1L, 1L, 0L))
  expect_identical(unname(.resolve_phi_tv_mask("crosslag")), c(0L, 1L, 1L, 0L))
  expect_identical(unname(.resolve_phi_tv_mask(c("phi11", "phi21"))), c(1L, 0L, 1L, 0L))
  expect_identical(names(.resolve_phi_tv_mask(TRUE)),
                   c("phi11", "phi12", "phi21", "phi22"))

  expect_error(.resolve_phi_tv_mask("nonsense"), "Invalid")
  expect_error(.resolve_phi_tv_mask(NA), "logical")
  expect_error(.resolve_phi_tv_mask(c(TRUE, FALSE)), "logical")
  expect_error(.resolve_phi_tv_mask(character(0)), "empty")
})

test_that("prepare_dcvar_data resolves a tv_phi selector to a coefficient mask", {
  df <- data.frame(time = 1:30, y1 = rnorm(30), y2 = rnorm(30))

  ar <- prepare_dcvar_data(df, vars = c("y1", "y2"), tv_phi = "ar")
  expect_identical(ar$tv_phi, 1L)
  expect_identical(ar$phi_tv_mask, c(1L, 0L, 0L, 1L))
  expect_identical(ar$tv_sigma, 0L)
  expect_identical(attr(ar, "phi_tv_mask"), .resolve_phi_tv_mask("ar"))

  cross <- prepare_dcvar_data(df, vars = c("y1", "y2"), tv_phi = "cross")
  expect_identical(cross$phi_tv_mask, c(0L, 1L, 1L, 0L))

  full <- prepare_dcvar_data(df, vars = c("y1", "y2"), tv_phi = TRUE)
  expect_identical(full$phi_tv_mask, c(1L, 1L, 1L, 1L))
})

test_that("prepare_dcvar_data with flags off is byte-identical to before", {
  df <- data.frame(time = 1:30, y1 = rnorm(30), y2 = rnorm(30))
  plain <- prepare_dcvar_data(df, vars = c("y1", "y2"))
  flagged_off <- prepare_dcvar_data(df, vars = c("y1", "y2"),
                                    tv_phi = FALSE, tv_sigma = FALSE)
  expect_identical(plain, flagged_off)
  expect_null(plain$tv_phi)
  expect_null(plain$family)
})

test_that("prepare_dcvar_data adds the TV fields when a flag is set", {
  df <- data.frame(time = 1:30, y1 = rnorm(30), y2 = rnorm(30))
  out <- prepare_dcvar_data(df, vars = c("y1", "y2"), tv_phi = TRUE, tv_sigma = TRUE)

  expect_identical(out$tv_phi, 1L)
  expect_identical(out$tv_sigma, 1L)
  expect_identical(out$family, c(1L, 1L))
  expect_identical(out$skew_direction, c(1, 1))
  expect_identical(out$tau_phi_prior, 0.025)
  expect_identical(out$tau_sigma_prior, 0.05)
  expect_identical(out$barrier_k, 8)
  expect_true(attr(out, "tv_phi"))
  expect_true(attr(out, "tv_sigma"))

  mixed <- prepare_dcvar_data(df, vars = c("y1", "y2"),
                              margins = c("normal", "exponential"),
                              skew_direction = c(1, -1),
                              tv_phi = TRUE)
  expect_identical(mixed$family, c(1L, 2L))
  expect_identical(mixed$skew_direction, c(1, -1))
})

test_that("prepare_dcvar_data validates the TV flags", {
  df <- data.frame(time = 1:30, y1 = rnorm(30), y2 = rnorm(30))
  expect_error(prepare_dcvar_data(df, c("y1", "y2"), tv_phi = NA), "tv_phi")
  expect_error(prepare_dcvar_data(df, c("y1", "y2"), tv_phi = c(TRUE, TRUE)), "tv_phi")
  expect_error(prepare_dcvar_data(df, c("y1", "y2"), tv_phi = TRUE,
                                  prior_tau_phi_rate = -1), "prior_tau_phi_rate")
})

test_that(".init_dcvar_tv_params sizes the walk containers by flag", {
  i_both <- .init_dcvar_tv_params(2, 40, "normal", tv_phi = TRUE, tv_sigma = TRUE)
  expect_length(i_both$tau_phi, 4)
  expect_identical(dim(i_both$phi_raw), c(39L, 4L))
  expect_length(i_both$tau_sigma, 2)
  expect_identical(dim(i_both$sigma_raw), c(39L, 2L))
  # Full union initialised regardless of margins
  expect_length(i_both$shape_gam, 2)
  expect_length(i_both$omega_raw, 39)

  i_phi <- .init_dcvar_tv_params(2, 40, "normal", tv_phi = TRUE, tv_sigma = FALSE)
  expect_length(i_phi$tau_sigma, 0)
  expect_identical(dim(i_phi$sigma_raw), c(0L, 2L))

  i_sig <- .init_dcvar_tv_params(2, 40, "normal", tv_phi = FALSE, tv_sigma = TRUE)
  expect_length(i_sig$tau_phi, 0)
  # phi_raw has zero active columns when no coefficient varies
  expect_identical(dim(i_sig$phi_raw), c(0L, 0L))

  # AR selector: 2 active coefficients -> 2 walk columns
  i_ar <- .init_dcvar_tv_params(2, 40, "normal", tv_phi = "ar", tv_sigma = FALSE)
  expect_length(i_ar$tau_phi, 2)
  expect_identical(dim(i_ar$phi_raw), c(39L, 2L))
})

# --- Smoke fit: structure, extractors, diagnostics name pin -----------------

test_that("dcvar(tv_phi, tv_sigma) fits and returns the subclass", {
  skip_if_no_rstan()

  fit <- get_dcvar_tv_fit()
  expect_s3_class(fit, "dcvar_tv_fit")
  expect_identical(class(fit), c("dcvar_tv_fit", "dcvar_fit", "dcvar_model_fit"))
  expect_identical(fit$model, "dcvar_tv")
  expect_true(fit$tv_phi)
  expect_true(fit$tv_sigma)
  expect_named(fit$priors, c("mu_sd", "phi_sd", "sigma_eps_rate", "sigma_omega_rate",
                             "rho_init_sd", "tau_phi_rate", "tau_sigma_rate"))
})

test_that("TV trajectory extractors return the documented structure", {
  skip_if_no_rstan()

  fit <- get_dcvar_tv_fit()
  n_eff <- fit$stan_data$n_time - 1L

  phi_df <- phi_trajectory(fit)
  expect_s3_class(phi_df, "data.frame")
  expect_equal(nrow(phi_df), 4 * n_eff)
  expect_true(all(c("time", "coefficient", "mean", "sd", "q2.5", "q97.5") %in% names(phi_df)))
  expect_identical(sort(unique(phi_df$coefficient)), c("phi11", "phi12", "phi21", "phi22"))

  sigma_df <- sigma_trajectory(fit)
  expect_equal(nrow(sigma_df), 2 * n_eff)
  expect_identical(sort(unique(sigma_df$variable)), sort(fit$vars))
  expect_true(all(sigma_df$mean > 0))

  rho_df <- rho_trajectory(fit)
  expect_equal(nrow(rho_df), n_eff)
  expect_true(all(abs(rho_df$mean) <= 1))
})

test_that("TV fit methods run and report the time-varying components", {
  skip_if_no_rstan()

  fit <- get_dcvar_tv_fit()

  out <- capture.output(print(fit))
  expect_true(any(grepl("Phi\\(t\\)", out)))

  s <- summary(fit)
  expect_s3_class(s, "dcvar_tv_summary")
  expect_false(is.null(s$phi_summary))
  expect_false(is.null(s$sigma_summary))
  expect_silent(invisible(capture.output(print(s))))

  co <- coef(fit)
  expect_true(all(c("mu", "Phi", "tau_phi", "tau_sigma", "sigma_omega") %in% names(co)))
  expect_length(co$tau_phi, 4)

  vp <- var_params(fit)
  expect_true(all(c("mu", "Phi", "tau_phi", "tau_sigma", "sigma_omega") %in% names(vp)))

  expect_s3_class(plot(fit, type = "phi"), "ggplot")
  expect_s3_class(plot(fit, type = "sigma"), "ggplot")
  expect_s3_class(plot(fit, type = "rho"), "ggplot")

  pit_df <- pit_values(fit)
  expect_true(all(pit_df$pit >= 0 & pit_df$pit <= 1, na.rm = TRUE))

  pred <- predict(fit)
  expect_true(all(pred$lower <= pred$upper))

  ll <- loo(fit)
  expect_s3_class(ll, "loo")
})

test_that("diagnostics monitoring list matches the sampled Stan parameters", {
  skip_if_no_rstan()

  fit <- get_dcvar_tv_fit()
  monitored <- .diagnostic_parameter_variables(fit)
  available <- posterior::variables(draws(fit))
  missing <- setdiff(monitored, available)
  expect_identical(missing, character(0))

  diag <- dcvar_diagnostics(fit)
  expect_true(is.finite(diag$max_rhat))
})

test_that("phi/sigma trajectories tile constants when a flag is off", {
  skip_if_no_rstan()

  fit <- get_dcvar_tv_sigma_only_fit()
  expect_false(fit$tv_phi)
  n_eff <- fit$stan_data$n_time - 1L

  phi_df <- phi_trajectory(fit)
  expect_equal(nrow(phi_df), 4 * n_eff)
  # Constant baseline tiled: one unique mean per coefficient
  per_coef <- tapply(phi_df$mean, phi_df$coefficient, function(x) length(unique(x)))
  expect_true(all(per_coef == 1L))

  sigma_df <- sigma_trajectory(fit)
  expect_equal(nrow(sigma_df), 2 * n_eff)
})

test_that("AR-only tv_phi fits, sizes tau_phi compactly, and pins diagnostics", {
  skip_if_no_rstan()

  fit <- get_dcvar_tv_ar_fit()
  expect_s3_class(fit, "dcvar_tv_fit")
  expect_true(fit$tv_phi)
  expect_false(fit$tv_sigma)
  expect_identical(unname(fit$phi_tv_mask), c(1L, 0L, 0L, 1L))
  expect_identical(names(fit$phi_tv_mask), c("phi11", "phi12", "phi21", "phi22"))

  # tau_phi is sized to the 2 active coefficients and labelled by name
  co <- coef(fit)
  expect_length(co$tau_phi, 2)
  expect_identical(names(co$tau_phi), c("tau_phi[phi11]", "tau_phi[phi22]"))

  vp <- var_params(fit)
  expect_identical(vp$tau_phi$variable, c("tau_phi[phi11]", "tau_phi[phi22]"))

  # The monitored diagnostic parameters must exist in the draws AND the draws
  # must not contain a tau_phi index beyond the active count -- pins the masked
  # tau_phi sizing in both directions (over- and under-counting).
  monitored <- .diagnostic_parameter_variables(fit)
  available <- posterior::variables(draws(fit))
  expect_identical(setdiff(monitored, available), character(0))
  expect_true(all(c("tau_phi[1]", "tau_phi[2]") %in% available))
  expect_false("tau_phi[3]" %in% available)

  # phi_trajectory still returns all four coefficients; the cross-lagged ones
  # are constant (deviation forced to zero), while the active AR coefficients
  # genuinely vary -- the latter ties the masked-in coefficients to actual
  # time variation (catches a mask that only leaves phi12/phi21 off).
  phi_df <- phi_trajectory(fit)
  expect_identical(sort(unique(phi_df$coefficient)), c("phi11", "phi12", "phi21", "phi22"))
  uniq <- tapply(phi_df$mean, phi_df$coefficient, function(x) length(unique(x)))
  expect_equal(unname(uniq[["phi12"]]), 1L)
  expect_equal(unname(uniq[["phi21"]]), 1L)
  # phi11 is simulated as a decreasing path, so its posterior-mean trajectory
  # must not collapse to a constant.
  expect_gt(unname(uniq[["phi11"]]), 1L)

  out <- capture.output(print(fit))
  expect_true(any(grepl("phi11", out)))

  expect_s3_class(plot(fit, type = "diagnostics"), "ggplot")
})

test_that("dcvar() warns about ignored TV prior args and exp/gamma scales", {
  df <- data.frame(time = 1:25, y1 = rnorm(25), y2 = rnorm(25))
  expect_warning(
    tryCatch(
      dcvar(df, c("y1", "y2"), tv_phi = FALSE, prior_tau_phi_rate = 0.5,
            chains = 1, iter_warmup = 2, iter_sampling = 2, refresh = 0,
            stan_file = NA_character_),
      error = function(e) NULL
    ),
    "prior_tau_phi_rate"
  )
})

test_that("dcvar_compare warns when mixing TV and non-TV fits", {
  skip_if_no_rstan()

  # The cached fits use different data lengths, so loo_compare() itself errors
  # after the warning under test has fired; swallow the error.
  w <- testthat::capture_warnings(
    try(dcvar_compare(tv = get_dcvar_tv_fit(), constant = get_constant_fit()),
        silent = TRUE)
  )
  expect_true(any(grepl("latent paths", w)))
})

test_that("soft-barrier exp margin with tv_sigma fits and exposes a scale path", {
  skip_if_no_rstan()

  fit <- get_dcvar_tv_exp_fit()
  expect_s3_class(fit, "dcvar_tv_fit")
  expect_true(fit$tv_sigma)
  expect_identical(unname(rep(fit$margins, length.out = 2)), c("normal", "exponential"))

  # sigma_t exists for both dims and is positive; the exp dim genuinely varies
  sigma_df <- sigma_trajectory(fit)
  expect_identical(sort(unique(sigma_df$variable)), sort(fit$vars))
  expect_true(all(sigma_df$mean > 0))
  exp_dim <- sigma_df[sigma_df$variable == fit$vars[2], ]
  expect_gt(length(unique(exp_dim$mean)), 1L)

  # Monitored diagnostics exist (the soft-barrier path keeps the union names)
  monitored <- .diagnostic_parameter_variables(fit)
  expect_identical(setdiff(monitored, posterior::variables(draws(fit))), character(0))

  ll <- loo(fit)
  expect_s3_class(ll, "loo")
})

# --- Gated recovery tests (DCVAR_SLOW_TESTS) ---------------------------------

test_that("TV fit recovers a phi12 step and a sigma ramp directionally", {
  skip_if_no_rstan()
  skip_if_not_slow()

  n <- 200
  sim <- simulate_dcvar(
    n, rho_constant(n, 0.4),
    phi_trajectory = list(
      rho_constant(n, 0.3),
      rho_step(n, rho_before = 0.5, rho_after = 0.0, breakpoint = 0.5),
      rho_constant(n, 0.1),
      rho_constant(n, 0.3)
    ),
    sigma_trajectory = cbind(seq(0.7, 1.8, length.out = n - 1), rep(1, n - 1)),
    seed = 7
  )
  fit <- dcvar(sim$Y_df, vars = c("y1", "y2"),
               tv_phi = TRUE, tv_sigma = TRUE,
               standardize = FALSE,
               chains = 2,
               iter_warmup = margin_iter_warmup,
               iter_sampling = margin_iter_sampling,
               refresh = 0, seed = 11)

  phi_df <- phi_trajectory(fit)
  phi12 <- phi_df[phi_df$coefficient == "phi12", ]
  n_eff <- nrow(phi12)
  expect_gt(mean(phi12$mean[1:50]), mean(phi12$mean[(n_eff - 49):n_eff]))

  sigma_df <- sigma_trajectory(fit)
  s1 <- sigma_df[sigma_df$variable == "y1", ]
  expect_gt(mean(s1$mean[(n_eff - 49):n_eff]), mean(s1$mean[1:50]))
})

test_that("tv_phi = 'ar' ties the fitted variation to the correct coefficient", {
  skip_if_no_rstan()
  skip_if_not_slow()

  # Only phi11 genuinely varies (a clear ramp); phi22 is constant. Under
  # tv_phi = "ar" both phi11 and phi22 carry a walk, so this test confirms the
  # FITTED phi11(t) path moves substantially more than phi22(t) -- i.e. the
  # active tau_phi/phi_dev columns map to the coefficient they are labelled
  # with. A Stan-idx vs R-name transposition (phi11 <-> phi22) would fail here.
  n <- 200
  sim <- simulate_dcvar(
    n, rho_constant(n, 0.4),
    phi_trajectory = list(
      rho_decreasing(n, 0.5, 0.05),  # phi11: strong ramp
      rho_constant(n, 0.1),          # phi12: constant (fixed by the mask)
      rho_constant(n, 0.1),          # phi21: constant (fixed by the mask)
      rho_constant(n, 0.3)           # phi22: constant
    ),
    seed = 13
  )
  fit <- dcvar(sim$Y_df, vars = c("y1", "y2"),
               tv_phi = "ar",
               standardize = FALSE,
               chains = 2,
               iter_warmup = margin_iter_warmup,
               iter_sampling = margin_iter_sampling,
               refresh = 0, seed = 17)

  phi_df <- phi_trajectory(fit)
  rng <- tapply(phi_df$mean, phi_df$coefficient, function(x) diff(range(x)))
  # phi11 (the truly varying coefficient) must move far more than phi22.
  expect_gt(rng[["phi11"]], rng[["phi22"]])
  # And phi11's recovered path must decrease (it was simulated decreasing).
  phi11 <- phi_df[phi_df$coefficient == "phi11", ]
  n_eff <- nrow(phi11)
  expect_gt(mean(phi11$mean[1:50]), mean(phi11$mean[(n_eff - 49):n_eff]))
})

test_that("soft-barrier recovers a time-varying exponential scale ramp", {
  skip_if_no_rstan()
  skip_if_not_slow()

  n <- 200
  scale_path <- seq(0.6, 1.8, length.out = n - 1)  # rising exp scale
  sim <- simulate_dcvar(
    n, rho_constant(n, 0.4),
    margins = c("normal", "exponential"),
    skew_direction = c(1, 1),
    sigma_trajectory = cbind(rep(1, n - 1), scale_path),
    tv_sigma_k = 8,
    seed = 23
  )
  fit <- dcvar(sim$Y_df, vars = c("y1", "y2"),
               margins = c("normal", "exponential"),
               skew_direction = c(1, 1),
               tv_sigma = TRUE, tv_sigma_k = 8,
               standardize = FALSE,
               chains = 2,
               iter_warmup = margin_iter_warmup,
               iter_sampling = margin_iter_sampling,
               refresh = 0, seed = 29)

  sigma_df <- sigma_trajectory(fit)
  exp_dim <- sigma_df[sigma_df$variable == "y2", ]
  n_eff <- nrow(exp_dim)
  # The recovered scale path must rise (matching the simulated ramp).
  expect_gt(mean(exp_dim$mean[(n_eff - 49):n_eff]), mean(exp_dim$mean[1:50]))
})
