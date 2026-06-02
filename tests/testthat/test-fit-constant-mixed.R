# Per-variable (mixed) margins for the constant copula model (Phase 1).

test_that("dcvar_constant() fits mixed normal + exponential margins", {
  skip_if_no_rstan()

  fit <- get_constant_mixed_fit()
  expect_s3_class(fit, "dcvar_constant_fit")
  expect_equal(fit$margins, c("normal", "exponential"))
  # The mixed Stan model receives a per-dimension family array.
  expect_equal(fit$stan_data$family, c(1L, 2L))
})

test_that("mixed fit emits only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_constant_mixed_fit_warnings(), "constant mixed")
})

test_that("coef() reports each dimension under its own family", {
  skip_if_no_rstan()

  co <- coef(get_constant_mixed_fit())
  expect_named(co, c("mu", "Phi", "sigma_eps", "sigma_exp", "rho"))
  # Normal dim reports sigma_eps[1]; exponential dim reports sigma_exp[2].
  expect_equal(names(co$sigma_eps), "sigma_eps[1]")
  expect_equal(names(co$sigma_exp), "sigma_exp[2]")
})

test_that("var_params() returns per-dimension scale groups for mixed fits", {
  skip_if_no_rstan()

  vp <- var_params(get_constant_mixed_fit())
  expect_true(all(c("mu", "Phi", "sigma_eps", "sigma_exp") %in% names(vp)))
  expect_equal(vp$sigma_eps$variable, "sigma_eps[1]")
  expect_equal(vp$sigma_exp$variable, "sigma_exp[2]")
})

test_that("summary() and print() run for mixed fits", {
  skip_if_no_rstan()

  fit <- get_constant_mixed_fit()
  expect_no_error(capture.output(print(fit)))
  s <- summary(fit)
  out <- capture.output(print(s))
  expect_true(any(grepl("sigma_eps", out)))
  expect_true(any(grepl("sigma_exp", out)))
})

test_that("pit_values() dispatches per dimension for mixed fits", {
  skip_if_no_rstan()

  pit_df <- pit_values(get_constant_mixed_fit())
  expect_true(all(c("time", "variable", "pit") %in% names(pit_df)))
  expect_setequal(unique(pit_df$variable), c("y1", "y2"))
  expect_true(all(pit_df$pit >= 0 & pit_df$pit <= 1))
})

test_that("predict() rejects mixed margins (non-normal)", {
  skip_if_no_rstan()

  expect_error(predict(get_constant_mixed_fit()), "normal margins")
})

test_that("mixed normal + exponential recovers the constant rho", {
  skip_if_no_rstan()

  rho_df <- rho_trajectory(get_constant_mixed_fit())
  # Standardisation shrinks the recovered correlation slightly; allow a
  # generous band around the simulated rho = 0.6.
  expect_true(rho_df$mean[1] > 0.35 && rho_df$mean[1] < 0.8)
})

test_that("identical-family vectors reuse the specialised model and collapse margins", {
  skip_if_no_rstan()

  sim <- simulate_dcvar(
    n_time = 40, rho_trajectory = rho_constant(40, 0.5), seed = 42
  )
  fit <- dcvar_constant(
    sim$Y_df, vars = c("y1", "y2"),
    margins = c("normal", "normal"),
    chains = 1, iter_warmup = 75, iter_sampling = 75,
    refresh = 0, seed = 123
  )
  # Homogeneous vector is treated exactly like the scalar form.
  expect_identical(fit$margins, "normal")
  expect_null(fit$stan_data$family)
  expect_named(coef(fit), c("mu", "Phi", "sigma_eps", "rho"))
})
