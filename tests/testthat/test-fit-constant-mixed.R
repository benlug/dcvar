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


# ---------------------------------------------------------------------------
# Coverage across every ordered family pair
#
# Exercises each marginal branch of constant_mixed.stan (normal / exponential /
# skew_normal / gamma) in combination, and checks the R extraction picks up the
# right per-dimension scale/shape variables. These are short smoke fits: they
# only need the Stan variables to exist, not full convergence, so warnings are
# suppressed and recovery is asserted separately below.
# ---------------------------------------------------------------------------

.simulate_mixed_pair <- function(fams, n = 120, rho = 0.6, seed = 7) {
  args <- list(n_time = n, rho_trajectory = rho_constant(n, rho),
               margins = fams, seed = seed)
  if (any(fams %in% c("exponential", "gamma"))) args$skew_direction <- c(1, 1)
  sp <- list()
  if (any(fams == "gamma")) sp$shape <- 2
  if (any(fams == "skew_normal")) sp$alpha <- c(3, 3)
  if (length(sp) > 0L) args$skew_params <- sp
  do.call(simulate_dcvar, args)
}

.fit_mixed_pair <- function(fams, seed = 7,
                            iter_warmup = 150, iter_sampling = 150,
                            chains = 1, adapt_delta = 0.95,
                            max_treedepth = 11) {
  sim <- .simulate_mixed_pair(fams, seed = seed)
  args <- list(
    sim$Y_df, vars = c("y1", "y2"), margins = fams,
    chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    refresh = 0, seed = 123
  )
  if (any(fams %in% c("exponential", "gamma"))) args$skew_direction <- c(1, 1)
  suppressWarnings(do.call(dcvar_constant, args))
}

mixed_family_pairs <- list(
  c("normal", "gamma"),
  c("normal", "skew_normal"),
  c("exponential", "gamma"),
  c("exponential", "skew_normal"),
  c("gamma", "skew_normal")
)

for (.fams in mixed_family_pairs) {
  local({
    fams <- .fams
    label <- paste(fams, collapse = " + ")
    test_that(sprintf("mixed %s fits and reports the right per-dim families", label), {
      skip_if_no_rstan()
      if (any(fams == "skew_normal")) skip_if_not_installed("sn")

      fit <- .fit_mixed_pair(fams)
      expect_s3_class(fit, "dcvar_constant_fit")
      expect_equal(fit$margins, fams)
      expect_equal(fit$stan_data$family, unname(.family_codes[fams]))

      # coef() reports each dimension under its own family's scale/shape group.
      co <- coef(fit)
      expect_true("rho" %in% names(co))
      expected_groups <- names(.mixed_margin_report_vars(fams))
      expect_true(all(expected_groups %in% names(co)))

      # var_params() exposes the same per-dimension groups.
      vp <- var_params(fit)
      expect_true(all(expected_groups %in% names(vp)))

      # PIT runs per dimension and stays in the unit interval.
      pit_df <- pit_values(fit)
      expect_true(all(pit_df$pit >= 0 & pit_df$pit <= 1))
    })
  })
}

test_that("mixed normal + gamma recovers rho and the VAR autoregression", {
  skip_if_no_rstan()

  fams <- c("normal", "gamma")
  sim <- .simulate_mixed_pair(fams, n = 150, rho = 0.6, seed = 13)
  fit <- suppressWarnings(dcvar_constant(
    sim$Y_df, vars = c("y1", "y2"), margins = fams,
    skew_direction = c(1, 1),
    chains = 2, iter_warmup = 500, iter_sampling = 500,
    adapt_delta = 0.999, max_treedepth = 14, refresh = 0, seed = 123
  ))

  rho_mean <- rho_trajectory(fit)$mean[1]
  expect_true(rho_mean > 0.3 && rho_mean < 0.8)

  # The VAR(1) diagonal (true 0.3 on each series) should be recovered loosely
  # even after standardisation.
  phi <- coef(fit)$Phi
  expect_true(phi[["Phi[1,1]"]] > 0 && phi[["Phi[1,1]"]] < 0.6)
  expect_true(phi[["Phi[2,2]"]] > 0 && phi[["Phi[2,2]"]] < 0.6)
})
