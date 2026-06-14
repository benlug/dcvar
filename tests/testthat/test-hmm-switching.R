# ============================================================================
# Full Markov-switching HMM (hmm_switching.stan): routing, resolvers, data prep,
# init shapes, positional compatibility, and (gated) MCMC behaviour.
# ============================================================================

# --- routing (fast, no fit) -------------------------------------------------

test_that("the switching engine is reachable and the public HMM path is unchanged", {
  sw_path <- dcvar_stan_path("hmm_switching")
  expect_true(file.exists(sw_path))
  expect_match(sw_path, "hmm_switching\\.stan$")

  # The public default HMM path is the legacy file; the engine is internal.
  expect_match(dcvar_stan_path("hmm"), "hmm_copula_model\\.stan$")
  expect_identical(.margin_stan_file("hmm_switching", "normal"), "hmm_switching.stan")

  # The switching predicate only recognises the bundled engine file (not NULL).
  expect_false(.uses_hmm_switching_stan_file(NULL))
  expect_true(.uses_hmm_switching_stan_file(sw_path))
  expect_false(.uses_hmm_switching_stan_file(dcvar_stan_path("hmm")))
})

test_that("dcvar_hmm keeps time_var positional and adds switch after skew_direction", {
  fmls <- names(formals(dcvar_hmm))
  # Positional contract: data, vars, K, time_var, standardize, margins, skew_direction
  expect_identical(fmls[1:7],
                   c("data", "vars", "K", "time_var", "standardize", "margins", "skew_direction"))
  expect_gt(which(fmls == "switch"), which(fmls == "skew_direction"))
})

# --- switch resolver --------------------------------------------------------

test_that(".resolve_switch_spec maps selectors and enforces the rho anchor", {
  rho_only <- .resolve_switch_spec("rho")
  expect_identical(rho_only$mu, 0L)
  expect_identical(unname(rho_only$phi_mask), c(0L, 0L, 0L, 0L))
  expect_identical(rho_only$rho, 1L)
  expect_identical(rho_only$margins, 0L)

  expect_identical(.resolve_switch_spec(c("rho", "mu"))$mu, 1L)
  expect_identical(unname(.resolve_switch_spec(c("rho", "ar"))$phi_mask), c(1L, 0L, 0L, 1L))
  expect_identical(unname(.resolve_switch_spec(c("rho", "cross"))$phi_mask), c(0L, 1L, 1L, 0L))
  expect_identical(unname(.resolve_switch_spec(c("rho", "phi"))$phi_mask), c(1L, 1L, 1L, 1L))
  expect_identical(.resolve_switch_spec(c("rho", "sigma"))$margins, 1L)

  full <- .resolve_switch_spec(TRUE)
  expect_identical(full$mu, 1L)
  expect_identical(unname(full$phi_mask), c(1L, 1L, 1L, 1L))
  expect_identical(full$margins, 1L)

  # rho is mandatory whenever anything else switches
  expect_error(.resolve_switch_spec("mu"), "rho")
  expect_error(.resolve_switch_spec(c("mu", "phi")), "rho")
  # degenerate / invalid selectors
  expect_error(.resolve_switch_spec(FALSE), "at least one")
  expect_error(.resolve_switch_spec("none"), "at least one")
  expect_error(.resolve_switch_spec(character(0)), "at least one")
  expect_error(.resolve_switch_spec(c("rho", "nonsense")), "Invalid")
})

# --- per-state margin configuration -----------------------------------------

test_that(".hmm_margin_config handles global, per-state, and collapse", {
  # Global scalar
  g <- .hmm_margin_config("normal", NULL, K = 2)
  expect_false(g$per_state)
  expect_identical(g$family, matrix(1L, 2, 2))
  expect_identical(g$margins_global, "normal")

  # Per-state with differing families
  ps <- .hmm_margin_config(list(c("normal", "normal"), c("exponential", "gamma")),
                           skew_direction = c(1, 1), K = 2)
  expect_true(ps$per_state)
  expect_true(ps$per_state_differ)
  # state 1 = (normal, normal) -> (1, 1); state 2 = (exponential, gamma) -> (2, 4)
  expect_identical(ps$family, matrix(c(1L, 1L, 2L, 4L), nrow = 2, byrow = TRUE))

  # All-identical per-state list collapses to the global form (legacy path)
  col <- .hmm_margin_config(list("normal", "normal"), NULL, K = 2)
  expect_false(col$per_state)
  expect_false(col$per_state_differ)

  col_exp <- .hmm_margin_config(list("exponential", "exponential"),
                                skew_direction = list(c(1, -1), c(1, -1)),
                                K = 2)
  expect_false(col_exp$per_state)
  expect_identical(col_exp$skew_global, c(1L, -1L))

  # Differing consulted skew orientation is itself state-specific and needs the engine.
  skew_diff <- .hmm_margin_config(list("exponential", "exponential"),
                                  skew_direction = list(c(1, 1), c(-1, -1)),
                                  K = 2)
  expect_true(skew_diff$per_state)
  expect_true(skew_diff$per_state_differ)

  # Length mismatch and missing skew error (naming the state)
  expect_error(.hmm_margin_config(list("normal"), NULL, K = 2), "length-2 list")
  expect_error(
    .hmm_margin_config(list(c("normal", "normal"), c("exponential", "normal")), NULL, K = 2),
    "State 2"
  )
})

test_that(".as_hmm_switching_stan_data builds the engine block; init has the right shapes", {
  df <- data.frame(time = 1:30, y1 = rnorm(30), y2 = rnorm(30))
  base <- prepare_hmm_data(df, c("y1", "y2"), K = 2)
  sw <- .resolve_switch_spec(c("rho", "mu", "phi"))
  cfg <- .hmm_margin_config(list(c("normal", "normal"), c("exponential", "exponential")),
                            skew_direction = c(1, 1), K = 2)
  sw$margins <- 1L
  aug <- .as_hmm_switching_stan_data(base, sw, cfg$family, cfg$skew, prior_phi_dev_sd = 0.5)
  expect_identical(dim(aug$family), c(2L, 2L))
  expect_identical(dim(aug$skew_direction), c(2L, 2L))
  expect_identical(aug$switch_mu, 1L)
  expect_identical(aug$switch_phi, 1L)
  expect_identical(aug$phi_switch_mask, c(1L, 1L, 1L, 1L))
  expect_identical(aug$switch_margins, 1L)
  expect_true(is.numeric(aug$prior_phi_dev_sd))

  # init: state-indexed switching groups; Phi_dev present (phi switches)
  i1 <- .init_hmm_switching_params(2, 2, sw)
  expect_identical(dim(i1$mu), c(2L, 2L))
  expect_identical(dim(i1$sigma_eps), c(2L, 2L))
  expect_identical(dim(i1$Phi_dev), c(2L, 2L, 2L))
  expect_length(i1$z_rho, 2L)

  # rho-only: shared groups length 1, no Phi_dev
  i0 <- .init_hmm_switching_params(2, 2, .resolve_switch_spec("rho"))
  expect_identical(dim(i0$mu), c(1L, 2L))
  expect_identical(dim(i0$sigma_eps), c(1L, 2L))
  expect_null(i0$Phi_dev)
})

test_that("the engine config / init helpers handle K > 2 with per-state families", {
  cfg <- .hmm_margin_config(
    list(c("normal", "normal"), c("exponential", "exponential"), c("gamma", "gamma")),
    skew_direction = c(1, 1), K = 3
  )
  expect_true(cfg$per_state)
  expect_identical(dim(cfg$family), c(3L, 2L))
  expect_identical(cfg$family[3, ], c(4L, 4L))

  sw <- .resolve_switch_spec(c("rho", "mu", "phi"))
  sw$margins <- 1L # the wrapper forces this for differing per-state families
  init <- .init_hmm_switching_params(2, 3, sw)
  expect_identical(dim(init$mu), c(3L, 2L))
  expect_identical(dim(init$sigma_eps), c(3L, 2L))
  expect_identical(dim(init$Phi_dev), c(3L, 2L, 2L))
  expect_length(init$z_rho, 3L)

  df <- data.frame(time = 1:30, y1 = rnorm(30), y2 = rnorm(30))
  base <- prepare_hmm_data(df, c("y1", "y2"), K = 3)
  aug <- .as_hmm_switching_stan_data(base, sw, cfg$family, cfg$skew, prior_phi_dev_sd = 0.5)
  expect_identical(dim(aug$family), c(3L, 2L))
})

# --- MCMC behaviour (gated) -------------------------------------------------

test_that("a state-specific HMM fit exposes per-state parameters and pins diagnostics", {
  skip_if_no_rstan()

  fit <- get_hmm_switching_fit()
  expect_s3_class(fit, "dcvar_hmm_fit")
  expect_true(isTRUE(fit$switching))

  sp <- hmm_state_params(fit)
  expect_identical(dim(sp$mu), c(2L, 2L))
  expect_length(sp$Phi, 2L)
  expect_length(sp$rho_state, 2L)
  expect_identical(dim(sp$families), c(2L, 2L))
  expect_true(all(seq_len(2) %in% sp$scales$state))

  # Name-pin: every monitored name exists in the draws.
  monitored <- .diagnostic_parameter_variables(fit)
  available <- posterior::variables(draws(fit))
  expect_identical(setdiff(monitored, available), character(0))

  # Per-state coef + var_params resolve (state-indexed Phi baseline + deviations).
  co <- coef(fit)
  expect_false(is.null(co$mu))
  expect_false(is.null(co$Phi))
  vp <- var_params(fit)
  expect_false(is.null(vp$Phi))

  # fitted() returns the gamma-weighted one-step prediction.
  ft <- fitted(fit)
  expect_equal(nrow(ft), fit$stan_data$n_time - 1L)
  expect_true(all(is.finite(ft[[2]])))

  # Mixture predictive / PIT / PPC are guarded for switching fits.
  expect_error(predict(fit), "state-specific")
  expect_error(pit_values(fit), "state-specific")
  expect_error(plot_phi(fit), "state-specific")
  expect_error(plot_ppc(fit), "state-specific")
  expect_error(plot(fit, type = "ppc"), "state-specific")
  expect_s3_class(plot(fit, type = "diagnostics"), "ggplot")
})

test_that("default dcvar_hmm is unchanged (legacy path, not the engine)", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  expect_false(isTRUE(fit$switching))
  co <- coef(fit)
  expect_false(is.null(co$sigma_eps))
  # Legacy downstream methods keep working.
  expect_s3_class(fitted(fit), "data.frame")
  expect_true(all(predict(fit)$lower <= predict(fit)$upper))
})

test_that("the switching engine fits under cmdstanr with per-state margins", {
  skip_if_no_cmdstanr_backend()
  skip_if_not_slow()

  sim <- simulate_dcvar(n_time = 45, rho_trajectory = rho_step(45, 0.8, 0.1), seed = 42)
  fit <- dcvar_hmm(
    sim$Y_df, vars = c("y1", "y2"), K = 2,
    switch = c("rho", "mu", "phi"),
    margins = list(c("normal", "normal"), c("exponential", "exponential")),
    skew_direction = c(1, 1),
    chains = 1, iter_warmup = 60, iter_sampling = 60, refresh = 0,
    seed = 13, backend = "cmdstanr"
  )
  expect_s3_class(fit, "dcvar_hmm_fit")
  expect_identical(fit$backend, "cmdstanr")
  expect_true(isTRUE(fit$switching))
  # Name-pin holds on the cmdstanr serialization too (state-indexed names).
  expect_identical(
    setdiff(.diagnostic_parameter_variables(fit), posterior::variables(draws(fit))),
    character(0)
  )
})

test_that("the switching engine reproduces the legacy HMM on a rho-only config", {
  skip_if_no_rstan()
  skip_if_not_slow()

  sim <- simulate_dcvar(
    n_time = 60,
    rho_trajectory = rho_step(60, rho_before = 0.8, rho_after = 0.1),
    seed = 42
  )
  common <- list(
    data = sim$Y_df, vars = c("y1", "y2"), K = 2,
    chains = 2, iter_warmup = 500, iter_sampling = 500,
    adapt_delta = 0.99, max_treedepth = 13, refresh = 0, seed = 321
  )

  fit_legacy <- do.call(dcvar_hmm, common)
  # Force the engine on the same rho-only config via the bundled switching file.
  fit_engine <- do.call(dcvar_hmm, c(common, list(stan_file = dcvar_stan_path("hmm_switching"))))
  expect_true(isTRUE(fit_engine$switching))

  cl <- coef(fit_legacy)
  ce <- coef(fit_engine)
  expect_equal(unname(cl$mu), unname(ce$mu), tolerance = 0.08)
  expect_equal(sort(unname(cl$rho_state)), sort(unname(ce$rho_state)), tolerance = 0.08)
  expect_equal(unname(cl$sigma_eps), unname(ce$sigma_eps), tolerance = 0.08)

  rl <- rho_trajectory(fit_legacy)$mean
  re <- rho_trajectory(fit_engine)$mean
  expect_equal(rl, re, tolerance = 0.08)
})
