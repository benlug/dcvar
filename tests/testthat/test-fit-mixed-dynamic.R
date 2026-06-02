# Per-variable (mixed) margins for the dynamic (DC-VAR) and HMM models (Phase 2).

# --- routing (fast, no fit) -------------------------------------------------

test_that("mixed margins route to the generic dynamic and HMM Stan models", {
  expect_equal(.margin_stan_file("dcvar", c("normal", "exponential")),
               "dcvar_mixed_ncp.stan")
  expect_equal(.margin_stan_file("hmm", c("normal", "gamma")),
               "hmm_mixed.stan")
  # Identical-family vectors still route to the specialised models.
  expect_equal(.margin_stan_file("dcvar", c("normal", "normal")),
               "dcvar_model_ncp.stan")
  expect_equal(.margin_stan_file("hmm", c("gamma", "gamma")),
               "hmm_GG.stan")
  # Still unsupported elsewhere.
  expect_error(.margin_stan_file("multilevel", c("normal", "exponential")),
               "constant.*dcvar.*hmm")
})

test_that("prepare_dcvar_data / prepare_hmm_data build the family array for mixed", {
  df <- data.frame(time = 1:40, y1 = rnorm(40), y2 = cumsum(rnorm(40)))
  sd_dc <- prepare_dcvar_data(df, vars = c("y1", "y2"),
                              margins = c("normal", "exponential"),
                              skew_direction = c(1, 1))
  expect_equal(sd_dc$family, c(1L, 2L))
  expect_false(is.null(sd_dc$sigma_eps_prior))
  expect_equal(sd_dc$skew_direction, c(1, 1))

  sd_hmm <- prepare_hmm_data(df, vars = c("y1", "y2"), K = 2,
                             margins = c("normal", "skew_normal"))
  expect_equal(sd_hmm$family, c(1L, 3L))
  # skew_direction defaults to +1 when no exp/gamma dimension is present.
  expect_equal(sd_hmm$skew_direction, c(1, 1))
})

# --- DC-VAR mixed fit -------------------------------------------------------

test_that("dcvar() fits mixed normal + exponential margins", {
  skip_if_no_rstan()

  fit <- get_dcvar_mixed_fit()
  expect_s3_class(fit, "dcvar_fit")
  expect_equal(fit$margins, c("normal", "exponential"))
  expect_equal(fit$stan_data$family, c(1L, 2L))
})

test_that("dcvar mixed fit emits only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_dcvar_mixed_fit_warnings(), "dcvar mixed")
})

test_that("dcvar mixed coef() and var_params() report per-dimension families", {
  skip_if_no_rstan()

  fit <- get_dcvar_mixed_fit()
  co <- coef(fit)
  expect_equal(names(co$sigma_eps), "sigma_eps[1]")
  expect_equal(names(co$sigma_exp), "sigma_exp[2]")

  vp <- var_params(fit)
  expect_true(all(c("sigma_eps", "sigma_exp", "sigma_omega") %in% names(vp)))
})

test_that("dcvar mixed recovers a declining rho trajectory", {
  skip_if_no_rstan()

  rt <- rho_trajectory(get_dcvar_mixed_fit())
  expect_equal(nrow(rt), 119)
  expect_true(all(rt$mean > -1 & rt$mean < 1))
  # True trajectory is decreasing: the early correlation should exceed the late.
  expect_true(mean(rt$mean[1:30]) > mean(rt$mean[90:119]))
})

test_that("dcvar mixed pit_values dispatches per dimension", {
  skip_if_no_rstan()

  pit_df <- pit_values(get_dcvar_mixed_fit())
  expect_setequal(unique(pit_df$variable), c("y1", "y2"))
  expect_true(all(pit_df$pit >= 0 & pit_df$pit <= 1))
})

# --- HMM mixed fit ----------------------------------------------------------

test_that("dcvar_hmm() fits mixed normal + exponential margins", {
  skip_if_no_rstan()

  fit <- get_hmm_mixed_fit()
  expect_s3_class(fit, "dcvar_hmm_fit")
  expect_equal(fit$margins, c("normal", "exponential"))
  expect_equal(fit$stan_data$family, c(1L, 2L))
})

test_that("hmm mixed fit emits only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_hmm_mixed_fit_warnings(), "hmm mixed")
})

test_that("hmm mixed coef() reports per-dimension families", {
  skip_if_no_rstan()

  co <- coef(get_hmm_mixed_fit())
  expect_equal(names(co$sigma_eps), "sigma_eps[1]")
  expect_equal(names(co$sigma_exp), "sigma_exp[2]")
})

test_that("hmm mixed recovers two separated correlation regimes", {
  skip_if_no_rstan()

  st <- hmm_states(get_hmm_mixed_fit())
  # Ordered states; the true regimes are rho ~ 0.1 and rho ~ 0.8.
  expect_true(st$rho_state$mean[1] < st$rho_state$mean[2])
  expect_true(st$rho_state$mean[1] < 0.5)
  expect_true(st$rho_state$mean[2] > 0.4)
  expect_equal(dim(st$gamma), c(119L, 2L))
})
