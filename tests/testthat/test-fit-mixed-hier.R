# Per-variable (mixed) margins for the multilevel and SEM models (Phase 3).

# --- multilevel routing / prep (fast) ---------------------------------------

test_that("multilevel mixed margins route and build a family array", {
  expect_equal(.margin_stan_file("multilevel", c("normal", "exponential")),
               "multilevel_mixed.stan")
  # Homogeneous normal/exponential keep the specialised models.
  expect_equal(.margin_stan_file("multilevel", c("normal", "normal")),
               "multilevel_copula_var.stan")
  expect_equal(.margin_stan_file("multilevel", c("exponential", "exponential")),
               "multilevel_EG.stan")

  dfm <- data.frame(id = rep(1:3, each = 12), time = rep(1:12, 3),
                    y1 = rnorm(36), y2 = rnorm(36))
  sd_ml <- prepare_multilevel_data(dfm, vars = c("y1", "y2"), id_var = "id",
                                   margins = c("normal", "gamma"),
                                   skew_direction = c(1, 1))
  expect_equal(sd_ml$family, c(1L, 4L))
  expect_equal(sd_ml$skew_direction, c(1, 1))
})

test_that("single-family multilevel still rejects gamma/skew_normal", {
  dfm <- data.frame(id = rep(1:3, each = 12), time = rep(1:12, 3),
                    y1 = rnorm(36), y2 = rnorm(36))
  expect_error(
    prepare_multilevel_data(dfm, vars = c("y1", "y2"), id_var = "id",
                            margins = "gamma", skew_direction = c(1, 1)),
    "only.*normal"
  )
})

# --- multilevel mixed fit ---------------------------------------------------

test_that("dcvar_multilevel() fits mixed normal + exponential margins", {
  skip_if_no_rstan()

  fit <- get_multilevel_mixed_fit()
  expect_s3_class(fit, "dcvar_multilevel_fit")
  expect_equal(fit$margins, c("normal", "exponential"))
  expect_equal(fit$stan_data$family, c(1L, 2L))
})

test_that("multilevel mixed fit emits only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_multilevel_mixed_fit_warnings(), "multilevel mixed")
})

test_that("multilevel mixed coef()/var_params() report per-dimension families", {
  skip_if_no_rstan()

  fit <- get_multilevel_mixed_fit()
  co <- coef(fit)
  expect_equal(names(co$sigma_eps), "sigma_eps[1]")
  expect_equal(names(co$sigma_exp), "sigma_exp[2]")

  vp <- var_params(fit)
  expect_true(all(c("phi_bar", "tau_phi", "sigma_eps", "sigma_exp", "rho") %in% names(vp)))
})

test_that("multilevel mixed recovers the global copula correlation", {
  skip_if_no_rstan()

  rho <- coef(get_multilevel_mixed_fit())$rho
  expect_true(rho > 0.2 && rho < 0.8)
})


# --- SEM routing / prep (fast) ----------------------------------------------

test_that("SEM mixed margins route and build a family array", {
  expect_equal(.margin_stan_file("sem", c("normal", "exponential")),
               "sem_mixed.stan")
  expect_equal(.margin_stan_file("sem_naive", c("normal", "gamma")),
               "sem_naive_mixed.stan")
  # Homogeneous normal/exponential keep the specialised models.
  expect_equal(.margin_stan_file("sem", c("normal", "normal")),
               "sem_copula_var.stan")

  df <- data.frame(time = 1:30, y1_1 = rnorm(30), y1_2 = rnorm(30),
                   y2_1 = rnorm(30), y2_2 = rnorm(30))
  indicators <- list(latent1 = c("y1_1", "y1_2"), latent2 = c("y2_1", "y2_2"))
  sd_sem <- prepare_sem_data(df, indicators = indicators, J = 2,
                             lambda = c(0.8, 0.8), sigma_e = 0.5,
                             margins = c("normal", "exponential"),
                             skew_direction = c(1, 1))
  expect_equal(sd_sem$family, c(1L, 2L))
})

# --- SEM mixed fit ----------------------------------------------------------

test_that("dcvar_sem() fits mixed margins (indicator method)", {
  skip_if_no_rstan()

  fit <- get_sem_mixed_fit()
  expect_s3_class(fit, "dcvar_sem_fit")
  expect_equal(fit$margins, c("normal", "exponential"))
  expect_equal(fit$stan_data$family, c(1L, 2L))

  co <- coef(fit)
  expect_equal(names(co$sigma_eps), "sigma_eps[1]")
  expect_equal(names(co$sigma_exp), "sigma_exp[2]")

  vp <- var_params(fit)
  expect_true(all(c("mu", "Phi", "sigma_eps", "sigma_exp", "rho") %in% names(vp)))
})

test_that("SEM mixed fit emits only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_sem_mixed_fit_warnings(), "sem mixed")
})

test_that("SEM mixed (indicator) recovers a positive copula correlation", {
  skip_if_no_rstan()

  rho <- coef(get_sem_mixed_fit())$rho
  expect_true(rho > 0.1 && rho < 0.9)
})

test_that("dcvar_sem() fits mixed margins (naive method)", {
  skip_if_no_rstan()

  fit <- get_sem_naive_mixed_fit()
  expect_s3_class(fit, "dcvar_sem_fit")
  expect_equal(fit$method, "naive")
  expect_equal(fit$margins, c("normal", "exponential"))
  expect_equal(fit$stan_data$family, c(1L, 2L))
  expect_true(all(c("sigma_eps", "sigma_exp") %in% names(coef(fit))))
})
