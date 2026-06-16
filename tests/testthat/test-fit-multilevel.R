# ============================================================================
# Tests for dcvar_multilevel() fit and methods
# ============================================================================

# --- Fit object structure ----------------------------------------------------

make_multilevel_mixed_engine_stub_fit <- function(margins = c("skew_normal", "gamma")) {
  margins <- match.arg(margins)
  N <- 2
  n_time <- 4
  variables <- c(
    paste0("phi_bar[", 1:4, "]"),
    paste0("tau_phi[", 1:4, "]"),
    paste0(
      "z_phi[",
      rep(seq_len(N), each = 4),
      ",",
      rep(seq_len(4), times = N),
      "]"
    ),
    paste0(
      "phi_unit[",
      rep(seq_len(N), each = 4),
      ",",
      rep(seq_len(4), times = N),
      "]"
    ),
    "rho",
    paste0(
      "log_lik[",
      rep(seq_len(N), each = n_time - 1),
      ",",
      rep(seq_len(n_time - 1), times = N),
      "]"
    )
  )
  if (identical(margins, "skew_normal")) {
    variables <- c(variables, "omega[1]", "omega[2]", "delta[1]", "delta[2]")
  } else {
    variables <- c(
      variables,
      "eta[1]", "eta[2]",
      "sigma_gam[1]", "sigma_gam[2]",
      "shape_gam[1]", "shape_gam[2]"
    )
  }

  draws <- array(
    rnorm(40 * length(variables), 0, 0.1),
    dim = c(40, 1, length(variables)),
    dimnames = list(NULL, NULL, variables)
  )
  draws[, , paste0("tau_phi[", 1:4, "]")] <- runif(40 * 4, 0.05, 0.2)
  draws[, , "rho"] <- runif(40, -0.4, 0.4)
  draws[, , grep("^log_lik\\[", variables, value = TRUE)] <-
    -abs(draws[, , grep("^log_lik\\[", variables, value = TRUE)])
  if (identical(margins, "skew_normal")) {
    draws[, , c("omega[1]", "omega[2]")] <- runif(40 * 2, 0.5, 1.5)
    draws[, , c("delta[1]", "delta[2]")] <- runif(40 * 2, -0.5, 0.5)
  } else {
    draws[, , c("sigma_gam[1]", "sigma_gam[2]")] <- runif(40 * 2, 0.5, 1.5)
    draws[, , c("shape_gam[1]", "shape_gam[2]")] <- runif(40 * 2, 0.5, 2.0)
  }

  stan_data <- structure(
    list(
      N = N,
      n_time = n_time,
      y = list(matrix(rnorm(n_time * 2), n_time, 2), matrix(rnorm(n_time * 2), n_time, 2))
    ),
    ids = c("a", "b"),
    time_values = 11:14,
    margins = margins,
    skew_direction = if (identical(margins, "gamma")) c(1, -1) else NULL,
    person_means = matrix(0, N, 2)
  )

  new_dcvar_multilevel_fit(
    fit = posterior::as_draws_array(draws),
    stan_data = stan_data,
    N = N,
    vars = c("y1", "y2"),
    centered = TRUE,
    person_means = attr(stan_data, "person_means"),
    margins = margins,
    skew_direction = attr(stan_data, "skew_direction"),
    backend = "rstan",
    priors = list(),
    meta = list(chains = 1, iter_warmup = 10, iter_sampling = 40)
  )
}

test_that("dcvar_multilevel returns correct class", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  expect_s3_class(fit, "dcvar_multilevel_fit")
  expect_s3_class(fit, "dcvar_model_fit")
})

test_that("dcvar_multilevel rejects center = FALSE with bundled Stan model", {
  df <- data.frame(
    id = rep(1:2, each = 4),
    time = rep(1:4, times = 2),
    y1 = rnorm(8),
    y2 = rnorm(8)
  )

  expect_error(
    dcvar_multilevel(df, vars = c("y1", "y2"), center = FALSE),
    "not supported by the bundled multilevel model"
  )
})

test_that("dcvar_multilevel rejects explicit bundled stan_file with center = FALSE", {
  df <- data.frame(
    id = rep(1:2, each = 4),
    time = rep(1:4, times = 2),
    y1 = rnorm(8),
    y2 = rnorm(8)
  )

  expect_error(
    dcvar_multilevel(
      df,
      vars = c("y1", "y2"),
      center = FALSE,
      stan_file = dcvar_stan_path("multilevel")
    ),
    "not supported by the bundled multilevel model"
  )
})

test_that("multilevel coef() returns expected structure", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  co <- coef(fit)

  expect_type(co, "list")
  expect_named(co, c("phi_bar", "tau_phi", "sigma", "rho"))
  expect_true(length(co$phi_bar) == 4)
  expect_true(length(co$tau_phi) == 4)
  expect_true(length(co$sigma) == 2)
  expect_true(length(co$rho) == 1)
})

test_that("multilevel summary() returns correct class", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  s <- summary(fit)

  expect_s3_class(s, "dcvar_multilevel_summary")
  expect_true(all(c("var_params", "random_effects", "diagnostics") %in% names(s)))
})

test_that("multilevel print() works", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  expect_output(print(fit), "multilevel|Multilevel")
})

# --- Extraction functions ----------------------------------------------------

test_that("multilevel rho_trajectory() returns expected structure", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  rho_df <- rho_trajectory(fit)

  expect_s3_class(rho_df, "data.frame")
  expect_true(all(c("time", "mean", "sd") %in% names(rho_df)))
  # Constant rho: all means should be identical
  expect_equal(length(unique(rho_df$mean)), 1)
})

test_that("multilevel var_params() returns expected structure", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  vp <- var_params(fit)

  expect_type(vp, "list")
  expect_named(vp, c("phi_bar", "tau_phi", "sigma", "rho"))
  expect_s3_class(vp$phi_bar, "data.frame")
  expect_true(all(c("variable", "mean", "sd", "q2.5", "q97.5") %in% names(vp$phi_bar)))
})

test_that("multilevel homogeneous skew_normal and gamma stubs report mixed-engine params", {
  skew_fit <- make_multilevel_mixed_engine_stub_fit("skew_normal")
  skew_co <- coef(skew_fit)
  skew_vp <- var_params(skew_fit)
  expect_named(skew_co, c("phi_bar", "tau_phi", "omega", "delta", "rho"))
  expect_equal(names(skew_co$omega), c("omega[1]", "omega[2]"))
  expect_true(all(c("omega", "delta") %in% names(skew_vp)))

  gamma_fit <- make_multilevel_mixed_engine_stub_fit("gamma")
  gamma_co <- coef(gamma_fit)
  gamma_vp <- var_params(gamma_fit)
  expect_named(gamma_co, c("phi_bar", "tau_phi", "sigma_gam", "shape_gam", "rho"))
  expect_equal(names(gamma_co$shape_gam), c("shape_gam[1]", "shape_gam[2]"))
  expect_true(all(c("sigma_gam", "shape_gam") %in% names(gamma_vp)))
})

test_that("dcvar_multilevel fits homogeneous skew_normal and gamma margins", {
  skip_if_no_rstan()
  skip_if_not_slow()

  for (margin in c("skew_normal", "gamma")) {
    skew_direction <- if (identical(margin, "gamma")) c(1, 1) else NULL
    sim <- simulate_dcvar_multilevel(
      N = 6, n_time = 40, rho = 0.5,
      margins = margin,
      skew_direction = skew_direction,
      skew_params = if (identical(margin, "skew_normal")) list(alpha = c(2, 2)) else list(shape = 2),
      seed = 600 + match(margin, c("skew_normal", "gamma"))
    )
    fit <- dcvar_multilevel(
      sim$data, vars = c("y1", "y2"), id_var = "id",
      margins = margin,
      skew_direction = skew_direction,
      chains = 2, iter_warmup = 300, iter_sampling = 300,
      adapt_delta = 0.95, max_treedepth = 13, refresh = 0, seed = 123
    )

    expect_equal(fit$margins, margin)
    expect_equal(fit$stan_data$family, rep(unname(.family_codes[[margin]]), 2))
    expect_true(coef(fit)$rho > -0.2)
  }
})

test_that("multilevel random_effects() returns expected structure", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  re <- random_effects(fit)

  expect_s3_class(re, "data.frame")
  expect_true(all(c("unit", "parameter", "mean", "sd", "q2.5", "q97.5") %in% names(re)))
  expect_true(all(re$parameter %in% c("phi11", "phi12", "phi21", "phi22")))
  expect_equal(nrow(re), fit$N * 4)
})

# --- Plotting ----------------------------------------------------------------

test_that("multilevel plot random_effects returns ggplot", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  p <- plot(fit, type = "random_effects")
  expect_s3_class(p, "ggplot")
})

test_that("multilevel plot diagnostics returns ggplot", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  p <- plot(fit, type = "diagnostics")
  expect_s3_class(p, "ggplot")
})

# --- Fitted values and prediction --------------------------------------------

test_that("multilevel fitted() returns unit-specific trajectories", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()

  fit_link <- fitted(fit, type = "link")
  fit_resp <- fitted(fit, type = "response")

  expect_s3_class(fit_link, "data.frame")
  expect_named(fit_link, c("unit", "time", "y1", "y2"))
  expect_equal(nrow(fit_link), fit$N * (fit$stan_data$n_time - 1))
  expect_equal(length(unique(fit_link$unit)), fit$N)
  expect_false(identical(fit_link[, fit$vars], fit_resp[, fit$vars]))
})

test_that("multilevel predict() returns unit-specific intervals", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()

  pred <- predict(fit)

  expect_s3_class(pred, "data.frame")
  expect_named(pred, c("unit", "time", "variable", "mean", "lower", "upper"))
  expect_equal(nrow(pred), fit$N * (fit$stan_data$n_time - 1) * length(fit$vars))
  expect_equal(sort(unique(pred$variable)), sort(fit$vars))
  expect_true(all(pred$lower <= pred$mean))
  expect_true(all(pred$upper >= pred$mean))
})

test_that("multilevel pit_values() errors informatively", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  expect_error(pit_values(fit), "not.*supported")
})

# --- Diagnostics -------------------------------------------------------------

test_that("multilevel diagnostics are finite", {
  skip_if_no_rstan()
  fit <- get_multilevel_fit()
  diag <- dcvar_diagnostics(fit)

  expect_lte(diag$n_divergent, 1)
  expect_equal(diag$n_max_treedepth, 0)
  expect_true(is.finite(diag$max_rhat))
  expect_true(diag$max_rhat < 1.25)
  expect_true(diag$min_ess_bulk > 8)
  expect_true(diag$min_ess_tail > 8)
})

test_that("multilevel fit cache emits only known diagnostic warnings", {
  skip_if_no_rstan()

  expect_known_fit_warnings(get_multilevel_fit_warnings(), "multilevel")
})
