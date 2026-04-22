test_that("plot_rho() returns ggplot for dcvar", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
  p <- plot_rho(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_rho() returns ggplot for hmm", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  p <- plot_rho(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_rho() returns ggplot for constant", {
  skip_if_no_rstan()

  fit <- get_constant_fit()
  p <- plot_rho(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_rho() rejects invalid interval levels", {
  fit <- structure(
    list(model = "dcvar", stan_data = list(n_time = 5)),
    class = "dcvar_fit"
  )
  expect_error(plot_rho(fit, ci_level = 0), "ci_level")
  expect_error(plot_rho(fit, ci_level = 0.9, inner_level = 0.9), "inner_level")
  expect_error(plot_rho(fit, ci_level = 0.9, inner_level = 1), "inner_level")
})

test_that("plot_phi() returns ggplot", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
  p <- plot_phi(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_diagnostics() works for dcvar", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
  p <- plot_diagnostics(fit)
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("plot_diagnostics() works for hmm", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  p <- plot_diagnostics(fit)
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("plot_diagnostics() works for constant", {
  skip_if_no_rstan()

  fit <- get_constant_fit()
  p <- plot_diagnostics(fit)
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("plot_hmm_states() returns ggplot with expected structure", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  p <- plot_hmm_states(fit)
  expect_s3_class(p, "ggplot")
  # Verify the plot has data
  expect_true(nrow(p$data) > 0)
})

test_that("plot_rho() uses preserved time values for true_rho overlays", {
  fit <- structure(
    list(model = "dcvar", stan_data = list(n_time = 5)),
    class = "dcvar_fit"
  )
  time_values <- seq.Date(
    as.Date("2021-01-01"),
    by = "day",
    length.out = fit$stan_data$n_time
  )
  fit$stan_data <- structure(fit$stan_data, time_values = time_values)

  testthat::local_mocked_bindings(
    rho_trajectory = function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
      time_axis <- .observed_time_values(object$stan_data, drop_first = TRUE)
      data.frame(
        time = time_axis,
        mean = rep(0.2, length(time_axis)),
        sd = rep(0.05, length(time_axis)),
        q2.5 = rep(0.1, length(time_axis)),
        q10 = rep(0.15, length(time_axis)),
        q50 = rep(0.2, length(time_axis)),
        q90 = rep(0.25, length(time_axis)),
        q97.5 = rep(0.3, length(time_axis))
      )
    },
    .package = "dcvar"
  )

  true_rho <- rep(0.2, fit$stan_data$n_time - 1)
  p <- plot_rho(fit, true_rho = true_rho)
  overlay_data <- p$layers[[length(p$layers)]]$data

  expect_equal(overlay_data$time, time_values[-1])
  expect_equal(overlay_data$rho, true_rho)
})

test_that("plot_latent_states() uses preserved time values for true_states overlays", {
  fit <- structure(
    list(
      vars = c("latent1", "latent2"),
      stan_data = list(n_time = 5)
    ),
    class = "dcvar_sem_fit"
  )
  time_values <- seq.Date(
    as.Date("2021-01-01"),
    by = "day",
    length.out = fit$stan_data$n_time
  )
  fit$stan_data <- structure(fit$stan_data, time_values = time_values)

  testthat::local_mocked_bindings(
    latent_states = function(object, probs = c(0.025, 0.5, 0.975), ...) {
      data.frame(
        time = rep(time_values, 2),
        variable = rep(object$vars, each = length(time_values)),
        mean = rep(0, 2 * length(time_values)),
        sd = rep(1, 2 * length(time_values)),
        q2.5 = rep(-1, 2 * length(time_values)),
        q50 = rep(0, 2 * length(time_values)),
        q97.5 = rep(1, 2 * length(time_values))
      )
    },
    .package = "dcvar"
  )

  true_states <- cbind(
    latent1 = seq_len(fit$stan_data$n_time),
    latent2 = seq_len(fit$stan_data$n_time) * 2
  )
  p <- plot_latent_states(fit, true_states = true_states)
  overlay_data <- p$layers[[length(p$layers)]]$data

  expect_equal(overlay_data$time, rep(time_values, 2))
  expect_equal(overlay_data$value, c(true_states[, 1], true_states[, 2]))
})

test_that("plot_ppc() returns ggplot for dcvar", {
  skip_if_no_rstan()

  fit <- get_dcvar_fit()
  expect_no_warning(p <- plot_ppc(fit))
  expect_s3_class(p, "ggplot")
  expect_false(all(is.na(p$data$cor)))
})

test_that("plot_ppc() returns ggplot for hmm", {
  skip_if_no_rstan()

  fit <- get_hmm_fit()
  expect_no_warning(p <- plot_ppc(fit))
  expect_s3_class(p, "ggplot")
  expect_false(all(is.na(p$data$cor)))
})

test_that("plot_ppc() returns ggplot for constant", {
  skip_if_no_rstan()

  fit <- get_constant_fit()
  expect_no_warning(p <- plot_ppc(fit))
  expect_s3_class(p, "ggplot")
  expect_false(all(is.na(p$data$cor)))
})

test_that("plot_ppc() works for exponential fits", {
  skip_if_no_rstan()

  fit <- get_dcvar_exponential_fit()
  expect_no_warning(p <- plot_ppc(fit))
  expect_s3_class(p, "ggplot")
  expect_false(all(is.na(p$data$cor)))
})

test_that("plot_ppc() rejects unsupported non-normal margins", {
  gamma_fit <- structure(list(margins = "gamma"), class = "dcvar_fit")
  skew_fit <- structure(list(margins = "skew_normal"), class = "dcvar_fit")

  expect_error(plot_ppc(gamma_fit), "not supported")
  expect_error(plot_ppc(skew_fit), "not supported")
})

test_that("plot_ppc() rejects unsupported fitted gamma and skew-normal objects", {
  skip_if_no_rstan()

  expect_error(plot_ppc(get_constant_gamma_fit()), "not supported")

  skip_if_not_installed("sn")
  expect_error(plot_ppc(get_constant_skew_normal_fit()), "not supported")
})

test_that("plot methods fail clearly when required Stan outputs are missing", {
  make_stub_draws <- function(variables) {
    posterior::as_draws_array(
      array(
        seq_len(4L * length(variables)),
        dim = c(4L, 1L, length(variables)),
        dimnames = list(NULL, NULL, variables)
      )
    )
  }

  phi_fit <- structure(
    list(
      fit = make_stub_draws(c("mu[1]", "mu[2]")),
      model = "dcvar",
      vars = c("y1", "y2"),
      margins = "normal",
      stan_data = list(n_time = 5, D = 2),
      backend = "rstan",
      meta = list(iter_sampling = 10, chains = 1)
    ),
    class = c("dcvar_fit", "dcvar_model_fit")
  )
  expect_error(plot_phi(phi_fit), "Custom Stan files must preserve")

  ppc_fit <- structure(
    list(
      fit = make_stub_draws(c("eps[1,1]", "eps[1,2]", "eps[2,1]", "eps[2,2]")),
      model = "dcvar",
      vars = c("y1", "y2"),
      margins = "normal",
      stan_data = list(n_time = 3, D = 2),
      backend = "rstan",
      meta = list(iter_sampling = 10, chains = 1)
    ),
    class = c("dcvar_fit", "dcvar_model_fit")
  )
  expect_error(plot_ppc(ppc_fit), "Custom Stan files must preserve")
})

test_that("plot_trajectories() returns ggplot", {
  p <- plot_trajectories(50)
  expect_s3_class(p, "ggplot")
  expect_true(nrow(p$data) > 0)
})

test_that("plot_trajectories() works with subset of scenarios", {
  p <- plot_trajectories(50, scenarios = c("constant", "decreasing"))
  expect_s3_class(p, "ggplot")
})

test_that("plot(constant_fit, type = 'rho') works", {
  skip_if_no_rstan()
  skip_if_not_installed("ggplot2")
  fit <- get_constant_fit()
  p <- plot(fit, type = "rho")
  expect_s3_class(p, "gg")
})
