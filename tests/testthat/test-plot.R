test_that("plot_rho() returns ggplot for dcvar", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  p <- plot_rho(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_rho() returns ggplot for hmm", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  p <- plot_rho(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_rho() returns ggplot for constant", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  p <- plot_rho(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_phi() returns ggplot", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  p <- plot_phi(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_diagnostics() works for dcvar", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  p <- plot_diagnostics(fit)
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("plot_diagnostics() works for hmm", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  p <- plot_diagnostics(fit)
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("plot_diagnostics() works for constant", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  p <- plot_diagnostics(fit)
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("plot_hmm_states() returns ggplot with expected structure", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  p <- plot_hmm_states(fit)
  expect_s3_class(p, "ggplot")
  # Verify the plot has data
  expect_true(nrow(p$data) > 0)
})

test_that("plot_ppc() returns ggplot for dcvar", {
  skip_if_no_cmdstanr()

  fit <- get_dcvar_fit()
  expect_no_warning(p <- plot_ppc(fit))
  expect_s3_class(p, "ggplot")
  expect_false(all(is.na(p$data$cor)))
})

test_that("plot_ppc() returns ggplot for hmm", {
  skip_if_no_cmdstanr()

  fit <- get_hmm_fit()
  expect_no_warning(p <- plot_ppc(fit))
  expect_s3_class(p, "ggplot")
  expect_false(all(is.na(p$data$cor)))
})

test_that("plot_ppc() returns ggplot for constant", {
  skip_if_no_cmdstanr()

  fit <- get_constant_fit()
  expect_no_warning(p <- plot_ppc(fit))
  expect_s3_class(p, "ggplot")
  expect_false(all(is.na(p$data$cor)))
})

test_that("plot_ppc() works for exponential fits", {
  skip_if_no_cmdstanr()

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
  skip_if_no_cmdstanr()

  expect_error(plot_ppc(get_constant_gamma_fit()), "not supported")

  skip_if_not_installed("sn")
  expect_error(plot_ppc(get_constant_skew_normal_fit()), "not supported")
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
  skip_if_no_cmdstanr()
  skip_if_not_installed("ggplot2")
  fit <- get_constant_fit()
  p <- plot(fit, type = "rho")
  expect_s3_class(p, "gg")
})
