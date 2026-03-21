test_that("dcvar_diagnostics scopes headline convergence to sampled parameters", {
  object <- structure(
    list(
      fit = structure(list(), class = "mock_fit"),
      stan_data = list(T = 5, D = 2),
      model = "dcvar",
      vars = c("y1", "y2"),
      standardized = TRUE,
      margins = "normal",
      backend = "rstan",
      priors = list(),
      meta = list()
    ),
    class = c("dcvar_fit", "dcvar_model_fit")
  )

  requested_variables <- NULL

  testthat::local_mocked_bindings(
    .fit_diagnostic_summary = function(fit, backend) {
      list(num_divergent = 1L, num_max_treedepth = 2L)
    },
    .fit_sampler_diagnostics = function(fit, backend) {
      array(
        c(0.9, 0.8),
        dim = c(1L, 2L, 1L),
        dimnames = list(NULL, NULL, "accept_stat__")
      )
    },
    .fit_summary = function(fit, variables = NULL, backend, ...) {
      requested_variables <<- variables
      data.frame(
        variable = variables,
        rhat = seq(1.01, by = 0.01, length.out = length(variables)),
        ess_bulk = seq(200, by = -1, length.out = length(variables)),
        ess_tail = seq(150, by = -1, length.out = length(variables))
      )
    },
    .package = "dcvar"
  )

  diag <- dcvar_diagnostics(object)

  expect_identical(
    requested_variables,
    c(
      "mu[1]", "mu[2]",
      "Phi[1,1]", "Phi[1,2]", "Phi[2,1]", "Phi[2,2]",
      "sigma_eps[1]", "sigma_eps[2]",
      "z_rho_init", "sigma_omega",
      "omega_raw[1]", "omega_raw[2]", "omega_raw[3]", "omega_raw[4]"
    )
  )
  expect_equal(diag$n_divergent, 1L)
  expect_equal(diag$n_max_treedepth, 2L)
  expect_equal(diag$max_rhat, max(seq(1.01, by = 0.01, length.out = length(requested_variables))))
  expect_equal(diag$min_ess_bulk, min(seq(200, by = -1, length.out = length(requested_variables))))
  expect_equal(diag$min_ess_tail, min(seq(150, by = -1, length.out = length(requested_variables))))
  expect_equal(diag$mean_accept_prob, 0.85)
})
