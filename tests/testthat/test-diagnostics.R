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

test_that("dcvar_diagnostics handles SEM fits without stan_data$D", {
  object <- structure(
    list(
      fit = structure(list(), class = "mock_fit"),
      stan_data = list(T = 3),
      model = "sem",
      backend = "rstan",
      priors = list(),
      meta = list()
    ),
    class = c("dcvar_sem_fit", "dcvar_model_fit")
  )

  requested_variables <- NULL

  testthat::local_mocked_bindings(
    .fit_diagnostic_summary = function(fit, backend) {
      list(num_divergent = 0L, num_max_treedepth = 0L)
    },
    .fit_sampler_diagnostics = function(fit, backend) {
      array(
        0.9,
        dim = c(1L, 1L, 1L),
        dimnames = list(NULL, NULL, "accept_stat__")
      )
    },
    .fit_summary = function(fit, variables = NULL, backend, ...) {
      requested_variables <<- variables
      data.frame(
        variable = variables,
        rhat = rep(1.01, length(variables)),
        ess_bulk = rep(200, length(variables)),
        ess_tail = rep(150, length(variables))
      )
    },
    .package = "dcvar"
  )

  diag <- dcvar_diagnostics(object)

  expect_identical(
    requested_variables,
    c(
      "mu[1]", "mu[2]",
      "phi11", "phi12", "phi21", "phi22",
      "sigma[1]", "sigma[2]",
      "rho_raw",
      "zeta[1,1]", "zeta[1,2]",
      "zeta[2,1]", "zeta[2,2]",
      "zeta[3,1]", "zeta[3,2]"
    )
  )
  expect_equal(diag$max_rhat, 1.01)
})

test_that("dcvar_diagnostics handles multilevel fits without stan_data$D", {
  object <- structure(
    list(
      fit = structure(list(), class = "mock_fit"),
      stan_data = list(T = 4),
      model = "multilevel",
      N = 2,
      backend = "rstan",
      priors = list(),
      meta = list()
    ),
    class = c("dcvar_multilevel_fit", "dcvar_model_fit")
  )

  requested_variables <- NULL

  testthat::local_mocked_bindings(
    .fit_diagnostic_summary = function(fit, backend) {
      list(num_divergent = 0L, num_max_treedepth = 0L)
    },
    .fit_sampler_diagnostics = function(fit, backend) {
      array(
        0.95,
        dim = c(1L, 1L, 1L),
        dimnames = list(NULL, NULL, "accept_stat__")
      )
    },
    .fit_summary = function(fit, variables = NULL, backend, ...) {
      requested_variables <<- variables
      data.frame(
        variable = variables,
        rhat = rep(1.02, length(variables)),
        ess_bulk = rep(220, length(variables)),
        ess_tail = rep(170, length(variables))
      )
    },
    .package = "dcvar"
  )

  diag <- dcvar_diagnostics(object)

  expect_identical(
    requested_variables,
    c(
      "phi_bar[1]", "phi_bar[2]", "phi_bar[3]", "phi_bar[4]",
      "tau_phi[1]", "tau_phi[2]", "tau_phi[3]", "tau_phi[4]",
      "z_phi[1,1]", "z_phi[1,2]", "z_phi[1,3]", "z_phi[1,4]",
      "z_phi[2,1]", "z_phi[2,2]", "z_phi[2,3]", "z_phi[2,4]",
      "sigma[1]", "sigma[2]",
      "rho"
    )
  )
  expect_equal(diag$max_rhat, 1.02)
})

test_that("dcvar_diagnostics errors clearly when D is missing for single-level models", {
  object <- structure(
    list(
      fit = structure(list(), class = "mock_fit"),
      stan_data = list(T = 5),
      model = "dcvar",
      margins = "normal",
      backend = "rstan",
      priors = list(),
      meta = list()
    ),
    class = c("dcvar_fit", "dcvar_model_fit")
  )

  expect_error(
    dcvar_diagnostics(object),
    "requires .*D.*positive integer.*stan_data"
  )
})
