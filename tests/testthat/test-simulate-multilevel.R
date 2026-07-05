test_that("simulate_dcvar_multilevel returns correct structure", {
  set.seed(42)
  sim <- simulate_dcvar_multilevel(N = 5, n_time = 20, rho = 0.5, seed = 1)

  expect_type(sim, "list")
  expect_named(sim, c("data", "true_params", "person_means"))
  expect_s3_class(sim$data, "data.frame")
  expect_true(all(c("id", "time", "y1", "y2") %in% names(sim$data)))
})

test_that("simulate_dcvar_multilevel data has correct dimensions", {
  N <- 5
  n_time_val <- 20
  sim <- simulate_dcvar_multilevel(N = N, n_time = n_time_val, rho = 0.3, seed = 2)

  expect_equal(nrow(sim$data), N * n_time_val)
  expect_equal(ncol(sim$data), 4)
  expect_equal(length(unique(sim$data$id)), N)
})

test_that("simulate_dcvar_multilevel true_params contains expected elements", {
  sim <- simulate_dcvar_multilevel(N = 3, n_time = 15, rho = 0.4, seed = 3)

  expect_named(sim$true_params,
               c("phi_bar", "tau_phi", "sigma", "rho", "margins",
                 "skew_direction", "skew_params", "Phi_mat", "Phi_list",
                 "Phi_population", "Phi_unit_paths"))
  expect_equal(sim$true_params$margins, "normal")
  expect_equal(sim$true_params$rho, 0.4)
  expect_equal(length(sim$true_params$phi_bar), 4)
  expect_equal(length(sim$true_params$tau_phi), 4)
  expect_equal(length(sim$true_params$sigma), 2)
  expect_equal(nrow(sim$true_params$Phi_mat), 3)
  expect_equal(ncol(sim$true_params$Phi_mat), 4)
  expect_equal(length(sim$true_params$Phi_list), 3)
})

test_that("simulate_dcvar_multilevel person_means has correct shape", {
  N <- 8
  sim <- simulate_dcvar_multilevel(N = N, n_time = 25, rho = 0.5, seed = 4)

  expect_equal(nrow(sim$person_means), N)
  expect_equal(ncol(sim$person_means), 2)
})

test_that("simulate_dcvar_multilevel works with different N and T values", {
  sim_small <- simulate_dcvar_multilevel(N = 2, n_time = 10, rho = 0.5, seed = 5)
  expect_equal(nrow(sim_small$data), 2 * 10)

  sim_large <- simulate_dcvar_multilevel(N = 10, n_time = 50, rho = 0.5, seed = 6)
  expect_equal(nrow(sim_large$data), 10 * 50)
})

test_that("simulate_dcvar_multilevel validates rho bounds", {
  expect_error(
    simulate_dcvar_multilevel(N = 2, n_time = 10, rho = 1.2, seed = 7),
    "must be a single finite numeric value in \\[-1, 1\\]"
  )
})

test_that("simulate_dcvar_multilevel validates vector lengths and center flag", {
  expect_error(
    simulate_dcvar_multilevel(tau_phi = c(0.1, 0.2, 0.3), seed = 1),
    "tau_phi"
  )
  expect_error(
    simulate_dcvar_multilevel(sigma = c(1, 1, 1), seed = 1),
    "sigma"
  )
  expect_error(
    simulate_dcvar_multilevel(center = NA, seed = 1),
    "center"
  )
})

test_that("simulate_dcvar_multilevel uses recorded sigma for skew-normal and gamma innovations", {
  skip_if_not_installed("sn")

  sim_sn <- simulate_dcvar_multilevel(
    N = 1,
    n_time = 8000,
    burnin = 0,
    phi_bar = c(0, 0, 0, 0),
    tau_phi = c(0, 0, 0, 0),
    margins = "skew_normal",
    sigma = c(2, 2.5),
    skew_params = list(alpha = c(0, 0)),
    rho = 0,
    center = FALSE,
    seed = 31
  )
  y_sn <- as.matrix(sim_sn$data[, c("y1", "y2")])
  expect_equal(sim_sn$true_params$sigma, c(2, 2.5))
  expect_equal(unname(apply(y_sn[-1, , drop = FALSE], 2, stats::sd)), sim_sn$true_params$sigma, tolerance = 0.06)

  sim_gam <- simulate_dcvar_multilevel(
    N = 1,
    n_time = 8000,
    burnin = 0,
    phi_bar = c(0, 0, 0, 0),
    tau_phi = c(0, 0, 0, 0),
    margins = "gamma",
    sigma = c(2, 2.5),
    skew_direction = c(1, 1),
    skew_params = list(shape = 2),
    rho = 0,
    center = FALSE,
    seed = 32
  )
  y_gam <- as.matrix(sim_gam$data[, c("y1", "y2")])
  expect_equal(sim_gam$true_params$sigma, c(2, 2.5))
  expect_equal(unname(apply(y_gam[-1, , drop = FALSE], 2, stats::sd)), sim_gam$true_params$sigma, tolerance = 0.08)
})

test_that("simulate_dcvar_multilevel preserves sampled nonstationary Phi matrices", {
  sim <- simulate_dcvar_multilevel(
    N = 1,
    n_time = 20,
    phi_bar = c(1.2, 0, 0, 1.1),
    tau_phi = c(0, 0, 0, 0),
    rho = 0.5,
    seed = 8
  )

  ev <- eigen(sim$true_params$Phi_list[[1]], only.values = TRUE)$values
  expect_true(any(Mod(ev) >= 1))
  expect_equal(
    sim$true_params$Phi_list[[1]],
    matrix(c(1.2, 0, 0, 1.1), 2, 2, byrow = TRUE)
  )
})

test_that("simulate_dcvar_multilevel is reproducible with seed", {
  s1 <- simulate_dcvar_multilevel(N = 3, n_time = 15, rho = 0.5, seed = 99)
  s2 <- simulate_dcvar_multilevel(N = 3, n_time = 15, rho = 0.5, seed = 99)
  expect_identical(s1$data, s2$data)
  expect_identical(s1$true_params$Phi_mat, s2$true_params$Phi_mat)
})

test_that("simulate_dcvar_multilevel centering works", {
  sim_centered <- simulate_dcvar_multilevel(N = 3, n_time = 30, rho = 0.5,
                                            center = TRUE, seed = 10)
  sim_uncentered <- simulate_dcvar_multilevel(N = 3, n_time = 30, rho = 0.5,
                                              center = FALSE, seed = 10)

  # Centered data should have near-zero person means
  for (i in unique(sim_centered$data$id)) {
    sub <- sim_centered$data[sim_centered$data$id == i, ]
    expect_true(abs(mean(sub$y1)) < 0.5)
    expect_true(abs(mean(sub$y2)) < 0.5)
  }

  # Uncentered data differs from centered data
  expect_false(identical(sim_centered$data$y1, sim_uncentered$data$y1))
})
