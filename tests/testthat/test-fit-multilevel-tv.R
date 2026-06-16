# ============================================================================
# Tests for time-varying multilevel support
# ============================================================================

make_multilevel_tv_stub_fit <- function(tv_phi = TRUE, tv_sigma = TRUE) {
  n_draw <- 30
  N <- 2L
  n_time <- 4L
  n_eff <- n_time - 1L
  phi_mask <- .resolve_phi_tv_mask(if (isTRUE(tv_phi)) "cross" else FALSE)
  n_phi_tv <- sum(phi_mask)

  variables <- c(
    paste0("phi_bar[", seq_len(4), "]"),
    paste0("tau_phi[", seq_len(4), "]"),
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
    "sigma_eps[1]", "sigma_eps[2]",
    paste0(
      "log_lik[",
      rep(seq_len(N), each = n_eff),
      ",",
      rep(seq_len(n_eff), times = N),
      "]"
    )
  )
  if (isTRUE(tv_phi)) {
    variables <- c(
      variables,
      paste0("tau_phi_tv[", seq_len(n_phi_tv), "]"),
      paste0(
        "phi_raw[",
        rep(seq_len(n_eff), each = n_phi_tv),
        ",",
        rep(seq_len(n_phi_tv), times = n_eff),
        "]"
      ),
      paste0(
        "phi_t[",
        rep(seq_len(n_eff), each = 4),
        ",",
        rep(seq_len(4), times = n_eff),
        "]"
      ),
      paste0(
        "phi_unit_t[",
        rep(seq_len(N), each = n_eff * 4),
        ",",
        rep(rep(seq_len(n_eff), each = 4), times = N),
        ",",
        rep(seq_len(4), times = N * n_eff),
        "]"
      )
    )
  }
  if (isTRUE(tv_sigma)) {
    variables <- c(
      variables,
      "tau_sigma[1]", "tau_sigma[2]",
      paste0(
        "sigma_raw[",
        rep(seq_len(n_eff), each = 2),
        ",",
        rep(seq_len(2), times = n_eff),
        "]"
      ),
      paste0(
        "sigma_t[",
        rep(seq_len(n_eff), each = 2),
        ",",
        rep(seq_len(2), times = n_eff),
        "]"
      )
    )
  }

  draws <- array(
    0,
    dim = c(n_draw, 1, length(variables)),
    dimnames = list(NULL, NULL, variables)
  )
  fill <- function(name, value) {
    draws[, , name] <<- rep(value, length.out = n_draw)
  }

  phi_bar <- c(0.15, 0.05, -0.02, 0.25)
  unit_phi <- rbind(
    c(0.20, 0.40, -0.10, 0.30),
    c(0.10, -0.20, 0.30, 0.40)
  )
  for (k in seq_len(4)) {
    fill(paste0("phi_bar[", k, "]"), phi_bar[k])
    fill(paste0("tau_phi[", k, "]"), 0.10 + 0.01 * k)
    for (i in seq_len(N)) {
      fill(paste0("z_phi[", i, ",", k, "]"), 0)
      fill(paste0("phi_unit[", i, ",", k, "]"), unit_phi[i, k])
    }
  }
  fill("rho", 0.35)
  fill("sigma_eps[1]", 0.9)
  fill("sigma_eps[2]", 1.1)
  for (i in seq_len(N)) {
    for (t in seq_len(n_eff)) {
      fill(paste0("log_lik[", i, ",", t, "]"), -1)
    }
  }

  if (isTRUE(tv_phi)) {
    for (k in seq_len(n_phi_tv)) fill(paste0("tau_phi_tv[", k, "]"), 0.03 + 0.01 * k)
    for (t in seq_len(n_eff)) {
      for (k in seq_len(n_phi_tv)) {
        fill(paste0("phi_raw[", t, ",", k, "]"), 0.01 * t)
      }
    }
    for (t in seq_len(n_eff)) {
      drift <- c(0.01, 0.02, -0.01, 0.015) * (t - 1L)
      for (k in seq_len(4)) {
        fill(paste0("phi_t[", t, ",", k, "]"), phi_bar[k] + drift[k])
        for (i in seq_len(N)) {
          fill(paste0("phi_unit_t[", i, ",", t, ",", k, "]"), unit_phi[i, k] + drift[k])
        }
      }
    }
  }

  if (isTRUE(tv_sigma)) {
    fill("tau_sigma[1]", 0.04)
    fill("tau_sigma[2]", 0.05)
    for (t in seq_len(n_eff)) {
      for (d in seq_len(2)) {
        fill(paste0("sigma_raw[", t, ",", d, "]"), 0.01 * d)
      }
    }
    for (t in seq_len(n_eff)) {
      fill(paste0("sigma_t[", t, ",1]"), 0.9 + 0.04 * t)
      fill(paste0("sigma_t[", t, ",2]"), 1.1 - 0.02 * t)
    }
  }

  y_list <- list(
    matrix(c(1, 2, 2, 1, 3, 0, 4, -1), ncol = 2, byrow = TRUE),
    matrix(c(-1, 1, 0, 2, 1, 3, 2, 4), ncol = 2, byrow = TRUE)
  )
  stan_data <- structure(
    list(
      N = N,
      n_time = n_time,
      y = y_list,
      family = c(1L, 1L),
      tv_phi = as.integer(isTRUE(tv_phi)),
      phi_tv_mask = unname(phi_mask),
      tv_sigma = as.integer(isTRUE(tv_sigma)),
      tau_phi_prior = 0.025,
      tau_sigma_prior = 0.05,
      barrier_k = 8
    ),
    ids = c("a", "b"),
    time_values = 11:14,
    margins = "normal",
    vars = c("y1", "y2"),
    person_means = matrix(c(10, 20, 30, 40), nrow = N, byrow = TRUE),
    tv_phi = isTRUE(tv_phi),
    phi_tv_mask = phi_mask,
    tv_sigma = isTRUE(tv_sigma)
  )

  new_dcvar_multilevel_tv_fit(
    fit = posterior::as_draws_array(draws),
    stan_data = stan_data,
    N = N,
    vars = c("y1", "y2"),
    centered = TRUE,
    person_means = attr(stan_data, "person_means"),
    margins = "normal",
    tv_phi = isTRUE(tv_phi),
    phi_tv_mask = phi_mask,
    tv_sigma = isTRUE(tv_sigma),
    backend = "rstan",
    priors = list(
      phi_bar_sd = 0.5,
      tau_phi_scale = 0.2,
      sigma_sd = 1,
      rho_sd = 0.5,
      tau_phi_rate = 0.025,
      tau_sigma_rate = 0.05
    ),
    meta = list(chains = 1, iter_warmup = 10, iter_sampling = n_draw)
  )
}

test_that("multilevel TV routing and prepared data expose the mixed engine contract", {
  expect_identical(.margin_stan_file("multilevel_tv", "normal"), "multilevel_tv_mixed.stan")
  expect_identical(.margin_stan_file("multilevel_tv", c("normal", "gamma")), "multilevel_tv_mixed.stan")
  expect_match(dcvar_stan_path("multilevel_tv"), "multilevel_tv_mixed\\.stan$")
  expect_true(file.exists(dcvar_stan_path("multilevel_tv")))

  sim <- simulate_dcvar_multilevel(
    N = 3,
    n_time = 6,
    burnin = 0,
    phi_trajectory = list(rep(0.3, 5), rep(0.1, 5), rep(-0.05, 5), seq(0.2, 0.4, length.out = 5)),
    sigma_trajectory = cbind(seq(0.8, 1.1, length.out = 5), rep(1, 5)),
    seed = 202
  )
  expect_equal(dim(sim$true_params$Phi_population), c(5L, 4L))
  expect_length(sim$true_params$Phi_unit_paths, 3L)
  expect_equal(dim(sim$true_params$sigma), c(5L, 2L))

  md <- prepare_multilevel_data(
    sim$data,
    vars = c("y1", "y2"),
    tv_phi = "cross",
    tv_sigma = TRUE
  )
  expect_identical(md$tv_phi, 1L)
  expect_identical(md$phi_tv_mask, c(0L, 1L, 1L, 0L))
  expect_identical(md$tv_sigma, 1L)
  expect_identical(md$family, c(1L, 1L))
  expect_identical(md$skew_direction, c(1, 1))
  expect_true(attr(md, "tv_phi"))
  expect_true(attr(md, "tv_sigma"))
})

test_that("multilevel TV initial values size active random-walk containers only", {
  init <- .init_multilevel_tv_params(2, 3, 6, "normal", tv_phi = "ar", tv_sigma = TRUE)
  expect_length(init$tau_phi_tv, 2)
  expect_identical(dim(init$phi_raw), c(5L, 2L))
  expect_length(init$tau_sigma, 2)
  expect_identical(dim(init$sigma_raw), c(5L, 2L))
  expect_identical(dim(init$z_phi), c(3L, 4L))

  constant <- .init_multilevel_tv_params(2, 3, 6, "normal", tv_phi = FALSE, tv_sigma = FALSE)
  expect_null(constant$tau_phi_tv)
  expect_null(constant$phi_raw)
  expect_null(constant$tau_sigma)
  expect_null(constant$sigma_raw)
})

test_that("multilevel TV fit object methods report population and unit trajectories", {
  fit <- make_multilevel_tv_stub_fit(tv_phi = TRUE, tv_sigma = TRUE)
  expect_s3_class(fit, "dcvar_multilevel_tv_fit")
  expect_identical(class(fit), c("dcvar_multilevel_tv_fit", "dcvar_multilevel_fit", "dcvar_model_fit"))
  expect_identical(fit$model, "multilevel_tv")
  expect_true(fit$tv_phi)
  expect_true(fit$tv_sigma)
  expect_identical(unname(fit$phi_tv_mask), c(0L, 1L, 1L, 0L))

  co <- coef(fit)
  expect_true(all(c("phi_bar", "tau_phi", "sigma_eps", "rho", "tau_phi_tv", "tau_sigma") %in% names(co)))
  expect_identical(names(co$tau_phi_tv), c("tau_phi_tv[phi12]", "tau_phi_tv[phi21]"))

  vp <- var_params(fit)
  expect_true(all(c("phi_bar", "tau_phi", "sigma_eps", "rho", "tau_phi_tv", "tau_sigma") %in% names(vp)))
  expect_identical(vp$tau_phi_tv$variable, c("tau_phi_tv[phi12]", "tau_phi_tv[phi21]"))

  rho_df <- rho_trajectory(fit)
  expect_equal(nrow(rho_df), 3L)
  expect_identical(rho_df$time, 12:14)
  expect_true(length(unique(rho_df$mean)) == 1L)

  phi_pop <- phi_trajectory(fit)
  expect_equal(nrow(phi_pop), 12L)
  expect_false("unit" %in% names(phi_pop))
  expect_identical(sort(unique(phi_pop$coefficient)), c("phi11", "phi12", "phi21", "phi22"))

  phi_units <- phi_trajectory(fit, unit = "all")
  expect_equal(nrow(phi_units), 24L)
  expect_identical(sort(unique(as.character(phi_units$unit))), c("a", "b"))

  sigma_df <- sigma_trajectory(fit)
  expect_equal(nrow(sigma_df), 6L)
  expect_true(all(sigma_df$mean > 0))

  dep_df <- dependence_summary(fit)
  expect_equal(nrow(dep_df), 3L)

  expect_output(print(fit), "TV Multilevel")
  s <- summary(fit)
  expect_s3_class(s, "dcvar_multilevel_tv_summary")
  expect_false(is.null(s$phi_summary))
  expect_false(is.null(s$sigma_summary))
  expect_silent(invisible(capture.output(print(s))))
})

test_that("multilevel TV fitted and predict use unit-time Phi and sigma paths", {
  fit <- make_multilevel_tv_stub_fit(tv_phi = TRUE, tv_sigma = TRUE)

  fv <- fitted(fit)
  expect_equal(nrow(fv), 6L)
  expect_identical(names(fv), c("unit", "time", "y1", "y2"))
  expect_equal(fv$y1[1], 1 * 0.20 + 2 * 0.40)
  expect_equal(fv$y2[1], 1 * -0.10 + 2 * 0.30)

  response <- fitted(fit, type = "response")
  expect_equal(response$y1[1], fv$y1[1] + 10)
  expect_equal(response$y2[1], fv$y2[1] + 20)

  pr <- predict(fit, ci_level = 0.8)
  expect_equal(nrow(pr), 12L)
  expect_true(all(c("unit", "time", "variable", "mean", "lower", "upper") %in% names(pr)))
  expect_true(all(pr$lower <= pr$mean))
  expect_true(all(pr$mean <= pr$upper))
})

test_that("multilevel TV phi_trajectory tiles baselines when tv_phi is off", {
  fit <- make_multilevel_tv_stub_fit(tv_phi = FALSE, tv_sigma = TRUE)
  expect_false(fit$tv_phi)

  phi_pop <- phi_trajectory(fit)
  expect_equal(nrow(phi_pop), 12L)
  per_coef <- tapply(phi_pop$mean, phi_pop$coefficient, function(x) length(unique(x)))
  expect_true(all(per_coef == 1L))

  phi_units <- phi_trajectory(fit, unit = c("a", "b"))
  expect_equal(nrow(phi_units), 24L)
  expect_identical(sort(unique(as.character(phi_units$unit))), c("a", "b"))
})
