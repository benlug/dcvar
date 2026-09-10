causal_effect_fixture <- function(model, parameters, weights = c(0.5, 0.5)) {
  columns <- list()
  for (name in names(parameters)) {
    x <- parameters[[name]]
    if (is.matrix(x)) {
      indices <- expand.grid(seq_len(nrow(x)), seq_len(ncol(x)))
      nms <- sprintf("%s[%d,%d]", name, indices[[1]], indices[[2]])
    } else nms <- sprintf("%s[%d]", name, seq_along(x))
    columns <- c(columns, stats::setNames(as.list(as.vector(x)), nms))
  }
  a <- array(rep(unlist(columns), each = 8), c(4, 2, length(columns)),
             dimnames = list(NULL, NULL, names(columns)))
  identity <- list(center = 0, scale = 1)
  structure(list(
    fit = posterior::as_draws_array(a), model = model, backend = "rstan",
    meta = list(target_weights = weights, roles = list(outcome = "Y"),
                scales = list(outcome = identity, covariate = identity))
  ), class = "dcvar_causal_sem_fit")
}

test_that("regression effects use weighted means and input units", {
  fit <- causal_effect_fixture("latent_covariate", list(
    alpha = c(0.2, 1.1), beta = c(0.3, 0.7), mu_v = c(-0.2, 0.8)
  ), c(0.7, 0.3))
  e <- causal_effects(fit)
  expect_equal(e$mean, c(0.94, 0.4))
  expect_equal(attr(e, "target_weights"), c(0.7, 0.3))
  expect_equal(dim(causal_effects(fit, summary = FALSE)), c(4, 2, 2))
  expect_equal(causal_effects(fit, target_weights = c(1, 0))$mean[1], 0.82)

  fit$model <- "latent_outcome"
  fit$meta$scales$covariate <- list(center = 4, scale = 2)
  fit$meta$scales$outcome <- list(center = 10, scale = 3)
  e <- causal_effects(fit, covariate_values = c(2, 4, 6))
  expect_equal(e$mean, c(2.82, 0.6, 1.5, 2.7, 3.9))
  expect_equal(causal_effects(fit, scale = "model")$mean, c(0.94, 0.4))
})

test_that("mediation uses the paper's adjusted means", {
  p <- list(alpha_m = c(0.1, 0.8), B = rbind(c(0.3, 0.2), c(0.7, -0.1)),
            alpha_y = c(-0.2, 0.5), D = rbind(c(0.1, 0.4), c(-0.2, 0.3)),
            d = c(0.5, 1.2), mu_q = rbind(c(0, 1), c(2, -1)))
  for (model in c("latent_mediator", "latent_mediator_baseline")) {
    fit <- causal_effect_fixture(model, p, c(0.25, 0.75))
    e <- causal_effects(fit)
    # Direct adjustment uses mediator mean 1.8. Total adjustment uses 0.45 and 1.9.
    expect_equal(e$mean, c(2.355, 1.56, 0.795))
    raw <- causal_effects(fit, summary = FALSE)
    expect_equal(as.numeric(raw[, , "ATE"]), as.numeric(raw[, , "ADE"] + raw[, , "AIE"]))
    expect_error(causal_effects(fit, covariate_values = 0), "regression models")
    reversed <- lapply(p, function(x) if (is.matrix(x)) x[2:1, , drop = FALSE] else x[2:1])
    swapped <- causal_effect_fixture(model, reversed, c(0.75, 0.25))
    expect_equal(causal_effects(swapped)$mean, -e$mean)
  }
})

test_that("effects retain posterior dependence between paths", {
  p <- list(alpha_m = c(0, 0), B = matrix(0, 2, 2), alpha_y = c(0, 0),
            D = matrix(0, 2, 2), d = c(0, 0), mu_q = matrix(0, 2, 2))
  fit <- causal_effect_fixture("latent_mediator", p)
  fit$fit[, , "alpha_m[2]"] <- rep(c(1, 3), 4)
  fit$fit[, , "d[2]"] <- rep(c(3, 1), 4)
  expect_equal(causal_effects(fit)$mean, c(3, 1.5, 1.5))
  # Multiplication of separate posterior means would give total effect 4.
  expect_false(isTRUE(all.equal(causal_effects(fit)$mean[1], 4)))
})

test_that("a change of covariate units preserves the outcome effect", {
  p <- list(alpha = c(0.2, 1.1), beta = c(0.3, 0.7), mu_v = c(-0.2, 0.8))
  original <- causal_effect_fixture("latent_covariate", p)
  changed <- p
  changed$alpha <- p$alpha - p$beta * 2 / 3
  changed$beta <- p$beta / 3
  changed$mu_v <- 2 + 3 * p$mu_v
  transformed <- causal_effect_fixture("latent_covariate", changed)
  expect_equal(causal_effects(transformed)$mean[1], causal_effects(original)$mean[1])
  expect_equal(causal_effects(transformed)$mean[2], causal_effects(original)$mean[2] / 3)
})

test_that("effect arguments and unsupported objects fail clearly", {
  fit <- causal_effect_fixture("latent_covariate", list(
    alpha = c(0, 0), beta = c(0, 0), mu_v = c(0, 0)))
  for (w in list(c(1, 1), c(-1, 2), NA_real_, c(0.1, 0.2, 0.7))) {
    expect_error(causal_effects(fit, target_weights = w), "target_weights")
  }
  expect_error(causal_effects(fit, summary = NA), "summary")
  expect_error(causal_effects(fit, level = 1), "level")
  expect_error(causal_effects(fit, covariate_values = Inf), "finite")
  expect_error(causal_effects(list()), "dcvar_causal_sem_fit")
  expect_error(predict(fit, newdata = data.frame(x = 1)), "New data")
  expect_error(predict(fit, type = "mediator"), "mediation model")
})

test_that("sampler validation runs before compilation", {
  expect_error(dcvar_causal_sem(chains = 0), "chains")
  expect_error(dcvar_causal_sem(seed = -1), "seed")
  expect_error(dcvar_causal_sem(refresh = -1), "refresh")
})

test_that("all causal Stan models are available to both parsers", {
  models <- c("latent_covariate", "latent_outcome", "latent_mediator", "latent_mediator_baseline")
  for (model in models) {
    path <- .causal_sem_stan_path(model)
    expect_true(file.exists(path))
    result <- rstan::stanc(file = path, isystem = dirname(path))
    expect_true(result$status)
  }
})
