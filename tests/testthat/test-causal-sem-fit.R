test_that("causal SEM models complete the public workflow", {
  skip_if_no_backend("rstan")
  models <- c("latent_covariate", "latent_outcome", "latent_mediator", "latent_mediator_baseline")
  for (model in models) {
    sim <- simulate_dcvar_causal_sem(n = c(30, 30), model = model, seed = 907)
    fit <- suppressWarnings(do.call(dcvar_causal_sem, c(sim$args, list(
      backend = "rstan", chains = 2, cores = 2, iter_warmup = 200,
      iter_sampling = 200, refresh = 0, seed = 1907))))
    expect_s3_class(fit, "dcvar_causal_sem_fit")
    expect_false(inherits(fit, "dcvar_model_fit"))
    expect_true(all(is.finite(coef(fit))))
    e <- causal_effects(fit)
    expect_equal(e$effect, names(sim$effects))
    expect_true(all(is.finite(e$mean)))
    expect_true(all(e$lower <= e$upper))
    expect_s3_class(summary(fit), "dcvar_causal_sem_summary")
    expect_output(print(fit), "Bayesian causal SEM")
    expect_s3_class(draws(fit), "draws_array")
    p <- predict(fit, type = "indicators")
    expect_equal(dim(p)[3], 3 * nrow(sim$data))
    expect_true(all(p %in% 1:3))
    expect_true(all(is.finite(predict(fit))))
    diag <- dcvar_diagnostics(fit)
    expect_true(is.finite(diag$max_rhat))
    expect_true(is.finite(diag$n_divergent))
    expect_false(any(grepl("rep\\[|^lambda\\[1\\]|^L_q\\[[12],1,", diag$parameters$parameter)))
  }
})

test_that("CmdStan fits retain their draws across serialization", {
  skip_if_no_cmdstanr_toolchain()
  sim <- simulate_dcvar_causal_sem(n = c(25, 25), model = "latent_outcome", seed = 3009)
  fit <- suppressWarnings(do.call(dcvar_causal_sem, c(sim$args, list(
    backend = "cmdstanr", chains = 2, cores = 2, iter_warmup = 150,
    iter_sampling = 150, refresh = 0, seed = 53009))))
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_equal(causal_effects(restored), causal_effects(fit))
  expect_equal(predict(restored), predict(fit))
})
