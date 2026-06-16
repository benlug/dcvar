test_that(".validate_margins accepts valid margin types", {
  expect_invisible(.validate_margins("normal"))
  expect_invisible(.validate_margins("exponential", c(1, -1)))
  expect_invisible(.validate_margins("skew_normal"))
  expect_invisible(.validate_margins("gamma", c(1, 1)))
})

test_that(".validate_copula accepts and rejects copula families", {
  expect_invisible(.validate_copula("gaussian"))
  expect_invisible(.validate_copula("clayton"))
  expect_error(.validate_copula("frank"), "must be one of")
})

test_that(".validate_margins rejects invalid margin types", {
  expect_error(.validate_margins("invalid"), "must be one of")
  expect_error(.validate_margins("student_t"), "must be one of")
})

test_that(".validate_margins requires skew_direction for exponential", {
  expect_error(
    .validate_margins("exponential"),
    "skew_direction.*required"
  )
})

test_that(".validate_margins requires skew_direction for gamma", {
  expect_error(
    .validate_margins("gamma"),
    "skew_direction.*required"
  )
})

test_that(".validate_margins validates skew_direction length and values", {
  expect_error(
    .validate_margins("exponential", c(1)),
    "length-2"
  )
  expect_error(
    .validate_margins("exponential", c(1, -1, 1)),
    "length-2"
  )
  expect_error(
    .validate_margins("exponential", c(1, 2)),
    "\\+1 or -1"
  )
  expect_error(
    .validate_margins("gamma", c(0, 1)),
    "\\+1 or -1"
  )
})

test_that(".margin_stan_suffix returns correct suffix for each margin type", {
  expect_equal(.margin_stan_suffix("normal"), "")
  expect_equal(.margin_stan_suffix("exponential"), "_EG")
  expect_equal(.margin_stan_suffix("skew_normal"), "_SNG")
  expect_equal(.margin_stan_suffix("gamma"), "_GG")
})

test_that(".margin_stan_file returns correct filename for normal margins", {
  expect_equal(.margin_stan_file("constant", "normal"),
               "constant_copula_var.stan")
  expect_equal(.margin_stan_file("dcvar", "normal"),
               "dcvar_model_ncp.stan")
  expect_equal(.margin_stan_file("hmm", "normal"),
               "hmm_copula_model.stan")
  expect_equal(.margin_stan_file("multilevel", "normal"),
               "multilevel_copula_var.stan")
  expect_equal(.margin_stan_file("sem_naive", "normal"),
               "sem_naive_NG.stan")
})

test_that(".margin_stan_file returns correct filename for non-normal margins", {
  expect_equal(.margin_stan_file("constant", "exponential"),
               "constant_EG.stan")
  expect_equal(.margin_stan_file("dcvar", "exponential"),
               "dcvar_EG_ncp.stan")
  expect_equal(.margin_stan_file("hmm", "exponential"),
               "hmm_EG.stan")
  expect_equal(.margin_stan_file("multilevel", "exponential"),
               "multilevel_EG.stan")
  expect_equal(.margin_stan_file("sem_naive", "exponential"),
               "sem_naive_EG.stan")

  expect_equal(.margin_stan_file("dcvar", "skew_normal"),
               "dcvar_SNG_ncp.stan")
  expect_equal(.margin_stan_file("hmm", "gamma"),
               "hmm_GG.stan")
  expect_equal(.margin_stan_file("constant", "normal", copula = "clayton"),
               "constant_NCl.stan")
})

test_that(".margin_cache_key returns expected cache key", {
  expect_equal(.margin_cache_key("constant", "normal"),
               "constant_model")
  expect_equal(.margin_cache_key("dcvar", "normal"),
               "dcvar_model")
  expect_equal(.margin_cache_key("hmm", "exponential"),
               "hmm_EG_model")
  expect_equal(.margin_cache_key("dcvar", "skew_normal"),
               "dcvar_SNG_model")
  expect_equal(.margin_cache_key("constant", "gamma"),
               "constant_GG_model")
  expect_equal(.margin_cache_key("constant", "normal", copula = "clayton"),
               "constant_clayton_model")
})


# ---------------------------------------------------------------------------
# Per-variable (mixed) margins
# ---------------------------------------------------------------------------

test_that(".family_codes maps families to the Stan dispatch integers", {
  expect_equal(unname(.family_codes[c("normal", "exponential", "skew_normal", "gamma")]),
               c(1L, 2L, 3L, 4L))
})

test_that(".normalize_margins_spec collapses identical vectors only", {
  expect_identical(.normalize_margins_spec(c("normal", "normal")), "normal")
  expect_identical(.normalize_margins_spec("normal"), "normal")
  expect_identical(.normalize_margins_spec(c("normal", "exponential")),
                   c("normal", "exponential"))
})

test_that(".is_mixed_margins detects genuinely mixed specs", {
  expect_false(.is_mixed_margins("normal"))
  expect_false(.is_mixed_margins(c("gamma", "gamma")))
  expect_true(.is_mixed_margins(c("normal", "exponential")))
})

test_that(".validate_margins accepts length-2 vectors", {
  expect_invisible(.validate_margins(c("normal", "exponential"), c(1, 1)))
  expect_invisible(.validate_margins(c("normal", "skew_normal")))
  expect_invisible(.validate_margins(c("normal", "normal")))
})

test_that(".validate_margins requires skew_direction when any dim is exp/gamma", {
  expect_error(
    .validate_margins(c("normal", "exponential")),
    "skew_direction.*required"
  )
  expect_error(
    .validate_margins(c("gamma", "skew_normal")),
    "skew_direction.*required"
  )
})

test_that(".validate_margins rejects invalid entries and wrong length", {
  expect_error(.validate_margins(c("normal", "student_t")), "must be one of")
  expect_error(.validate_margins(c("normal", "exponential", "gamma")),
               "length 1 or 2")
})

test_that(".margin_stan_file routes mixed margins to the generic model", {
  expect_equal(.margin_stan_file("constant", c("normal", "exponential")),
               "constant_mixed.stan")
  expect_equal(.margin_stan_file("constant", c("skew_normal", "gamma")),
               "constant_mixed.stan")
  # Identical-family vectors keep routing to the specialised models.
  expect_equal(.margin_stan_file("constant", c("normal", "normal")),
               "constant_copula_var.stan")
  expect_equal(.margin_stan_file("constant", c("gamma", "gamma")),
               "constant_GG.stan")
})

test_that(".margin_stan_file routes mixed margins for every supported model", {
  expect_equal(.margin_stan_file("multilevel", c("normal", "exponential")),
               "multilevel_mixed.stan")
  expect_equal(.margin_stan_file("multilevel", "skew_normal"),
               "multilevel_mixed.stan")
  expect_equal(.margin_stan_file("multilevel", "gamma"),
               "multilevel_mixed.stan")
  expect_equal(.margin_stan_file("sem", c("normal", "exponential")),
               "sem_mixed.stan")
  expect_equal(.margin_stan_file("sem", "skew_normal"),
               "sem_mixed.stan")
  expect_equal(.margin_stan_file("sem", "gamma"),
               "sem_mixed.stan")
  expect_equal(.margin_stan_file("sem_naive", c("normal", "exponential")),
               "sem_naive_mixed.stan")
  expect_equal(.margin_stan_file("sem_naive", "skew_normal"),
               "sem_naive_mixed.stan")
  expect_equal(.margin_stan_file("sem_naive", "gamma"),
               "sem_naive_mixed.stan")
  # Clayton mixed is constant-only; the Gaussian copula covers the rest.
  expect_equal(
    .margin_stan_file("constant", c("normal", "gamma"), copula = "clayton"),
    "constant_mixed_clayton.stan"
  )
  expect_error(
    .margin_stan_file("dcvar", c("normal", "exponential"), copula = "clayton"),
    "Clayton.*constant"
  )
})

test_that(".margin_cache_key encodes the full family vector for mixed fits", {
  expect_equal(.margin_cache_key("constant", c("normal", "exponential")),
               "constant_mixed12_model")
  expect_equal(.margin_cache_key("constant", c("exponential", "gamma")),
               "constant_mixed24_model")
  expect_equal(.margin_cache_key("constant", c("skew_normal", "normal")),
               "constant_mixed31_model")
  expect_equal(.margin_cache_key("multilevel", "skew_normal"),
               "multilevel_mixed33_model")
  expect_equal(.margin_cache_key("sem", "gamma"),
               "sem_mixed44_model")
  expect_equal(.margin_cache_key("sem_naive", "gamma"),
               "sem_naive_mixed44_model")
})

test_that(".mixed_margin_report_vars restricts each group to its own dims", {
  specs <- .mixed_margin_report_vars(c("normal", "exponential"))
  expect_equal(specs$sigma_eps, "sigma_eps[1]")
  expect_equal(specs$sigma_exp, "sigma_exp[2]")
  expect_null(specs$omega)
})

test_that("routing helpers validate family names before dispatching", {
  # Two distinct but invalid families must not silently route to the mixed model.
  expect_error(.margin_stan_file("constant", c("not_a_margin", "also_bad")),
               "must be one of")
  expect_error(.margin_cache_key("constant", c("bad", "worse")),
               "must be one of")
  expect_error(dcvar_stan_path("constant", margins = c("not_a_margin", "x")),
               "must be one of")
  # The exported path helper still resolves a valid mixed request.
  expect_match(dcvar_stan_path("constant", margins = c("normal", "exponential")),
               "constant_mixed\\.stan$")
})

test_that("homogeneous SEM gamma/skew_normal use mixed-engine Stan data", {
  df <- data.frame(time = 1:30, y1_1 = rnorm(30), y1_2 = rnorm(30),
                   y2_1 = rnorm(30), y2_2 = rnorm(30))
  indicators <- list(latent1 = c("y1_1", "y1_2"), latent2 = c("y2_1", "y2_2"))
  gamma_data <- prepare_sem_data(
    df, indicators = indicators, J = 2,
    lambda = c(0.8, 0.8), sigma_e = 0.5,
    margins = "gamma", skew_direction = c(1, 1)
  )
  skew_data <- prepare_sem_data(
    df, indicators = indicators, J = 2,
    lambda = c(0.8, 0.8), sigma_e = 0.5,
    margins = "skew_normal"
  )

  expect_equal(gamma_data$family, c(4L, 4L))
  expect_equal(gamma_data$skew_direction, c(1, 1))
  expect_equal(attr(gamma_data, "margins"), "gamma")
  expect_equal(skew_data$family, c(3L, 3L))
  expect_equal(skew_data$skew_direction, c(1, 1))
  expect_equal(attr(skew_data, "margins"), "skew_normal")
})
