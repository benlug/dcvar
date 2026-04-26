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
