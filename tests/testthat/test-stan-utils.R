test_that("dcvar_stan_path returns valid file paths", {
  for (m in c("dcvar", "hmm", "constant", "multilevel", "sem")) {
    path <- dcvar_stan_path(m)
    expect_true(file.exists(path), info = paste("model:", m))
    expect_true(grepl("\\.stan$", path), info = paste("model:", m))
  }

  sem_exp <- dcvar_stan_path("sem", margins = "exponential")
  expect_true(file.exists(sem_exp))
  expect_match(basename(sem_exp), "^sem_EG\\.stan$")
})

test_that("dcvar_stan_path rejects invalid model names", {
  expect_error(dcvar_stan_path("invalid"))
})

test_that(".resolve_backend() resolves auto to rstan", {
  expect_equal(.resolve_backend("auto"), "rstan")
})

test_that(".compiled_exe_path uses the correct platform suffix", {
  path <- .compiled_exe_path(tempdir(), "dcvar-model")

  if (.Platform$OS.type == "windows") {
    expect_match(path, "\\.exe$")
  } else {
    expect_false(grepl("\\.exe$", path, ignore.case = TRUE))
  }
})

test_that(".compile_model_backend() rejects missing Stan source files", {
  missing_file <- file.path(tempdir(), "does-not-exist.stan")
  expect_error(
    .compile_model_backend(missing_file, backend = "rstan"),
    "does not exist"
  )
})

test_that(".stan_include_paths() keeps custom and bundled include roots", {
  custom_dir <- tempfile("custom-stan-dir")
  dir.create(custom_dir)
  custom_file <- file.path(custom_dir, "model.stan")
  package_dir <- system.file("stan", package = "dcvar")
  custom_dir_norm <- normalizePath(custom_dir, winslash = "/", mustWork = TRUE)
  package_dir_norm <- normalizePath(package_dir, winslash = "/", mustWork = TRUE)

  include_paths <- .stan_include_paths(custom_file, package_dir)

  expect_equal(include_paths[[1]], custom_dir_norm)
  expect_true(package_dir_norm %in% include_paths)
  expect_equal(length(include_paths), length(unique(include_paths)))
})

test_that(".stan_cache_fingerprint changes when recursively included files change", {
  root <- tempfile("stan-cache-root")
  dir.create(file.path(root, "functions"), recursive = TRUE)

  model_file <- file.path(root, "model.stan")
  shared_file <- file.path(root, "functions", "shared.stan")
  nested_file <- file.path(root, "functions", "nested.stan")

  writeLines("#include functions/shared.stan", model_file)
  writeLines("#include functions/nested.stan", shared_file)
  writeLines("functions { real nested() { return 1; } }", nested_file)

  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  fingerprint_1 <- .stan_cache_fingerprint(model_file, c(root, system.file("stan", package = "dcvar")))
  writeLines("functions { real nested() { return 2; } }", nested_file)
  fingerprint_2 <- .stan_cache_fingerprint(model_file, c(root, system.file("stan", package = "dcvar")))

  expect_false(identical(fingerprint_1, fingerprint_2))
})

test_that(".stan_cache_fingerprint respects include search path order", {
  model_root <- tempfile("stan-cache-model")
  extra_root <- tempfile("stan-cache-extra")
  package_root <- tempfile("stan-cache-package")

  dir.create(model_root)
  dir.create(file.path(extra_root, "functions"), recursive = TRUE)
  dir.create(file.path(package_root, "functions"), recursive = TRUE)

  model_file <- file.path(model_root, "model.stan")
  writeLines("#include functions/shared.stan", model_file)
  writeLines("functions { real shared() { return 1; } }", file.path(extra_root, "functions", "shared.stan"))
  writeLines("functions { real shared() { return 2; } }", file.path(package_root, "functions", "shared.stan"))

  on.exit(unlink(c(model_root, extra_root, package_root), recursive = TRUE), add = TRUE)

  fingerprint_1 <- .stan_cache_fingerprint(model_file, c(extra_root, package_root))
  fingerprint_2 <- .stan_cache_fingerprint(model_file, c(package_root, extra_root))

  expect_false(identical(fingerprint_1, fingerprint_2))
})

test_that(".fit_summary() assigns stable names to unnamed summary functions", {
  draws <- posterior::as_draws_array(
    array(
      rnorm(40),
      dim = c(10, 1, 4),
      dimnames = list(NULL, NULL, c("phi_unit[1,1]", "phi_unit[1,2]", "phi_unit[1,3]", "phi_unit[1,4]"))
    )
  )

  summ <- .fit_summary(
    draws,
    variables = NULL,
    backend = "rstan",
    required = paste0("phi_unit[1,", 1:4, "]"),
    required_type = "exact",
    context = "test",
    output_type = "parameter",
    mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )

  expect_true(all(c("variable", "mean", "sd", "q2.5", "q97.5") %in% names(summ)))
})

test_that("all bundled top-level Stan files translate successfully", {
  skip_if_no_cmdstanr_toolchain()

  stanc <- file.path(cmdstanr::cmdstan_path(), "bin", "stanc")
  stan_dir <- system.file("stan", package = "dcvar")
  stan_files <- list.files(stan_dir, pattern = "\\.stan$", full.names = TRUE, recursive = TRUE)
  stan_files <- stan_files[!grepl("/functions/", stan_files, fixed = TRUE)]

  for (stan_file in stan_files) {
    output_file <- tempfile(
      pattern = paste0(tools::file_path_sans_ext(basename(stan_file)), "_"),
      fileext = ".hpp"
    )
    on.exit(unlink(output_file), add = TRUE)

    result <- system2(
      stanc,
      c("--include-paths", stan_dir, "--output", output_file, normalizePath(stan_file)),
      stdout = TRUE,
      stderr = TRUE
    )
    status <- attr(result, "status")
    if (is.null(status)) status <- 0L
    expect_equal(status, 0L, info = paste(stan_file, paste(result, collapse = "\n")))
    expect_true(file.exists(output_file), info = paste("missing output for", stan_file))
  }
})
