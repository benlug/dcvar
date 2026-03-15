test_that("dcvar_stan_path returns valid file paths", {
  for (m in c("dcvar", "hmm", "constant", "multilevel", "sem")) {
    path <- dcvar_stan_path(m)
    expect_true(file.exists(path), info = paste("model:", m))
    expect_true(grepl("\\.stan$", path), info = paste("model:", m))
  }
})

test_that("dcvar_stan_path rejects invalid model names", {
  expect_error(dcvar_stan_path("invalid"))
})

test_that(".compiled_exe_path uses the correct platform suffix", {
  path <- .compiled_exe_path(tempdir(), "dcvar-model")

  if (.Platform$OS.type == "windows") {
    expect_match(path, "\\.exe$")
  } else {
    expect_false(grepl("\\.exe$", path, ignore.case = TRUE))
  }
})

test_that("all bundled top-level Stan files translate successfully", {
  skip_if_no_cmdstanr()

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
