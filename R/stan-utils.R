#' Get path to bundled Stan model file
#'
#' Returns the file path to a Stan model file included with the package.
#'
#' @param model Character string: `"dcVar"`, `"hmm"`, `"constant"`,
#'   `"multilevel"`, or `"sem"`.
#' @param margins Character string: margin type (`"normal"`, `"exponential"`,
#'   `"skew_normal"`, `"gamma"`). Default: `"normal"`.
#'
#' @return File path to the Stan model file.
#' @export
#'
#' @examples
#' dcVar_stan_path("dcVar")
#' dcVar_stan_path("constant", margins = "exponential")
dcVar_stan_path <- function(model = c("dcVar", "hmm", "constant", "multilevel", "sem"),
                            margins = "normal") {
  model <- match.arg(model)

  if (model %in% c("multilevel", "sem")) {
    file <- switch(model,
      multilevel = "multilevel_copula_var.stan",
      sem = "sem_copula_var.stan"
    )
  } else {
    file <- .margin_stan_file(model, margins)
  }

  system.file("stan", file, package = "dcVar", mustWork = TRUE)
}

#' Validate common MCMC sampling arguments
#'
#' @param chains,iter_warmup,iter_sampling,adapt_delta,max_treedepth
#'   Sampling arguments to validate.
#' @return Invisible `TRUE` if all checks pass.
#' @noRd
.validate_sampling_args <- function(chains, iter_warmup, iter_sampling,
                                    adapt_delta, max_treedepth) {
  if (!is.numeric(chains) || chains < 1) {
    cli_abort("{.arg chains} must be a positive integer, got {.val {chains}}.")
  }
  if (!is.numeric(iter_warmup) || iter_warmup < 1) {
    cli_abort("{.arg iter_warmup} must be a positive integer, got {.val {iter_warmup}}.")
  }
  if (!is.numeric(iter_sampling) || iter_sampling < 1) {
    cli_abort("{.arg iter_sampling} must be a positive integer, got {.val {iter_sampling}}.")
  }
  if (!is.numeric(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1) {
    cli_abort("{.arg adapt_delta} must be in (0, 1), got {.val {adapt_delta}}.")
  }
  if (!is.numeric(max_treedepth) || max_treedepth < 1) {
    cli_abort("{.arg max_treedepth} must be a positive integer, got {.val {max_treedepth}}.")
  }
  invisible(TRUE)
}


#' Check cmdstanr availability
#'
#' @return Invisible `TRUE` if cmdstanr and CmdStan are available.
#' @noRd
.check_cmdstanr <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    cli_abort(c(
      "Package {.pkg cmdstanr} is required but not installed.",
      "i" = "Install with: {.code install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))}",
      "i" = "Then run: {.code cmdstanr::install_cmdstan()}"
    ))
  }
  cmdstan_path <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) NULL
  )
  if (is.null(cmdstan_path)) {
    cli_abort(c(
      "CmdStan is not installed.",
      "i" = "Run: {.code cmdstanr::install_cmdstan()}"
    ))
  }
  invisible(TRUE)
}

#' Build a platform-appropriate cached executable path
#'
#' @param cache_dir Directory used for compiled model caching.
#' @param cache_key Basename for the cached executable.
#' @return Absolute path to the cached executable.
#' @noRd
.compiled_exe_path <- function(cache_dir, cache_key) {
  exe_file <- file.path(cache_dir, cache_key)
  if (.Platform$OS.type == "windows" && !grepl("\\.exe$", exe_file, ignore.case = TRUE)) {
    exe_file <- paste0(exe_file, ".exe")
  }
  exe_file
}

#' Compile a Stan model with caching
#'
#' @param model_type Character string: `"dcVar"`, `"hmm"`, `"constant"`,
#'   `"multilevel"`, or `"sem"`.
#' @param margins Character: margin type (default: `"normal"`).
#' @param stan_file Optional path to a custom Stan file. If `NULL`, uses the
#'   bundled model.
#' @param force_recompile Force recompilation even if cached.
#' @param quiet Suppress compilation messages.
#'
#' @return A `CmdStanModel` object.
#' @noRd
.compile_model <- function(model_type = c("dcVar", "hmm", "constant", "multilevel", "sem"),
                           margins = "normal",
                           stan_file = NULL,
                           force_recompile = FALSE,
                           quiet = FALSE) {
  .check_cmdstanr()
  model_type <- match.arg(model_type)

  if (is.null(stan_file)) {
    stan_file <- dcVar_stan_path(model_type, margins = margins)
  }

  cache_dir <- file.path(tempdir(), "dcVar-cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Include the specific Stan source in the cache key so custom files cannot
  # reuse an executable built from a different model.
  stan_hash <- unname(tools::md5sum(stan_file))
  if (length(stan_hash) != 1 || is.na(stan_hash) || !nzchar(stan_hash)) {
    cli_abort("Failed to hash Stan source file {.file {stan_file}} for caching.")
  }
  cache_key <- paste0(.margin_cache_key(model_type, margins), "_", stan_hash)
  exe_file <- .compiled_exe_path(cache_dir, cache_key)

  # Include path for shared Stan functions (e.g., gaussian_copula.stan)
  include_paths <- system.file("stan", package = "dcVar", mustWork = TRUE)

  cmdstanr::cmdstan_model(
    stan_file = stan_file,
    exe_file = exe_file,
    include_paths = include_paths,
    force_recompile = force_recompile,
    quiet = quiet
  )
}
