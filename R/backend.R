# ============================================================================
# Backend Adapter Layer
#
# Provides a unified interface for rstan and cmdstanr backends.
# rstan is the default (in Imports); cmdstanr is an optional fallback
# (in Suggests). All backend-specific code is isolated here.
# ============================================================================

#' @importFrom rstan sampling stan_model
NULL

.rstan_model_cache <- new.env(parent = emptyenv())


# -- Backend resolution -------------------------------------------------------

#' Resolve which backend to use
#' @param backend Character: `"auto"`, `"rstan"`, or `"cmdstanr"`.
#' @return Resolved backend string.
#' @noRd
.resolve_backend <- function(backend = c("auto", "rstan", "cmdstanr")) {
  backend <- match.arg(backend)
  if (backend == "auto") {
    return("rstan")
  }

  if (backend == "cmdstanr" && !.backend_available("cmdstanr")) {
    cli_abort(c(
      "Backend {.val cmdstanr} requested but not available.",
      "i" = "Install with: {.code install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))}",
      "i" = "Then run: {.code cmdstanr::install_cmdstan()}"
    ))
  }

  backend
}

#' Check whether a backend is available
#' @param backend Character: `"rstan"` or `"cmdstanr"`.
#' @return Logical.
#' @noRd
.backend_available <- function(backend) {
  switch(backend,
    rstan = requireNamespace("rstan", quietly = TRUE),
    cmdstanr = {
      requireNamespace("cmdstanr", quietly = TRUE) &&
        !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))
    },
    FALSE
  )
}


# -- Core resolution ----------------------------------------------------------

#' Resolve a safe default core count
#' @param chains Number of chains requested.
#' @return Positive integer core count.
#' @noRd
.default_cores <- function(chains) {
  cores_opt <- getOption("mc.cores")
  if (is.numeric(cores_opt) && length(cores_opt) == 1L &&
      is.finite(cores_opt) && cores_opt >= 1) {
    cores <- as.integer(cores_opt)
  } else {
    detected <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (!is.numeric(detected) || length(detected) != 1L ||
        !is.finite(detected) || detected < 1) {
      cores <- 1L
    } else {
      cores <- as.integer(detected)
    }
  }

  limit_cores <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (limit_cores %in% c("true", "t", "1")) {
    cores <- min(cores, 2L)
  }

  min(cores, as.integer(chains))
}


#' Validate and normalize the requested core count
#' @param cores Requested core count or `NULL`.
#' @param chains Number of chains requested.
#' @return Positive integer core count.
#' @noRd
.normalize_cores <- function(cores, chains) {
  if (is.null(cores)) {
    return(.default_cores(chains))
  }

  if (!is.numeric(cores) || length(cores) != 1L ||
      !is.finite(cores) || cores < 1) {
    cli_abort("{.arg cores} must be a positive integer or {.code NULL}, got {.val {cores}}.")
  }

  min(as.integer(cores), as.integer(chains))
}


# -- Model compilation --------------------------------------------------------

#' Compile a Stan model via the resolved backend
#'
#' @param stan_file Path to the `.stan` file.
#' @param backend Character: `"rstan"` or `"cmdstanr"`.
#' @param cache_dir Writable directory for caching compiled models.
#' @param force_recompile Logical; bypass cache.
#' @param quiet Logical; suppress compiler output.
#' @return A compiled model object (class depends on backend).
#' @noRd
.compile_model_backend <- function(stan_file,
                                   backend,
                                   cache_dir = file.path(tempdir(), "dcvar-cache"),
                                   force_recompile = FALSE,
                                   quiet = FALSE) {
  if (!file.exists(stan_file)) {
    cli_abort("Stan source file {.file {stan_file}} does not exist.")
  }

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  include_path <- .stan_include_paths(
    stan_file,
    system.file("stan", package = "dcvar", mustWork = TRUE)
  )
  stan_hash <- .stan_cache_fingerprint(stan_file, include_path)

  switch(backend,
    rstan = .compile_rstan(stan_file, include_path, cache_dir, stan_hash,
                           force_recompile, quiet),
    cmdstanr = .compile_cmdstanr(stan_file, include_path, cache_dir, stan_hash,
                                 force_recompile, quiet)
  )
}

#' Internal: resolve Stan include paths for bundled and custom models
#' @noRd
.stan_include_paths <- function(stan_file, package_include_path) {
  unique(vapply(
    c(dirname(stan_file), package_include_path),
    normalizePath,
    character(1),
    winslash = "/",
    mustWork = TRUE
  ))
}


#' Internal: collect direct `#include` targets from a Stan file
#' @noRd
.stan_direct_includes <- function(stan_file) {
  lines <- readLines(stan_file, warn = FALSE)
  include_lines <- grepl("^\\s*#include\\s+", lines)
  if (!any(include_lines)) {
    return(character())
  }

  includes <- sub("^\\s*#include\\s+", "", lines[include_lines])
  includes <- sub("\\s*(//.*)?$", "", includes)
  includes <- trimws(includes)
  includes <- gsub('^[\"<]|[\">]$', "", includes)
  includes[nzchar(includes)]
}


#' Internal: resolve a Stan include against the active include search path
#' @noRd
.stan_resolve_include <- function(include_name, current_dir, include_path) {
  if (grepl("^(/|[A-Za-z]:[\\\\/])", include_name)) {
    if (file.exists(include_name)) {
      return(normalizePath(include_name, winslash = "/", mustWork = TRUE))
    }
    return(NA_character_)
  }

  candidates <- c(
    file.path(current_dir, include_name),
    file.path(include_path, include_name)
  )

  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  NA_character_
}


#' Internal: fingerprint a Stan file and all recursively included dependencies
#' @noRd
.stan_cache_fingerprint <- function(stan_file, include_path) {
  include_path <- unique(normalizePath(include_path, winslash = "/", mustWork = TRUE))

  visited <- character()
  fingerprint <- character()

  add_file <- function(file_path) {
    file_path <- normalizePath(file_path, winslash = "/", mustWork = TRUE)
    if (file_path %in% visited) {
      return(invisible(NULL))
    }

    visited <<- c(visited, file_path)
    fingerprint <<- c(fingerprint, unname(tools::md5sum(file_path)))

    for (include_name in .stan_direct_includes(file_path)) {
      resolved <- .stan_resolve_include(include_name, dirname(file_path), include_path)
      if (is.na(resolved)) {
        cli_abort(c(
          "Stan include {.val {include_name}} could not be resolved while computing the compile cache key.",
          "i" = "Checked {.file {dirname(file_path)}} and the configured include paths."
        ))
      }
      add_file(resolved)
    }

    invisible(NULL)
  }

  add_file(stan_file)

  signature_file <- tempfile("dcvar-stan-cache-signature-")
  on.exit(unlink(signature_file), add = TRUE)
  writeLines(
    c(
      paste0("include_paths:", paste(include_path, collapse = "|")),
      fingerprint
    ),
    signature_file,
    useBytes = TRUE
  )

  unname(tools::md5sum(signature_file))
}


#' @noRd
.compile_rstan <- function(stan_file, include_path, cache_dir, stan_hash,
                           force_recompile, quiet) {
  cache_key <- paste0("rstan_", stan_hash)

  if (!force_recompile && exists(cache_key, envir = .rstan_model_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .rstan_model_cache, inherits = FALSE))
  }

  model <- rstan::stan_model(
    file = stan_file,
    isystem = include_path,
    auto_write = FALSE,
    verbose = !quiet
  )

  assign(cache_key, model, envir = .rstan_model_cache)

  model
}


#' @noRd
.compile_cmdstanr <- function(stan_file, include_path, cache_dir, stan_hash,
                              force_recompile, quiet) {
  .check_cmdstanr()

  cache_key <- paste0(stan_hash, "_cmdstanr")
  exe_file <- .compiled_exe_path(cache_dir, cache_key)

  cmdstanr::cmdstan_model(
    stan_file = stan_file,
    exe_file = exe_file,
    include_paths = include_path,
    force_recompile = force_recompile,
    quiet = quiet
  )
}


# -- Sampling -----------------------------------------------------------------

#' Run MCMC sampling via the resolved backend
#'
#' @param compiled_model A compiled model object.
#' @param stan_data Named list of data for Stan.
#' @param backend Character: `"rstan"` or `"cmdstanr"`.
#' @param chains,iter_warmup,iter_sampling,adapt_delta,max_treedepth
#'   Standard sampler arguments.
#' @param seed Random seed (or `NULL`).
#' @param cores Number of parallel chains.
#' @param init Init specification (function or list).
#' @param refresh Print frequency.
#' @param ... Additional backend-specific arguments.
#' @return A fitted model object (class depends on backend).
#' @noRd
.sample_model <- function(compiled_model, stan_data, backend,
                          chains, iter_warmup, iter_sampling,
                          adapt_delta, max_treedepth,
                          seed, cores, init, refresh, ...) {
  switch(backend,
    rstan = .sample_rstan(compiled_model, stan_data,
                          chains, iter_warmup, iter_sampling,
                          adapt_delta, max_treedepth,
                          seed, cores, init, refresh, ...),
    cmdstanr = .sample_cmdstanr(compiled_model, stan_data,
                                chains, iter_warmup, iter_sampling,
                                adapt_delta, max_treedepth,
                                seed, cores, init, refresh, ...)
  )
}


#' @noRd
.sample_rstan <- function(compiled_model, stan_data,
                          chains, iter_warmup, iter_sampling,
                          adapt_delta, max_treedepth,
                          seed, cores, init, refresh, ...) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)

  rstan::sampling(
    object = compiled_model,
    data = stan_data,
    chains = chains,
    warmup = iter_warmup,
    iter = iter_warmup + iter_sampling,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    seed = seed,
    cores = cores,
    init = init,
    refresh = refresh,
    ...
  )
}


#' @noRd
.sample_cmdstanr <- function(compiled_model, stan_data,
                             chains, iter_warmup, iter_sampling,
                             adapt_delta, max_treedepth,
                             seed, cores, init, refresh, ...) {
  compiled_model$sample(
    data = stan_data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    parallel_chains = cores,
    init = init,
    refresh = refresh,
    ...
  )
}


# -- Post-fit extraction ------------------------------------------------------

#' Extract posterior draws in a `posterior`-compatible format
#'
#' @param fit Raw backend fit object.
#' @param variables Character vector of parameter names, or `NULL` for all.
#' @param format One of `"draws_array"`, `"draws_matrix"`, `"draws_df"`.
#' @param backend Character: `"rstan"` or `"cmdstanr"`.
#' @return A `posterior::draws_*` object.
#' @noRd
.fit_draws <- function(fit, variables = NULL,
                       format = "draws_array", backend,
                       required = NULL,
                       required_type = c("exact", "pattern"),
                       context = "This method",
                       output_type = "Stan output") {
  required_type <- match.arg(required_type)
  d <- if (backend == "cmdstanr") {
    fit$draws(format = "draws_array")
  } else {
    posterior::as_draws_array(fit)
  }

  .validate_required_stan_outputs(
    d,
    required = required,
    required_type = required_type,
    context = context,
    output_type = output_type
  )

  if (!is.null(variables)) {
    d <- tryCatch(
      posterior::subset_draws(d, variable = variables),
      error = function(e) {
        cli_abort(c(
          "Stan output required by this method is missing from the fitted model.",
          "i" = "Requested variable{?s}: {.val {variables}}.",
          "i" = "Custom Stan files must preserve the expected parameter and generated-quantity names.",
          "x" = conditionMessage(e)
        ))
      }
    )
  }

  switch(format,
    draws_array = d,
    draws_matrix = posterior::as_draws_matrix(d),
    draws_df = posterior::as_draws_df(d),
    d
  )
}


#' Internal: assign stable names to common summary functions
#' @noRd
.normalize_summary_dots <- function(dots) {
  if (length(dots) == 0L) {
    return(dots)
  }

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }

  for (i in seq_along(dots)) {
    if (nzchar(dot_names[[i]])) {
      next
    }

    dot_i <- dots[[i]]
    if (!is.function(dot_i)) {
      next
    }

    if (identical(dot_i, base::mean)) {
      dot_names[[i]] <- "mean"
    } else if (identical(dot_i, stats::sd)) {
      dot_names[[i]] <- "sd"
    } else if (identical(dot_i, posterior::rhat)) {
      dot_names[[i]] <- "rhat"
    } else if (identical(dot_i, posterior::ess_bulk)) {
      dot_names[[i]] <- "ess_bulk"
    } else if (identical(dot_i, posterior::ess_tail)) {
      dot_names[[i]] <- "ess_tail"
    }
  }

  names(dots) <- dot_names
  dots
}


#' Summarise posterior draws
#'
#' When called without extra summary functions in `...`, returns the standard
#' summary (mean, sd, q2.5, q97.5, rhat, ess_bulk, ess_tail).  When `...`
#' contains summary functions (e.g. `mean, sd, ~quantile2(...)`) they are
#' forwarded to [posterior::summarise_draws()].
#'
#' @param fit Raw backend fit object.
#' @param variables Character vector or `NULL`.
#' @param backend Character.
#' @param ... Optional summary functions forwarded to
#'   [posterior::summarise_draws()].
#' @return A tibble.
#' @noRd
.fit_summary <- function(fit, variables = NULL, backend, ...) {
  dots <- list(...)
  required <- dots$required %||% NULL
  required_type <- dots$required_type %||% c("exact", "pattern")
  context <- dots$context %||% "This method"
  output_type <- dots$output_type %||% "Stan output"

  d <- .fit_draws(
    fit,
    variables = variables,
    format = "draws_array",
    backend = backend,
    required = required,
    required_type = required_type,
    context = context,
    output_type = output_type
  )

  dots$required <- NULL
  dots$required_type <- NULL
  dots$context <- NULL
  dots$output_type <- NULL
  if (length(dots) > 0) {
    dots <- .normalize_summary_dots(dots)
    suppressWarnings(do.call(posterior::summarise_draws, c(list(d), dots)))
  } else {
    suppressWarnings(
      posterior::summarise_draws(
        d,
        mean = mean,
        sd = stats::sd,
        ~posterior::quantile2(.x, probs = c(0.025, 0.975)),
        rhat = posterior::rhat,
        ess_bulk = posterior::ess_bulk,
        ess_tail = posterior::ess_tail
      )
    )
  }
}


#' Extract sampler diagnostics as a 3-D array
#'
#' Returns an array with dimensions `[iteration, chain, parameter]`.
#'
#' @param fit Raw backend fit object.
#' @param backend Character.
#' @return A 3-D numeric array.
#' @noRd
.fit_sampler_diagnostics <- function(fit, backend) {
  if (backend == "cmdstanr") {
    return(fit$sampler_diagnostics())
  }

  if (inherits(fit, "draws")) {
    return(array(
      numeric(0),
      dim = c(0L, 0L, 0L),
      dimnames = list(NULL, NULL, character())
    ))
  }

  # rstan: get_sampler_params returns a list of data frames (one per chain)
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  n_chains <- length(sampler_params)
  if (n_chains == 0L) return(array(numeric(0), dim = c(0L, 0L, 0L)))

  n_iter <- nrow(sampler_params[[1L]])
  param_names <- colnames(sampler_params[[1L]])

  result <- array(NA_real_,
                  dim = c(n_iter, n_chains, length(param_names)),
                  dimnames = list(NULL, NULL, param_names))
  for (ch in seq_len(n_chains)) {
    result[, ch, ] <- as.matrix(sampler_params[[ch]])
  }
  result
}


#' Summarise divergences and treedepth issues
#'
#' @param fit Raw backend fit object.
#' @param backend Character.
#' @return A named list with `num_divergent` and `num_max_treedepth`.
#' @noRd
.fit_diagnostic_summary <- function(fit, backend) {
  if (backend == "cmdstanr") {
    return(fit$diagnostic_summary(quiet = TRUE))
  }

  # rstan path
  sampler_diags <- .fit_sampler_diagnostics(fit, "rstan")
  diag_names <- dimnames(sampler_diags)[[3L]]

  num_divergent <- 0L
  if ("divergent__" %in% diag_names) {
    num_divergent <- sum(sampler_diags[, , "divergent__"])
  }

  num_max_treedepth <- 0L
  if ("treedepth__" %in% diag_names) {
    # rstan stores max_treedepth in the stan_args slot
    max_td <- tryCatch(
      as.numeric(fit@stan_args[[1L]]$control$max_treedepth),
      error = function(e) 10
    )
    if (is.na(max_td)) max_td <- 10
    num_max_treedepth <- sum(sampler_diags[, , "treedepth__"] >= max_td)
  }

  list(
    num_divergent = as.integer(num_divergent),
    num_max_treedepth = as.integer(num_max_treedepth)
  )
}
