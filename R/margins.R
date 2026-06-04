# ============================================================================
# Margin Type Infrastructure
# ============================================================================

#' Valid margin types
#' @noRd
.valid_margins <- c("normal", "exponential", "skew_normal", "gamma")

#' Per-dimension margin family codes passed to the generic mixed Stan model
#'
#' These integer codes label each dimension's marginal family for
#' `constant_mixed.stan` (and future generic models). The order must match the
#' `family[i] == k` dispatch in the Stan code.
#' @noRd
.family_codes <- c(normal = 1L, exponential = 2L, skew_normal = 3L, gamma = 4L)

#' Collapse an all-identical margins vector to a single string
#'
#' Per-dimension margins are accepted as a length-1 (recycled) or length-D
#' character vector. When every entry is identical the specification is
#' equivalent to the scalar form, so it is collapsed back to a single string.
#' This keeps the homogeneous code paths (specialised Stan models, scalar
#' `margins` attribute, all existing extractors) completely unchanged. The
#' result is therefore always either length-1 (homogeneous) or a length-D vector
#' with at least two distinct families (genuinely mixed).
#' @noRd
.normalize_margins_spec <- function(margins) {
  if (length(margins) > 1L && length(unique(margins)) == 1L) {
    margins[[1L]]
  } else {
    margins
  }
}

#' Is this a genuinely mixed (per-variable) margin specification?
#'
#' Returns `TRUE` only when the margins vector contains more than one distinct
#' family, i.e. the fit must route to the generic mixed Stan model.
#' @noRd
.is_mixed_margins <- function(margins) {
  length(margins) > 1L && length(unique(margins)) > 1L
}

#' Sampled per-dimension margin parameters to monitor for a mixed fit
#'
#' Returns the indexed names of the *sampled* parameter each dimension uses for
#' its family (the union's unused parameters merely sample from their priors and
#' need not be monitored for convergence). Used by the diagnostics machinery.
#' @noRd
.mixed_diagnostic_margin_vars <- function(margins) {
  unlist(lapply(seq_along(margins), function(i) {
    switch(margins[[i]],
      normal = paste0("sigma_eps[", i, "]"),
      exponential = paste0("eta[", i, "]"),
      skew_normal = c(paste0("omega[", i, "]"), paste0("delta[", i, "]")),
      gamma = c(paste0("eta[", i, "]"), paste0("shape_gam[", i, "]"))
    )
  }))
}

#' Per-dimension scale/shape variables to display in trace plots for a mixed fit
#'
#' Returns the indexed names of the interpretable scale/shape variable each
#' dimension uses (e.g. `sigma_exp[d]` for an exponential dimension). Used by the
#' diagnostic plots.
#' @noRd
.mixed_plot_margin_vars <- function(margins) {
  unlist(lapply(seq_along(margins), function(i) {
    switch(margins[[i]],
      normal = paste0("sigma_eps[", i, "]"),
      exponential = paste0("sigma_exp[", i, "]"),
      skew_normal = c(paste0("omega[", i, "]"), paste0("delta[", i, "]")),
      gamma = c(paste0("sigma_gam[", i, "]"), paste0("shape_gam[", i, "]"))
    )
  }))
}

#' Per-dimension scale/shape variables to report for a mixed fit
#'
#' Maps a mixed margins vector to the indexed Stan variable names that should be
#' reported for each dimension, restricted to the dimensions of each family.
#' For example `c("normal", "exponential")` yields
#' `list(sigma_eps = "sigma_eps[1]", sigma_exp = "sigma_exp[2]")`.
#' @noRd
.mixed_margin_report_vars <- function(margins) {
  out <- list()
  add <- function(name, dims) {
    if (length(dims) > 0L) out[[name]] <<- paste0(name, "[", dims, "]")
  }
  add("sigma_eps", which(margins == "normal"))
  add("sigma_exp", which(margins == "exponential"))
  add("omega", which(margins == "skew_normal"))
  add("delta", which(margins == "skew_normal"))
  add("sigma_gam", which(margins == "gamma"))
  add("shape_gam", which(margins == "gamma"))
  out
}

#' Valid copula families
#' @noRd
.valid_copulas <- c("gaussian", "clayton")

#' Validate copula specification
#'
#' @param copula Character string: one of `"gaussian"` or `"clayton"`.
#' @return Invisible TRUE if valid.
#' @noRd
.validate_copula <- function(copula) {
  if (!is.character(copula) || length(copula) != 1) {
    cli_abort("{.arg copula} must be a single character string.")
  }
  if (!copula %in% .valid_copulas) {
    cli_abort(
      "{.arg copula} must be one of {.val {(.valid_copulas)}}, got {.val {copula}}."
    )
  }

  invisible(TRUE)
}

#' Validate margin family names (length and membership only)
#'
#' Checks that `margins` is a length-1 or length-2 character vector of known
#' families, without requiring `skew_direction`. Used by routing/path helpers
#' (such as [dcvar_stan_path()]) where the marginal scale parameters are
#' irrelevant to file lookup.
#' @noRd
.validate_margin_families <- function(margins) {
  if (!is.character(margins) || length(margins) < 1L || length(margins) > 2L) {
    cli_abort("{.arg margins} must be a character vector of length 1 or 2.")
  }
  invalid <- setdiff(margins, .valid_margins)
  if (length(invalid) > 0L) {
    cli_abort(
      "{.arg margins} must be one of {.val {(.valid_margins)}}, got {.val {invalid}}."
    )
  }
  invisible(TRUE)
}

#' Validate margin specification
#'
#' Accepts either a single margin string (applied to both variables) or a
#' length-2 vector specifying a per-variable (mixed) margin for the bivariate
#' model. `skew_direction` is required whenever *any* dimension uses an
#' exponential or gamma margin.
#'
#' @param margins Character vector of length 1 or 2; each entry one of
#'   "normal", "exponential", "skew_normal", "gamma".
#' @param skew_direction Length-2 integer vector of +1/-1. Required when any
#'   dimension uses an exponential or gamma margin.
#' @return Invisible TRUE if valid.
#' @noRd
.validate_margins <- function(margins, skew_direction = NULL) {
  .validate_margin_families(margins)

  needs_skew <- intersect(margins, c("exponential", "gamma"))
  if (length(needs_skew) > 0L) {
    if (is.null(skew_direction)) {
      cli_abort(
        "{.arg skew_direction} is required for {.val {needs_skew}} margins."
      )
    }
    if (length(skew_direction) != 2 || !all(skew_direction %in% c(-1, 1))) {
      cli_abort(
        "{.arg skew_direction} must be a length-2 vector of +1 or -1."
      )
    }
  }

  invisible(TRUE)
}


#' Validate SEM margin specification
#'
#' SEM currently supports the normal and exponential latent innovation
#' parameterizations only.
#'
#' @param margins Character string: one of `"normal"` or `"exponential"`.
#' @param skew_direction Length-2 integer vector of +1/-1 for exponential
#'   margins.
#' @return Invisible TRUE if valid.
#' @noRd
.validate_sem_margins <- function(margins, skew_direction = NULL) {
  margins <- .normalize_margins_spec(margins)
  .validate_margins(margins, skew_direction)

  # Homogeneous SEM fits use the specialised normal/exponential models; mixed
  # (per-variable) fits use the generic sem_mixed / sem_naive_mixed models,
  # which support all four families per dimension.
  if (!.is_mixed_margins(margins) && !all(margins %in% c("normal", "exponential"))) {
    cli_abort(c(
      "Single-family {.fun dcvar_sem} supports only {.val {c('normal', 'exponential')}} margins.",
      "i" = "Use a per-variable {.arg margins} vector (e.g. {.code c('normal', 'gamma')}) for other families."
    ))
  }

  invisible(TRUE)
}


#' Get Stan model suffix for a given margin type
#'
#' Maps margin names to Stan file suffixes used in model file naming.
#'
#' @param margins Character: margin type.
#' @return Character: suffix for Stan file (e.g., "_EG", "_SNG", "_GG"),
#'   or empty string for normal margins.
#' @noRd
.margin_stan_suffix <- function(margins) {
  switch(margins,
    normal = "",
    exponential = "_EG",
    skew_normal = "_SNG",
    gamma = "_GG",
    cli_abort("Unknown margin type: {.val {margins}}")
  )
}


#' Get Stan file name for a given model type and margin
#'
#' @param model_type Character: model family key.
#' @param margins Character: margin type.
#' @param copula Character: copula family.
#' @return Character: Stan file name (without path).
#' @noRd
.margin_stan_file <- function(model_type, margins, copula = "gaussian") {
  .validate_copula(copula)
  .validate_margin_families(margins)

  if (.is_mixed_margins(margins)) {
    if (identical(copula, "clayton")) {
      if (!identical(model_type, "constant")) {
        cli_abort(c(
          "Mixed margins with the Clayton copula are currently implemented only for the constant model.",
          "i" = "Use {.fun dcvar_constant} with {.code copula = \"clayton\"}."
        ))
      }
      return("constant_mixed_clayton.stan")
    }
    mixed_file <- switch(model_type,
      constant = "constant_mixed.stan",
      dcvar = "dcvar_mixed_ncp.stan",
      hmm = "hmm_mixed.stan",
      multilevel = "multilevel_mixed.stan",
      sem = "sem_mixed.stan",
      sem_naive = "sem_naive_mixed.stan",
      cli_abort(c(
        "Mixed per-variable margins are not implemented for the {.val {model_type}} model.",
        "i" = "Supported: {.val constant}, {.val dcvar}, {.val hmm}, {.val multilevel}, {.val sem}, {.val sem_naive}."
      ))
    )
    return(mixed_file)
  }
  margins <- margins[[1L]]

  if (identical(copula, "clayton")) {
    if (identical(model_type, "constant") && identical(margins, "normal")) {
      return("constant_NCl.stan")
    }
    cli_abort(c(
      "The Clayton copula is currently implemented only for the constant model with normal margins.",
      "i" = "Use {.code model_type = 'constant'}, {.code margins = 'normal'}, and {.code copula = 'clayton'}."
    ))
  }

  base_files <- c(
    constant = "constant_copula_var",
    dcvar = "dcvar_model_ncp",
    hmm = "hmm_copula_model",
    sem = "sem_copula_var"
  )
  if (identical(model_type, "multilevel")) {
    if (margins == "normal") {
      return("multilevel_copula_var.stan")
    }
    if (margins == "exponential") {
      return("multilevel_EG.stan")
    }
    cli_abort(
      "Multilevel Stan models currently support only {.val {c('normal', 'exponential')}} margins, got {.val {margins}}."
    )
  }
  if (identical(model_type, "sem")) {
    if (margins == "normal") {
      return("sem_copula_var.stan")
    }
    if (margins == "exponential") {
      return("sem_EG.stan")
    }
    cli_abort(
      "SEM Stan models currently support only {.val {c('normal', 'exponential')}} margins, got {.val {margins}}."
    )
  }
  if (identical(model_type, "sem_naive")) {
    if (margins == "normal") {
      return("sem_naive_NG.stan")
    }
    if (margins == "exponential") {
      return("sem_naive_EG.stan")
    }
    cli_abort(
      "Naive SEM Stan models currently support only {.val {c('normal', 'exponential')}} margins, got {.val {margins}}."
    )
  }
  if (margins == "normal") {
    paste0(base_files[model_type], ".stan")
  } else {
    suffix <- .margin_stan_suffix(margins)
    base_short <- c(
      constant = "constant",
      dcvar = "dcvar",
      hmm = "hmm"
    )
    # dcvar non-normal Stan files have _ncp suffix (e.g., dcvar_EG_ncp.stan)
    ncp <- if (model_type == "dcvar") "_ncp" else ""
    paste0(base_short[model_type], suffix, ncp, ".stan")
  }
}


#' Get the cache key for a compiled model
#'
#' Includes margin type to prevent cache collisions between different
#' margin specifications of the same base model.
#'
#' @param model_type Character: model family key.
#' @param margins Character: margin type.
#' @param copula Character: copula family.
#' @return Character: cache key for the compiled model.
#' @noRd
.margin_cache_key <- function(model_type, margins, copula = "gaussian") {
  .validate_copula(copula)
  .validate_margin_families(margins)
  if (.is_mixed_margins(margins)) {
    suffix <- paste0("_mixed", paste(.family_codes[margins], collapse = ""))
    if (!identical(copula, "gaussian")) {
      suffix <- paste0(suffix, "_", copula)
    }
    return(paste0(model_type, suffix, "_model"))
  }
  margins <- margins[[1L]]
  suffix <- .margin_stan_suffix(margins)
  if (!identical(copula, "gaussian")) {
    suffix <- paste0(suffix, "_", copula)
  }
  paste0(model_type, suffix, "_model")
}
