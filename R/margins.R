# ============================================================================
# Margin Type Infrastructure
# ============================================================================

#' Valid margin types
#' @noRd
.valid_margins <- c("normal", "exponential", "skew_normal", "gamma")

#' Skew-normal residual SD as a fraction of omega
#' @noRd
.skew_normal_residual_sd_factor <- function(delta) {
  sqrt(1 - 2 * delta^2 / pi)
}

#' Per-dimension margin family codes passed to the generic mixed Stan model
#'
#' These integer codes label each dimension's marginal family for the generic
#' mixed Stan models (`constant_mixed`, `constant_mixed_clayton`,
#' `dcvar_mixed_ncp`, `hmm_mixed`, `multilevel_mixed`, `multilevel_tv_mixed`,
#' `sem_mixed`, `sem_tv_mixed`, and `sem_naive_mixed`). The order must match the
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

#' Does this SEM/multilevel margin spec use the generic mixed engine?
#'
#' Multilevel and SEM have specialised Stan files for homogeneous normal and
#' exponential margins. Homogeneous skew-normal and gamma are served by the
#' generic mixed-margin engines, using the same family code on both dimensions.
#' @noRd
.uses_sem_multilevel_mixed_engine <- function(model_type, margins) {
  if (!model_type %in% c("multilevel", "sem", "sem_naive")) {
    return(FALSE)
  }
  if (.is_mixed_margins(margins)) {
    return(TRUE)
  }

  length(margins) >= 1L &&
    length(unique(margins)) == 1L &&
    margins[[1L]] %in% c("skew_normal", "gamma")
}

#' Margin vector passed to a SEM/multilevel mixed engine
#' @noRd
.sem_multilevel_engine_margins <- function(model_type, margins, D = 2L) {
  if (!.uses_sem_multilevel_mixed_engine(model_type, margins)) {
    return(margins)
  }
  if (.is_mixed_margins(margins)) {
    return(margins)
  }

  rep(margins[[1L]], D)
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

#' Reportable per-state, per-dimension scale/shape variables for the switching HMM
#'
#' Like [.mixed_margin_report_vars()] but lifted to the engine's 2D-indexed
#' `[m, d]` layout: maps each margin config `m` (state when margins switch, else
#' the shared config) and dimension to the interpretable scale/shape variable for
#' its family (e.g. `sigma_exp[m, d]` for an exponential dimension).
#' @noRd
.hmm_switching_report_vars <- function(margins_char, Mrg_K) {
  D <- ncol(margins_char)
  groups <- list()
  add <- function(nm, m, d) {
    groups[[nm]] <<- c(groups[[nm]], paste0(nm, "[", m, ",", d, "]"))
  }
  for (m in seq_len(Mrg_K)) {
    for (d in seq_len(D)) {
      switch(margins_char[m, d],
        normal = add("sigma_eps", m, d),
        exponential = add("sigma_exp", m, d),
        skew_normal = {
          add("omega", m, d)
          add("delta", m, d)
        },
        gamma = {
          add("sigma_gam", m, d)
          add("shape_gam", m, d)
        }
      )
    }
  }
  groups
}

#' Sampled per-state, per-dimension margin parameters for the switching HMM
#'
#' The switching engine declares the margin union as `array[Mrg_K] vector[D]`, so
#' the sampled names are 2D-indexed `[m, d]`. For each margin config `m` (state
#' when margins switch, otherwise the single shared config) and dimension `d`,
#' returns the indexed name of the parameter that config/dimension actually uses
#' for its family. `margins_char` is the K x D character family matrix; only the
#' first `Mrg_K` rows are consulted.
#' @noRd
.hmm_switching_diagnostic_margin_vars <- function(margins_char, Mrg_K) {
  D <- ncol(margins_char)
  out <- character()
  for (m in seq_len(Mrg_K)) {
    for (d in seq_len(D)) {
      v <- switch(margins_char[m, d],
        normal = paste0("sigma_eps[", m, ",", d, "]"),
        exponential = paste0("eta[", m, ",", d, "]"),
        skew_normal = c(paste0("omega[", m, ",", d, "]"), paste0("delta[", m, ",", d, "]")),
        gamma = c(paste0("eta[", m, ",", d, "]"), paste0("shape_gam[", m, ",", d, "]"))
      )
      out <- c(out, v)
    }
  }
  out
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
#' SEM supports all four latent innovation margin families. Homogeneous normal
#' and exponential fits use specialised Stan files; homogeneous skew-normal,
#' homogeneous gamma, and per-variable specs use the generic mixed engines.
#'
#' @param margins Character vector of margin families.
#' @param skew_direction Length-2 integer vector of +1/-1 for exponential or
#'   gamma margins.
#' @return Invisible TRUE if valid.
#' @noRd
.validate_sem_margins <- function(margins, skew_direction = NULL) {
  margins <- .normalize_margins_spec(margins)
  .validate_margins(margins, skew_direction)

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

  # The unified dynamic engine is a single generic file: per-dimension family
  # codes handle homogeneous and mixed margins alike.
  if (identical(model_type, "dcvar_dynamic")) {
    if (!identical(copula, "gaussian")) {
      cli_abort(c(
        "Time-varying VAR components are currently implemented only for the Gaussian copula.",
        "i" = "Use {.code copula = \"gaussian\"} (the default) with {.arg tv_phi} / {.arg tv_sigma}."
      ))
    }
    return("dcvar_dynamic.stan")
  }

  # The Markov-switching HMM engine is a single generic file: per-state, per-
  # dimension family codes handle homogeneous, mixed, and per-state margins alike.
  if (identical(model_type, "hmm_switching")) {
    if (!identical(copula, "gaussian")) {
      cli_abort(c(
        "State-specific HMM switching is currently implemented only for the Gaussian copula.",
        "i" = "Use {.code copula = \"gaussian\"} (the default)."
      ))
    }
    return("hmm_switching.stan")
  }

  # The public TV path remains the legacy-compatible Stan file, matching the
  # exported prepare_dcvar_data(tv_*) output. The dcvar() wrapper routes bundled
  # TV fits to dcvar_dynamic.stan internally.
  if (identical(model_type, "dcvar_tv")) {
    if (!identical(copula, "gaussian")) {
      cli_abort(c(
        "Time-varying VAR components are currently implemented only for the Gaussian copula.",
        "i" = "Use {.code copula = \"gaussian\"} (the default) with {.arg tv_phi} / {.arg tv_sigma}."
      ))
    }
    return("dcvar_tv_mixed.stan")
  }

  if (model_type %in% c("sem_tv", "multilevel_tv")) {
    if (!identical(copula, "gaussian")) {
      cli_abort(c(
        "Time-varying SEM and multilevel components are currently implemented only for the Gaussian copula.",
        "i" = "Use {.code copula = \"gaussian\"} (the default) with {.arg tv_phi} / {.arg tv_sigma}."
      ))
    }
    return(switch(model_type,
      sem_tv = "sem_tv_mixed.stan",
      multilevel_tv = "multilevel_tv_mixed.stan"
    ))
  }

  if (.is_mixed_margins(margins) ||
      .uses_sem_multilevel_mixed_engine(model_type, margins)) {
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
    if (margins %in% c("skew_normal", "gamma")) {
      return("multilevel_mixed.stan")
    }
    cli_abort(
      "Unknown multilevel margin family {.val {margins}}."
    )
  }
  if (identical(model_type, "sem")) {
    if (margins == "normal") {
      return("sem_copula_var.stan")
    }
    if (margins == "exponential") {
      return("sem_EG.stan")
    }
    if (margins %in% c("skew_normal", "gamma")) {
      return("sem_mixed.stan")
    }
    cli_abort(
      "Unknown SEM margin family {.val {margins}}."
    )
  }
  if (identical(model_type, "sem_naive")) {
    if (margins == "normal") {
      return("sem_naive_NG.stan")
    }
    if (margins == "exponential") {
      return("sem_naive_EG.stan")
    }
    if (margins %in% c("skew_normal", "gamma")) {
      return("sem_naive_mixed.stan")
    }
    cli_abort(
      "Unknown naive SEM margin family {.val {margins}}."
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


#' Get a descriptive model key (margin/copula-labelled)
#'
#' Builds a human-readable model key that encodes the margin/copula family. This
#' is a labelling helper only: the actual compiled-model cache is keyed on a
#' content hash of the Stan source and its recursive includes (see
#' [.stan_cache_fingerprint()]), so the same `dcvar_dynamic.stan` engine is
#' compiled once and reused across configurations regardless of this key.
#'
#' @param model_type Character: model family key.
#' @param margins Character: margin type.
#' @param copula Character: copula family.
#' @return Character: cache key for the compiled model.
#' @noRd
.margin_cache_key <- function(model_type, margins, copula = "gaussian") {
  .validate_copula(copula)
  .validate_margin_families(margins)
  if (identical(model_type, "dcvar_dynamic")) {
    margins_vec <- rep(margins, length.out = 2L)
    return(paste0("dcvar_dynamic", paste(.family_codes[margins_vec], collapse = ""), "_model"))
  }
  if (identical(model_type, "hmm_switching")) {
    margins_vec <- rep(margins, length.out = 2L)
    return(paste0("hmm_switching", paste(.family_codes[margins_vec], collapse = ""), "_model"))
  }
  if (identical(model_type, "dcvar_tv")) {
    margins_vec <- rep(margins, length.out = 2L)
    return(paste0("dcvar_tv_mixed", paste(.family_codes[margins_vec], collapse = ""), "_model"))
  }
  if (.is_mixed_margins(margins) ||
      .uses_sem_multilevel_mixed_engine(model_type, margins)) {
    margins_vec <- .sem_multilevel_engine_margins(model_type, margins)
    suffix <- paste0("_mixed", paste(.family_codes[margins_vec], collapse = ""))
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
