# ============================================================================
# Margin Type Infrastructure
# ============================================================================

#' Valid margin types
#' @noRd
.valid_margins <- c("normal", "exponential", "skew_normal", "gamma")

#' Validate margin specification
#'
#' @param margins Character string: one of "normal", "exponential",
#'   "skew_normal", "gamma".
#' @param skew_direction Length-2 integer vector of +1/-1. Required for
#'   exponential and gamma margins.
#' @return Invisible TRUE if valid.
#' @noRd
.validate_margins <- function(margins, skew_direction = NULL) {
  if (!is.character(margins) || length(margins) != 1) {
    cli_abort("{.arg margins} must be a single character string.")
  }
  if (!margins %in% .valid_margins) {
    cli_abort(
      "{.arg margins} must be one of {.val {(.valid_margins)}}, got {.val {margins}}."
    )
  }

  if (margins %in% c("exponential", "gamma")) {
    if (is.null(skew_direction)) {
      cli_abort(
        "{.arg skew_direction} is required for {.val {margins}} margins."
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
  .validate_margins(margins, skew_direction)

  if (!margins %in% c("normal", "exponential")) {
    cli_abort(
      "{.arg margins} for {.fun dcvar_sem} must be one of {.val {c('normal', 'exponential')}}, got {.val {margins}}."
    )
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
#' @param model_type Character: "constant", "dcvar", or "hmm".
#' @param margins Character: margin type.
#' @return Character: Stan file name (without path).
#' @noRd
.margin_stan_file <- function(model_type, margins) {
  base_files <- c(
    constant = "constant_copula_var",
    dcvar = "dcvar_model_ncp",
    hmm = "hmm_copula_model",
    sem = "sem_copula_var"
  )
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
#' @param model_type Character: "constant", "dcvar", or "hmm".
#' @param margins Character: margin type.
#' @return Character: cache key for the compiled model.
#' @noRd
.margin_cache_key <- function(model_type, margins) {
  suffix <- .margin_stan_suffix(margins)
  paste0(model_type, suffix, "_model")
}
