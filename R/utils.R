# ============================================================================
# Internal Utility Functions
# ============================================================================

#' Safe correlation that returns NA for constant or too-short vectors
#'
#' @param x,y Numeric vectors.
#' @return Correlation coefficient, or `NA_real_` if either vector has zero
#'   variance or fewer than two observations (where `sd()` returns `NA` and
#'   `cor()` errors).
#' @noRd
.safe_cor <- function(x, y) {
  if (length(x) < 2L || length(y) < 2L) {
    return(NA_real_)
  }
  if (sd(x) < 1e-10 || sd(y) < 1e-10) {
    return(NA_real_)
  }

  cor(x, y)
}

#' Validate a credible/prediction interval level
#'
#' @param level Numeric scalar interval level.
#' @param arg_name Name of the argument for error messages.
#' @return Invisibly returns `level`.
#' @noRd
.validate_interval_level <- function(level, arg_name = "ci_level") {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    cli_abort("{.arg {arg_name}} must be a single numeric value strictly between 0 and 1.")
  }

  invisible(level)
}

#' Build a regex that matches a Stan output group by base name
#'
#' @param name Stan variable name or indexed element.
#' @return Regex matching the base output name with or without indices.
#' @noRd
.stan_output_group_pattern <- function(name) {
  base_name <- sub("\\[.*$", "", name)
  escaped <- gsub("([][{}()+*^$.|?\\\\])", "\\\\\\1", base_name)
  paste0("^", escaped, "(\\[|$)")
}

#' Internal: extract required coefficient summaries
#' @noRd
.extract_required_coef <- function(summ, pattern, label = NULL, context = NULL) {
  rows <- grep(pattern, summ$variable)
  if (length(rows) == 0) {
    if (is.null(label)) label <- pattern
    if (is.null(context)) context <- "Coefficient extraction"
    cli_abort(c(
      "{context} requires Stan output that is not present in the fitted model.",
      "i" = "Missing coefficient group: {.val {label}}.",
      "i" = "Custom Stan files must preserve the expected parameter names."
    ))
  }

  setNames(summ$mean[rows], summ$variable[rows])
}

#' Internal: validate that a fitted model exposes required Stan outputs
#' @noRd
.validate_required_stan_outputs <- function(draws, required = NULL,
                                            required_type = c("exact", "pattern"),
                                            context = "This method",
                                            output_type = "Stan output") {
  if (is.null(required) || length(required) == 0L) {
    return(invisible(TRUE))
  }

  required_type <- match.arg(required_type)
  vars <- posterior::variables(draws)

  missing <- switch(required_type,
    exact = setdiff(required, vars),
    pattern = required[!vapply(required, function(pattern) {
      any(grepl(pattern, vars))
    }, logical(1))]
  )

  if (length(missing) > 0) {
    output_label <- switch(output_type,
      "generated quantity" = "generated quantities",
      "generated quantities" = "generated quantities",
      "parameter group" = "parameter groups",
      "parameter" = "parameters",
      "parameters" = "parameters",
      "Stan output" = "Stan output",
      output_type
    )

    cli_abort(c(
      "{context} requires Stan output that is not present in the fitted model.",
      "i" = "Missing {output_label}: {.val {missing}}.",
      "i" = "Custom Stan files must preserve the expected parameter and generated-quantity names."
    ))
  }

  invisible(TRUE)
}


# ============================================================================
# Shared Initialization Helpers
# ============================================================================

# Default VAR initialization values shared across all model types.
# These are moderate starting values chosen for standardised data:
#   - Phi diagonal 0.25: mild autoregression, well within stationary region
#   - Phi off-diagonal jitter SD 0.05: small cross-lag starting values
#   - sigma_eps in [0.8, 1.2]: near unit variance (appropriate for z-scored data)
#   - mu near zero: centred around zero (appropriate for z-scored data)

#' Generate default VAR initialization values
#'
#' @param D Number of variables.
#' @param margins Margin type.
#' @return A named list with `mu`, `Phi`, and margin-specific scale params.
#' @noRd
.init_var_params <- function(D, margins = "normal") {
  base <- list(
    mu = rnorm(D, 0, 0.1),
    Phi = diag(0.25, D) + matrix(rnorm(D^2, 0, 0.05), D, D)
  )

  # Mixed (per-variable) margins use the generic Stan model, whose parameter
  # block is the union of all per-family parameters. Initialise the whole union
  # so every sampled parameter (including those that only sample from their
  # prior on a given dimension) starts from a sensible value. shape_gam is a
  # per-dimension vector in the generic model (unlike the shared scalar in the
  # specialised gamma model).
  if (.is_mixed_margins(margins)) {
    return(c(base, list(
      sigma_eps = runif(D, 0.8, 1.2),
      eta = rnorm(D, 0, 0.3),
      omega = runif(D, 0.5, 1.5),
      delta = runif(D, -0.3, 0.3),
      shape_gam = runif(D, 0.5, 2.0)
    )))
  }

  switch(margins[[1L]],
    normal = c(base, list(sigma_eps = runif(D, 0.8, 1.2))),
    exponential = c(base, list(eta = rnorm(D, 0, 0.3))),
    skew_normal = c(base, list(
      omega = runif(D, 0.5, 1.5),
      delta = runif(D, -0.3, 0.3)
    )),
    gamma = c(base, list(
      eta = rnorm(D, 0, 0.3),
      shape_gam = runif(1, 0.5, 2.0)
    )),
    c(base, list(sigma_eps = runif(D, 0.8, 1.2)))
  )
}


#' Generate default DC-VAR initialization values
#'
#' @param D Number of variables.
#' @param T_obs Number of time points.
#' @param margins Margin type.
#' @return A named list with VAR params plus `sigma_omega`, `z_rho_init`,
#'   `omega_raw`.
#' @noRd
.init_dcvar_params <- function(D, T_obs, margins = "normal") {
  c(
    .init_var_params(D, margins),
    list(
      sigma_omega = runif(1, 0.05, 0.15),
      z_rho_init = rnorm(1, 0, 0.1),
      omega_raw = rnorm(T_obs - 1, 0, 0.1)
    )
  )
}


#' Row-major VAR(1) coefficient labels (phi11, phi12, phi21, phi22)
#' @noRd
.phi_coef_names <- c("phi11", "phi12", "phi21", "phi22")

#' Resolve a `tv_phi` specification to a length-4 coefficient mask
#'
#' Accepts a logical scalar (`TRUE` = all four coefficients vary, `FALSE` =
#' none) or a character selector: `"ar"` (the autoregressive coefficients
#' phi11, phi22), `"cross"` / `"crosslag"` (the cross-lagged coefficients
#' phi12, phi21), or specific coefficient names from
#' `c("phi11", "phi12", "phi21", "phi22")`. Returns a named length-4 integer
#' vector in row-major order.
#' @noRd
.resolve_phi_tv_mask <- function(tv_phi) {
  mask <- stats::setNames(integer(4), .phi_coef_names)

  # Idempotent: an already-resolved length-4 0/1 mask is returned as-is
  # (so callers can pass either the user spec or the resolved mask).
  if (is.numeric(tv_phi) && !is.logical(tv_phi) &&
      length(tv_phi) == 4L && all(tv_phi %in% c(0, 1))) {
    mask[] <- as.integer(tv_phi)
    return(mask)
  }

  if (is.logical(tv_phi)) {
    if (length(tv_phi) != 1L || is.na(tv_phi)) {
      cli_abort("{.arg tv_phi} must be a single non-missing logical, or a character selector.")
    }
    mask[] <- as.integer(tv_phi)
    return(mask)
  }

  if (is.character(tv_phi)) {
    if (length(tv_phi) == 0L) {
      cli_abort("{.arg tv_phi} must not be an empty character vector; use {.code FALSE} for no time variation.")
    }
    selected <- unlist(lapply(tv_phi, function(token) {
      switch(token,
        ar = c("phi11", "phi22"),
        cross = ,
        crosslag = ,
        cl = c("phi12", "phi21"),
        token
      )
    }), use.names = FALSE)
    invalid <- setdiff(selected, .phi_coef_names)
    if (length(invalid) > 0L) {
      cli_abort(c(
        "Invalid {.arg tv_phi} selector{?s}: {.val {invalid}}.",
        "i" = "Use {.val TRUE}/{.val FALSE}, {.val ar}, {.val cross}, or names from {.val {(.phi_coef_names)}}."
      ))
    }
    mask[.phi_coef_names %in% selected] <- 1L
    return(mask)
  }

  cli_abort("{.arg tv_phi} must be a logical scalar or a character selector ({.val ar}, {.val cross}, or coefficient names).")
}

#' Resolve a `switch` specification for the Markov-switching HMM
#'
#' Maps the public `switch` argument of [dcvar_hmm()] to canonical flags. Accepts
#' `TRUE` (every component switches) or a character vector drawn from
#' `c("rho", "mu", "phi", "ar", "cross", "sigma")` plus explicit Phi coefficient
#' names. `rho` is mandatory (it is the `ordered[K]` label-switching anchor): a
#' selector that switches other components without `rho` is rejected, and a
#' selector that switches nothing is rejected (use [dcvar_constant()]). Phi
#' granularity reuses [.resolve_phi_tv_mask()].
#'
#' @param switch The user's `switch` argument.
#' @return A named list with `mu` (0/1), `phi_mask` (named length-4 integer),
#'   `rho` (always 1L), and `margins` (0/1, requesting state-specific scales).
#' @noRd
.resolve_switch_spec <- function(switch) {
  if (isTRUE(switch)) {
    return(list(mu = 1L, phi_mask = .resolve_phi_tv_mask(TRUE), rho = 1L, margins = 1L))
  }
  if (isFALSE(switch) ||
      (is.character(switch) && (length(switch) == 0L || identical(switch, "none")))) {
    cli_abort(c(
      "{.arg switch} must select at least one switching component.",
      "i" = "A Markov-switching model with no state-specific parameter is not identified; use {.fun dcvar_constant} or {.fun dcvar_hmm} with {.code switch = \"rho\"}."
    ))
  }
  if (!is.character(switch)) {
    cli_abort("{.arg switch} must be {.val TRUE} or a character selector (e.g. {.val rho}, {.val mu}, {.val phi}).")
  }

  valid <- c("rho", "mu", "sigma", "margins", "phi",
             "ar", "cross", "crosslag", "cl", .phi_coef_names)
  invalid <- setdiff(switch, valid)
  if (length(invalid) > 0L) {
    cli_abort(c(
      "Invalid {.arg switch} selector{?s}: {.val {invalid}}.",
      "i" = "Use {.val rho}, {.val mu}, {.val phi} (or {.val ar}/{.val cross}/coefficient names), {.val sigma}, or {.val TRUE}."
    ))
  }

  mu <- as.integer("mu" %in% switch)
  margins <- as.integer(any(c("sigma", "margins") %in% switch))
  phi_tokens <- intersect(switch, c("phi", "ar", "cross", "crosslag", "cl", .phi_coef_names))
  phi_mask <- if (length(phi_tokens) == 0L) {
    stats::setNames(integer(4), .phi_coef_names)
  } else if ("phi" %in% phi_tokens) {
    .resolve_phi_tv_mask(TRUE)
  } else {
    .resolve_phi_tv_mask(phi_tokens)
  }
  any_phi <- sum(phi_mask) > 0L

  if (!("rho" %in% switch)) {
    cli_abort(c(
      "{.arg switch} must include {.val rho}: the copula correlation is the label-switching anchor.",
      "i" = "State-specific {.val mu}/{.val phi}/{.val sigma} without an ordered {.code rho} are not identified."
    ))
  }

  list(mu = mu, phi_mask = phi_mask, rho = 1L, margins = margins)
}

#' Generate default TV-VAR initialization values
#'
#' The generic time-varying model declares the full mixed-margin parameter
#' union plus conditionally sized random-walk containers. rstan requires every
#' supplied init entry to match the declared dimensions exactly, so the
#' deviation entries are zero-length when their flag is off.
#'
#' @param D Number of variables.
#' @param T_obs Number of time points.
#' @param margins Margin type (any spec; the union is always initialised).
#' @param tv_phi Logical scalar or character selector (see
#'   [.resolve_phi_tv_mask()]); `tv_sigma` is a logical flag.
#' @param tv_sigma Logical; time-varying residual scales.
#' @return A named list of initial values matching the conditionally sized
#'   Stan parameters (zero-length where a component is off).
#' @noRd
.init_dcvar_tv_params <- function(D, T_obs, margins = "normal",
                                  tv_phi = FALSE, tv_sigma = FALSE) {
  n_eff <- T_obs - 1L
  phi_mask <- .resolve_phi_tv_mask(tv_phi)
  n_phi_tv <- sum(phi_mask)
  any_phi <- n_phi_tv > 0L

  base <- list(
    mu = rnorm(D, 0, 0.1),
    Phi = diag(0.25, D) + matrix(rnorm(D^2, 0, 0.05), D, D),
    # Full mixed-margin union (the generic model declares every family's
    # parameters regardless of the requested margins)
    sigma_eps = runif(D, 0.8, 1.2),
    eta = rnorm(D, 0, 0.3),
    omega = runif(D, 0.5, 1.5),
    delta = runif(D, -0.3, 0.3),
    shape_gam = runif(D, 0.5, 2.0),
    sigma_omega = runif(1, 0.05, 0.15),
    z_rho_init = rnorm(1, 0, 0.1),
    omega_raw = rnorm(n_eff, 0, 0.1)
  )

  # Only supply inits for the ACTIVE components. Inactive components are
  # zero-sized Stan parameters; omitting them lets the sampler initialise them
  # trivially and avoids passing a zero-extent matrix, which cmdstanr
  # serialises to a shapeless empty array (CmdStan then aborts with a
  # "declared (0, D); found (0)" dimension mismatch). Partial inits are
  # accepted by both backends.
  tv <- list()
  if (any_phi) {
    tv$tau_phi <- runif(n_phi_tv, 0.02, 0.08)
    tv$phi_raw <- matrix(rnorm(n_eff * n_phi_tv, 0, 0.1), n_eff, n_phi_tv)
  }
  if (tv_sigma) {
    tv$tau_sigma <- runif(D, 0.02, 0.08)
    tv$sigma_raw <- matrix(rnorm(n_eff * D, 0, 0.1), n_eff, D)
  }

  c(base, tv)
}


#' Generate default initialization values for the unified dynamic engine
#'
#' Inits for `dcvar_dynamic.stan`, which always carries residual drift, the full
#' mixed-margin parameter union, the covariate effects, and the conditionally
#' sized time-varying random-walk containers. As with [.init_dcvar_tv_params()],
#' only the ACTIVE components are supplied: inactive ones (e.g. `beta` when
#' `P = 0`, the TV deviation containers when their flag is off) are omitted so a
#' zero-extent matrix is never serialised (cmdstanr aborts on those). Length-1
#' vectors are wrapped with `array(..., dim = n)` so they keep their declared
#' shape. The sampled intercept is `z_rho_init` (the engine exposes `beta_0` as
#' a transformed-parameter alias, so it is not an init).
#'
#' @param D Number of variables.
#' @param T_obs Number of time points.
#' @param margins Margin spec (any; the union is always initialised).
#' @param tv_phi Logical scalar or character selector (see
#'   [.resolve_phi_tv_mask()]).
#' @param tv_sigma Logical; time-varying residual scales.
#' @param P Number of covariates (0 omits `beta`).
#' @param zero_init_eta Logical; whether `eta[1]` is fixed at zero (sets the
#'   `omega_raw` length to `(T_obs - 1) - zero_init_eta`).
#' @return A named list of initial values matching the engine's conditionally
#'   sized parameters.
#' @noRd
.init_dcvar_dynamic_params <- function(D, T_obs, margins = "normal",
                                       tv_phi = FALSE, tv_sigma = FALSE,
                                       P = 0L, zero_init_eta = FALSE) {
  n_eff <- T_obs - 1L
  phi_mask <- .resolve_phi_tv_mask(tv_phi)
  n_phi_tv <- sum(phi_mask)
  any_phi <- n_phi_tv > 0L
  n_omega <- n_eff - as.integer(isTRUE(zero_init_eta))

  base <- list(
    mu = rnorm(D, 0, 0.1),
    Phi = diag(0.25, D) + matrix(rnorm(D^2, 0, 0.05), D, D),
    # Full mixed-margin union (the engine declares every family's parameters).
    sigma_eps = runif(D, 0.8, 1.2),
    eta = rnorm(D, 0, 0.3),
    omega = runif(D, 0.5, 1.5),
    delta = runif(D, -0.3, 0.3),
    shape_gam = runif(D, 0.5, 2.0),
    sigma_omega = runif(1, 0.05, 0.15),
    z_rho_init = rnorm(1, 0, 0.1),
    omega_raw = array(rnorm(n_omega, 0, 0.1), dim = n_omega)
  )

  extra <- list()
  if (P > 0L) {
    extra$beta <- array(rnorm(P, 0, 0.05), dim = P)
  }
  if (any_phi) {
    extra$tau_phi <- runif(n_phi_tv, 0.02, 0.08)
    extra$phi_raw <- matrix(rnorm(n_eff * n_phi_tv, 0, 0.1), n_eff, n_phi_tv)
  }
  if (tv_sigma) {
    extra$tau_sigma <- runif(D, 0.02, 0.08)
    extra$sigma_raw <- matrix(rnorm(n_eff * D, 0, 0.1), n_eff, D)
  }

  c(base, extra)
}


#' Generate default covariate DC-VAR initialization values
#'
#' @param D Number of variables.
#' @param T_obs Number of time points.
#' @param P Number of covariates.
#' @param drift Logical; include residual random-walk drift.
#' @param zero_init_eta Logical; whether `eta[1]` is fixed at zero.
#' @return A named list with VAR params, covariate effects, and optionally
#'   residual drift parameters.
#' @noRd
.init_dcvar_covariate_params <- function(D, T_obs, P, drift = TRUE,
                                         zero_init_eta = TRUE) {
  base <- c(
    .init_var_params(D, "normal"),
    list(
      beta_0 = rnorm(1, 0, 0.1),
      beta = array(rnorm(P, 0, 0.05), dim = P)
    )
  )

  if (!isTRUE(drift)) {
    return(base)
  }

  n_omega <- (T_obs - 1L) - as.integer(isTRUE(zero_init_eta))
  c(
    base,
    list(
      sigma_omega = runif(1, 0.05, 0.15),
      omega_raw = array(rnorm(n_omega, 0, 0.1), dim = n_omega)
    )
  )
}


#' Generate default HMM initialization values
#'
#' @param D Number of variables.
#' @param K Number of hidden states.
#' @param margins Margin type.
#' @return A named list with VAR params plus `z_rho`, `pi0`, `A`.
#' @noRd
.init_hmm_params <- function(D, K, margins = "normal") {
  # Ordered z_rho: evenly spaced with jitter, sorted to satisfy ordering constraint
  z_rho_init <- sort(seq(0.2, 0.8, length.out = K) + rnorm(K, 0, 0.1))

  # Near-identity transition matrix: ~95% self-transition with small jitter
  A_init <- matrix(0.05 / (K - 1), K, K)
  diag(A_init) <- 0.95
  A_init <- A_init + matrix(runif(K * K, -0.01, 0.01), K, K)
  A_init <- pmax(A_init, 0.01)
  A_init <- A_init / rowSums(A_init)

  c(
    .init_var_params(D, margins),
    list(
      z_rho = z_rho_init,
      pi0 = rep(1 / K, K),
      A = A_init
    )
  )
}


#' Generate default initialization values for the Markov-switching HMM engine
#'
#' Inits for `hmm_switching.stan`. State-indexed groups are sized to the engine's
#' conditional extent (`K` when switching, `1` when shared); the full mixed-margin
#' union is initialised per margin config; `Phi_dev` (the per-state coefficient
#' deviations) is supplied only when any Phi coefficient switches, so a zero-extent
#' array is never serialised (cmdstanr hazard). `z_rho` is the ordered anchor; the
#' transition matrix uses a near-identity sticky init.
#'
#' @param D Number of variables.
#' @param K Number of hidden states.
#' @param switch_spec Resolved spec from [.resolve_switch_spec()].
#' @return A named list of initial values matching the engine's parameters.
#' @noRd
.init_hmm_switching_params <- function(D, K, switch_spec) {
  Mu_K <- if (isTRUE(switch_spec$mu == 1L)) K else 1L
  Mrg_K <- if (isTRUE(switch_spec$margins == 1L)) K else 1L
  any_phi <- sum(switch_spec$phi_mask) > 0L

  z_rho_init <- sort(seq(0.2, 0.8, length.out = K) + rnorm(K, 0, 0.1))

  A_init <- matrix(0.05 / (K - 1), K, K)
  diag(A_init) <- 0.95
  A_init <- A_init + matrix(runif(K * K, -0.01, 0.01), K, K)
  A_init <- pmax(A_init, 0.01)
  A_init <- A_init / rowSums(A_init)

  init <- list(
    # array[Mu_K] vector[D]  -> (Mu_K x D)
    mu = array(rnorm(Mu_K * D, 0, 0.1), dim = c(Mu_K, D)),
    Phi_base = diag(0.25, D) + matrix(rnorm(D^2, 0, 0.05), D, D),
    z_rho = z_rho_init,
    pi0 = rep(1 / K, K),
    A = A_init,
    # Full mixed-margin union per config (array[Mrg_K] vector[D] -> Mrg_K x D).
    sigma_eps = array(runif(Mrg_K * D, 0.8, 1.2), dim = c(Mrg_K, D)),
    eta = array(rnorm(Mrg_K * D, 0, 0.3), dim = c(Mrg_K, D)),
    omega = array(runif(Mrg_K * D, 0.5, 1.5), dim = c(Mrg_K, D)),
    delta = array(runif(Mrg_K * D, -0.3, 0.3), dim = c(Mrg_K, D)),
    shape_gam = array(runif(Mrg_K * D, 0.5, 2.0), dim = c(Mrg_K, D))
  )
  # Per-state Phi deviations only when a coefficient switches (else zero-extent).
  if (any_phi) {
    init$Phi_dev <- array(rnorm(K * D * D, 0, 0.05), dim = c(K, D, D))
  }
  init
}


#' Generate default constant copula initialization values
#'
#' @param D Number of variables.
#' @param margins Margin type.
#' @param copula Copula family.
#' @return A named list with VAR params plus `z_rho` or `theta`.
#' @noRd
.init_constant_params <- function(D, margins = "normal", copula = "gaussian") {
  dep_init <- if (identical(copula, "clayton")) {
    list(theta = runif(1, 0.2, 2.5))
  } else {
    list(z_rho = rnorm(1, 0, 0.3))
  }

  c(.init_var_params(D, margins), dep_init)
}


# ============================================================================
# Shared Coefficient Extraction Helper
# ============================================================================

#' Generate default multilevel initialization values
#'
#' @param D Number of variables.
#' @param N Number of units.
#' @param margins Margin type.
#' @return A named list with multilevel params.
#' @noRd
.init_multilevel_params <- function(D, N, margins = "normal") {
  base <- list(
    phi_bar = rnorm(4, 0, 0.1),
    tau_phi = runif(4, 0.05, 0.15),
    z_phi = matrix(rnorm(N * 4, 0, 0.5), N, 4),
    rho = runif(1, -0.3, 0.3)
  )
  if (.uses_sem_multilevel_mixed_engine("multilevel", margins)) {
    # Generic mixed multilevel model: union of per-dimension margin parameters.
    return(c(base, list(
      sigma_eps = runif(D, 0.8, 1.2),
      eta = rnorm(D, 0, 0.3),
      omega = runif(D, 0.5, 1.5),
      delta = runif(D, -0.3, 0.3),
      shape_gam = runif(D, 0.5, 2.0)
    )))
  }
  if (identical(margins[[1L]], "exponential")) {
    c(base, list(eta = rep(log(0.2), D)))
  } else {
    c(base, list(sigma = runif(D, 0.8, 1.2)))
  }
}


#' Generate default TV multilevel initialization values
#' @noRd
.init_multilevel_tv_params <- function(D, N, T_obs, margins = "normal",
                                       tv_phi = FALSE, tv_sigma = FALSE) {
  n_eff <- T_obs - 1L
  phi_mask <- .resolve_phi_tv_mask(tv_phi)
  n_phi_tv <- sum(phi_mask)

  init <- list(
    phi_bar = rnorm(4, 0, 0.1),
    tau_phi = runif(4, 0.05, 0.15),
    z_phi = matrix(rnorm(N * 4, 0, 0.5), N, 4),
    rho = runif(1, -0.3, 0.3),
    # Generic TV engine union of margin parameters.
    sigma_eps = runif(D, 0.8, 1.2),
    eta = rnorm(D, 0, 0.3),
    omega = runif(D, 0.5, 1.5),
    delta = runif(D, -0.3, 0.3),
    shape_gam = runif(D, 0.5, 2.0)
  )
  if (n_phi_tv > 0L) {
    init$tau_phi_tv <- runif(n_phi_tv, 0.02, 0.08)
    init$phi_raw <- matrix(rnorm(n_eff * n_phi_tv, 0, 0.1), n_eff, n_phi_tv)
  }
  if (isTRUE(tv_sigma)) {
    init$tau_sigma <- runif(D, 0.02, 0.08)
    init$sigma_raw <- matrix(rnorm(n_eff * D, 0, 0.1), n_eff, D)
  }
  init
}


#' Generate default SEM initialization values
#'
#' @param T_obs Number of time points.
#' @param margins Margin type. Homogeneous skew-normal/gamma and per-variable
#'   specs use the mixed engine, in which case the union of per-family parameters
#'   is initialised.
#' @return A named list with SEM params.
#' @noRd
.init_sem_params <- function(T_obs, margins = "normal") {
  init <- list(
    mu = rnorm(2, 0, 0.1),
    phi11 = runif(1, 0.1, 0.5),
    phi12 = runif(1, -0.2, 0.2),
    phi21 = runif(1, -0.2, 0.2),
    phi22 = runif(1, 0.1, 0.5),
    rho_raw = rnorm(1, 0, 0.3),
    zeta = matrix(rnorm(T_obs * 2, 0, 0.5), T_obs, 2)
  )

  if (.uses_sem_multilevel_mixed_engine("sem", margins)) {
    init$sigma_eps <- runif(2, 0.8, 1.2)
    init$eta <- rnorm(2, 0, 0.3)
    init$omega <- runif(2, 0.5, 1.5)
    init$delta <- runif(2, -0.3, 0.3)
    init$shape_gam <- runif(2, 0.5, 2.0)
  } else if (identical(margins[[1L]], "exponential")) {
    init$eta <- rnorm(2, 0, 0.2)
  } else {
    init$sigma <- runif(2, 0.5, 1.5)
  }

  init
}


#' Generate default TV SEM initialization values
#' @noRd
.init_sem_tv_params <- function(T_obs, margins = "normal",
                                tv_phi = FALSE, tv_sigma = FALSE) {
  n_eff <- T_obs - 1L
  phi_mask <- .resolve_phi_tv_mask(tv_phi)
  n_phi_tv <- sum(phi_mask)

  init <- list(
    mu = rnorm(2, 0, 0.1),
    Phi = diag(c(runif(1, 0.1, 0.5), runif(1, 0.1, 0.5))) +
      matrix(c(0, runif(1, -0.2, 0.2), runif(1, -0.2, 0.2), 0), 2, 2),
    z_rho_init = rnorm(1, 0, 0.3),
    sigma_omega = runif(1, 0.05, 0.15),
    omega_raw = array(rnorm(n_eff, 0, 0.1), dim = n_eff),
    zeta = matrix(rnorm(T_obs * 2, 0, 0.5), T_obs, 2),
    # Generic TV engine union of margin parameters.
    sigma_eps = runif(2, 0.8, 1.2),
    eta = rnorm(2, 0, 0.3),
    omega = runif(2, 0.5, 1.5),
    delta = runif(2, -0.3, 0.3),
    shape_gam = runif(2, 0.5, 2.0)
  )
  if (n_phi_tv > 0L) {
    init$tau_phi <- runif(n_phi_tv, 0.02, 0.08)
    init$phi_raw <- matrix(rnorm(n_eff * n_phi_tv, 0, 0.1), n_eff, n_phi_tv)
  }
  if (isTRUE(tv_sigma)) {
    init$tau_sigma <- runif(2, 0.02, 0.08)
    init$sigma_raw <- matrix(rnorm(n_eff * 2, 0, 0.1), n_eff, 2)
  }
  init
}

#' Generate default naive SEM initialization values
#'
#' @param y T x 2 matrix of row-mean factor scores.
#' @param margins Margin type. Homogeneous skew-normal/gamma and per-variable
#'   specs use the mixed engine, in which case the union of per-family parameters
#'   is initialised.
#' @return A named list with observed-score VAR params.
#' @noRd
.init_sem_naive_params <- function(y, margins = "normal") {
  y <- as.matrix(y)
  ylag <- rbind(c(0, 0), y[-nrow(y), , drop = FALSE])

  fit1 <- try(stats::lm(y[, 1] ~ ylag[, 1] + ylag[, 2]), silent = TRUE)
  fit2 <- try(stats::lm(y[, 2] ~ ylag[, 1] + ylag[, 2]), silent = TRUE)

  coef_or <- function(fit, idx, fallback) {
    if (inherits(fit, "try-error")) {
      return(fallback)
    }
    val <- stats::coef(fit)[idx]
    if (!is.finite(val)) fallback else val
  }
  clip_phi <- function(x) max(min(x, 0.9), -0.9)

  init <- list(
    mu = c(0, 0),
    phi11 = clip_phi(coef_or(fit1, 2, 0.55)),
    phi12 = clip_phi(coef_or(fit1, 3, 0.10)),
    phi21 = clip_phi(coef_or(fit2, 2, 0.10)),
    phi22 = clip_phi(coef_or(fit2, 3, 0.25)),
    rho_raw = rnorm(1, 0, 0.5)
  )

  if (.uses_sem_multilevel_mixed_engine("sem_naive", margins)) {
    init$sigma_eps <- runif(2, 0.8, 1.2)
    init$eta <- rnorm(2, 0, 0.3)
    init$omega <- runif(2, 0.5, 1.5)
    init$delta <- runif(2, -0.3, 0.3)
    init$shape_gam <- runif(2, 0.5, 2.0)
  } else if (identical(margins[[1L]], "exponential")) {
    init$eta <- rep(log(0.2), 2)
  } else {
    init$sigma <- runif(2, 0.8, 1.2)
  }

  init
}


#' Extract margin-specific coefficients from a summary
#'
#' @param summ A summary data frame from `object$fit$summary()`.
#' @param margins Character: margin type.
#' @return A named list of margin-specific coefficient vectors.
#' @noRd
.extract_margin_coefs <- function(summ, margins) {
  if (length(margins) > 1L) {
    # Report each dimension's own family scale/shape (restricted to that
    # family's dimensions), e.g. sigma_eps[1] for a normal dim and
    # sigma_exp[2] for an exponential dim. This also covers homogeneous
    # SEM/multilevel gamma/skew-normal fits routed through mixed engines.
    specs <- .mixed_margin_report_vars(margins)
    return(lapply(specs, function(vars) {
      rows <- match(vars, summ$variable)
      rows <- rows[!is.na(rows)]
      if (length(rows) == 0L) return(NULL)
      setNames(summ$mean[rows], summ$variable[rows])
    }))
  }
  switch(margins[[1L]],
    normal = list(sigma_eps = .extract_coef(summ, "^sigma_eps\\[")),
    exponential = list(sigma_exp = .extract_coef(summ, "^sigma_exp\\[")),
    skew_normal = list(
      omega = .extract_coef(summ, "^omega\\["),
      delta = .extract_coef(summ, "^delta\\[")
    ),
    gamma = list(
      sigma_gam = .extract_coef(summ, "^sigma_gam\\["),
      shape_gam = .extract_coef(summ, "^shape_gam$")
    ),
    list(sigma_eps = .extract_coef(summ, "^sigma_eps\\["))
  )
}


#' Extract posterior means from a fit summary by regex pattern
#'
#' Shared helper used by all `coef.*` methods to avoid duplicated grep logic.
#'
#' @param summ A summary data frame from `object$fit$summary()`.
#' @param pattern Regex pattern to match against `summ$variable`.
#' @return A named numeric vector of posterior means, or `NULL` if no match.
#' @noRd
.extract_coef <- function(summ, pattern) {
  rows <- grep(pattern, summ$variable)
  if (length(rows) == 0) return(NULL)
  setNames(summ$mean[rows], summ$variable[rows])
}


# ============================================================================
# Shared Print Helpers
# ============================================================================

#' Print margin-specific scale parameters from var_params
#'
#' Encapsulates the if/else logic for printing the appropriate scale
#' parameters depending on which margin was used (normal, exponential,
#' skew_normal, gamma).
#'
#' @param var_params A named list (the `var_params` element of a summary object).
#' @return Called for its side effect (printing); returns `NULL` invisibly.
#' @noRd
.print_margin_params <- function(var_params) {
  cols <- c("variable", "mean", "q2.5", "q97.5")
  # Independent checks (not else-if) so mixed fits print every family's scale
  # parameters; single-family fits still print exactly one group.
  for (nm in c("sigma_eps", "sigma_exp", "omega", "delta", "sigma_gam", "shape_gam")) {
    if (!is.null(var_params[[nm]])) {
      cat(sprintf("\n  %s:\n", nm))
      print(var_params[[nm]][, cols], row.names = FALSE)
    }
  }
  invisible(NULL)
}


#' Print shared header for all fit objects
#' @noRd
.print_fit_header <- function(x, model_label) {
  cat(model_label, "\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("Variables: %s\n", paste(x$vars, collapse = ", ")))
}


#' Print shared MCMC info and diagnostics footer
#' @noRd
.print_fit_footer <- function(x) {
  diag <- dcvar_diagnostics(x)
  cat(sprintf("Chains: %d | Warmup: %d | Sampling: %d\n",
              x$meta$chains, x$meta$iter_warmup, x$meta$iter_sampling))
  cat(sprintf("Divergences: %d | Max Rhat: %.3f\n",
              diag$n_divergent, diag$max_rhat))
}
