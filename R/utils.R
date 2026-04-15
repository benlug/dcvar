# ============================================================================
# Internal Utility Functions
# ============================================================================

#' Safe correlation that returns NA for constant vectors
#'
#' @param x,y Numeric vectors.
#' @return Correlation coefficient, or `NA_real_` if either vector has zero
#'   variance.
#' @noRd
.safe_cor <- function(x, y) {
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

  switch(margins,
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


#' Generate default constant copula initialization values
#'
#' @param D Number of variables.
#' @param margins Margin type.
#' @return A named list with VAR params plus `z_rho`.
#' @noRd
.init_constant_params <- function(D, margins = "normal") {
  c(
    .init_var_params(D, margins),
    list(z_rho = rnorm(1, 0, 0.3))
  )
}


# ============================================================================
# Shared Coefficient Extraction Helper
# ============================================================================

#' Generate default multilevel initialization values
#'
#' @param D Number of variables.
#' @param N Number of units.
#' @return A named list with multilevel params.
#' @noRd
.init_multilevel_params <- function(D, N) {
  list(
    phi_bar = rnorm(4, 0, 0.1),
    tau_phi = runif(4, 0.05, 0.15),
    z_phi = matrix(rnorm(N * 4, 0, 0.5), N, 4),
    sigma = runif(D, 0.8, 1.2),
    rho = runif(1, -0.3, 0.3)
  )
}


#' Generate default SEM initialization values
#'
#' @param T_obs Number of time points.
#' @return A named list with SEM params.
#' @noRd
.init_sem_params <- function(T_obs) {
  list(
    mu = rnorm(2, 0, 0.1),
    phi11 = runif(1, 0.1, 0.5),
    phi12 = runif(1, -0.2, 0.2),
    phi21 = runif(1, -0.2, 0.2),
    phi22 = runif(1, 0.1, 0.5),
    sigma = runif(2, 0.5, 1.5),
    rho_raw = rnorm(1, 0, 0.3),
    zeta = matrix(rnorm(T_obs * 2, 0, 0.5), T_obs, 2)
  )
}


#' Extract margin-specific coefficients from a summary
#'
#' @param summ A summary data frame from `object$fit$summary()`.
#' @param margins Character: margin type.
#' @return A named list of margin-specific coefficient vectors.
#' @noRd
.extract_margin_coefs <- function(summ, margins) {
  switch(margins,
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
  if (!is.null(var_params$sigma_eps)) {
    cat("\n  sigma_eps:\n")
    print(var_params$sigma_eps[, cols], row.names = FALSE)
  } else if (!is.null(var_params$sigma_exp)) {
    cat("\n  sigma_exp:\n")
    print(var_params$sigma_exp[, cols], row.names = FALSE)
  } else if (!is.null(var_params$omega)) {
    cat("\n  omega:\n")
    print(var_params$omega[, cols], row.names = FALSE)
    cat("\n  delta:\n")
    print(var_params$delta[, cols], row.names = FALSE)
  } else if (!is.null(var_params$sigma_gam)) {
    cat("\n  sigma_gam:\n")
    print(var_params$sigma_gam[, cols], row.names = FALSE)
    cat("\n  shape_gam:\n")
    print(var_params$shape_gam[, cols], row.names = FALSE)
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
