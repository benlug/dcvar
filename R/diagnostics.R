# ============================================================================
# MCMC Diagnostics
# ============================================================================

#' Extract MCMC diagnostics
#'
#' Returns a summary of sampling diagnostics including divergences,
#' tree depth warnings, Rhat, and effective sample size.
#' The convergence headline is computed from sampled parameters only and
#' excludes generated quantities and deterministic transformed outputs.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments (unused).
#'
#' @return A named list with:
#'   - `n_divergent`: total number of divergent transitions
#'   - `n_max_treedepth`: transitions hitting max tree depth
#'   - `max_rhat`: worst (highest) Rhat across sampled parameters
#'   - `min_ess_bulk`: smallest bulk ESS among sampled parameters
#'   - `min_ess_tail`: smallest tail ESS among sampled parameters
#'   - `mean_accept_prob`: mean acceptance probability
#' @export
dcvar_diagnostics <- function(object, ...) {
  UseMethod("dcvar_diagnostics")
}

#' @rdname dcvar_diagnostics
#' @export
dcvar_diagnostics.default <- function(object, ...) {
  cli_abort("{.fun dcvar_diagnostics} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' Internal: build the exact sampled-parameter set for diagnostics
#' @noRd
.diagnostic_parameter_variables <- function(object) {
  margins <- object$margins %||% "normal"
  model <- object$model

  mu_vars <- paste0("mu[", seq_len(object$stan_data$D), "]")
  phi_vars <- paste0(
    "Phi[",
    rep(seq_len(object$stan_data$D), each = object$stan_data$D),
    ",",
    rep(seq_len(object$stan_data$D), times = object$stan_data$D),
    "]"
  )
  margin_vars <- switch(margins,
    normal = paste0("sigma_eps[", seq_len(object$stan_data$D), "]"),
    exponential = paste0("eta[", seq_len(object$stan_data$D), "]"),
    skew_normal = c(
      paste0("omega[", seq_len(object$stan_data$D), "]"),
      paste0("delta[", seq_len(object$stan_data$D), "]")
    ),
    gamma = c(
      paste0("eta[", seq_len(object$stan_data$D), "]"),
      "shape_gam"
    ),
    paste0("sigma_eps[", seq_len(object$stan_data$D), "]")
  )

  switch(model,
    constant = c(mu_vars, phi_vars, margin_vars, "z_rho"),
    dcvar = c(
      mu_vars,
      phi_vars,
      margin_vars,
      "z_rho_init",
      "sigma_omega",
      paste0("omega_raw[", seq_len(object$stan_data$T - 1), "]")
    ),
    hmm = c(
      mu_vars,
      phi_vars,
      margin_vars,
      paste0("z_rho[", seq_len(object$K), "]"),
      paste0("pi0[", seq_len(object$K), "]"),
      paste0(
        "A[",
        rep(seq_len(object$K), each = object$K),
        ",",
        rep(seq_len(object$K), times = object$K),
        "]"
      )
    ),
    multilevel = c(
      paste0("phi_bar[", seq_len(4), "]"),
      paste0("tau_phi[", seq_len(4), "]"),
      paste0(
        "z_phi[",
        rep(seq_len(object$N), each = 4),
        ",",
        rep(seq_len(4), times = object$N),
        "]"
      ),
      paste0("sigma[", seq_len(2), "]"),
      "rho"
    ),
    sem = c(
      "mu[1]", "mu[2]",
      "phi11", "phi12", "phi21", "phi22",
      paste0("sigma[", seq_len(2), "]"),
      "rho_raw",
      paste0(
        "zeta[",
        rep(seq_len(object$stan_data$T), each = 2),
        ",",
        rep(seq_len(2), times = object$stan_data$T),
        "]"
      )
    ),
    c(mu_vars, phi_vars, margin_vars)
  )
}

#' Internal: extract common diagnostics from a CmdStanMCMC fit
#' @noRd
.sampling_diagnostics_from_fit <- function(fit, backend = "rstan", object = NULL) {
  vars <- if (is.null(object)) NULL else .diagnostic_parameter_variables(object)
  diag_summ <- .fit_diagnostic_summary(fit, backend)
  summ <- suppressWarnings(.fit_summary(fit, variables = vars, backend = backend))
  sampler_diags <- .fit_sampler_diagnostics(fit, backend)
  accept_stat <- as.numeric(sampler_diags[, , "accept_stat__", drop = TRUE])
  accept_stat <- accept_stat[is.finite(accept_stat)]

  rhat <- summ$rhat[is.finite(summ$rhat)]
  ess_bulk <- summ$ess_bulk[is.finite(summ$ess_bulk)]
  ess_tail <- summ$ess_tail[is.finite(summ$ess_tail)]

  list(
    n_divergent = sum(diag_summ$num_divergent),
    n_max_treedepth = sum(diag_summ$num_max_treedepth),
    max_rhat = if (length(rhat) > 0) max(rhat) else NA_real_,
    min_ess_bulk = if (length(ess_bulk) > 0) min(ess_bulk) else NA_real_,
    min_ess_tail = if (length(ess_tail) > 0) min(ess_tail) else NA_real_,
    mean_accept_prob = if (length(accept_stat) > 0) mean(accept_stat) else NA_real_
  )
}

#' Internal: report post-sampling diagnostics to the user
#' @noRd
.report_sampling_outcome <- function(fit, model_label, chains = NA_integer_,
                                     rhat_threshold = 1.10, backend = "rstan") {
  diag <- .sampling_diagnostics_from_fit(fit, backend = backend)
  issues <- character()

  if (diag$n_divergent > 0) {
    issues <- c(
      issues,
      sprintf("%d divergent transition%s.", diag$n_divergent,
              if (diag$n_divergent == 1) "" else "s")
    )
  }
  if (diag$n_max_treedepth > 0) {
    issues <- c(
      issues,
      sprintf("%d transition%s hit the maximum treedepth.", diag$n_max_treedepth,
              if (diag$n_max_treedepth == 1) "" else "s")
    )
  }
  if (isTRUE(chains > 1) && is.finite(diag$max_rhat) && diag$max_rhat > rhat_threshold) {
    issues <- c(
      issues,
      sprintf("Max R-hat is %.3f.", diag$max_rhat)
    )
  }

  if (length(issues) == 0) {
    cli_alert_success("{model_label} sampling complete.")
  } else {
    cli_warn(c(
      "{model_label} sampling finished with diagnostic issues.",
      setNames(issues, rep("!", length(issues))),
      "i" = "Inspect {.fun dcvar_diagnostics} before using this fit for inference."
    ))
  }

  invisible(diag)
}

#' @rdname dcvar_diagnostics
#' @export
dcvar_diagnostics.dcvar_model_fit <- function(object, ...) {
  .sampling_diagnostics_from_fit(object$fit, backend = object$backend %||% "rstan", object = object)
}
