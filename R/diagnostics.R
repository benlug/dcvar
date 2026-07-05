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

#' Internal: validate model dimensions used by diagnostics
#' @noRd
.diagnostic_positive_int <- function(value, name, location) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value != as.integer(value)) {
    cli_abort(
      "{.fun dcvar_diagnostics} requires {.field {name}} to be a positive integer in {.field {location}}."
    )
  }

  as.integer(value)
}

#' Internal: build the exact sampled-parameter set for diagnostics
#' @noRd
.diagnostic_parameter_variables <- function(object) {
  margins <- object$margins %||% "normal"
  model <- object$model

  if (identical(model, "multilevel_tv")) {
    N <- .diagnostic_positive_int(object$N, "N", "object")
    n_time_obs <- .diagnostic_positive_int(object$stan_data$n_time, "n_time", "stan_data")
    n_eff <- n_time_obs - 1L
    margin_vars <- .mixed_diagnostic_margin_vars(rep(margins, length.out = 2L))
    vars <- c(
      paste0("phi_bar[", seq_len(4), "]"),
      paste0("tau_phi[", seq_len(4), "]"),
      paste0(
        "z_phi[",
        rep(seq_len(N), each = 4),
        ",",
        rep(seq_len(4), times = N),
        "]"
      ),
      margin_vars,
      "rho"
    )
    if (isTRUE(object$tv_phi)) {
      n_phi_tv <- sum(object$phi_tv_mask %||% .resolve_phi_tv_mask(TRUE))
      vars <- c(
        vars,
        paste0("tau_phi_tv[", seq_len(n_phi_tv), "]"),
        paste0(
          "phi_raw[",
          rep(seq_len(n_eff), each = n_phi_tv), ",", rep(seq_len(n_phi_tv), times = n_eff),
          "]"
        )
      )
    }
    if (isTRUE(object$tv_sigma)) {
      vars <- c(
        vars,
        paste0("tau_sigma[", seq_len(2), "]"),
        paste0(
          "sigma_raw[",
          rep(seq_len(n_eff), each = 2), ",", rep(seq_len(2), times = n_eff),
          "]"
        )
      )
    }
    return(vars)
  }

  if (identical(model, "multilevel")) {
    N <- .diagnostic_positive_int(object$N, "N", "object")
    margin_vars <- if (.uses_sem_multilevel_mixed_engine("multilevel", margins)) {
      .mixed_diagnostic_margin_vars(
        .sem_multilevel_engine_margins("multilevel", margins, 2L)
      )
    } else if (identical(margins, "exponential")) {
      paste0("eta[", seq_len(2), "]")
    } else {
      paste0("sigma[", seq_len(2), "]")
    }
    return(c(
      paste0("phi_bar[", seq_len(4), "]"),
      paste0("tau_phi[", seq_len(4), "]"),
      paste0(
        "z_phi[",
        rep(seq_len(N), each = 4),
        ",",
        rep(seq_len(4), times = N),
        "]"
      ),
      margin_vars,
      "rho"
    ))
  }

  if (identical(model, "sem_tv")) {
    n_time_obs <- .diagnostic_positive_int(object$stan_data$n_time, "n_time", "stan_data")
    n_eff <- n_time_obs - 1L
    margin_vars <- .mixed_diagnostic_margin_vars(rep(margins, length.out = 2L))
    vars <- c(
      "mu[1]", "mu[2]",
      "Phi[1,1]", "Phi[1,2]", "Phi[2,1]", "Phi[2,2]",
      margin_vars,
      "z_rho_init",
      "sigma_omega",
      paste0("omega_raw[", seq_len(n_eff), "]"),
      paste0(
        "zeta[",
        rep(seq_len(n_time_obs), each = 2),
        ",",
        rep(seq_len(2), times = n_time_obs),
        "]"
      )
    )
    if (isTRUE(object$tv_phi)) {
      n_phi_tv <- sum(object$phi_tv_mask %||% .resolve_phi_tv_mask(TRUE))
      vars <- c(
        vars,
        paste0("tau_phi[", seq_len(n_phi_tv), "]"),
        paste0(
          "phi_raw[",
          rep(seq_len(n_eff), each = n_phi_tv), ",", rep(seq_len(n_phi_tv), times = n_eff),
          "]"
        )
      )
    }
    if (isTRUE(object$tv_sigma)) {
      vars <- c(
        vars,
        paste0("tau_sigma[", seq_len(2), "]"),
        paste0(
          "sigma_raw[",
          rep(seq_len(n_eff), each = 2), ",", rep(seq_len(2), times = n_eff),
          "]"
        )
      )
    }
    return(vars)
  }

  if (identical(model, "sem")) {
    n_time_obs <- .diagnostic_positive_int(object$stan_data$n_time, "n_time", "stan_data")
    method <- object$method %||% "indicator"
    margin_vars <- if (.uses_sem_multilevel_mixed_engine("sem", margins)) {
      .mixed_diagnostic_margin_vars(
        .sem_multilevel_engine_margins("sem", margins, 2L)
      )
    } else if (identical(margins, "exponential")) {
      paste0("eta[", seq_len(2), "]")
    } else {
      paste0("sigma[", seq_len(2), "]")
    }
    vars <- c(
      "mu[1]", "mu[2]",
      "phi11", "phi12", "phi21", "phi22",
      margin_vars,
      "rho_raw"
    )
    if (identical(method, "naive")) {
      return(vars)
    }
    return(c(
      vars,
      paste0(
        "zeta[",
        rep(seq_len(n_time_obs), each = 2),
        ",",
        rep(seq_len(2), times = n_time_obs),
        "]"
      )
    ))
  }

  D <- .diagnostic_positive_int(object$stan_data$D, "D", "stan_data")
  mu_vars <- paste0("mu[", seq_len(D), "]")
  phi_vars <- paste0(
    "Phi[",
    rep(seq_len(D), each = D),
    ",",
    rep(seq_len(D), times = D),
    "]"
  )
  margin_vars <- if (.is_mixed_margins(margins)) {
    # Generic mixed model: check the sampled parameter each dimension actually
    # uses for its own family (the union's unused parameters merely sample from
    # their priors and need not be monitored).
    .mixed_diagnostic_margin_vars(margins)
  } else {
    switch(margins[[1L]],
      normal = paste0("sigma_eps[", seq_len(D), "]"),
      exponential = paste0("eta[", seq_len(D), "]"),
      skew_normal = c(
        paste0("omega[", seq_len(D), "]"),
        paste0("delta[", seq_len(D), "]")
      ),
      gamma = c(
        paste0("eta[", seq_len(D), "]"),
        "shape_gam"
      ),
      paste0("sigma_eps[", seq_len(D), "]")
    )
  }

  if (identical(model, "constant")) {
    copula <- object$copula %||% "gaussian"
    dependence_var <- if (identical(copula, "clayton")) "theta" else "z_rho"
    return(c(mu_vars, phi_vars, margin_vars, dependence_var))
  }

  if (identical(model, "dcvar_tv")) {
    n_time_obs <- .diagnostic_positive_int(object$stan_data$n_time, "n_time", "stan_data")
    n_eff <- n_time_obs - 1L
    # The generic TV model always uses the mixed-union parameter block, so the
    # per-dimension sampled margin parameters follow the mixed convention even
    # for homogeneous margin requests.
    tv_margin_vars <- .mixed_diagnostic_margin_vars(rep(margins, length.out = D))
    vars <- c(
      mu_vars,
      phi_vars,
      tv_margin_vars,
      "z_rho_init",
      "sigma_omega",
      paste0("omega_raw[", seq_len(n_eff), "]")
    )
    if (isTRUE(object$tv_phi)) {
      # tau_phi / phi_raw are sized to the number of active coefficients.
      n_phi_tv <- sum(object$phi_tv_mask %||% .resolve_phi_tv_mask(TRUE))
      vars <- c(
        vars,
        paste0("tau_phi[", seq_len(n_phi_tv), "]"),
        paste0(
          "phi_raw[",
          rep(seq_len(n_eff), each = n_phi_tv), ",", rep(seq_len(n_phi_tv), times = n_eff),
          "]"
        )
      )
    }
    if (isTRUE(object$tv_sigma)) {
      vars <- c(
        vars,
        paste0("tau_sigma[", seq_len(D), "]"),
        paste0(
          "sigma_raw[",
          rep(seq_len(n_eff), each = D), ",", rep(seq_len(D), times = n_eff),
          "]"
        )
      )
    }
    return(vars)
  }

  if (identical(model, "dcvar")) {
    n_time_obs <- .diagnostic_positive_int(object$stan_data$n_time, "n_time", "stan_data")
    return(c(
      mu_vars,
      phi_vars,
      margin_vars,
      "z_rho_init",
      "sigma_omega",
      paste0("omega_raw[", seq_len(n_time_obs - 1L), "]")
    ))
  }

  if (identical(model, "dcvar_covariate") || identical(model, "dcvar_covariate_nodrift")) {
    P <- .diagnostic_positive_int(object$stan_data$P, "P", "stan_data")
    intercept_var <- if (identical(model, "dcvar_covariate") &&
                         isTRUE(object$dynamic_engine)) {
      "z_rho_init"
    } else {
      "beta_0"
    }
    vars <- c(
      mu_vars,
      phi_vars,
      margin_vars,
      intercept_var,
      paste0("beta[", seq_len(P), "]")
    )
    if (identical(model, "dcvar_covariate")) {
      n_time_obs <- .diagnostic_positive_int(object$stan_data$n_time, "n_time", "stan_data")
      n_omega <- (n_time_obs - 1L) - as.integer(isTRUE(object$zero_init_eta))
      if (n_omega > 0L) {
        vars <- c(vars, "sigma_omega", paste0("omega_raw[", seq_len(n_omega), "]"))
      } else {
        vars <- c(vars, "sigma_omega")
      }
    }
    return(vars)
  }

  if (identical(model, "hmm")) {
    K <- .diagnostic_positive_int(object$K, "K", "object")
    hmm_tail <- c(
      paste0("z_rho[", seq_len(K), "]"),
      paste0("pi0[", seq_len(K), "]"),
      paste0(
        "A[",
        rep(seq_len(K), each = K),
        ",",
        rep(seq_len(K), times = K),
        "]"
      )
    )

    if (isTRUE(object$switching)) {
      # The Markov-switching engine samples state-indexed mu / Phi_base+Phi_dev /
      # margin-union, so the monitored names follow the array layout.
      sw <- object$switch
      Mu_K <- if (isTRUE(sw$mu == 1L)) K else 1L
      Mrg_K <- if (isTRUE(sw$margins == 1L)) K else 1L
      mu_sw <- paste0(
        "mu[",
        rep(seq_len(Mu_K), each = D), ",", rep(seq_len(D), times = Mu_K),
        "]"
      )
      phi_sw <- paste0(
        "Phi_base[",
        rep(seq_len(D), each = D), ",", rep(seq_len(D), times = D),
        "]"
      )
      if (any(sw$phi_mask > 0L)) {
        ij <- list(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(2L, 2L))
        active <- which(sw$phi_mask > 0L)
        dev <- unlist(lapply(seq_len(K), function(k) {
          vapply(active, function(a) {
            sprintf("Phi_dev[%d,%d,%d]", k, ij[[a]][1], ij[[a]][2])
          }, character(1))
        }))
        phi_sw <- c(phi_sw, dev)
      }
      margin_sw <- .hmm_switching_diagnostic_margin_vars(object$margins_matrix, Mrg_K)
      return(c(mu_sw, phi_sw, margin_sw, hmm_tail))
    }

    return(c(
      mu_vars,
      phi_vars,
      margin_vars,
      hmm_tail
    ))
  }

  c(mu_vars, phi_vars, margin_vars)
}

#' Internal: extract common diagnostics from a CmdStanMCMC fit
#' @noRd
.sampling_diagnostics_from_fit <- function(fit, backend = "rstan", object = NULL) {
  vars <- if (is.null(object)) NULL else .diagnostic_parameter_variables(object)
  diag_summ <- .fit_diagnostic_summary(fit, backend)
  summ <- suppressWarnings(.fit_summary(fit, variables = vars, backend = backend))
  sampler_diags <- .fit_sampler_diagnostics(fit, backend)
  diag_names <- dimnames(sampler_diags)[[3L]]
  accept_stat <- if ("accept_stat__" %in% diag_names) {
    as.numeric(sampler_diags[, , "accept_stat__", drop = TRUE])
  } else {
    numeric()
  }
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
                                     rhat_threshold = 1.10, backend = "rstan",
                                     object = NULL) {
  diag <- .sampling_diagnostics_from_fit(fit, backend = backend, object = object)
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
