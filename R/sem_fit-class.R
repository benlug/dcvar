# ============================================================================
# S3 Class: dcvar_sem_fit
# ============================================================================

#' Construct a dcvar_sem_fit object
#' @noRd
new_dcvar_sem_fit <- function(fit, stan_data, vars, J, lambda, sigma_e,
                               indicators, margins = "normal",
                               method = "indicator",
                               skew_direction = NULL, backend = "rstan",
                               priors, meta) {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "sem",
      vars = vars,
      standardized = FALSE,
      margins = margins,
      method = method,
      J = J,
      lambda = lambda,
      sigma_e = sigma_e,
      indicators = indicators,
      skew_direction = skew_direction,
      backend = backend,
      priors = priors,
      meta = meta
    ),
    class = c("dcvar_sem_fit", "dcvar_model_fit")
  )
}


#' Construct a dcvar_sem_tv_fit object
#' @noRd
new_dcvar_sem_tv_fit <- function(fit, stan_data, vars, J, lambda, sigma_e,
                                 indicators, margins = "normal",
                                 method = "indicator",
                                 skew_direction = NULL,
                                 tv_phi = FALSE, phi_tv_mask = NULL,
                                 tv_sigma = FALSE,
                                 backend = "rstan", priors, meta) {
  phi_tv_mask <- phi_tv_mask %||% .resolve_phi_tv_mask(tv_phi)
  out <- new_dcvar_sem_fit(
    fit = fit,
    stan_data = stan_data,
    vars = vars,
    J = J,
    lambda = lambda,
    sigma_e = sigma_e,
    indicators = indicators,
    margins = margins,
    method = method,
    skew_direction = skew_direction,
    backend = backend,
    priors = priors,
    meta = meta
  )
  out$model <- "sem_tv"
  out$tv_phi <- isTRUE(sum(phi_tv_mask) > 0L)
  out$phi_tv_mask <- phi_tv_mask
  out$tv_sigma <- isTRUE(tv_sigma)
  class(out) <- c("dcvar_sem_tv_fit", "dcvar_sem_fit", "dcvar_model_fit")
  out
}


#' S3 methods for dcvar_sem_fit objects
#'
#' @param x,object A `dcvar_sem_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcvar_sem_fit-methods
NULL

#' @describeIn dcvar_sem_fit-methods Print a concise overview.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_sem_fit <- function(x, ...) {
  title <- if (identical(x$method %||% "indicator", "naive")) {
    "Naive SEM Copula VAR Model Fit"
  } else {
    "SEM Copula VAR Model Fit"
  }
  .print_fit_header(x, title)
  cat(sprintf("n_time = %d, J = %d indicators per latent\n", x$stan_data$n_time, x$J))
  .print_fit_footer(x)

  cat(sprintf("rho: %.3f\n", coef(x)$rho[[1]]))
  invisible(x)
}


#' @describeIn dcvar_sem_fit-methods Produce a detailed summary.
#' @return A `dcvar_sem_summary` object (a list).
#' @export
summary.dcvar_sem_fit <- function(object, ...) {
  vp <- var_params(object)
  diag <- dcvar_diagnostics(object)

  out <- list(
    model = "sem",
    method = object$method %||% "indicator",
    n_time = object$stan_data$n_time,
    J = object$J,
    lambda = object$lambda,
    sigma_e = object$sigma_e,
    var_params = vp,
    diagnostics = diag
  )
  class(out) <- "dcvar_sem_summary"
  out
}


#' @describeIn dcvar_sem_fit-methods Print a concise overview of a TV SEM fit.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_sem_tv_fit <- function(x, ...) {
  .print_fit_header(x, "TV SEM Copula VAR Model Fit")
  cat(sprintf("n_time = %d, J = %d indicators per latent\n", x$stan_data$n_time, x$J))
  cat(sprintf("time-varying components: %s\n", .tv_components_label(x)))
  .print_fit_footer(x)

  rho_df <- rho_trajectory(x)
  cat(sprintf("rho range: [%.3f, %.3f] (posterior mean)\n",
              min(rho_df$mean), max(rho_df$mean)))
  invisible(x)
}


#' @describeIn dcvar_sem_fit-methods Produce a detailed TV SEM summary.
#' @param probs Numeric vector of quantile probabilities for trajectories.
#' @return A `dcvar_sem_tv_summary` object (a list).
#' @export
summary.dcvar_sem_tv_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  vp <- var_params(object)
  diag <- dcvar_diagnostics(object)
  rho_df <- rho_trajectory(object, probs = probs)

  out <- list(
    model = "sem_tv",
    method = object$method %||% "indicator",
    n_time = object$stan_data$n_time,
    J = object$J,
    lambda = object$lambda,
    sigma_e = object$sigma_e,
    components = .tv_components_label(object),
    rho_trajectory = rho_df,
    rho_summary = data.frame(
      statistic = c("start", "end", "delta", "min", "max", "mean"),
      value = c(
        rho_df$mean[1],
        rho_df$mean[nrow(rho_df)],
        rho_df$mean[nrow(rho_df)] - rho_df$mean[1],
        min(rho_df$mean),
        max(rho_df$mean),
        mean(rho_df$mean)
      )
    ),
    phi_summary = if (isTRUE(object$tv_phi)) {
      .tv_trajectory_summary_block(phi_trajectory(object, probs = probs), "coefficient")
    } else {
      NULL
    },
    sigma_summary = if (isTRUE(object$tv_sigma)) {
      .tv_trajectory_summary_block(sigma_trajectory(object, probs = probs), "variable")
    } else {
      NULL
    },
    var_params = vp,
    diagnostics = diag
  )
  class(out) <- c("dcvar_sem_tv_summary", "dcvar_sem_summary")
  out
}


#' Print a dcvar_sem_summary object
#' @param x A `dcvar_sem_summary` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_sem_summary <- function(x, ...) {
  if (identical(x$method %||% "indicator", "naive")) {
    cat("Naive SEM Copula VAR Model Summary\n")
  } else {
    cat("SEM Copula VAR Model Summary\n")
  }
  cat(strrep("=", 50), "\n")
  cat(sprintf("n_time = %d, J = %d indicators per latent\n", x$n_time, x$J))
  if (!identical(x$method %||% "indicator", "naive")) {
    cat(sprintf("Fixed lambda: %s\n", paste(round(x$lambda, 3), collapse = ", ")))
    cat(sprintf("Fixed sigma_e: %.3f\n", x$sigma_e))
  }
  cat("\n")

  cat("Latent VAR Parameters:\n")
  if (!is.null(x$var_params$mu)) {
    cat("  mu:\n")
    print(x$var_params$mu[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  if (!is.null(x$var_params$Phi)) {
    cat("\n  Phi:\n")
    print(x$var_params$Phi[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  if (!is.null(x$var_params$sigma)) {
    cat("\n  sigma:\n")
    print(x$var_params$sigma[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  } else {
    .print_margin_params(x$var_params)
  }
  if (!is.null(x$var_params$rho)) {
    cat("\n  rho:\n")
    print(x$var_params$rho[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }

  cat("\nDiagnostics:\n")
  cat(sprintf("  Divergences: %d\n", x$diagnostics$n_divergent))
  cat(sprintf("  Max Rhat: %.3f\n", x$diagnostics$max_rhat))
  cat(sprintf("  Min ESS bulk: %.0f\n", x$diagnostics$min_ess_bulk))

  invisible(x)
}


#' Print a dcvar_sem_tv_summary object
#' @param x A `dcvar_sem_tv_summary` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_sem_tv_summary <- function(x, ...) {
  cat("TV SEM Copula VAR Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("n_time = %d, J = %d indicators per latent\n", x$n_time, x$J))
  cat(sprintf("time-varying components: %s\n\n", x$components))

  cat("Rho Trajectory:\n")
  print(x$rho_summary, row.names = FALSE)
  if (!is.null(x$phi_summary)) {
    cat("\nPhi(t) Trajectories (posterior mean):\n")
    print(x$phi_summary, row.names = FALSE)
  }
  if (!is.null(x$sigma_summary)) {
    cat("\nsigma(t) Trajectories (posterior mean):\n")
    print(x$sigma_summary, row.names = FALSE)
  }

  cat("\nLatent VAR Parameters:\n")
  cat("  mu:\n")
  print(x$var_params$mu[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  cat("\n  Phi (baseline):\n")
  print(x$var_params$Phi[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  if (!is.null(x$var_params$tau_phi)) {
    cat("\n  tau_phi (walk SDs):\n")
    print(x$var_params$tau_phi[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  if (!is.null(x$var_params$tau_sigma)) {
    cat("\n  tau_sigma (log-scale walk SDs):\n")
    print(x$var_params$tau_sigma[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  .print_margin_params(x$var_params)
  if (!is.null(x$var_params$sigma_omega)) {
    cat("\n  sigma_omega (rho random-walk SD):\n")
    print(x$var_params$sigma_omega[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }

  cat("\nDiagnostics:\n")
  cat(sprintf("  Divergences: %d\n", x$diagnostics$n_divergent))
  cat(sprintf("  Max Rhat: %.3f\n", x$diagnostics$max_rhat))
  cat(sprintf("  Min ESS bulk: %.0f\n", x$diagnostics$min_ess_bulk))

  invisible(x)
}


#' @describeIn dcvar_sem_fit-methods Extract posterior means of
#'   latent VAR coefficients.
#' @return A named list of posterior means.
#' @export
coef.dcvar_sem_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  result <- list(
    mu = .extract_required_coef(summ, "^mu\\[", "mu", "coef.dcvar_sem_fit()"),
    Phi = .extract_required_coef(summ, "^Phi\\[", "Phi", "coef.dcvar_sem_fit()"),
    rho = .extract_required_coef(summ, "^rho$", "rho", "coef.dcvar_sem_fit()")
  )
  margins <- object$margins %||% "normal"
  margin_coefs <- if (identical(margins, "normal")) {
    list(sigma = .extract_required_coef(summ, "^sigma\\[", "sigma", "coef.dcvar_sem_fit()"))
  } else if (.uses_sem_multilevel_mixed_engine("sem", margins)) {
    .extract_margin_coefs(
      summ,
      .sem_multilevel_engine_margins("sem", margins, 2L)
    )
  } else {
    .extract_margin_coefs(summ, margins)
  }
  c(result[1:2], margin_coefs, result[3])
}


#' @export
coef.dcvar_sem_tv_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  margins_vec <- rep(object$margins %||% "normal", length.out = 2L)
  result <- list(
    mu = .extract_required_coef(summ, "^mu\\[", "mu", "coef.dcvar_sem_tv_fit()"),
    Phi = .extract_required_coef(summ, "^Phi\\[", "Phi", "coef.dcvar_sem_tv_fit()")
  )
  margin_coefs <- .extract_margin_coefs(summ, margins_vec)
  extras <- list(
    sigma_omega = .extract_required_coef(summ, "^sigma_omega$", "sigma_omega", "coef.dcvar_sem_tv_fit()"),
    rho = .extract_required_coef(summ, "^rho\\[", "rho", "coef.dcvar_sem_tv_fit()")
  )
  if (isTRUE(object$tv_phi)) {
    tau_phi <- .extract_required_coef(summ, "^tau_phi\\[", "tau_phi", "coef.dcvar_sem_tv_fit()")
    names(tau_phi) <- paste0("tau_phi[", .tv_active_phi_coefs(object), "]")
    extras$tau_phi <- tau_phi
  }
  if (isTRUE(object$tv_sigma)) {
    extras$tau_sigma <- .extract_required_coef(summ, "^tau_sigma\\[", "tau_sigma", "coef.dcvar_sem_tv_fit()")
  }
  c(result, margin_coefs, extras)
}


#' @describeIn dcvar_sem_fit-methods Dispatch to a plot type.
#' @param type Character; one of `"latent_states"`, `"rho"`, `"diagnostics"`.
#' @return A ggplot object.
#' @export
plot.dcvar_sem_fit <- function(x,
                                type = c("latent_states", "rho", "diagnostics"),
                                ...) {
  type <- match.arg(type)
  if (identical(type, "latent_states") && identical(x$method %||% "indicator", "naive")) {
    cli_abort("Latent-state plots are not available for naive SEM fits because no latent measurement model is estimated.")
  }
  switch(type,
    latent_states = plot_latent_states(x, ...),
    rho = plot_rho(x, ...),
    diagnostics = plot_diagnostics(x, ...)
  )
}
