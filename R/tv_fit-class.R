# ============================================================================
# S3 Class: dcvar_tv_fit (time-varying DC-VAR; subclass of dcvar_fit)
# ============================================================================

#' Construct a dcvar_tv_fit object
#'
#' Subclass of `dcvar_fit`: rho-based methods (`rho_trajectory()`,
#' `plot_rho()`, `loo()`, `dependence_summary()`) inherit unchanged, while
#' Phi/sigma-shaped methods dispatch to the trajectory-aware overrides.
#' @noRd
new_dcvar_tv_fit <- function(fit, stan_data, vars, standardized,
                             margins = "normal", skew_direction = NULL,
                             tv_phi = FALSE, phi_tv_mask = NULL, tv_sigma = FALSE,
                             backend = "rstan", priors, meta) {
  if (is.null(phi_tv_mask)) {
    phi_tv_mask <- .resolve_phi_tv_mask(tv_phi)
  }
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "dcvar_tv",
      vars = vars,
      standardized = standardized,
      margins = margins,
      skew_direction = skew_direction,
      tv_phi = tv_phi,
      phi_tv_mask = phi_tv_mask,
      tv_sigma = tv_sigma,
      backend = backend,
      priors = priors,
      meta = meta
    ),
    class = c("dcvar_tv_fit", "dcvar_fit", "dcvar_model_fit")
  )
}

#' Internal: active (time-varying) coefficient names of a TV fit
#' @noRd
.tv_active_phi_coefs <- function(object) {
  mask <- object$phi_tv_mask %||% .resolve_phi_tv_mask(isTRUE(object$tv_phi))
  names(mask)[mask == 1L]
}


#' S3 methods for dcvar_tv_fit objects
#'
#' Methods for fits from [dcvar()] with `tv_phi = TRUE` and/or
#' `tv_sigma = TRUE` (time-varying VAR coefficients and residual scales).
#'
#' @param x,object A `dcvar_tv_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcvar_tv_fit-methods
NULL

#' Internal: human-readable list of time-varying components
#' @noRd
.tv_components_label <- function(x) {
  phi_label <- if (isTRUE(x$tv_phi)) {
    active <- .tv_active_phi_coefs(x)
    if (length(active) == 4L) "Phi(t)" else paste0("Phi(t): ", paste(active, collapse = ", "))
  }
  paste(
    c("rho(t)", phi_label, if (isTRUE(x$tv_sigma)) "sigma(t)"),
    collapse = ", "
  )
}

#' @describeIn dcvar_tv_fit-methods Print a concise overview of the TV DC-VAR
#'   fit.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_tv_fit <- function(x, ...) {
  .print_fit_header(x, "TV DC-VAR Model Fit")
  cat(sprintf("n_time = %d, D = %d\n", x$stan_data$n_time, x$stan_data$D))
  cat(sprintf("time-varying components: %s\n", .tv_components_label(x)))
  .print_fit_footer(x)

  rho_df <- rho_trajectory(x)
  cat(sprintf("rho range: [%.3f, %.3f] (posterior mean)\n",
              min(rho_df$mean), max(rho_df$mean)))
  if (isTRUE(x$tv_phi)) {
    phi_df <- phi_trajectory(x, probs = 0.5)
    rng <- range(phi_df$mean)
    cat(sprintf("Phi(t) coefficient range: [%.3f, %.3f] (posterior mean)\n",
                rng[1], rng[2]))
  }
  invisible(x)
}


#' @describeIn dcvar_tv_fit-methods Produce a detailed summary including the
#'   rho, Phi, and sigma trajectories, constant parameters, and diagnostics.
#' @param probs Numeric vector of quantile probabilities for the trajectories
#'   (default: `c(0.025, 0.5, 0.975)`).
#' @return A `dcvar_tv_summary` object (a list).
#' @export
summary.dcvar_tv_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  rho_df <- rho_trajectory(object, probs = probs)
  vp <- var_params(object)
  diag <- dcvar_diagnostics(object)

  trajectory_block <- function(df, group_col) {
    do.call(rbind, lapply(split(df, df[[group_col]]), function(g) {
      data.frame(
        group = g[[group_col]][1],
        start = g$mean[1],
        end = g$mean[nrow(g)],
        delta = g$mean[nrow(g)] - g$mean[1],
        min = min(g$mean),
        max = max(g$mean)
      )
    }))
  }

  phi_summary <- if (isTRUE(object$tv_phi)) {
    trajectory_block(phi_trajectory(object, probs = probs), "coefficient")
  } else {
    NULL
  }
  sigma_summary <- if (isTRUE(object$tv_sigma)) {
    trajectory_block(sigma_trajectory(object, probs = probs), "variable")
  } else {
    NULL
  }

  out <- list(
    model = "dcvar_tv",
    n_time = object$stan_data$n_time,
    D = object$stan_data$D,
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
    phi_summary = phi_summary,
    sigma_summary = sigma_summary,
    var_params = vp,
    diagnostics = diag
  )
  class(out) <- "dcvar_tv_summary"
  out
}


#' Print a dcvar_tv_summary object
#'
#' @param x A `dcvar_tv_summary` object as returned by
#'   [summary.dcvar_tv_fit()].
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_tv_summary <- function(x, ...) {
  cat("TV DC-VAR Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("n_time = %d, D = %d\n", x$n_time, x$D))
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

  cat("\nConstant Parameters:\n")
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
    cat("\n  sigma_omega:\n")
    print(x$var_params$sigma_omega[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }

  cat("\nDiagnostics:\n")
  cat(sprintf("  Divergences: %d\n", x$diagnostics$n_divergent))
  cat(sprintf("  Max Rhat: %.3f\n", x$diagnostics$max_rhat))
  cat(sprintf("  Min ESS bulk: %.0f\n", x$diagnostics$min_ess_bulk))

  invisible(x)
}


#' @describeIn dcvar_tv_fit-methods Extract posterior means of the constant
#'   model coefficients (baselines, walk innovation SDs, margin shape
#'   parameters). Use [phi_trajectory()] and [sigma_trajectory()] for the
#'   time-varying paths.
#' @return A named list of posterior means.
#' @export
coef.dcvar_tv_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  result <- list(
    mu = .extract_required_coef(summ, "^mu\\[", "mu", "coef.dcvar_tv_fit()"),
    Phi = .extract_required_coef(summ, "^Phi\\[", "Phi", "coef.dcvar_tv_fit()")
  )

  # The generic TV model reports per-dimension margin parameters (union
  # block), so use the per-dimension report regardless of homogeneity.
  margins_vec <- rep(object$margins %||% "normal", length.out = object$stan_data$D)
  specs <- .mixed_margin_report_vars(margins_vec)
  margin_coefs <- lapply(specs, function(vars) {
    rows <- match(vars, summ$variable)
    rows <- rows[!is.na(rows)]
    if (length(rows) == 0L) return(NULL)
    setNames(summ$mean[rows], summ$variable[rows])
  })
  result <- c(result, Filter(Negate(is.null), margin_coefs))

  if (isTRUE(object$tv_phi)) {
    tau_phi <- .extract_required_coef(summ, "^tau_phi\\[", "tau_phi", "coef.dcvar_tv_fit()")
    # Stan indexes tau_phi by active position; relabel by coefficient name.
    names(tau_phi) <- paste0("tau_phi[", .tv_active_phi_coefs(object), "]")
    result$tau_phi <- tau_phi
  }
  if (isTRUE(object$tv_sigma)) {
    result$tau_sigma <- .extract_required_coef(summ, "^tau_sigma\\[", "tau_sigma", "coef.dcvar_tv_fit()")
  }
  result$sigma_omega <- .extract_required_coef(summ, "^sigma_omega$", "sigma_omega", "coef.dcvar_tv_fit()")
  result
}


#' @describeIn dcvar_tv_fit-methods Dispatch to a plot type: `"rho"`,
#'   `"phi"` (coefficient trajectories), `"sigma"` (scale trajectories),
#'   `"diagnostics"`, or `"pit"`.
#' @param type Character; one of `"rho"`, `"phi"`, `"sigma"`,
#'   `"diagnostics"`, or `"pit"`.
#' @return A ggplot object.
#' @export
plot.dcvar_tv_fit <- function(x, type = c("rho", "phi", "sigma", "diagnostics", "pit"), ...) {
  type <- match.arg(type)
  switch(type,
    rho = plot_rho(x, ...),
    phi = plot_phi_trajectory(x, ...),
    sigma = plot_sigma_trajectory(x, ...),
    diagnostics = plot_diagnostics(x, ...),
    pit = plot_pit(x, ...)
  )
}
