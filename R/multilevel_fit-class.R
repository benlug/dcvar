# ============================================================================
# S3 Class: dcvar_multilevel_fit
# ============================================================================

#' Construct a dcvar_multilevel_fit object
#' @noRd
new_dcvar_multilevel_fit <- function(fit, stan_data, N, vars, centered,
                                     person_means, priors, meta,
                                     standardized = FALSE,
                                     margins = "normal",
                                     skew_direction = NULL,
                                     backend = "rstan") {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = "multilevel",
      N = N,
      vars = vars,
      standardized = standardized,
      margins = margins,
      skew_direction = skew_direction,
      centered = centered,
      person_means = person_means,
      backend = backend,
      priors = priors,
      meta = meta
    ),
    class = c("dcvar_multilevel_fit", "dcvar_model_fit")
  )
}


#' Construct a dcvar_multilevel_tv_fit object
#' @noRd
new_dcvar_multilevel_tv_fit <- function(fit, stan_data, N, vars, centered,
                                        person_means, priors, meta,
                                        standardized = FALSE,
                                        margins = "normal",
                                        skew_direction = NULL,
                                        tv_phi = FALSE, phi_tv_mask = NULL,
                                        tv_sigma = FALSE,
                                        backend = "rstan") {
  phi_tv_mask <- phi_tv_mask %||% .resolve_phi_tv_mask(tv_phi)
  out <- new_dcvar_multilevel_fit(
    fit = fit,
    stan_data = stan_data,
    N = N,
    vars = vars,
    centered = centered,
    person_means = person_means,
    priors = priors,
    meta = meta,
    standardized = standardized,
    margins = margins,
    skew_direction = skew_direction,
    backend = backend
  )
  out$model <- "multilevel_tv"
  out$tv_phi <- isTRUE(sum(phi_tv_mask) > 0L)
  out$phi_tv_mask <- phi_tv_mask
  out$tv_sigma <- isTRUE(tv_sigma)
  class(out) <- c("dcvar_multilevel_tv_fit", "dcvar_multilevel_fit", "dcvar_model_fit")
  out
}


#' S3 methods for dcvar_multilevel_fit objects
#'
#' @param x,object A `dcvar_multilevel_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @name dcvar_multilevel_fit-methods
NULL

#' @describeIn dcvar_multilevel_fit-methods Print a concise overview.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_multilevel_fit <- function(x, ...) {
  .print_fit_header(x, "Multilevel Copula VAR Model Fit")
  cat(sprintf("N = %d units, n_time = %d per unit\n", x$N, x$stan_data$n_time))
  .print_fit_footer(x)

  cat(sprintf("rho (global): %.3f\n", coef(x)$rho[[1]]))
  invisible(x)
}


#' @describeIn dcvar_multilevel_fit-methods Produce a detailed summary.
#' @return A `dcvar_multilevel_summary` object (a list).
#' @export
summary.dcvar_multilevel_fit <- function(object, ...) {
  vp <- var_params(object)
  diag <- dcvar_diagnostics(object)
  re <- random_effects(object)

  out <- list(
    model = "multilevel",
    N = object$N,
    n_time = object$stan_data$n_time,
    var_params = vp,
    random_effects = re,
    diagnostics = diag
  )
  class(out) <- "dcvar_multilevel_summary"
  out
}


#' Internal: human-readable list of time-varying multilevel components
#' @noRd
.multilevel_tv_components_label <- function(x) {
  phi_label <- if (isTRUE(x$tv_phi)) {
    active <- .tv_active_phi_coefs(x)
    if (length(active) == 4L) "shared Phi(t) drift" else paste0("shared Phi(t) drift: ", paste(active, collapse = ", "))
  }
  paste(
    c(phi_label, if (isTRUE(x$tv_sigma)) "sigma(t)"),
    collapse = ", "
  )
}


#' @describeIn dcvar_multilevel_fit-methods Print a concise overview of a TV
#'   multilevel fit.
#' @return Invisibly returns `x`.
#' @export
print.dcvar_multilevel_tv_fit <- function(x, ...) {
  .print_fit_header(x, "TV Multilevel Copula VAR Model Fit")
  cat(sprintf("N = %d units, n_time = %d per unit\n", x$N, x$stan_data$n_time))
  cat(sprintf("time-varying components: %s\n", .multilevel_tv_components_label(x)))
  .print_fit_footer(x)

  cat(sprintf("rho (global): %.3f\n", coef(x)$rho[[1]]))
  invisible(x)
}


#' @describeIn dcvar_multilevel_fit-methods Produce a detailed TV multilevel
#'   summary.
#' @param probs Numeric vector of quantile probabilities for trajectories.
#' @return A `dcvar_multilevel_tv_summary` object (a list).
#' @export
summary.dcvar_multilevel_tv_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  vp <- var_params(object)
  diag <- dcvar_diagnostics(object)
  re <- random_effects(object)
  rho_df <- rho_trajectory(object, probs = probs)
  rho_quantiles <- grep("^q", names(rho_df), value = TRUE)

  out <- list(
    model = "multilevel_tv",
    N = object$N,
    n_time = object$stan_data$n_time,
    components = .multilevel_tv_components_label(object),
    rho_trajectory = rho_df,
    rho_summary = data.frame(
      statistic = c("mean", "sd", rho_quantiles),
      value = c(rho_df$mean[1], rho_df$sd[1], unname(as.numeric(rho_df[1, rho_quantiles, drop = TRUE])))
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
    random_effects = re,
    diagnostics = diag
  )
  class(out) <- c("dcvar_multilevel_tv_summary", "dcvar_multilevel_summary")
  out
}


#' Print a dcvar_multilevel_summary object
#' @param x A `dcvar_multilevel_summary` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_multilevel_summary <- function(x, ...) {
  cat("Multilevel Copula VAR Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("N = %d units, n_time = %d per unit\n\n", x$N, x$n_time))

  cat("Population-Level Parameters:\n")
  if (!is.null(x$var_params$phi_bar)) {
    cat("  phi_bar:\n")
    print(x$var_params$phi_bar[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  if (!is.null(x$var_params$tau_phi)) {
    cat("\n  tau_phi:\n")
    print(x$var_params$tau_phi[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  if (!is.null(x$var_params$sigma)) {
    cat("\n  sigma:\n")
    print(x$var_params$sigma[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  } else {
    # Exponential and per-variable (mixed) margins report scale/shape under
    # family-specific names (sigma_exp / omega / delta / sigma_gam / shape_gam).
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


#' Print a dcvar_multilevel_tv_summary object
#' @param x A `dcvar_multilevel_tv_summary` object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_multilevel_tv_summary <- function(x, ...) {
  cat("TV Multilevel Copula VAR Model Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("N = %d units, n_time = %d per unit\n", x$N, x$n_time))
  cat(sprintf("time-varying components: %s\n\n", x$components))

  cat("Rho (global):\n")
  print(x$rho_summary, row.names = FALSE)
  if (!is.null(x$phi_summary)) {
    cat("\nPopulation Phi(t) Drift (posterior mean):\n")
    print(x$phi_summary, row.names = FALSE)
  }
  if (!is.null(x$sigma_summary)) {
    cat("\nsigma(t) Trajectories (posterior mean):\n")
    print(x$sigma_summary, row.names = FALSE)
  }

  cat("\nPopulation-Level Parameters:\n")
  cat("  phi_bar:\n")
  print(x$var_params$phi_bar[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  cat("\n  tau_phi:\n")
  print(x$var_params$tau_phi[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  if (!is.null(x$var_params$tau_phi_tv)) {
    cat("\n  tau_phi_tv (shared drift walk SDs):\n")
    print(x$var_params$tau_phi_tv[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  if (!is.null(x$var_params$tau_sigma)) {
    cat("\n  tau_sigma (log-scale walk SDs):\n")
    print(x$var_params$tau_sigma[, c("variable", "mean", "q2.5", "q97.5")], row.names = FALSE)
  }
  .print_margin_params(x$var_params)

  cat("\nDiagnostics:\n")
  cat(sprintf("  Divergences: %d\n", x$diagnostics$n_divergent))
  cat(sprintf("  Max Rhat: %.3f\n", x$diagnostics$max_rhat))
  cat(sprintf("  Min ESS bulk: %.0f\n", x$diagnostics$min_ess_bulk))

  invisible(x)
}


#' @describeIn dcvar_multilevel_fit-methods Extract posterior means of
#'   population-level coefficients.
#' @details
#' Unlike single-level models (where `coef()` returns `$Phi`), the multilevel
#' model returns hierarchical parameters:
#' \describe{
#'   \item{`phi_bar`}{Population-mean VAR coefficients (analogous to `Phi`
#'     in single-level models, vectorised as phi11, phi12, phi21, phi22).}
#'   \item{`tau_phi`}{Between-unit SD of VAR coefficients.}
#'   \item{scale parameters}{`sigma` (innovation SDs) for normal margins,
#'     `sigma_exp` for exponential margins, or per-family scale/shape parameters
#'     (e.g. `sigma_eps`, `sigma_gam`, `shape_gam`) for per-variable,
#'     skew-normal, and gamma mixed-engine margins.}
#'   \item{`rho`}{Copula correlation (constant across units).}
#' }
#' Use [random_effects()] to obtain unit-specific VAR coefficients.
#' @return A named list of posterior means.
#' @export
coef.dcvar_multilevel_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  margins <- object$margins %||% "normal"
  scale_coef <- if (.uses_sem_multilevel_mixed_engine("multilevel", margins)) {
    .extract_margin_coefs(
      summ,
      .sem_multilevel_engine_margins("multilevel", margins, 2L)
    )
  } else if (identical(margins, "exponential")) {
    list(sigma_exp = .extract_required_coef(summ, "^sigma_exp\\[", "sigma_exp", "coef.dcvar_multilevel_fit()"))
  } else {
    list(sigma = .extract_required_coef(summ, "^sigma\\[", "sigma", "coef.dcvar_multilevel_fit()"))
  }

  c(list(
    phi_bar = .extract_required_coef(summ, "^phi_bar\\[", "phi_bar", "coef.dcvar_multilevel_fit()"),
    tau_phi = .extract_required_coef(summ, "^tau_phi\\[", "tau_phi", "coef.dcvar_multilevel_fit()")
  ), scale_coef, list(
    rho = .extract_required_coef(summ, "^rho$", "rho", "coef.dcvar_multilevel_fit()")
  ))
}


#' @export
coef.dcvar_multilevel_tv_fit <- function(object, ...) {
  summ <- .fit_summary(object$fit, backend = object$backend)
  margins_vec <- rep(object$margins %||% "normal", length.out = 2L)
  margin_coefs <- .extract_margin_coefs(summ, margins_vec)
  result <- list(
    phi_bar = .extract_required_coef(summ, "^phi_bar\\[", "phi_bar", "coef.dcvar_multilevel_tv_fit()"),
    tau_phi = .extract_required_coef(summ, "^tau_phi\\[", "tau_phi", "coef.dcvar_multilevel_tv_fit()")
  )
  extras <- list(
    rho = .extract_required_coef(summ, "^rho$", "rho", "coef.dcvar_multilevel_tv_fit()")
  )
  if (isTRUE(object$tv_phi)) {
    tau_phi_tv <- .extract_required_coef(summ, "^tau_phi_tv\\[", "tau_phi_tv", "coef.dcvar_multilevel_tv_fit()")
    names(tau_phi_tv) <- paste0("tau_phi_tv[", .tv_active_phi_coefs(object), "]")
    extras$tau_phi_tv <- tau_phi_tv
  }
  if (isTRUE(object$tv_sigma)) {
    extras$tau_sigma <- .extract_required_coef(summ, "^tau_sigma\\[", "tau_sigma", "coef.dcvar_multilevel_tv_fit()")
  }
  c(result, margin_coefs, extras)
}


#' @describeIn dcvar_multilevel_fit-methods Dispatch to a plot type.
#' @param type Character; one of `"random_effects"`, `"diagnostics"`.
#' @return A ggplot object.
#' @export
plot.dcvar_multilevel_fit <- function(x,
                                      type = c("random_effects", "diagnostics"),
                                      ...) {
  type <- match.arg(type)
  switch(type,
    random_effects = plot_random_effects(x, ...),
    diagnostics = plot_diagnostics(x, ...)
  )
}
