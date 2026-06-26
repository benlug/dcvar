# ============================================================================
# Applicability check for the flexible-margin copula-VAR
# ============================================================================

#' Check whether the flexible-margin copula-VAR is appropriate for the data
#'
#' Screens a fitted constant copula-VAR for the floor/ceiling pile-up pathology.
#' When a non-trivial point mass of observations sits exactly at a scale bound
#' (a floor or ceiling effect, common for raw single-item slider or Likert
#' responses), a flexible skewed margin (skew-normal, gamma, exponential) can
#' reproduce that pile-up only by concentrating the innovations at the bound,
#' which is achieved when the autoregressive prediction is flat. The margin then
#' absorbs the marginal shape and the VAR coefficients collapse toward zero, even
#' when the series is clearly autocorrelated. The decision whether the approach
#' suits a dataset cannot be made from the response format (a visual analog scale
#' looks continuous yet can have a heavy floor pile-up); it must be made from the
#' realized data and the fit.
#'
#' Three signals are combined:
#' \enumerate{
#'   \item \strong{Boundary atom}: the fraction of observations exactly at each
#'     series minimum and maximum, counted only when at least two observations
#'     coincide at the bound (a lone extreme value is not a point mass). A
#'     continuous margin cannot reproduce a point mass, so a non-trivial atom is
#'     a warning sign.
#'   \item \strong{Dynamics collapse}: the fitted self-lags against an OLS VAR(1)
#'     anchor (and, when supplied, a normal-margin reference fit). A
#'     flexible-margin fit whose self-lags collapse toward zero while an anchor
#'     shows clear autoregression (in either direction) is the signature of the
#'     pathology.
#'   \item \strong{Boundary and convergence flags}: a skew-normal slant estimated
#'     at its boundary, divergent transitions, or elevated split-Rhat.
#' }
#'
#' @param fit A `dcvar_constant_fit`.
#' @param reference Optional normal-margin `dcvar_constant_fit` on the same
#'   variables and number of observations. When supplied, the dynamics-collapse
#'   signal compares the fitted self-lags against this reference fit's self-lags
#'   in addition to the OLS VAR(1) anchor (which is always used). If `NULL`
#'   (default), only the OLS VAR(1) anchor is used.
#' @param atom_tol Minimum fraction at a bound to flag a point mass (default `0.05`).
#' @param delta_tol Magnitude of the skew-normal slant treated as a boundary value
#'   (default `0.98`).
#' @param phi_collapse_tol Self-lag magnitude treated as collapsed (default `0.10`).
#' @param ols_clear_tol Anchor self-lag magnitude treated as clear
#'   autoregression when detecting a dynamics collapse (default `0.20`).
#' @param rhat_tol Split-Rhat threshold (default `1.01`).
#'
#' @return An object of class `dcvar_applicability`: a list with the `verdict`
#'   (`"suitable"`, `"caution"`, or `"unsuitable"`), a `recommendation`, the
#'   triggered `reasons`, the variable names (`vars`) and `margins`, the
#'   boundary-atom fractions (`atom`, `atom_min`, `atom_max`), the fitted and
#'   OLS self-lags (`fit_self_lag`, `ols_self_lag`, plus `reference_self_lag`
#'   when a reference is supplied), the skew-normal slant means (`delta`, or
#'   `NULL`), the convergence `diagnostics`, and the `thresholds` used. Has a
#'   `print` method.
#'
#' @seealso [dcvar_constant()], [dcvar_diagnostics()], [var_params()].
#'
#' @examples
#' \donttest{
#' sim <- simulate_dcvar(
#'   n_time = 60,
#'   rho_trajectory = rho_constant(60, rho = 0.4),
#'   margins = "skew_normal",
#'   skew_params = list(alpha = c(4, 4)),
#'   seed = 1
#' )
#' fit <- dcvar_constant(
#'   sim$Y_df, vars = c("y1", "y2"), margins = "skew_normal",
#'   chains = 2, iter_warmup = 500, iter_sampling = 500, refresh = 0, seed = 1
#' )
#' applicability_check(fit)
#' }
#' @export
applicability_check <- function(fit, reference = NULL, atom_tol = 0.05,
                                delta_tol = 0.98, phi_collapse_tol = 0.10,
                                ols_clear_tol = 0.20, rhat_tol = 1.01) {
  if (!inherits(fit, "dcvar_constant_fit")) {
    cli_abort("{.fun applicability_check} expects a {.cls dcvar_constant_fit}, got {.cls {class(fit)[[1]]}}.")
  }
  Y <- fit$stan_data$Y
  if (is.null(Y)) cli_abort("Fitted object does not carry the data matrix {.field stan_data$Y}.")
  Y <- as.matrix(Y)
  if (anyNA(Y)) {
    cli_abort("{.field stan_data$Y} contains missing values; cannot run {.fun applicability_check}.")
  }
  n <- nrow(Y); D <- ncol(Y)
  vars <- fit$vars %||% paste0("y", seq_len(D))
  margins <- fit$margins %||% "normal"
  # margins is stored as a single string when homogeneous; expand to length D so
  # per-dimension indexing (e.g. which dims are skew-normal) is well defined.
  margins_full <- if (length(margins) == 1L) rep(margins, D) else margins
  is_flexible <- any(margins %in% c("skew_normal", "gamma", "exponential"))

  if (!is.null(reference)) {
    if (!inherits(reference, "dcvar_constant_fit")) {
      cli_abort("{.arg reference} must be a {.cls dcvar_constant_fit}, got {.cls {class(reference)[[1]]}}.")
    }
    Yref <- reference$stan_data$Y
    if (is.null(Yref)) {
      cli_abort("{.arg reference} does not carry the data matrix {.field stan_data$Y}.")
    }
    Yref <- as.matrix(Yref)
    if (ncol(Yref) != D) {
      cli_abort("{.arg reference} must be fitted on the same number of variables as {.arg fit} ({D}).")
    }
    if (nrow(Yref) != n) {
      cli_abort("{.arg reference} must be fitted on the same number of observations as {.arg fit} ({n}).")
    }
    ref_vars <- reference$vars %||% vars
    if (length(ref_vars) != D || !identical(as.character(ref_vars), as.character(vars))) {
      cli_abort("{.arg reference} must be fitted on the same variables as {.arg fit}, in the same order.")
    }
    ref_margins <- reference$margins %||% "normal"
    if (!all(ref_margins == "normal")) {
      cli_abort("{.arg reference} must be a normal-margin fit ({.arg margins = 'normal'}).")
    }
  }

  ## var_params() runs a full posterior re-summary, so compute it once and reuse.
  vp <- var_params(fit)

  ## (1) boundary atom: fraction exactly at each series min / max, but only when
  ## at least two observations tie at the bound (a single extreme is not a mass).
  atom_frac <- function(col, extreme) {
    k <- sum(col == extreme)
    if (k >= 2L) k / n else 0
  }
  atom_min <- vapply(seq_len(D), function(d) atom_frac(Y[, d], min(Y[, d])), numeric(1))
  atom_max <- vapply(seq_len(D), function(d) atom_frac(Y[, d], max(Y[, d])), numeric(1))
  atom <- pmax(atom_min, atom_max)
  atom_present <- any(atom >= atom_tol)

  ## (2) dynamics: OLS VAR(1) anchor (and optional reference) vs fitted self-lags
  zv <- vapply(seq_len(D), function(d) isTRUE(stats::sd(Y[, d]) == 0), logical(1))
  Ylag <- Y[-n, , drop = FALSE]; Ycur <- Y[-1, , drop = FALSE]
  ols_self <- vapply(seq_len(D), function(d) {
    unname(stats::coef(stats::lm(Ycur[, d] ~ Ylag))[1L + d])
  }, numeric(1))
  self_from <- function(phi_df) vapply(seq_len(D), function(d) {
    row <- phi_df[phi_df$variable == sprintf("Phi[%d,%d]", d, d), , drop = FALSE]
    if (nrow(row) == 0) NA_real_ else row$mean[[1]]
  }, numeric(1))
  fit_self <- self_from(vp$Phi)
  ref_self <- if (is.null(reference)) NULL else self_from(var_params(reference)$Phi)
  fit_self_bad <- !is.finite(fit_self)
  ols_self_bad <- !is.finite(ols_self)
  ref_self_bad <- if (is.null(ref_self)) rep(FALSE, D) else !is.finite(ref_self)
  collapsed_vs <- function(anchor) {
    any(is.finite(fit_self) & is.finite(anchor) &
          abs(fit_self) < phi_collapse_tol &
          abs(anchor) > ols_clear_tol, na.rm = TRUE)
  }
  collapse <- collapsed_vs(ols_self) || (!is.null(ref_self) && collapsed_vs(ref_self))

  ## (3) slant at boundary (any skew-normal dimension) and convergence
  delta_vals <- NULL; delta_rail <- FALSE; delta_bad <- FALSE; delta_names <- NULL
  if (any(margins == "skew_normal")) {
    vpd <- vp$delta
    if (!is.null(vpd) && length(vpd$mean) > 0L) {
      delta_vals <- vpd$mean
      delta_bad <- any(!is.finite(delta_vals))
      delta_rail <- any(abs(delta_vals) >= delta_tol, na.rm = TRUE)
      ## $delta carries one row per skew-normal dimension (all D in the
      ## homogeneous case, the skew-normal subset in the mixed case).
      skew_idx <- which(margins_full == "skew_normal")
      delta_names <- if (length(skew_idx) == length(delta_vals)) {
        vars[skew_idx]
      } else {
        vars[seq_along(delta_vals)]
      }
    }
  }
  diag <- dcvar_diagnostics(fit)
  conv_bad <- diag$n_divergent > 0 ||
    (is.finite(diag$max_rhat) && diag$max_rhat > rhat_tol)

  ## reasons (every triggered signal, regardless of the headline verdict)
  reasons <- character(0)
  if (atom_present)
    reasons <- c(reasons, sprintf("boundary atom present (up to %.0f%% of observations at a scale bound)", 100 * max(atom)))
  if (collapse)
    reasons <- c(reasons, "fitted self-lags collapse toward zero while the VAR(1) anchor shows clear autoregression")
  if (any(fit_self_bad))
    reasons <- c(reasons, sprintf("fitted self-lag summary undefined for series (%s)", paste(vars[fit_self_bad], collapse = ", ")))
  if (any(ols_self_bad))
    reasons <- c(reasons, sprintf("OLS VAR(1) anchor self-lag undefined for series (%s)", paste(vars[ols_self_bad], collapse = ", ")))
  if (any(ref_self_bad))
    reasons <- c(reasons, sprintf("reference self-lag undefined for series (%s)", paste(vars[ref_self_bad], collapse = ", ")))
  if (delta_rail)
    reasons <- c(reasons, sprintf("skew-normal slant at its boundary (delta = %s)", paste(sprintf("%.2f", delta_vals), collapse = ", ")))
  if (delta_bad)
    reasons <- c(reasons, sprintf("skew-normal slant summary is non-finite (delta = %s)", paste(sprintf("%.2f", delta_vals), collapse = ", ")))
  if (diag$n_divergent > 0)
    reasons <- c(reasons, sprintf("%d divergent transition%s", diag$n_divergent, if (diag$n_divergent == 1) "" else "s"))
  if (is.finite(diag$max_rhat) && diag$max_rhat > rhat_tol)
    reasons <- c(reasons, sprintf("max split-Rhat %.3f", diag$max_rhat))
  if (any(zv))
    reasons <- c(reasons, sprintf("zero-variance series (%s): self-lag undefined", paste(vars[zv], collapse = ", ")))

  ## verdict: convergence problems are a first-class signal and block "suitable"
  undefined_anchor <- any(fit_self_bad) || any(ols_self_bad) || any(ref_self_bad) || delta_bad
  verdict <- if (!atom_present && !collapse && !delta_rail && !conv_bad && !undefined_anchor) {
    "suitable"
  } else if (is_flexible && atom_present && (collapse || delta_rail || diag$n_divergent > 0)) {
    "unsuitable"
  } else {
    "caution"
  }
  recommendation <- switch(
    verdict,
    suitable = "No boundary atom detected and the sampler converged. The flexible-margin copula-VAR is appropriate; the asymmetry is a continuous feature of the innovations, separable from the dynamics.",
    unsuitable = paste0(
      "The flexible margin absorbs a boundary pile-up and collapses the dynamics. ",
      "Read the autoregressive effects from a normal margin (Gaussian quasi-maximum-likelihood recovers the conditional mean, and hence Phi, when the conditional mean is linear and only the innovation shape is distorted). ",
      "Note, however, that a genuine boundary atom usually reflects censoring, under which the conditional mean is nonlinear and the normal-margin Phi can itself be attenuated toward zero (the amount depends on the censoring severity); a censored latent VAR is then the principled remedy."
    ),
    caution = if (!is_flexible && atom_present && !collapse) {
      "A boundary atom is present. This normal-margin fit recovers the dynamics, but the margin does not model the atom, and a flexible skewed margin would risk collapsing the dynamics on these data. Prefer the normal margin for the dynamics, or a censored latent VAR if the boundedness itself is of interest (with a genuine atom the normal-margin Phi can itself be attenuated toward zero, depending on the censoring severity)."
    } else {
      sig <- c(
        if (atom_present) "a boundary atom is present",
        if (collapse) "the fitted self-lags collapse relative to the VAR(1) anchor",
        if (delta_rail) "the skew-normal slant sits at its boundary",
        if (undefined_anchor) "one or more summary anchors are undefined",
        if (conv_bad) "the sampler shows convergence problems (divergences or elevated Rhat)"
      )
      sprintf(
        "Mixed signals (%s). Inspect the fit before trusting the flexible-margin dynamics: check the slant, divergences, and the posterior-predictive fraction at the bound, and consider reading the dynamics from a normal margin or modelling the boundedness with a censored latent VAR.",
        paste(sig, collapse = "; ")
      )
    }
  )

  structure(
    list(
      verdict = verdict,
      recommendation = recommendation,
      reasons = reasons,
      vars = vars,
      margins = margins,
      atom = stats::setNames(round(atom, 3), vars),
      atom_min = stats::setNames(round(atom_min, 3), vars),
      atom_max = stats::setNames(round(atom_max, 3), vars),
      fit_self_lag = stats::setNames(round(fit_self, 3), vars),
      ols_self_lag = stats::setNames(round(ols_self, 3), vars),
      reference_self_lag = if (is.null(ref_self)) NULL else stats::setNames(round(ref_self, 3), vars),
      delta = if (is.null(delta_vals)) NULL else stats::setNames(round(delta_vals, 3), delta_names),
      diagnostics = diag,
      thresholds = list(atom_tol = atom_tol, delta_tol = delta_tol,
                        phi_collapse_tol = phi_collapse_tol,
                        ols_clear_tol = ols_clear_tol, rhat_tol = rhat_tol)
    ),
    class = "dcvar_applicability"
  )
}

#' Print a flexible-margin applicability check
#'
#' @param x A `dcvar_applicability` object, as returned by
#'   [applicability_check()].
#' @param ... Additional arguments (unused).
#' @return Invisibly returns `x`.
#' @export
print.dcvar_applicability <- function(x, ...) {
  status <- switch(x$verdict,
    suitable   = "SUITABLE",
    caution    = "CAUTION",
    unsuitable = "UNSUITABLE")
  cat(sprintf("Applicability check for the flexible-margin copula-VAR: %s\n", status))
  cat(sprintf("  margins: %s\n", paste(x$margins, collapse = ", ")))
  cat(sprintf("  boundary atom (max fraction at a bound): %s\n",
              paste(sprintf("%s %.0f%%", names(x$atom), 100 * x$atom), collapse = " | ")))
  cat(sprintf("  self-lags fitted vs OLS: %s\n",
              paste(sprintf("%s %.2f/%.2f", names(x$fit_self_lag), x$fit_self_lag, x$ols_self_lag), collapse = " | ")))
  if (!is.null(x$reference_self_lag))
    cat(sprintf("  self-lags fitted vs reference: %s\n",
                paste(sprintf("%s %.2f/%.2f", names(x$fit_self_lag), x$fit_self_lag, x$reference_self_lag), collapse = " | ")))
  if (!is.null(x$delta))
    cat(sprintf("  skew-normal slant delta: %s\n",
                paste(sprintf("%s %.2f", names(x$delta), x$delta), collapse = " | ")))
  cat(sprintf("  divergences: %d | max Rhat: %.3f\n", x$diagnostics$n_divergent, x$diagnostics$max_rhat))
  if (length(x$reasons))
    for (r in x$reasons) cat(sprintf("  - %s\n", r))
  cat(sprintf("  recommendation: %s\n", x$recommendation))
  invisible(x)
}
