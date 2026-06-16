# ============================================================================
# Extraction Functions
# ============================================================================

#' Extract the rho trajectory with credible intervals
#'
#' Returns a data frame with the posterior mean, SD, and quantiles of the
#' time-varying correlation at each time point.
#'
#' @param object A fitted model object (`dcvar_fit`, `dcvar_covariate_fit`,
#'   `dcvar_hmm_fit`, `dcvar_constant_fit`, `dcvar_multilevel_fit`, or
#'   `dcvar_sem_fit`).
#' @param probs Numeric vector of quantile probabilities (default: `c(0.025, 0.1, 0.5, 0.9, 0.975)`).
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `time`, `mean`, `sd`, and one column per
#'   quantile (e.g., `q2.5`, `q10`, `q50`, `q90`, `q97.5`). For
#'   `dcvar_constant_fit` objects, the constant rho is expanded to all `n_time - 1`
#'   time points for consistency with the time-varying models.
#'
#' @seealso [plot_rho()] to visualise the trajectory,
#'   [interpret_rho_trajectory()] for a text-based summary,
#'   [var_params()] for VAR parameter extraction.
#' @export
rho_trajectory <- function(object, ...) {
  UseMethod("rho_trajectory")
}

#' @rdname rho_trajectory
#' @export
rho_trajectory.default <- function(object, ...) {
  cli_abort("{.fun rho_trajectory} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' Extract a unified dependence summary
#'
#' Returns posterior summaries for Kendall's tau, using the fitted copula
#' family to transform the model-specific dependence parameter. For Gaussian
#' copulas, `tau = 2 / pi * asin(rho)`. For Clayton copulas,
#' `tau = theta / (theta + 2)`.
#'
#' @param object A fitted model object.
#' @param probs Numeric vector of quantile probabilities.
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `time`, `mean`, `sd`, and one column per
#'   requested quantile.
#' @export
dependence_summary <- function(object, ...) {
  UseMethod("dependence_summary")
}

#' @rdname dependence_summary
#' @export
dependence_summary.default <- function(object, ...) {
  cli_abort("{.fun dependence_summary} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' Internal helper: recover observed time values from Stan data
#' @noRd
.observed_time_values <- function(stan_data, drop_first = FALSE) {
  time_values <- attr(stan_data, "time_values")
  if (is.null(time_values)) {
    time_values <- seq_len(stan_data$n_time)
  }

  if (drop_first) {
    time_values[-1]
  } else {
    time_values
  }
}

#' Internal helper: summarise time-varying rho draws into a data frame
#' @noRd
.summarise_rho_draws <- function(rho_draws, probs, time_values = NULL) {
  if (!is.numeric(probs) || !all(probs >= 0 & probs <= 1)) {
    cli_abort("{.arg probs} must be numeric values in [0, 1].")
  }

  n_time_eff <- ncol(rho_draws)
  if (is.null(time_values)) {
    time_values <- 2:(n_time_eff + 1)
  }
  if (length(time_values) != n_time_eff) {
    cli_abort("Time axis length does not match rho draw length.")
  }

  rho_summary <- data.frame(
    time = time_values,
    mean = colMeans(rho_draws),
    sd = apply(rho_draws, 2, sd)
  )

  quantiles <- apply(rho_draws, 2, quantile, probs = probs)
  if (is.null(dim(quantiles))) {
    quantiles <- matrix(quantiles, nrow = 1)
  }
  for (i in seq_along(probs)) {
    rho_summary[[paste0("q", probs[i] * 100)]] <- quantiles[i, ]
  }

  rho_summary
}

#' @rdname rho_trajectory
#' @export
rho_trajectory.dcvar_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  rho_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "rho", backend = object$backend,
    required = .stan_output_group_pattern("rho"),
    required_type = "pattern",
    context = "rho_trajectory.dcvar_fit()",
    output_type = "transformed parameter group"
  ))
  .summarise_rho_draws(rho_draws, probs, .observed_time_values(object$stan_data, drop_first = TRUE))
}

#' @rdname dependence_summary
#' @export
dependence_summary.dcvar_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  rho_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "rho", backend = object$backend,
    required = .stan_output_group_pattern("rho"),
    required_type = "pattern",
    context = "dependence_summary.dcvar_fit()",
    output_type = "transformed parameter group"
  ))
  tau_draws <- 2 / pi * asin(rho_draws)
  .summarise_rho_draws(tau_draws, probs, .observed_time_values(object$stan_data, drop_first = TRUE))
}

#' @rdname rho_trajectory
#' @export
rho_trajectory.dcvar_covariate_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  rho_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "rho", backend = object$backend,
    required = .stan_output_group_pattern("rho"),
    required_type = "pattern",
    context = "rho_trajectory.dcvar_covariate_fit()",
    output_type = "transformed parameter group"
  ))
  .summarise_rho_draws(rho_draws, probs, .observed_time_values(object$stan_data, drop_first = TRUE))
}

#' @rdname dependence_summary
#' @export
dependence_summary.dcvar_covariate_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  rho_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "rho", backend = object$backend,
    required = .stan_output_group_pattern("rho"),
    required_type = "pattern",
    context = "dependence_summary.dcvar_covariate_fit()",
    output_type = "transformed parameter group"
  ))
  tau_draws <- 2 / pi * asin(rho_draws)
  .summarise_rho_draws(tau_draws, probs, .observed_time_values(object$stan_data, drop_first = TRUE))
}

#' @rdname rho_trajectory
#' @export
rho_trajectory.dcvar_hmm_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  rho_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "rho_hmm", backend = object$backend,
    required = .stan_output_group_pattern("rho_hmm"),
    required_type = "pattern",
    context = "rho_trajectory.dcvar_hmm_fit()",
    output_type = "generated quantity"
  ))
  .summarise_rho_draws(rho_draws, probs, .observed_time_values(object$stan_data, drop_first = TRUE))
}

#' @rdname dependence_summary
#' @export
dependence_summary.dcvar_hmm_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  # The HMM dependence at time t is a state mixture, so Kendall's tau must be
  # averaged on the tau scale: sum_k gamma[t, k] * (2/pi) asin(rho_state[k]).
  # Applying asin to the gamma-weighted average rho (rho_hmm) understates tau
  # whenever the smoothed state probabilities are mixed.
  gamma_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "gamma", backend = object$backend,
    required = .stan_output_group_pattern("gamma"),
    required_type = "pattern",
    context = "dependence_summary.dcvar_hmm_fit()",
    output_type = "generated quantity"
  ))
  rho_state_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "rho_state", backend = object$backend,
    required = .stan_output_group_pattern("rho_state"),
    required_type = "pattern",
    context = "dependence_summary.dcvar_hmm_fit()",
    output_type = "transformed parameter group"
  ))

  n_time_eff <- object$stan_data$n_time - 1L
  K <- ncol(rho_state_draws)
  tau_state <- 2 / pi * asin(rho_state_draws)

  tau_draws <- matrix(0, nrow(gamma_draws), n_time_eff)
  for (k in seq_len(K)) {
    gamma_cols <- paste0("gamma[", seq_len(n_time_eff), ",", k, "]")
    tau_draws <- tau_draws + as.matrix(gamma_draws[, gamma_cols, drop = FALSE]) * as.numeric(tau_state[, k])
  }
  .summarise_rho_draws(tau_draws, probs, .observed_time_values(object$stan_data, drop_first = TRUE))
}

#' Internal helper: constant rho trajectory (shared by constant and multilevel)
#' @noRd
.rho_trajectory_constant_impl <- function(object, probs) {
  if (!is.numeric(probs) || !all(probs >= 0 & probs <= 1)) {
    cli_abort("{.arg probs} must be numeric values in [0, 1].")
  }

  rho_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "rho", backend = object$backend,
    required = "rho",
    context = ".rho_trajectory_constant_impl()",
    output_type = "transformed parameter"
  ))

  n_time_obs <- object$stan_data$n_time
  n_time_eff <- n_time_obs - 1L

  rho_mean <- mean(rho_draws[, 1])
  rho_sd <- sd(rho_draws[, 1])
  quants <- quantile(rho_draws[, 1], probs = probs)
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)

  rho_summary <- data.frame(
    time = time_values,
    mean = rep(rho_mean, n_time_eff),
    sd = rep(rho_sd, n_time_eff)
  )

  for (i in seq_along(probs)) {
    rho_summary[[paste0("q", probs[i] * 100)]] <- rep(quants[i], n_time_eff)
  }

  rho_summary
}

#' Internal helper: constant Kendall's tau summary.
#' @noRd
.dependence_summary_constant_impl <- function(object, probs) {
  if (!is.numeric(probs) || !all(probs >= 0 & probs <= 1)) {
    cli_abort("{.arg probs} must be numeric values in [0, 1].")
  }

  copula <- object$copula %||% "gaussian"
  if (identical(copula, "clayton")) {
    theta_draws <- posterior::as_draws_matrix(.fit_draws(
      object$fit, "theta", backend = object$backend,
      required = "theta",
      context = ".dependence_summary_constant_impl()",
      output_type = "parameter"
    ))
    dep_draws <- theta_draws[, 1] / (theta_draws[, 1] + 2)
  } else {
    rho_draws <- posterior::as_draws_matrix(.fit_draws(
      object$fit, "rho", backend = object$backend,
      required = "rho",
      context = ".dependence_summary_constant_impl()",
      output_type = "transformed parameter"
    ))
    dep_draws <- 2 / pi * asin(rho_draws[, 1])
  }

  n_time_obs <- object$stan_data$n_time
  n_time_eff <- n_time_obs - 1L
  quants <- quantile(dep_draws, probs = probs)
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)

  dep_summary <- data.frame(
    time = time_values,
    mean = rep(mean(dep_draws), n_time_eff),
    sd = rep(sd(dep_draws), n_time_eff)
  )

  for (i in seq_along(probs)) {
    dep_summary[[paste0("q", probs[i] * 100)]] <- rep(quants[i], n_time_eff)
  }

  dep_summary
}

#' @rdname rho_trajectory
#' @export
rho_trajectory.dcvar_constant_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  if (identical(object$copula %||% "gaussian", "clayton")) {
    cli_abort("Clayton fits do not have a Gaussian copula rho. Use {.fun dependence_summary} for Kendall's tau.")
  }
  .rho_trajectory_constant_impl(object, probs)
}

#' @rdname dependence_summary
#' @export
dependence_summary.dcvar_constant_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  .dependence_summary_constant_impl(object, probs)
}


#' Extract VAR(1) parameter summaries
#'
#' Returns posterior summaries for the VAR parameters: intercepts (mu),
#' coefficients (Phi), innovation SDs (sigma_eps), and sigma_omega (DC-VAR only).
#'
#' @param object A fitted model object.
#' @param ... Additional arguments (unused).
#'
#' @return A named list of data frames with columns `variable`, `mean`, `sd`,
#'   `q2.5`, `q97.5`.
#' @export
var_params <- function(object, ...) {
  UseMethod("var_params")
}

#' @rdname var_params
#' @export
var_params.default <- function(object, ...) {
  cli_abort("{.fun var_params} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' @rdname var_params
#' @export
var_params.dcvar_model_fit <- function(object, ...) {
  margins <- object$margins %||% "normal"
  copula <- object$copula %||% "gaussian"
  mixed <- .is_mixed_margins(margins)
  required_patterns <- c("^mu\\[", "^Phi\\[")
  if (mixed) {
    margin_groups <- names(.mixed_margin_report_vars(margins))
    required_patterns <- c(required_patterns, paste0("^", margin_groups, "\\["))
  } else {
    switch(margins[[1L]],
      normal = {
        required_patterns <- c(required_patterns, "^sigma_eps\\[")
      },
      exponential = {
        required_patterns <- c(required_patterns, "^sigma_exp\\[")
      },
      skew_normal = {
        required_patterns <- c(required_patterns, "^omega\\[", "^delta\\[")
      },
      gamma = {
        required_patterns <- c(required_patterns, "^sigma_gam\\[", "^shape_gam$")
      },
      {
        required_patterns <- c(required_patterns, "^sigma_eps\\[")
      }
    )
  }
  has_sigma_omega <- identical(object$model, "dcvar") ||
    identical(object$model, "dcvar_covariate")
  if (has_sigma_omega) {
    required_patterns <- c(required_patterns, "^sigma_omega$")
  }
  if (identical(copula, "clayton")) {
    required_patterns <- c(required_patterns, "^theta$")
  }

  summ <- .fit_summary(
    object$fit, variables = NULL, backend = object$backend,
    required = required_patterns,
    required_type = "pattern",
    context = "var_params.dcvar_model_fit()",
    output_type = "parameter group",
    mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )
  extract_param <- function(pattern) {
    rows <- grep(pattern, summ$variable)
    data.frame(
      variable = summ$variable[rows],
      mean = summ$mean[rows],
      sd = summ$sd[rows],
      q2.5 = summ$q2.5[rows],
      q97.5 = summ$q97.5[rows]
    )
  }

  result <- list(
    mu = extract_param("^mu\\["),
    Phi = extract_param("^Phi\\[")
  )

  # Margin-specific scale parameters
  if (mixed) {
    # Report each dimension under its own family, restricted to that family's
    # dimensions (e.g. sigma_eps[1] for a normal dim, sigma_exp[2] for an
    # exponential dim).
    extract_vars <- function(vars) {
      rows <- match(vars, summ$variable)
      rows <- rows[!is.na(rows)]
      data.frame(
        variable = summ$variable[rows],
        mean = summ$mean[rows],
        sd = summ$sd[rows],
        q2.5 = summ$q2.5[rows],
        q97.5 = summ$q97.5[rows]
      )
    }
    specs <- .mixed_margin_report_vars(margins)
    for (nm in names(specs)) {
      result[[nm]] <- extract_vars(specs[[nm]])
    }
  } else {
    switch(margins[[1L]],
      normal = {
        result$sigma_eps <- extract_param("^sigma_eps\\[")
      },
      exponential = {
        result$sigma_exp <- extract_param("^sigma_exp\\[")
      },
      skew_normal = {
        result$omega <- extract_param("^omega\\[")
        result$delta <- extract_param("^delta\\[")
      },
      gamma = {
        result$sigma_gam <- extract_param("^sigma_gam\\[")
        result$shape_gam <- extract_param("^shape_gam$")
      },
      {
        result$sigma_eps <- extract_param("^sigma_eps\\[")
      }
    )
  }

  # sigma_omega is present in random-walk DC-VAR variants only.
  so <- if (has_sigma_omega) extract_param("^sigma_omega$") else NULL
  if (!is.null(so)) result$sigma_omega <- so
  if (identical(copula, "clayton")) {
    result$theta <- extract_param("^theta$")
  }

  result
}


#' Extract covariate effect summaries
#'
#' Returns posterior summaries for the Fisher-z intercept and covariate
#' effects. The residual random-walk innovation scale `sigma_omega` is
#' reported separately by [var_params()] for `drift = TRUE` fits.
#'
#' @param object A `dcvar_covariate_fit` object.
#' @param probs Numeric vector of quantile probabilities (default:
#'   `c(0.025, 0.5, 0.975)`).
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with one row per effect and columns `term`, `variable`,
#'   `mean`, `sd`, and one column per requested quantile.
#' @export
covariate_effects <- function(object, ...) {
  UseMethod("covariate_effects")
}

#' @rdname covariate_effects
#' @export
covariate_effects.default <- function(object, ...) {
  cli_abort("{.fun covariate_effects} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' @rdname covariate_effects
#' @export
covariate_effects.dcvar_covariate_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  if (!is.numeric(probs) || !all(probs >= 0 & probs <= 1)) {
    cli_abort("{.arg probs} must be numeric values in [0, 1].")
  }

  summ <- .fit_summary(
    object$fit,
    variables = NULL,
    backend = object$backend,
    required = c("^beta_0$", "^beta\\["),
    required_type = "pattern",
    context = "covariate_effects.dcvar_covariate_fit()",
    output_type = "parameter group",
    mean,
    sd,
    ~posterior::quantile2(.x, probs = probs)
  )

  rows <- c(
    grep("^beta_0$", summ$variable),
    grep("^beta\\[", summ$variable)
  )
  out <- data.frame(summ[rows, , drop = FALSE], row.names = NULL)
  out$term <- out$variable
  beta_rows <- grep("^beta\\[", out$variable)
  if (length(beta_rows) > 0L) {
    out$term[beta_rows] <- object$covariates
  }
  out$term[out$variable == "beta_0"] <- "(Intercept)"

  out[, c("term", setdiff(names(out), "term")), drop = FALSE]
}


#' Extract HMM state information
#'
#' Returns state posteriors, Viterbi path, state-specific rho values,
#' and the transition matrix from an HMM copula fit.
#'
#' @param object A `dcvar_hmm_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @return A named list with:
#'   - `gamma`: T_eff x K matrix of posterior state probabilities
#'   - `viterbi`: integer vector, the Viterbi (MAP) state path decoded from
#'     the posterior-mean emission log-likelihoods, transition matrix, and
#'     initial state probabilities (a plug-in estimator, not a full posterior
#'     summary of the path)
#'   - `rho_state`: list with `mean`, `lower`, `upper` for each state
#'   - `A`: K x K posterior mean transition matrix
#'   - `rho_hmm`: posterior-averaged rho trajectory
#' @export
hmm_states <- function(object, ...) {
  UseMethod("hmm_states")
}

#' Internal: Viterbi decoding for fixed log-scale HMM quantities
#' @noRd
.viterbi_path <- function(obs_ll, log_A, log_pi0) {
  n_time <- nrow(obs_ll)
  K <- ncol(obs_ll)

  log_delta <- matrix(-Inf, n_time, K)
  back_ptr <- matrix(NA_integer_, n_time, K)
  log_delta[1, ] <- log_pi0 + obs_ll[1, ]
  if (n_time >= 2) {
    for (t in 2:n_time) {
      for (k in seq_len(K)) {
        cand <- log_delta[t - 1, ] + log_A[, k]
        j <- which.max(cand)
        back_ptr[t, k] <- j
        log_delta[t, k] <- cand[j] + obs_ll[t, k]
      }
    }
  }

  path <- integer(n_time)
  path[n_time] <- which.max(log_delta[n_time, ])
  if (n_time >= 2) {
    for (t in seq(n_time - 1L, 1L)) {
      path[t] <- back_ptr[t + 1L, path[t + 1L]]
    }
  }
  path
}

#' @rdname hmm_states
#' @export
hmm_states.default <- function(object, ...) {
  cli_abort("{.fun hmm_states} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' @rdname hmm_states
#' @export
hmm_states.dcvar_hmm_fit <- function(object, ...) {
  K <- object$K
  n_time_eff <- object$stan_data$n_time - 1L

  .safe_draws <- function(var_name) {
    tryCatch(
      posterior::as_draws_matrix(.fit_draws(
        object$fit, var_name, backend = object$backend,
        required = .stan_output_group_pattern(var_name),
        required_type = "pattern",
        context = "hmm_states.dcvar_hmm_fit()",
        output_type = "Stan output group"
      )),
      error = function(e) {
        cli_abort("Failed to extract {.val {var_name}} draws: {e$message}")
      }
    )
  }

  # State posteriors: gamma[t, k]
  gamma_draws <- .safe_draws("gamma")

  # Validate gamma columns exist upfront
  expected_gamma_cols <- paste0("gamma[", rep(seq_len(n_time_eff), K), ",", rep(seq_len(K), each = n_time_eff), "]")
  missing_gamma <- setdiff(expected_gamma_cols, colnames(gamma_draws))
  if (length(missing_gamma) > 0) {
    cli_abort("Expected gamma columns missing from draws: {.val {head(missing_gamma, 5)}}")
  }

  gamma_mean <- matrix(NA_real_, n_time_eff, K)
  for (k in 1:K) {
    cols <- paste0("gamma[", seq_len(n_time_eff), ",", k, "]")
    gamma_mean[, k] <- colMeans(gamma_draws[, cols, drop = FALSE])
  }

  # State-specific rho
  rho_state_draws <- .safe_draws("rho_state")
  rho_state_mean <- colMeans(rho_state_draws)
  rho_state_lower <- apply(rho_state_draws, 2, quantile, 0.025)
  rho_state_upper <- apply(rho_state_draws, 2, quantile, 0.975)

  # Transition matrix
  A_draws <- .safe_draws("A")

  # Validate A columns exist upfront
  expected_A_cols <- paste0("A[", rep(seq_len(K), K), ",", rep(seq_len(K), each = K), "]")
  missing_A <- setdiff(expected_A_cols, colnames(A_draws))
  if (length(missing_A) > 0) {
    cli_abort("Expected transition matrix columns missing from draws: {.val {head(missing_A, 5)}}")
  }

  A_mean <- matrix(NA_real_, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      col <- paste0("A[", i, ",", j, "]")
      A_mean[i, j] <- mean(A_draws[, col])
    }
  }

  # MAP state sequence: Viterbi decoding on posterior-mean log quantities.
  # (The previous most-frequent-sampled-path estimator degenerates to an
  # arbitrary single draw's path when parameter uncertainty makes nearly all
  # sampled paths unique.)
  obs_ll_draws <- .safe_draws("obs_ll")
  obs_ll_mean <- matrix(NA_real_, n_time_eff, K)
  for (k in 1:K) {
    cols <- paste0("obs_ll[", seq_len(n_time_eff), ",", k, "]")
    obs_ll_mean[, k] <- colMeans(obs_ll_draws[, cols, drop = FALSE])
  }
  pi0_draws <- .safe_draws("pi0")
  pi0_mean <- colMeans(pi0_draws[, paste0("pi0[", seq_len(K), "]"), drop = FALSE])
  viterbi_map <- .viterbi_path(obs_ll_mean, log(A_mean), log(pi0_mean))

  # Posterior-averaged rho
  rho_hmm_draws <- .safe_draws("rho_hmm")
  rho_hmm_mean <- colMeans(rho_hmm_draws)

  list(
    gamma = gamma_mean,
    viterbi = viterbi_map,
    rho_state = list(
      mean = rho_state_mean,
      lower = rho_state_lower,
      upper = rho_state_upper
    ),
    A = A_mean,
    rho_hmm = rho_hmm_mean
  )
}


#' @rdname var_params
#' @details
#' For a Markov-switching HMM fit (`dcvar_hmm()` with `switch` beyond `"rho"` or
#' per-state `margins`), the VAR parameters are state-indexed: `mu[k, d]`, the Phi
#' baseline `Phi_base[i, j]` plus per-state deviations `Phi_dev[k, i, j]` on the
#' switching coefficients, and per-state margin scales. A non-switching HMM fit
#' delegates to the shared method. See [hmm_state_params()] for a per-state
#' assembled view.
#' @export
var_params.dcvar_hmm_fit <- function(object, ...) {
  if (!isTRUE(object$switching)) {
    return(NextMethod())
  }

  sw <- object$switch
  K <- object$K
  Mrg_K <- if (isTRUE(sw$margins == 1L)) K else 1L
  report <- .hmm_switching_report_vars(object$margins_matrix, Mrg_K)

  required_patterns <- c("^mu\\[", "^Phi_base\\[", "^rho_state\\[",
                         paste0("^", names(report), "\\["))
  if (any(sw$phi_mask > 0L)) {
    required_patterns <- c(required_patterns, "^Phi_dev\\[")
  }

  summ <- .fit_summary(
    object$fit, variables = NULL, backend = object$backend,
    required = required_patterns,
    required_type = "pattern",
    context = "var_params.dcvar_hmm_fit()",
    output_type = "parameter group",
    mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )
  extract_param <- function(pattern) {
    rows <- grep(pattern, summ$variable)
    data.frame(
      variable = summ$variable[rows],
      mean = summ$mean[rows],
      sd = summ$sd[rows],
      q2.5 = summ$q2.5[rows],
      q97.5 = summ$q97.5[rows]
    )
  }
  extract_vars <- function(vars) {
    rows <- match(vars, summ$variable)
    rows <- rows[!is.na(rows)]
    data.frame(
      variable = summ$variable[rows],
      mean = summ$mean[rows],
      sd = summ$sd[rows],
      q2.5 = summ$q2.5[rows],
      q97.5 = summ$q97.5[rows]
    )
  }

  result <- list(
    mu = extract_param("^mu\\["),
    Phi = extract_param("^Phi_base\\[")
  )
  if (any(sw$phi_mask > 0L)) {
    result$Phi_dev <- extract_param("^Phi_dev\\[")
  }
  for (nm in names(report)) {
    result[[nm]] <- extract_vars(report[[nm]])
  }
  result$rho_state <- extract_param("^rho_state\\[")
  result
}


#' Extract per-state parameters from a Markov-switching HMM fit
#'
#' Assembles the state-specific VAR parameters of a switching `dcvar_hmm()` fit
#' into a per-state view: the intercepts, the effective VAR(1) coefficient matrix
#' (`Phi_base` plus the per-state deviations on the switching coefficients), the
#' margin scale/shape parameters, and the per-state marginal family labels. For
#' components that are NOT state-specific the shared estimate is tiled across all
#' states, so the return shape does not depend on which components switch.
#'
#' @param object A `dcvar_hmm_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @return A list with `families` (a K x D character matrix), `rho_state`
#'   (length-K), `mu` (a K x D matrix of posterior means), `Phi` (a length-K list
#'   of D x D posterior-mean coefficient matrices), and `scales` (a data frame of
#'   per-state margin scale/shape posterior means; shared scales are tiled across
#'   states).
#' @seealso [hmm_states()], [var_params()]
#' @export
hmm_state_params <- function(object, ...) {
  UseMethod("hmm_state_params")
}

#' @rdname hmm_state_params
#' @export
hmm_state_params.default <- function(object, ...) {
  cli_abort("{.fun hmm_state_params} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' @rdname hmm_state_params
#' @export
hmm_state_params.dcvar_hmm_fit <- function(object, ...) {
  if (!isTRUE(object$switching)) {
    cli_abort(c(
      "{.fun hmm_state_params} is for state-specific HMM fits.",
      "i" = "This fit only switches the copula correlation; use {.fun hmm_states} and {.fun coef}."
    ))
  }
  K <- object$K
  D <- object$stan_data$D
  sw <- object$switch
  margins_char <- object$margins_matrix
  summ <- .fit_summary(object$fit, backend = object$backend)
  pm <- function(name) {
    row <- match(name, summ$variable)
    if (is.na(row)) NA_real_ else summ$mean[row]
  }

  # Per-state intercepts (tiled when mu is global).
  mu <- matrix(NA_real_, nrow = K, ncol = D)
  for (k in seq_len(K)) {
    mk <- if (isTRUE(sw$mu == 1L)) k else 1L
    for (d in seq_len(D)) mu[k, d] <- pm(sprintf("mu[%d,%d]", mk, d))
  }

  # Effective per-state Phi = baseline + masked per-state deviation.
  phi_base <- matrix(NA_real_, nrow = D, ncol = D)
  for (i in seq_len(D)) for (j in seq_len(D)) phi_base[i, j] <- pm(sprintf("Phi_base[%d,%d]", i, j))
  ij <- list(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(2L, 2L))
  active <- which(sw$phi_mask > 0L)
  Phi <- lapply(seq_len(K), function(k) {
    P <- phi_base
    for (a in active) {
      i <- ij[[a]][1]
      j <- ij[[a]][2]
      P[i, j] <- P[i, j] + pm(sprintf("Phi_dev[%d,%d,%d]", k, i, j))
    }
    P
  })

  rho_state <- vapply(seq_len(K), function(k) pm(sprintf("rho_state[%d]", k)), numeric(1))

  # Per-state margin scale/shape posterior means. Shared margin configs are
  # tiled across states so the returned shape is flag-agnostic.
  Mrg_K <- if (isTRUE(sw$margins == 1L)) K else 1L
  scale_rows <- list()
  row_idx <- 1L
  for (k in seq_len(K)) {
    mk <- if (Mrg_K == K) k else 1L
    for (d in seq_len(D)) {
      fam <- margins_char[k, d]
      vars <- switch(fam,
        normal = c(sigma_eps = sprintf("sigma_eps[%d,%d]", mk, d)),
        exponential = c(sigma_exp = sprintf("sigma_exp[%d,%d]", mk, d)),
        skew_normal = c(
          omega = sprintf("omega[%d,%d]", mk, d),
          delta = sprintf("delta[%d,%d]", mk, d)
        ),
        gamma = c(
          sigma_gam = sprintf("sigma_gam[%d,%d]", mk, d),
          shape_gam = sprintf("shape_gam[%d,%d]", mk, d)
        )
      )
      for (param in names(vars)) {
        scale_rows[[row_idx]] <- data.frame(
          state = k,
          variable = object$vars[d],
          dimension = d,
          family = fam,
          parameter = param,
          stan_variable = unname(vars[[param]]),
          mean = pm(unname(vars[[param]])),
          stringsAsFactors = FALSE
        )
        row_idx <- row_idx + 1L
      }
    }
  }
  scales <- do.call(rbind, scale_rows)
  rownames(scales) <- NULL

  list(
    families = margins_char,
    rho_state = rho_state,
    mu = mu,
    Phi = Phi,
    scales = scales
  )
}


#' Extract random effects from a multilevel fit
#'
#' Returns posterior summaries for unit-specific VAR coefficients.
#'
#' @param object A `dcvar_multilevel_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `unit`, `parameter`, `mean`, `sd`,
#'   `q2.5`, `q97.5`.
#' @export
random_effects <- function(object, ...) {
  UseMethod("random_effects")
}

#' @rdname random_effects
#' @export
random_effects.default <- function(object, ...) {
  cli_abort("{.fun random_effects} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' @rdname random_effects
#' @export
random_effects.dcvar_multilevel_fit <- function(object, ...) {
  N <- object$N
  unit_ids <- attr(object$stan_data, "ids")
  if (is.null(unit_ids)) {
    unit_ids <- seq_len(N)
  }
  summ <- .fit_summary(
    object$fit, variables = NULL, backend = object$backend,
    required = paste0(
      "phi_unit[",
      rep(seq_len(N), each = 4),
      ",",
      rep(seq_len(4), times = N),
      "]"
    ),
    required_type = "exact",
    context = "random_effects.dcvar_multilevel_fit()",
    output_type = "parameter",
    mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )

  param_names <- c("phi11", "phi12", "phi21", "phi22")
  results <- vector("list", N * 4)
  idx <- 1

  for (i in seq_len(N)) {
    for (k in 1:4) {
      var_name <- paste0("phi_unit[", i, ",", k, "]")
      row <- which(summ$variable == var_name)
      if (length(row) != 1L) {
        cli_abort(c(
          "{.fun random_effects} requires unit-specific VAR coefficients that are not present in the fitted model.",
          "i" = "Missing parameter: {.val {var_name}}.",
          "i" = "Custom Stan files must preserve the expected parameter names."
        ))
      }
      results[[idx]] <- data.frame(
        unit = unit_ids[i],
        parameter = param_names[k],
        mean = summ$mean[row],
        sd = summ$sd[row],
        q2.5 = summ$q2.5[row],
        q97.5 = summ$q97.5[row]
      )
      idx <- idx + 1
    }
  }

  do.call(rbind, results)
}


#' @rdname var_params
#' @details
#' For multilevel models, returns population-level parameters `phi_bar`
#' (mean VAR coefficients), `tau_phi` (between-unit SDs), `sigma`
#' (innovation SDs), and `rho` (copula correlation). These correspond to
#' `Phi`, `sigma_eps`, and `rho` in single-level models.
#' @export
var_params.dcvar_multilevel_fit <- function(object, ...) {
  margins <- object$margins %||% "normal"
  mixed <- .uses_sem_multilevel_mixed_engine("multilevel", margins)
  engine_margins <- .sem_multilevel_engine_margins("multilevel", margins, 2L)
  scale_pattern <- if (mixed) {
    paste0("^", names(.mixed_margin_report_vars(engine_margins)), "\\[")
  } else if (identical(margins, "exponential")) {
    "^sigma_exp\\["
  } else {
    "^sigma\\["
  }
  summ <- .fit_summary(
    object$fit, variables = NULL, backend = object$backend,
    required = c("^phi_bar\\[", "^tau_phi\\[", scale_pattern, "^rho$"),
    required_type = "pattern",
    context = "var_params.dcvar_multilevel_fit()",
    output_type = "parameter group",
    mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )
  extract_param <- function(pattern) {
    rows <- grep(pattern, summ$variable)
    data.frame(
      variable = summ$variable[rows],
      mean = summ$mean[rows],
      sd = summ$sd[rows],
      q2.5 = summ$q2.5[rows],
      q97.5 = summ$q97.5[rows]
    )
  }

  scale_param <- if (mixed) {
    specs <- .mixed_margin_report_vars(engine_margins)
    lapply(specs, function(vars) {
      rows <- match(vars, summ$variable)
      rows <- rows[!is.na(rows)]
      data.frame(
        variable = summ$variable[rows],
        mean = summ$mean[rows],
        sd = summ$sd[rows],
        q2.5 = summ$q2.5[rows],
        q97.5 = summ$q97.5[rows]
      )
    })
  } else if (identical(margins, "exponential")) {
    list(sigma_exp = extract_param("^sigma_exp\\["))
  } else {
    list(sigma = extract_param("^sigma\\["))
  }

  c(list(
    phi_bar = extract_param("^phi_bar\\["),
    tau_phi = extract_param("^tau_phi\\[")
  ), scale_param, list(
    rho = extract_param("^rho$")
  ))
}


#' @rdname rho_trajectory
#' @export
rho_trajectory.dcvar_multilevel_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  .rho_trajectory_constant_impl(object, probs)
}

#' @rdname dependence_summary
#' @export
dependence_summary.dcvar_multilevel_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  .dependence_summary_constant_impl(object, probs)
}


#' Extract latent states from a SEM fit
#'
#' Returns posterior summaries for the estimated latent states at each
#' time point.
#'
#' @param object A `dcvar_sem_fit` object.
#' @param probs Numeric vector of quantile probabilities.
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `time`, `variable`, `mean`, `sd`,
#'   and quantile columns.
#' @export
latent_states <- function(object, ...) {
  UseMethod("latent_states")
}

#' @rdname latent_states
#' @export
latent_states.default <- function(object, ...) {
  cli_abort("{.fun latent_states} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' @rdname latent_states
#' @export
latent_states.dcvar_sem_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  if (identical(object$method %||% "indicator", "naive")) {
    cli_abort("{.fun latent_states} is not defined for naive SEM fits because no latent measurement model is estimated.")
  }
  if (!is.numeric(probs) || !all(probs >= 0 & probs <= 1)) {
    cli_abort("{.arg probs} must be numeric values in [0, 1].")
  }

  state_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "state", backend = object$backend,
    required = .stan_output_group_pattern("state"),
    required_type = "pattern",
    context = "latent_states.dcvar_sem_fit()",
    output_type = "transformed parameter group"
  ))
  n_time_obs <- object$stan_data$n_time
  vars <- object$vars
  time_values <- .observed_time_values(object$stan_data)

  results <- vector("list", 2)
  for (d in 1:2) {
    cols <- paste0("state[", seq_len(n_time_obs), ",", d, "]")
    draws_d <- state_draws[, cols, drop = FALSE]

    df <- data.frame(
      time = time_values,
      variable = vars[d],
      mean = colMeans(draws_d),
      sd = apply(draws_d, 2, sd)
    )

    quants <- apply(draws_d, 2, quantile, probs = probs)
    if (is.null(dim(quants))) {
      quants <- matrix(quants, nrow = 1)
    }
    for (i in seq_along(probs)) {
      df[[paste0("q", probs[i] * 100)]] <- quants[i, ]
    }

    results[[d]] <- df
  }

  do.call(rbind, results)
}


#' @rdname var_params
#' @export
var_params.dcvar_sem_fit <- function(object, ...) {
  margins <- object$margins %||% "normal"
  mixed <- .uses_sem_multilevel_mixed_engine("sem", margins)
  engine_margins <- .sem_multilevel_engine_margins("sem", margins, 2L)
  required_patterns <- c("^mu\\[", "^Phi\\[", "^rho$")
  if (mixed) {
    required_patterns <- c(required_patterns,
                           paste0("^", names(.mixed_margin_report_vars(engine_margins)), "\\["))
  } else if (identical(margins, "normal")) {
    required_patterns <- c(required_patterns, "^sigma\\[")
  } else if (identical(margins, "exponential")) {
    required_patterns <- c(required_patterns, "^sigma_exp\\[")
  } else {
    cli_abort("Unsupported SEM margin type in {.fun var_params}: {.val {margins}}.")
  }

  summ <- .fit_summary(
    object$fit, variables = NULL, backend = object$backend,
    required = required_patterns,
    required_type = "pattern",
    context = "var_params.dcvar_sem_fit()",
    output_type = "parameter group",
    mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )
  extract_param <- function(pattern) {
    rows <- grep(pattern, summ$variable)
    data.frame(
      variable = summ$variable[rows],
      mean = summ$mean[rows],
      sd = summ$sd[rows],
      q2.5 = summ$q2.5[rows],
      q97.5 = summ$q97.5[rows]
    )
  }

  result <- list(
    mu = extract_param("^mu\\["),
    Phi = extract_param("^Phi\\["),
    rho = extract_param("^rho$")
  )
  scale_param <- if (mixed) {
    specs <- .mixed_margin_report_vars(engine_margins)
    lapply(specs, function(vars) {
      rows <- match(vars, summ$variable)
      rows <- rows[!is.na(rows)]
      data.frame(
        variable = summ$variable[rows],
        mean = summ$mean[rows],
        sd = summ$sd[rows],
        q2.5 = summ$q2.5[rows],
        q97.5 = summ$q97.5[rows]
      )
    })
  } else if (identical(margins, "normal")) {
    list(sigma = extract_param("^sigma\\["))
  } else {
    list(sigma_exp = extract_param("^sigma_exp\\["))
  }

  c(result[1:2], scale_param, result[3])
}


#' @rdname rho_trajectory
#' @export
rho_trajectory.dcvar_sem_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  # SEM has constant rho, expand to all n_time - 1 time points
  rho_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "rho", backend = object$backend,
    required = "rho",
    context = "rho_trajectory.dcvar_sem_fit()",
    output_type = "transformed parameter"
  ))
  n_time_obs <- object$stan_data$n_time
  n_time_eff <- n_time_obs - 1L

  rho_mean <- mean(rho_draws[, 1])
  rho_sd <- sd(rho_draws[, 1])
  quants <- quantile(rho_draws[, 1], probs = probs)
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)

  rho_summary <- data.frame(
    time = time_values,
    mean = rep(rho_mean, n_time_eff),
    sd = rep(rho_sd, n_time_eff)
  )

  for (i in seq_along(probs)) {
    rho_summary[[paste0("q", probs[i] * 100)]] <- rep(quants[i], n_time_eff)
  }

  rho_summary
}

#' @rdname dependence_summary
#' @export
dependence_summary.dcvar_sem_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  .dependence_summary_constant_impl(object, probs)
}


#' Extract posterior draws
#'
#' Extract posterior draws from a fitted model.
#'
#' @param object A fitted model object.
#' @param variable Character vector of parameter names. `NULL` returns all.
#' @param format Draw format: `"draws_array"`, `"draws_matrix"`, or
#'   `"draws_df"` (default: `"draws_array"`).
#' @param ... Additional arguments (unused).
#'
#' @return A posterior draws object.
#' @export
draws <- function(object, ...) {
  UseMethod("draws")
}

#' @rdname draws
#' @export
draws.default <- function(object, ...) {
  cli_abort("{.fun draws} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' @rdname draws
#' @export
draws.dcvar_model_fit <- function(object, variable = NULL, format = "draws_array", ...) {
  format <- match.arg(format, c("draws_array", "draws_matrix", "draws_df"))
  .fit_draws(object$fit, variables = variable, format = format, backend = object$backend)
}


# ============================================================================
# Time-varying trajectory extractors (dcvar_tv_fit)
# ============================================================================

#' Extract the time-varying VAR coefficient trajectories
#'
#' Returns posterior summaries of the four VAR(1) coefficient paths
#' phi11(t), phi12(t), phi21(t), phi22(t) from a time-varying DC-VAR fit.
#' For fits with `tv_phi = FALSE` the constant baseline coefficients are
#' tiled over time so the return shape does not depend on the flag.
#'
#' @param object A fitted model object.
#' @param probs Numeric vector of quantile probabilities
#'   (default: `c(0.025, 0.1, 0.5, 0.9, 0.975)`).
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `time`, `coefficient` (one of
#'   `"phi11"`, `"phi12"`, `"phi21"`, `"phi22"`), `mean`, `sd`, and one
#'   column per requested quantile.
#' @seealso [sigma_trajectory()], [rho_trajectory()], [plot_phi_trajectory()]
#' @export
phi_trajectory <- function(object, ...) {
  UseMethod("phi_trajectory")
}

#' @rdname phi_trajectory
#' @export
phi_trajectory.default <- function(object, ...) {
  cli_abort("{.fun phi_trajectory} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' Internal: coefficient labels in the row-major Stan column order
#' @noRd
.phi_tv_labels <- c("phi11", "phi12", "phi21", "phi22")

#' Internal: matching baseline Phi draw column for each phi_t column
#' @noRd
.phi_tv_baseline_cols <- c("Phi[1,1]", "Phi[1,2]", "Phi[2,1]", "Phi[2,2]")

#' @rdname phi_trajectory
#' @export
phi_trajectory.dcvar_tv_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  n_time_eff <- object$stan_data$n_time - 1L
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)

  if (isTRUE(object$tv_phi)) {
    phi_draws <- posterior::as_draws_matrix(.fit_draws(
      object$fit, "phi_t", backend = object$backend,
      required = .stan_output_group_pattern("phi_t"),
      required_type = "pattern",
      context = "phi_trajectory.dcvar_tv_fit()",
      output_type = "generated quantity"
    ))
    out <- lapply(1:4, function(k) {
      cols <- paste0("phi_t[", seq_len(n_time_eff), ",", k, "]")
      df <- .summarise_rho_draws(as.matrix(phi_draws[, cols, drop = FALSE]), probs, time_values)
      df$coefficient <- .phi_tv_labels[k]
      df
    })
  } else {
    # Constant Phi: tile the baseline so the contract is flag-agnostic.
    base_draws <- posterior::as_draws_matrix(.fit_draws(
      object$fit, "Phi", backend = object$backend,
      required = .stan_output_group_pattern("Phi"),
      required_type = "pattern",
      context = "phi_trajectory.dcvar_tv_fit()",
      output_type = "parameter group"
    ))
    out <- lapply(1:4, function(k) {
      one <- .summarise_rho_draws(
        as.matrix(base_draws[, .phi_tv_baseline_cols[k], drop = FALSE]),
        probs, time_values[1]
      )
      df <- one[rep(1L, n_time_eff), , drop = FALSE]
      df$time <- time_values
      rownames(df) <- NULL
      df$coefficient <- .phi_tv_labels[k]
      df
    })
  }

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[, c("time", "coefficient", setdiff(names(out), c("time", "coefficient")))]
}


#' Extract the time-varying residual scale trajectories
#'
#' Returns posterior summaries of the per-variable residual scale paths from
#' a time-varying DC-VAR fit. The reported value is each margin family's
#' natural scale: the innovation SD for normal dimensions, the direct
#' parameterization scale `omega` for skew-normal dimensions, and the
#' (time-constant) `sigma_exp` / `sigma_gam` for exponential and gamma
#' dimensions. For fits with `tv_sigma = FALSE` the constant baselines are
#' tiled over time so the return shape does not depend on the flag.
#'
#' @param object A fitted model object.
#' @param probs Numeric vector of quantile probabilities
#'   (default: `c(0.025, 0.1, 0.5, 0.9, 0.975)`).
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `time`, `variable`, `mean`, `sd`, and
#'   one column per requested quantile.
#' @seealso [phi_trajectory()], [rho_trajectory()], [plot_sigma_trajectory()]
#' @export
sigma_trajectory <- function(object, ...) {
  UseMethod("sigma_trajectory")
}

#' @rdname sigma_trajectory
#' @export
sigma_trajectory.default <- function(object, ...) {
  cli_abort("{.fun sigma_trajectory} is not defined for objects of class {.cls {class(object)[[1]]}}.")
}

#' Internal: baseline scale draw column for one dimension of a TV fit
#' @noRd
.sigma_tv_baseline_col <- function(family, d) {
  switch(family,
    normal = paste0("sigma_eps[", d, "]"),
    exponential = paste0("sigma_exp[", d, "]"),
    skew_normal = paste0("omega[", d, "]"),
    gamma = paste0("sigma_gam[", d, "]")
  )
}

#' @rdname sigma_trajectory
#' @export
sigma_trajectory.dcvar_tv_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  D <- object$stan_data$D
  n_time_eff <- object$stan_data$n_time - 1L
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)
  margins_vec <- rep(object$margins %||% "normal", length.out = D)

  if (isTRUE(object$tv_sigma)) {
    sigma_draws <- posterior::as_draws_matrix(.fit_draws(
      object$fit, "sigma_t", backend = object$backend,
      required = .stan_output_group_pattern("sigma_t"),
      required_type = "pattern",
      context = "sigma_trajectory.dcvar_tv_fit()",
      output_type = "generated quantity"
    ))
    out <- lapply(seq_len(D), function(d) {
      cols <- paste0("sigma_t[", seq_len(n_time_eff), ",", d, "]")
      df <- .summarise_rho_draws(as.matrix(sigma_draws[, cols, drop = FALSE]), probs, time_values)
      df$variable <- object$vars[d]
      df
    })
  } else {
    out <- lapply(seq_len(D), function(d) {
      col <- .sigma_tv_baseline_col(margins_vec[d], d)
      base_draws <- posterior::as_draws_matrix(.fit_draws(
        object$fit, sub("\\[.*$", "", col), backend = object$backend,
        required = paste0("^", sub("\\[.*$", "\\\\[", col)),
        required_type = "pattern",
        context = "sigma_trajectory.dcvar_tv_fit()",
        output_type = "Stan output group"
      ))
      one <- .summarise_rho_draws(
        as.matrix(base_draws[, col, drop = FALSE]),
        probs, time_values[1]
      )
      df <- one[rep(1L, n_time_eff), , drop = FALSE]
      df$time <- time_values
      rownames(df) <- NULL
      df$variable <- object$vars[d]
      df
    })
  }

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[, c("time", "variable", setdiff(names(out), c("time", "variable")))]
}


#' @rdname var_params
#' @export
var_params.dcvar_tv_fit <- function(object, ...) {
  D <- object$stan_data$D
  margins_vec <- rep(object$margins %||% "normal", length.out = D)
  margin_groups <- names(.mixed_margin_report_vars(margins_vec))

  required_patterns <- c("^mu\\[", "^Phi\\[", "^sigma_omega$",
                         paste0("^", margin_groups, "\\["))
  if (isTRUE(object$tv_phi)) {
    required_patterns <- c(required_patterns, "^tau_phi\\[")
  }
  if (isTRUE(object$tv_sigma)) {
    required_patterns <- c(required_patterns, "^tau_sigma\\[")
  }

  summ <- .fit_summary(
    object$fit, variables = NULL, backend = object$backend,
    required = required_patterns,
    required_type = "pattern",
    context = "var_params.dcvar_tv_fit()",
    output_type = "parameter group",
    mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )
  extract_param <- function(pattern) {
    rows <- grep(pattern, summ$variable)
    data.frame(
      variable = summ$variable[rows],
      mean = summ$mean[rows],
      sd = summ$sd[rows],
      q2.5 = summ$q2.5[rows],
      q97.5 = summ$q97.5[rows]
    )
  }

  result <- list(
    mu = extract_param("^mu\\["),
    Phi = extract_param("^Phi\\[")
  )
  for (group in margin_groups) {
    result[[group]] <- extract_param(paste0("^", group, "\\["))
  }
  if (isTRUE(object$tv_phi)) {
    tau_phi <- extract_param("^tau_phi\\[")
    # Stan indexes tau_phi by active position; relabel by coefficient name.
    tau_phi$variable <- paste0("tau_phi[", .tv_active_phi_coefs(object), "]")
    result$tau_phi <- tau_phi
  }
  if (isTRUE(object$tv_sigma)) {
    result$tau_sigma <- extract_param("^tau_sigma\\[")
  }
  result$sigma_omega <- extract_param("^sigma_omega$")
  result
}


# ============================================================================
# Time-varying SEM and multilevel extractors
# ============================================================================

#' @rdname rho_trajectory
#' @export
rho_trajectory.dcvar_sem_tv_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  rho_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "rho", backend = object$backend,
    required = .stan_output_group_pattern("rho"),
    required_type = "pattern",
    context = "rho_trajectory.dcvar_sem_tv_fit()",
    output_type = "transformed parameter group"
  ))
  .summarise_rho_draws(rho_draws, probs, .observed_time_values(object$stan_data, drop_first = TRUE))
}

#' @rdname dependence_summary
#' @export
dependence_summary.dcvar_sem_tv_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  rho_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "rho", backend = object$backend,
    required = .stan_output_group_pattern("rho"),
    required_type = "pattern",
    context = "dependence_summary.dcvar_sem_tv_fit()",
    output_type = "transformed parameter group"
  ))
  .summarise_rho_draws(2 / pi * asin(rho_draws), probs,
                       .observed_time_values(object$stan_data, drop_first = TRUE))
}

#' @rdname rho_trajectory
#' @export
rho_trajectory.dcvar_multilevel_tv_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  .rho_trajectory_constant_impl(object, probs)
}

#' @rdname dependence_summary
#' @export
dependence_summary.dcvar_multilevel_tv_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  .dependence_summary_constant_impl(object, probs)
}

.phi_trajectory_tv_matrix <- function(object, probs, context, draw_name = "phi_t") {
  n_time_eff <- object$stan_data$n_time - 1L
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)
  phi_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, draw_name, backend = object$backend,
    required = .stan_output_group_pattern(draw_name),
    required_type = "pattern",
    context = context,
    output_type = "generated quantity"
  ))
  out <- lapply(1:4, function(k) {
    cols <- paste0(draw_name, "[", seq_len(n_time_eff), ",", k, "]")
    df <- .summarise_rho_draws(as.matrix(phi_draws[, cols, drop = FALSE]), probs, time_values)
    df$coefficient <- .phi_tv_labels[k]
    df
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[, c("time", "coefficient", setdiff(names(out), c("time", "coefficient")))]
}

.phi_trajectory_tile_baseline <- function(draws, cols, labels, probs, time_values) {
  n_time_eff <- length(time_values)
  out <- lapply(seq_along(cols), function(k) {
    one <- .summarise_rho_draws(
      as.matrix(draws[, cols[k], drop = FALSE]),
      probs,
      time_values[1]
    )
    df <- one[rep(1L, n_time_eff), , drop = FALSE]
    df$time <- time_values
    rownames(df) <- NULL
    df$coefficient <- labels[k]
    df
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[, c("time", "coefficient", setdiff(names(out), c("time", "coefficient")))]
}

#' @rdname phi_trajectory
#' @export
phi_trajectory.dcvar_sem_tv_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  if (isTRUE(object$tv_phi)) {
    return(.phi_trajectory_tv_matrix(object, probs, "phi_trajectory.dcvar_sem_tv_fit()"))
  }
  base_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "Phi", backend = object$backend,
    required = .stan_output_group_pattern("Phi"),
    required_type = "pattern",
    context = "phi_trajectory.dcvar_sem_tv_fit()",
    output_type = "parameter group"
  ))
  .phi_trajectory_tile_baseline(
    base_draws, .phi_tv_baseline_cols, .phi_tv_labels, probs,
    .observed_time_values(object$stan_data, drop_first = TRUE)
  )
}

#' @rdname phi_trajectory
#' @export
phi_trajectory.dcvar_multilevel_tv_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975),
                                                   unit = NULL, ...) {
  n_time_eff <- object$stan_data$n_time - 1L
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)
  unit_ids <- attr(object$stan_data, "ids") %||% seq_len(object$N)

  if (is.null(unit)) {
    if (isTRUE(object$tv_phi)) {
      return(.phi_trajectory_tv_matrix(object, probs, "phi_trajectory.dcvar_multilevel_tv_fit()"))
    }
    base_draws <- posterior::as_draws_matrix(.fit_draws(
      object$fit, "phi_bar", backend = object$backend,
      required = .stan_output_group_pattern("phi_bar"),
      required_type = "pattern",
      context = "phi_trajectory.dcvar_multilevel_tv_fit()",
      output_type = "parameter group"
    ))
    return(.phi_trajectory_tile_baseline(
      base_draws, paste0("phi_bar[", 1:4, "]"), .phi_tv_labels, probs, time_values
    ))
  }

  selected_units <- if (identical(unit, "all")) unit_ids else unit
  selected_idx <- match(as.character(selected_units), as.character(unit_ids))
  if (anyNA(selected_idx)) {
    cli_abort("Unknown {.arg unit} value{?s}: {.val {selected_units[is.na(selected_idx)]}}.")
  }

  if (isTRUE(object$tv_phi)) {
    unit_draws <- posterior::as_draws_matrix(.fit_draws(
      object$fit, "phi_unit_t", backend = object$backend,
      required = .stan_output_group_pattern("phi_unit_t"),
      required_type = "pattern",
      context = "phi_trajectory.dcvar_multilevel_tv_fit()",
      output_type = "generated quantity"
    ))
    out <- lapply(seq_along(selected_idx), function(ii) {
      i <- selected_idx[ii]
      rows <- lapply(1:4, function(k) {
        cols <- paste0("phi_unit_t[", i, ",", seq_len(n_time_eff), ",", k, "]")
        df <- .summarise_rho_draws(as.matrix(unit_draws[, cols, drop = FALSE]), probs, time_values)
        df$unit <- unit_ids[i]
        df$coefficient <- .phi_tv_labels[k]
        df
      })
      do.call(rbind, rows)
    })
  } else {
    base_draws <- posterior::as_draws_matrix(.fit_draws(
      object$fit, "phi_unit", backend = object$backend,
      required = .stan_output_group_pattern("phi_unit"),
      required_type = "pattern",
      context = "phi_trajectory.dcvar_multilevel_tv_fit()",
      output_type = "transformed parameter group"
    ))
    out <- lapply(seq_along(selected_idx), function(ii) {
      i <- selected_idx[ii]
      df <- .phi_trajectory_tile_baseline(
        base_draws, paste0("phi_unit[", i, ",", 1:4, "]"),
        .phi_tv_labels, probs, time_values
      )
      df$unit <- unit_ids[i]
      df
    })
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[, c("unit", "time", "coefficient", setdiff(names(out), c("unit", "time", "coefficient")))]
}

.sigma_trajectory_tv_generic <- function(object, probs, context) {
  D <- 2L
  n_time_eff <- object$stan_data$n_time - 1L
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)
  margins_vec <- rep(object$margins %||% "normal", length.out = D)

  if (isTRUE(object$tv_sigma)) {
    sigma_draws <- posterior::as_draws_matrix(.fit_draws(
      object$fit, "sigma_t", backend = object$backend,
      required = .stan_output_group_pattern("sigma_t"),
      required_type = "pattern",
      context = context,
      output_type = "generated quantity"
    ))
    out <- lapply(seq_len(D), function(d) {
      cols <- paste0("sigma_t[", seq_len(n_time_eff), ",", d, "]")
      df <- .summarise_rho_draws(as.matrix(sigma_draws[, cols, drop = FALSE]), probs, time_values)
      df$variable <- object$vars[d]
      df
    })
  } else {
    out <- lapply(seq_len(D), function(d) {
      col <- .sigma_tv_baseline_col(margins_vec[d], d)
      base_draws <- posterior::as_draws_matrix(.fit_draws(
        object$fit, sub("\\[.*$", "", col), backend = object$backend,
        required = paste0("^", sub("\\[.*$", "\\\\[", col)),
        required_type = "pattern",
        context = context,
        output_type = "Stan output group"
      ))
      one <- .summarise_rho_draws(
        as.matrix(base_draws[, col, drop = FALSE]),
        probs,
        time_values[1]
      )
      df <- one[rep(1L, n_time_eff), , drop = FALSE]
      df$time <- time_values
      rownames(df) <- NULL
      df$variable <- object$vars[d]
      df
    })
  }

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[, c("time", "variable", setdiff(names(out), c("time", "variable")))]
}

#' @rdname sigma_trajectory
#' @export
sigma_trajectory.dcvar_sem_tv_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  .sigma_trajectory_tv_generic(object, probs, "sigma_trajectory.dcvar_sem_tv_fit()")
}

#' @rdname sigma_trajectory
#' @export
sigma_trajectory.dcvar_multilevel_tv_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  .sigma_trajectory_tv_generic(object, probs, "sigma_trajectory.dcvar_multilevel_tv_fit()")
}

.tv_margin_groups <- function(object, D = 2L) {
  names(.mixed_margin_report_vars(rep(object$margins %||% "normal", length.out = D)))
}

#' @rdname var_params
#' @export
var_params.dcvar_sem_tv_fit <- function(object, ...) {
  margin_groups <- .tv_margin_groups(object, 2L)
  required_patterns <- c("^mu\\[", "^Phi\\[", "^rho\\[", "^sigma_omega$",
                         paste0("^", margin_groups, "\\["))
  if (isTRUE(object$tv_phi)) required_patterns <- c(required_patterns, "^tau_phi\\[")
  if (isTRUE(object$tv_sigma)) required_patterns <- c(required_patterns, "^tau_sigma\\[")

  summ <- .fit_summary(
    object$fit, variables = NULL, backend = object$backend,
    required = required_patterns,
    required_type = "pattern",
    context = "var_params.dcvar_sem_tv_fit()",
    output_type = "parameter group",
    mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )
  extract_param <- function(pattern) {
    rows <- grep(pattern, summ$variable)
    data.frame(
      variable = summ$variable[rows],
      mean = summ$mean[rows],
      sd = summ$sd[rows],
      q2.5 = summ$q2.5[rows],
      q97.5 = summ$q97.5[rows]
    )
  }

  result <- list(mu = extract_param("^mu\\["), Phi = extract_param("^Phi\\["))
  for (group in margin_groups) result[[group]] <- extract_param(paste0("^", group, "\\["))
  result$sigma_omega <- extract_param("^sigma_omega$")
  if (isTRUE(object$tv_phi)) {
    tau_phi <- extract_param("^tau_phi\\[")
    tau_phi$variable <- paste0("tau_phi[", .tv_active_phi_coefs(object), "]")
    result$tau_phi <- tau_phi
  }
  if (isTRUE(object$tv_sigma)) result$tau_sigma <- extract_param("^tau_sigma\\[")
  result$rho <- extract_param("^rho\\[")
  result
}

#' @rdname var_params
#' @export
var_params.dcvar_multilevel_tv_fit <- function(object, ...) {
  margin_groups <- .tv_margin_groups(object, 2L)
  required_patterns <- c("^phi_bar\\[", "^tau_phi\\[", "^rho$",
                         paste0("^", margin_groups, "\\["))
  if (isTRUE(object$tv_phi)) required_patterns <- c(required_patterns, "^tau_phi_tv\\[")
  if (isTRUE(object$tv_sigma)) required_patterns <- c(required_patterns, "^tau_sigma\\[")

  summ <- .fit_summary(
    object$fit, variables = NULL, backend = object$backend,
    required = required_patterns,
    required_type = "pattern",
    context = "var_params.dcvar_multilevel_tv_fit()",
    output_type = "parameter group",
    mean, sd,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )
  extract_param <- function(pattern) {
    rows <- grep(pattern, summ$variable)
    data.frame(
      variable = summ$variable[rows],
      mean = summ$mean[rows],
      sd = summ$sd[rows],
      q2.5 = summ$q2.5[rows],
      q97.5 = summ$q97.5[rows]
    )
  }

  result <- list(
    phi_bar = extract_param("^phi_bar\\["),
    tau_phi = extract_param("^tau_phi\\[")
  )
  for (group in margin_groups) result[[group]] <- extract_param(paste0("^", group, "\\["))
  result$rho <- extract_param("^rho$")
  if (isTRUE(object$tv_phi)) {
    tau_phi_tv <- extract_param("^tau_phi_tv\\[")
    tau_phi_tv$variable <- paste0("tau_phi_tv[", .tv_active_phi_coefs(object), "]")
    result$tau_phi_tv <- tau_phi_tv
  }
  if (isTRUE(object$tv_sigma)) result$tau_sigma <- extract_param("^tau_sigma\\[")
  result
}
