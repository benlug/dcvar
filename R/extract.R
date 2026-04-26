# ============================================================================
# Extraction Functions
# ============================================================================

#' Extract the rho trajectory with credible intervals
#'
#' Returns a data frame with the posterior mean, SD, and quantiles of the
#' time-varying correlation at each time point.
#'
#' @param object A fitted model object (`dcvar_fit`, `dcvar_covariate_fit`,
#'   `dcvar_hmm_fit`, or `dcvar_constant_fit`).
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

#' @rdname rho_trajectory
#' @export
rho_trajectory.dcvar_constant_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  .rho_trajectory_constant_impl(object, probs)
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
  required_patterns <- c("^mu\\[", "^Phi\\[")
  switch(margins,
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
  has_sigma_omega <- identical(object$model, "dcvar") ||
    identical(object$model, "dcvar_covariate")
  if (has_sigma_omega) {
    required_patterns <- c(required_patterns, "^sigma_omega$")
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
  switch(margins,
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

  # sigma_omega is present in random-walk DC-VAR variants only.
  so <- if (has_sigma_omega) extract_param("^sigma_omega$") else NULL
  if (!is.null(so)) result$sigma_omega <- so

  result
}


#' Extract covariate effect summaries
#'
#' Returns posterior summaries for the Fisher-z intercept, covariate effects,
#' and, when present, the residual random-walk innovation scale.
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

  required_patterns <- c("^beta_0$", "^beta\\[")
  if (isTRUE(object$drift)) {
    required_patterns <- c(required_patterns, "^sigma_omega$")
  }

  summ <- .fit_summary(
    object$fit,
    variables = NULL,
    backend = object$backend,
    required = required_patterns,
    required_type = "pattern",
    context = "covariate_effects.dcvar_covariate_fit()",
    output_type = "parameter group",
    mean,
    sd,
    ~posterior::quantile2(.x, probs = probs)
  )

  rows <- c(
    grep("^beta_0$", summ$variable),
    grep("^beta\\[", summ$variable),
    grep("^sigma_omega$", summ$variable)
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
#'   - `viterbi`: integer vector of MAP state sequence
#'   - `rho_state`: list with `mean`, `lower`, `upper` for each state
#'   - `A`: K x K posterior mean transition matrix
#'   - `rho_hmm`: posterior-averaged rho trajectory
#' @export
hmm_states <- function(object, ...) {
  UseMethod("hmm_states")
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

  # Use the most frequent complete Viterbi path across draws so the returned
  # sequence is always a valid joint path.
  viterbi_draws <- .safe_draws("viterbi_state")
  viterbi_paths <- apply(viterbi_draws, 1, paste, collapse = ",")
  viterbi_mode <- as.integer(strsplit(names(sort(table(viterbi_paths), decreasing = TRUE))[1], ",", fixed = TRUE)[[1]])

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

  # Posterior-averaged rho
  rho_hmm_draws <- .safe_draws("rho_hmm")
  rho_hmm_mean <- colMeans(rho_hmm_draws)

  list(
    gamma = gamma_mean,
    viterbi = viterbi_mode,
    rho_state = list(
      mean = rho_state_mean,
      lower = rho_state_lower,
      upper = rho_state_upper
    ),
    A = A_mean,
    rho_hmm = rho_hmm_mean
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
  summ <- .fit_summary(
    object$fit, variables = NULL, backend = object$backend,
    required = c("^phi_bar\\[", "^tau_phi\\[", "^sigma\\[", "^rho$"),
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

  list(
    phi_bar = extract_param("^phi_bar\\["),
    tau_phi = extract_param("^tau_phi\\["),
    sigma = extract_param("^sigma\\["),
    rho = extract_param("^rho$")
  )
}


#' @rdname rho_trajectory
#' @export
rho_trajectory.dcvar_multilevel_fit <- function(object, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  .rho_trajectory_constant_impl(object, probs)
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
  required_patterns <- c("^mu\\[", "^Phi\\[", "^rho$")
  if (identical(margins, "normal")) {
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
  if (identical(margins, "normal")) {
    result <- c(result[1:2], list(sigma = extract_param("^sigma\\[")), result[3])
  } else {
    result <- c(result[1:2], list(sigma_exp = extract_param("^sigma_exp\\[")), result[3])
  }

  result
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
  .fit_draws(object$fit, variables = variable, format = format, backend = object$backend)
}
