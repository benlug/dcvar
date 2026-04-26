# ============================================================================
# predict() and fitted() methods
# ============================================================================

#' Fitted values from a copula VAR model
#'
#' Returns the one-step-ahead fitted values (posterior mean of y_hat) from the
#' VAR(1) component: `y_hat[t] = mu + Phi * (y[t-1] - mu)`.
#'
#' If the model was fit with `standardize = TRUE` (the default), fitted values
#' are on the standardized (z-scored) scale by default. Use `type = "response"`
#' to back-transform to the original data scale.
#'
#' `fitted()` and `predict()` are implemented for the three core single-level
#' models plus the multilevel and SEM variants. For multilevel fits, the
#' methods return unit-specific trajectories. For SEM fits, `type = "link"`
#' returns latent-state summaries and `type = "response"` returns observed
#' indicator-scale summaries.
#'
#' @param object A fitted model object.
#' @param type Character; `"link"` (default) returns values on the model's
#'   internal scale (standardized if applicable), `"response"` back-transforms
#'   to the original data scale.
#' @param ... Additional arguments (unused).
#'
#' @return A data frame of posterior-mean fitted values. Single-level fits
#'   return columns `time` plus one column per modeled variable. Multilevel
#'   fits additionally include `unit`. SEM fits return either latent-state
#'   columns (`type = "link"`) or observed-indicator columns
#'   (`type = "response"`).
#' @export
fitted.dcvar_model_fit <- function(object, type = c("link", "response"), ...) {
  type <- match.arg(type)
  D <- object$stan_data$D
  n_time_obs <- object$stan_data$n_time
  Y <- object$stan_data$Y
  mu_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "mu", backend = object$backend,
    required = .stan_output_group_pattern("mu"),
    required_type = "pattern",
    context = "fitted.dcvar_model_fit()",
    output_type = "parameter group"
  ))
  Phi_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "Phi", backend = object$backend,
    required = .stan_output_group_pattern("Phi"),
    required_type = "pattern",
    context = "fitted.dcvar_model_fit()",
    output_type = "parameter group"
  ))

  mu_cols <- paste0("mu[", seq_len(D), "]")
  phi_cols <- unlist(lapply(seq_len(D), function(i) paste0("Phi[", i, ",", seq_len(D), "]")))
  mu_mat <- mu_draws[, mu_cols, drop = FALSE]
  phi_mat <- Phi_draws[, phi_cols, drop = FALSE]

  y_hat <- matrix(NA_real_, n_time_obs - 1L, D)
  for (time_index in 2:n_time_obs) {
    centered_prev <- sweep(mu_mat, 2, Y[time_index - 1L, ], FUN = function(mu_val, y_prev) y_prev - mu_val)
    y_hat_draws <- matrix(NA_real_, nrow(mu_mat), D)

    for (d in seq_len(D)) {
      eta <- mu_mat[, d]
      for (j in seq_len(D)) {
        eta <- eta + phi_mat[, (d - 1) * D + j] * centered_prev[, j]
      }
      y_hat_draws[, d] <- eta
    }

    y_hat[time_index - 1L, ] <- colMeans(y_hat_draws)
  }

  if (type == "response" && isTRUE(object$standardized)) {
    Y_means <- attr(object$stan_data, "Y_means")
    Y_sds <- attr(object$stan_data, "Y_sds")
    if (!is.null(Y_means) && !is.null(Y_sds)) {
      for (d in seq_len(D)) {
        y_hat[, d] <- y_hat[, d] * Y_sds[d] + Y_means[d]
      }
    }
  }

  out <- data.frame(time = .observed_time_values(object$stan_data, drop_first = TRUE), y_hat)
  names(out) <- c("time", object$vars)
  out
}


#' Internal: extract posterior mean unit-specific Phi coefficients
#' @noRd
.multilevel_phi_means <- function(object) {
  re <- random_effects(object)
  unit_ids <- attr(object$stan_data, "ids")
  if (is.null(unit_ids)) {
    unit_ids <- seq_len(object$N)
  }
  unit_ids <- as.character(unit_ids)

  param_names <- c("phi11", "phi12", "phi21", "phi22")
  phi_means <- matrix(
    NA_real_,
    nrow = object$N,
    ncol = length(param_names),
    dimnames = list(unit_ids, param_names)
  )

  for (i in seq_len(object$N)) {
    for (p in param_names) {
      row <- which(as.character(re$unit) == unit_ids[i] & re$parameter == p)
      if (length(row) != 1L) {
        cli_abort(c(
          "{.fun fitted} and {.fun predict} require unit-specific VAR coefficients that are not present in the fitted model.",
          "i" = "Missing random effect: {.val {p}} for unit {.val {unit_ids[i]}}.",
          "i" = "Custom Stan files must preserve the expected parameter names."
        ))
      }
      phi_means[i, p] <- re$mean[row]
    }
  }

  phi_means
}


#' Internal: compute multilevel one-step-ahead fitted values
#' @noRd
.multilevel_fitted_values <- function(object, response = FALSE) {
  phi_means <- .multilevel_phi_means(object)
  y_list <- object$stan_data$y
  unit_ids <- attr(object$stan_data, "ids")
  if (is.null(unit_ids)) {
    unit_ids <- seq_len(object$N)
  }
  unit_ids <- as.character(unit_ids)

  person_means <- object$person_means %||% attr(object$stan_data, "person_means")
  response_shift <- isTRUE(object$centered) &&
    !is.null(person_means) &&
    all(is.finite(person_means))

  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)
  out <- vector("list", object$N)

  for (i in seq_len(object$N)) {
    phi <- phi_means[i, ]
    Phi <- matrix(
      c(phi[["phi11"]], phi[["phi12"]], phi[["phi21"]], phi[["phi22"]]),
      nrow = 2,
      byrow = FALSE
    )

    y_i <- y_list[[i]]
    T_eff <- nrow(y_i) - 1L
    fitted_mat <- matrix(NA_real_, nrow = T_eff, ncol = 2)
    for (t in 2:nrow(y_i)) {
      fitted_mat[t - 1L, ] <- as.numeric(y_i[t - 1L, , drop = FALSE] %*% Phi)
    }

    if (response && response_shift) {
      fitted_mat <- sweep(fitted_mat, 2, person_means[i, ], "+")
    }

    df <- data.frame(
      unit = rep(unit_ids[i], T_eff),
      time = time_values,
      fitted_mat,
      check.names = FALSE
    )
    names(df)[3:4] <- object$vars
    out[[i]] <- df
  }

  do.call(rbind, out)
}


#' @rdname fitted.dcvar_model_fit
#' @export
fitted.dcvar_multilevel_fit <- function(object, type = c("link", "response"), ...) {
  type <- match.arg(type)
  .multilevel_fitted_values(object, response = identical(type, "response"))
}


#' Internal: build SEM latent state summaries as a list
#' @noRd
.sem_state_summaries <- function(object) {
  states <- latent_states(object)
  time_values <- .observed_time_values(object$stan_data)
  latent_names <- object$vars

  if (!all(latent_names %in% states$variable)) {
    cli_abort("Latent state summaries are missing expected variable names.")
  }

  out <- vector("list", length(latent_names))
  names(out) <- latent_names
  for (i in seq_along(latent_names)) {
    df <- states[states$variable == latent_names[i], , drop = FALSE]
    df <- df[match(time_values, df$time), , drop = FALSE]
    if (nrow(df) != length(time_values) || anyNA(df$time)) {
      cli_abort(c(
        "{.fun fitted} and {.fun predict} require latent state summaries that are not present in the fitted model.",
        "i" = "Missing latent-state times for {.val {latent_names[i]}}.",
        "i" = "Custom Stan files must preserve the expected parameter and generated-quantity names."
      ))
    }
    out[[i]] <- df
  }

  out
}


#' Internal: observed indicator means from SEM latent states
#' @noRd
.sem_indicator_means <- function(object, state_summaries) {
  indicator_groups <- object$indicators
  latent_names <- names(indicator_groups)
  if (is.null(latent_names) || any(latent_names == "")) {
    latent_names <- object$vars
  }

  time_values <- .observed_time_values(object$stan_data)
  out <- data.frame(time = time_values)
  for (i in seq_along(latent_names)) {
    latent_df <- state_summaries[[latent_names[i]]]
    latent_mean <- latent_df$mean
    for (j in seq_len(object$J)) {
      out[[indicator_groups[[i]][j]]] <- object$lambda[j] * latent_mean
    }
  }

  out
}


#' Internal: SEM indicator prediction intervals
#' @noRd
.sem_indicator_predictions <- function(object, state_summaries, ci_level) {
  indicator_groups <- object$indicators
  latent_names <- names(indicator_groups)
  if (is.null(latent_names) || any(latent_names == "")) {
    latent_names <- object$vars
  }

  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
  rows <- vector("list", length(latent_names) * object$J)
  idx <- 1L

  for (i in seq_along(latent_names)) {
    latent_df <- state_summaries[[latent_names[i]]]
    for (j in seq_len(object$J)) {
      mean <- object$lambda[j] * latent_df$mean
      sd <- sqrt((object$lambda[j] * latent_df$sd)^2 + object$sigma_e^2)
      rows[[idx]] <- data.frame(
        time = latent_df$time,
        variable = indicator_groups[[i]][j],
        mean = mean,
        lower = mean - z_crit * sd,
        upper = mean + z_crit * sd
      )
      idx <- idx + 1L
    }
  }

  do.call(rbind, rows)
}


#' Internal: SEM latent-state predictions
#' @noRd
.sem_latent_predictions <- function(object, state_summaries, ci_level) {
  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
  rows <- vector("list", length(state_summaries))

  for (i in seq_along(state_summaries)) {
    latent_df <- state_summaries[[i]]
    rows[[i]] <- data.frame(
      time = latent_df$time,
      variable = latent_df$variable,
      mean = latent_df$mean,
      lower = latent_df$mean - z_crit * latent_df$sd,
      upper = latent_df$mean + z_crit * latent_df$sd
    )
  }

  do.call(rbind, rows)
}

#' Internal: fitted values for naive SEM score models.
#' @noRd
.sem_naive_fitted_values <- function(object) {
  y <- object$stan_data$y
  n_time_obs <- object$stan_data$n_time

  summ <- .fit_summary(
    object$fit,
    variables = NULL,
    backend = object$backend,
    required = c("^mu\\[", "^Phi\\["),
    required_type = "pattern",
    context = ".sem_naive_fitted_values()",
    output_type = "parameter group",
    mean
  )
  get_mean <- function(name) {
    row <- which(summ$variable == name)
    if (length(row) != 1L) {
      cli_abort("Naive SEM fitted values require Stan output {.val {name}}.")
    }
    summ$mean[row]
  }

  mu <- c(get_mean("mu[1]"), get_mean("mu[2]"))
  Phi <- matrix(
    c(
      get_mean("Phi[1,1]"), get_mean("Phi[2,1]"),
      get_mean("Phi[1,2]"), get_mean("Phi[2,2]")
    ),
    nrow = 2
  )
  fitted_mat <- matrix(NA_real_, nrow = n_time_obs - 1L, ncol = 2)
  for (t in 2:n_time_obs) {
    fitted_mat[t - 1L, ] <- as.numeric(mu + y[t - 1L, , drop = FALSE] %*% t(Phi))
  }

  out <- data.frame(
    time = .observed_time_values(object$stan_data, drop_first = TRUE),
    fitted_mat
  )
  names(out) <- c("time", object$vars)
  out
}


#' @rdname fitted.dcvar_model_fit
#' @export
fitted.dcvar_sem_fit <- function(object, type = c("link", "response"), ...) {
  type <- match.arg(type)
  if (identical(object$method %||% "indicator", "naive")) {
    return(.sem_naive_fitted_values(object))
  }
  state_summaries <- .sem_state_summaries(object)

  if (identical(type, "link")) {
    out <- data.frame(time = .observed_time_values(object$stan_data))
    for (nm in names(state_summaries)) {
      out[[nm]] <- state_summaries[[nm]]$mean
    }
    return(out)
  }

  .sem_indicator_means(object, state_summaries)
}


#' One-step-ahead predictions from a copula VAR model
#'
#' Returns point predictions and **marginal** prediction intervals by combining
#' the VAR(1) fitted values with the estimated innovation SDs. Intervals are
#' computed per-variable using a normal approximation and do not account for
#' the copula dependence structure between variables.
#'
#' `predict()` is implemented for the three core single-level models plus the
#' multilevel and SEM variants. For multilevel fits, the methods return
#' unit-specific trajectories. For SEM fits, `type = "link"` returns latent
#' states and `type = "response"` returns observed indicator predictions.
#'
#' @param object A fitted model object.
#' @param type Character; `"link"` (default) returns values on the model's
#'   internal scale (standardized if applicable), `"response"` back-transforms
#'   to the original data scale.
#' @param ci_level Prediction interval level (default: 0.95).
#' @param ... Additional arguments (unused).
#'
#' @return A data frame of marginal prediction intervals at the specified
#'   level. Single-level and SEM fits return columns `time`, `variable`,
#'   `mean`, `lower`, `upper`. Multilevel fits additionally include `unit`.
#' @export
predict.dcvar_model_fit <- function(object, type = c("link", "response"),
                                    ci_level = 0.95, ...) {
  type <- match.arg(type)
  .validate_interval_level(ci_level, arg_name = "ci_level")
  margins <- object$margins %||% "normal"
  if (margins != "normal") {
    cli_abort("Prediction intervals are currently only supported for normal margins, not {.val {margins}}.")
  }

  fit_df <- stats::fitted(object, type = "link")
  sigma_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "sigma_eps", backend = object$backend,
    required = .stan_output_group_pattern("sigma_eps"),
    required_type = "pattern",
    context = "predict.dcvar_model_fit()",
    output_type = "parameter group"
  ))
  sigma_cols <- paste0("sigma_eps[", seq_len(object$stan_data$D), "]")
  sigma_eps <- colMeans(sigma_draws[, sigma_cols, drop = FALSE])

  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)

  D <- object$stan_data$D
  vars <- object$vars

  rows <- list()
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)
  for (d in seq_len(D)) {
    rows[[d]] <- data.frame(
      time = time_values,
      variable = vars[d],
      mean = fit_df[[vars[d]]],
      lower = fit_df[[vars[d]]] - z_crit * sigma_eps[d],
      upper = fit_df[[vars[d]]] + z_crit * sigma_eps[d]
    )
  }

  out <- do.call(rbind, rows)

  if (type == "response" && isTRUE(object$standardized)) {
    Y_means <- attr(object$stan_data, "Y_means")
    Y_sds <- attr(object$stan_data, "Y_sds")
    if (!is.null(Y_means) && !is.null(Y_sds)) {
      for (d in seq_len(D)) {
        idx <- out$variable == vars[d]
        out$mean[idx] <- out$mean[idx] * Y_sds[d] + Y_means[d]
        out$lower[idx] <- out$lower[idx] * Y_sds[d] + Y_means[d]
        out$upper[idx] <- out$upper[idx] * Y_sds[d] + Y_means[d]
      }
    }
  }

  out
}

#' @rdname predict.dcvar_model_fit
#' @export
predict.dcvar_multilevel_fit <- function(object, type = c("link", "response"),
                                         ci_level = 0.95, ...) {
  type <- match.arg(type)
  .validate_interval_level(ci_level, arg_name = "ci_level")
  margins <- object$margins %||% "normal"
  if (!identical(margins, "normal")) {
    cli_abort("Prediction intervals are currently only supported for normal multilevel margins, not {.val {margins}}.")
  }

  fit_df <- .multilevel_fitted_values(object, response = identical(type, "response"))
  sigma <- coef(object)$sigma
  sigma <- as.numeric(sigma)
  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)

  rows <- vector("list", length(object$vars))
  for (d in seq_along(object$vars)) {
    mean <- fit_df[[object$vars[d]]]
    rows[[d]] <- data.frame(
      unit = fit_df$unit,
      time = fit_df$time,
      variable = object$vars[d],
      mean = mean,
      lower = mean - z_crit * sigma[d],
      upper = mean + z_crit * sigma[d],
      check.names = FALSE
    )
  }

  do.call(rbind, rows)
}

#' @rdname predict.dcvar_model_fit
#' @export
predict.dcvar_sem_fit <- function(object, type = c("link", "response"),
                                  ci_level = 0.95, ...) {
  type <- match.arg(type)
  .validate_interval_level(ci_level, arg_name = "ci_level")
  if (identical(object$method %||% "indicator", "naive")) {
    margins <- object$margins %||% "normal"
    if (!identical(margins, "normal")) {
      cli_abort("Prediction intervals are currently only supported for normal naive SEM margins, not {.val {margins}}.")
    }
    fit_df <- .sem_naive_fitted_values(object)
    sigma <- as.numeric(coef(object)$sigma)
    z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
    rows <- vector("list", length(object$vars))
    for (d in seq_along(object$vars)) {
      mean <- fit_df[[object$vars[d]]]
      rows[[d]] <- data.frame(
        time = fit_df$time,
        variable = object$vars[d],
        mean = mean,
        lower = mean - z_crit * sigma[d],
        upper = mean + z_crit * sigma[d]
      )
    }
    return(do.call(rbind, rows))
  }

  state_summaries <- .sem_state_summaries(object)
  if (identical(type, "link")) {
    .sem_latent_predictions(object, state_summaries, ci_level)
  } else {
    .sem_indicator_predictions(object, state_summaries, ci_level)
  }
}
