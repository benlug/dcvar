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
#' @param object A fitted model object (`dcVar_fit`, `dcVar_hmm_fit`, or
#'   `dcVar_constant_fit`).
#' @param type Character; `"link"` (default) returns values on the model's
#'   internal scale (standardized if applicable), `"response"` back-transforms
#'   to the original data scale.
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `time`, and one column per variable
#'   containing the posterior-mean fitted values.
#' @export
fitted.dcVar_model_fit <- function(object, type = c("link", "response"), ...) {
  type <- match.arg(type)
  D <- object$stan_data$D
  T_obs <- object$stan_data$T
  Y <- object$stan_data$Y
  mu_draws <- posterior::as_draws_matrix(object$fit$draws("mu"))
  Phi_draws <- posterior::as_draws_matrix(object$fit$draws("Phi"))

  mu_cols <- paste0("mu[", seq_len(D), "]")
  phi_cols <- unlist(lapply(seq_len(D), function(i) paste0("Phi[", i, ",", seq_len(D), "]")))
  mu_mat <- mu_draws[, mu_cols, drop = FALSE]
  phi_mat <- Phi_draws[, phi_cols, drop = FALSE]

  y_hat <- matrix(NA_real_, T_obs - 1, D)
  for (t in 2:T_obs) {
    centered_prev <- sweep(mu_mat, 2, Y[t - 1, ], FUN = function(mu_val, y_prev) y_prev - mu_val)
    y_hat_draws <- matrix(NA_real_, nrow(mu_mat), D)

    for (d in seq_len(D)) {
      eta <- mu_mat[, d]
      for (j in seq_len(D)) {
        eta <- eta + phi_mat[, (d - 1) * D + j] * centered_prev[, j]
      }
      y_hat_draws[, d] <- eta
    }

    y_hat[t - 1, ] <- colMeans(y_hat_draws)
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


#' @rdname fitted.dcVar_model_fit
#' @export
fitted.dcVar_multilevel_fit <- function(object, ...) {
  cli_abort(c(
    "{.fun fitted} is not yet implemented for multilevel models.",
    "i" = "Use {.fun var_params} and {.fun random_effects} to inspect parameters."
  ))
}

#' @rdname fitted.dcVar_model_fit
#' @export
fitted.dcVar_sem_fit <- function(object, ...) {
  cli_abort(c(
    "{.fun fitted} is not yet implemented for SEM models.",
    "i" = "Use {.fun latent_states} to extract estimated latent trajectories."
  ))
}


#' One-step-ahead predictions from a copula VAR model
#'
#' Returns point predictions and **marginal** prediction intervals by combining
#' the VAR(1) fitted values with the estimated innovation SDs. Intervals are
#' computed per-variable using a normal approximation and do not account for
#' the copula dependence structure between variables.
#'
#' @param object A fitted model object.
#' @param type Character; `"link"` (default) returns values on the model's
#'   internal scale (standardized if applicable), `"response"` back-transforms
#'   to the original data scale.
#' @param ci_level Prediction interval level (default: 0.95).
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `time`, `variable`, `mean`, `lower`,
#'   `upper` (marginal prediction interval at the specified level).
#' @export
predict.dcVar_model_fit <- function(object, type = c("link", "response"),
                                    ci_level = 0.95, ...) {
  type <- match.arg(type)
  margins <- object$margins %||% "normal"
  if (margins != "normal") {
    cli_abort("Prediction intervals are currently only supported for normal margins, not {.val {margins}}.")
  }

  fit_df <- stats::fitted(object, type = "link")
  sigma_draws <- posterior::as_draws_matrix(object$fit$draws("sigma_eps"))
  sigma_cols <- paste0("sigma_eps[", seq_len(object$stan_data$D), "]")
  sigma_eps <- colMeans(sigma_draws[, sigma_cols, drop = FALSE])

  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)

  D <- object$stan_data$D
  T_obs <- object$stan_data$T
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

#' @rdname predict.dcVar_model_fit
#' @export
predict.dcVar_multilevel_fit <- function(object, ...) {
  cli_abort(c(
    "{.fun predict} is not yet implemented for multilevel models.",
    "i" = "Use {.fun var_params} and {.fun random_effects} to inspect parameters."
  ))
}

#' @rdname predict.dcVar_model_fit
#' @export
predict.dcVar_sem_fit <- function(object, ...) {
  cli_abort(c(
    "{.fun predict} is not yet implemented for SEM models.",
    "i" = "Use {.fun latent_states} to extract estimated latent trajectories."
  ))
}
