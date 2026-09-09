#' Extract causal effects from a Bayesian SEM
#'
#' Computes effects from each joint posterior draw. The default target assigns
#' equal weight to both groups. The fit records any other chosen weights.
#'
#' @param object A `dcvar_causal_sem_fit` object.
#' @param summary Return a summary if `TRUE`. Otherwise return posterior draws.
#' @param level Probability covered by the central credible interval.
#' @param target_weights Optional weights for the control and treated groups,
#'   in that order. Two nonnegative values must sum to one. `NULL` uses the fit.
#' @param covariate_values Optional finite values for conditional treatment
#'   effects in the latent covariate or latent outcome model.
#' @param scale `"response"` returns effects in the input outcome units and
#'   takes covariate values in input units. `"model"` uses the internal scale.
#'   A latent variable always uses its identified factor scale.
#' @param ... Additional arguments for methods.
#'
#' @details
#' For the first two models, ATE is the intercept difference plus the slope
#' difference times the weighted covariate mean. The interaction is the slope
#' difference. For mediation, the direct effect adjusts both groups at the
#' same weighted baseline and mediator means. The total effect integrates the
#' mediator regression at the weighted baseline mean. The indirect effect is
#' ATE minus ADE. These effects use the definitions in the SEM study.
#'
#' @return A data frame with effect summaries, or a `posterior::draws_array`
#'   when `summary = FALSE`. Attributes record the scale and target weights.
#' @export
causal_effects <- function(object, ...) {
  UseMethod("causal_effects")
}

#' @rdname causal_effects
#' @export
causal_effects.default <- function(object, ...) {
  cli_abort("{.fun causal_effects} requires a {.cls dcvar_causal_sem_fit} object.")
}

#' @rdname causal_effects
#' @export
causal_effects.dcvar_causal_sem_fit <- function(object, summary = TRUE,
                                               level = 0.95,
                                               target_weights = NULL,
                                               covariate_values = NULL,
                                               scale = c("response", "model"),
                                               ...) {
  scale <- match.arg(scale)
  .validate_interval_level(level, "level")
  if (!is.logical(summary) || length(summary) != 1L || is.na(summary)) {
    cli_abort("{.arg summary} must be TRUE or FALSE.")
  }
  weights <- .causal_sem_effect_weights(target_weights %||% object$meta$target_weights)
  d <- .causal_sem_effect_draws(object, weights, covariate_values, scale)
  out <- if (summary) .causal_sem_summarise(d, level, "effect") else d
  attr(out, "target_weights") <- weights
  attr(out, "scale") <- scale
  attr(out, "outcome") <- object$meta$roles$outcome
  attr(out, "outcome_units") <- if (object$model == "latent_outcome") {
    "identified factor scale"
  } else if (scale == "response") "input units" else "internal units"
  if (!is.null(covariate_values)) attr(out, "covariate_values") <- covariate_values
  out
}

#' Check target weights
#' @noRd
.causal_sem_effect_weights <- function(weights) {
  .causal_sem_validate_weights(weights)
}

#' Compute effects without breaking the chain structure
#' @noRd
.causal_sem_effect_draws <- function(object, weights, covariate_values, scale) {
  regression <- object$model %in% c("latent_covariate", "latent_outcome")
  variables <- if (regression) c("alpha", "beta", "mu_v") else {
    c("alpha_m", "B", "alpha_y", "D", "d", "mu_q")
  }
  a <- draws(object, variable = variables)
  dims <- dim(a)
  get_draw <- function(name) {
    if (!name %in% dimnames(a)[[3L]]) {
      cli_abort("The fit does not contain required parameter {.val {name}}.")
    }
    matrix(a[, , name], nrow = dims[1L], ncol = dims[2L])
  }
  if (regression) {
    delta_a <- get_draw("alpha[2]") - get_draw("alpha[1]")
    delta_b <- get_draw("beta[2]") - get_draw("beta[1]")
    mean_v <- weights[1] * get_draw("mu_v[1]") + weights[2] * get_draw("mu_v[2]")
    effects <- list(ATE = delta_a + delta_b * mean_v, interaction = delta_b)
    if (!is.null(covariate_values)) {
      if (!is.numeric(covariate_values) || !length(covariate_values) ||
          any(!is.finite(covariate_values))) {
        cli_abort("{.arg covariate_values} must contain finite numeric values.")
      }
      x <- covariate_values
      if (scale == "response") {
        tr <- object$meta$scales$covariate
        x <- (x - tr$center) / tr$scale
      }
      for (i in seq_along(x)) effects[[paste0("conditional[", i, "]")]] <- delta_a + delta_b * x[i]
    }
  } else {
    if (!is.null(covariate_values)) {
      cli_abort("{.arg covariate_values} is only available for the two regression models.")
    }
    mu <- lapply(seq_len(2L), function(j) {
      weights[1] * get_draw(sprintf("mu_q[1,%d]", j)) +
        weights[2] * get_draw(sprintf("mu_q[2,%d]", j))
    })
    mean_m <- list()
    adjusted_m <- list()
    for (g in seq_len(2L)) {
      mean_m[[g]] <- get_draw(sprintf("alpha_m[%d]", g))
      adjusted_m[[g]] <- mean_m[[g]]
      for (j in seq_len(2L)) {
        slope <- get_draw(sprintf("B[%d,%d]", g, j))
        mean_m[[g]] <- mean_m[[g]] + slope * get_draw(sprintf("mu_q[%d,%d]", g, j))
        adjusted_m[[g]] <- adjusted_m[[g]] + slope * mu[[j]]
      }
    }
    pooled_m <- weights[1] * mean_m[[1]] + weights[2] * mean_m[[2]]
    total <- direct <- list()
    for (g in seq_len(2L)) {
      base <- get_draw(sprintf("alpha_y[%d]", g))
      for (j in seq_len(2L)) base <- base + get_draw(sprintf("D[%d,%d]", g, j)) * mu[[j]]
      d <- get_draw(sprintf("d[%d]", g))
      total[[g]] <- base + d * adjusted_m[[g]]
      direct[[g]] <- base + d * pooled_m
    }
    ate <- total[[2]] - total[[1]]
    ade <- direct[[2]] - direct[[1]]
    effects <- list(ATE = ate, ADE = ade, AIE = ate - ade)
  }
  if (scale == "response") {
    y_scale <- object$meta$scales$outcome$scale
    effects <- lapply(effects, function(x) x * y_scale)
    if (regression) effects$interaction <- effects$interaction / object$meta$scales$covariate$scale
  }
  out <- array(unlist(effects, use.names = FALSE),
               dim = c(dims[1L], dims[2L], length(effects)),
               dimnames = list(iteration = dimnames(a)[[1L]],
                               chain = dimnames(a)[[2L]], variable = names(effects)))
  posterior::as_draws_array(out)
}

#' Summarise draws with intervals and chain diagnostics
#' @noRd
.causal_sem_summarise <- function(d, level = 0.95, label = "parameter") {
  lo <- (1 - level) / 2
  out <- suppressWarnings(posterior::summarise_draws(
    d, mean = mean, median = stats::median, sd = stats::sd,
    lower = function(x) stats::quantile(x, lo, names = FALSE),
    upper = function(x) stats::quantile(x, 1 - lo, names = FALSE),
    probability_positive = function(x) mean(x > 0),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail, mcse_mean = posterior::mcse_mean
  ))
  out <- as.data.frame(out)
  names(out)[names(out) == "variable"] <- label
  out
}
