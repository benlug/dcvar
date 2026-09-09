#' Methods for a Bayesian causal SEM fit
#'
#' @param x,object A `dcvar_causal_sem_fit` object.
#' @param level Probability covered by the central credible interval.
#' @param ... Additional arguments, currently unused.
#' @name dcvar_causal_sem_fit-methods
NULL

#' @describeIn dcvar_causal_sem_fit-methods Print the model and its effects.
#' @return `print()` returns the input invisibly. `summary()` returns a list
#'   with parameter summaries, effects, and diagnostics. `coef()` returns
#'   posterior means on the internal model scale.
#' @export
print.dcvar_causal_sem_fit <- function(x, ...) {
  cat("Bayesian causal SEM\n")
  cat(sprintf("Model: %s\nPersons: %d (control: %d, treated: %d)\n",
              x$model, x$stan_data$N, x$meta$n_group[1], x$meta$n_group[2]))
  cat(sprintf("Backend: %s\nTarget weights: %.3f, %.3f\n",
              x$backend, x$meta$target_weights[1], x$meta$target_weights[2]))
  if (isTRUE(x$stan_data$prior_only == 1L)) cat("Draws use the prior only.\n")
  print(causal_effects(x)[c("effect", "mean", "lower", "upper")], row.names = FALSE)
  invisible(x)
}

#' @describeIn dcvar_causal_sem_fit-methods Summarise parameters and effects.
#' @export
summary.dcvar_causal_sem_fit <- function(object, level = 0.95, ...) {
  .validate_interval_level(level, "level")
  structure(list(
    model = object$model, n_group = object$meta$n_group,
    parameters = .causal_sem_summarise(draws(object, variable = .causal_sem_parameter_roots(object$model)), level),
    effects = causal_effects(object, level = level),
    diagnostics = dcvar_diagnostics(object), meta = object$meta
  ), class = "dcvar_causal_sem_summary")
}

#' @describeIn dcvar_causal_sem_fit-methods Print a fit summary.
#' @export
print.dcvar_causal_sem_summary <- function(x, ...) {
  cat(sprintf("Bayesian causal SEM: %s\n", x$model))
  print(x$effects, row.names = FALSE)
  cat("Parameter summaries use the internal model scale.\n")
  print(x$parameters, row.names = FALSE)
  print(x$diagnostics[c("n_divergent", "n_max_treedepth", "max_rhat",
                        "min_ess_bulk", "min_ess_tail")])
  invisible(x)
}

#' @describeIn dcvar_causal_sem_fit-methods Extract parameter means.
#' @export
coef.dcvar_causal_sem_fit <- function(object, ...) {
  d <- posterior::as_draws_matrix(draws(object, variable = .causal_sem_parameter_roots(object$model)))
  colMeans(d)
}

#' @rdname draws
#' @export
draws.dcvar_causal_sem_fit <- function(object, variable = NULL,
                                       format = "draws_array", ...) {
  format <- match.arg(format, c("draws_array", "draws_matrix", "draws_df"))
  .fit_draws(object$fit, variables = variable, format = format, backend = object$backend)
}

#' List interpretable parameter groups
#' @noRd
.causal_sem_parameter_roots <- function(model) {
  common <- c("lambda", "threshold_1", "threshold_2", "threshold_3", "item_sd")
  c(common, if (model %in% c("latent_covariate", "latent_outcome")) {
    c("alpha", "beta", "mu_v", "sigma_v", "sigma_y")
  } else {
    c("alpha_m", "B", "alpha_y", "D", "d", "mu_q", "sigma_q",
      "rho_q", "sigma_m", "sigma_y")
  })
}

#' @rdname dcvar_diagnostics
#' @export
dcvar_diagnostics.dcvar_causal_sem_fit <- function(object, ...) {
  d <- draws(object)
  roots <- sub("\\[.*$", "", posterior::variables(d))
  keep <- roots %in% .causal_sem_sampled_roots(object$model)
  keep <- keep & (roots != "L_q" | grepl("^L_q\\[[12],2,1\\]$", posterior::variables(d)))
  s <- .causal_sem_summarise(posterior::subset_draws(d, variable = posterior::variables(d)[keep]))
  available_max <- function(x) if (any(!is.na(x))) max(x[!is.na(x)]) else NA_real_
  available_min <- function(x) if (any(!is.na(x))) min(x[!is.na(x)]) else NA_real_
  raw <- .fit_sampler_diagnostics(object$fit, object$backend)
  nms <- dimnames(raw)[[3L]]
  energy <- if ("energy__" %in% nms) {
    e <- matrix(raw[, , "energy__"], ncol = dim(raw)[2L])
    apply(e, 2L, function(x) {
      if (length(x) < 3L || any(!is.finite(x)) || stats::var(x) == 0) return(NA_real_)
      mean(diff(x)^2) / stats::var(x)
    })
  } else NA_real_
  counts <- .fit_diagnostic_summary(object$fit, object$backend)
  effects <- causal_effects(object)
  effects$mcse_ratio <- ifelse(is.finite(effects$sd) & effects$sd > 0,
                               effects$mcse_mean / effects$sd, NA_real_)
  list(n_divergent = sum(counts$num_divergent),
       n_max_treedepth = sum(counts$num_max_treedepth),
       max_rhat = available_max(c(s$rhat, effects$rhat)),
       min_ess_bulk = available_min(c(s$ess_bulk, effects$ess_bulk)),
       min_ess_tail = available_min(c(s$ess_tail, effects$ess_tail)),
       max_parameter_rhat = available_max(s$rhat),
       min_parameter_ess_bulk = available_min(s$ess_bulk),
       min_parameter_ess_tail = available_min(s$ess_tail),
       max_effect_rhat = available_max(effects$rhat),
       min_effect_ess_bulk = available_min(effects$ess_bulk),
       min_effect_ess_tail = available_min(effects$ess_tail),
       max_effect_mcse_ratio = available_max(effects$mcse_ratio),
       high_effect_mcse = any(effects$mcse_ratio > 0.05, na.rm = TRUE),
       incomplete_diagnostics = any(!is.finite(c(s$rhat, s$ess_bulk,
         s$ess_tail, effects$rhat, effects$ess_bulk, effects$ess_tail,
         effects$mcse_ratio, energy))),
       mean_accept_prob = if ("accept_stat__" %in% nms) mean(raw[, , "accept_stat__"]) else NA_real_,
       ebfmi = energy, low_ebfmi = any(energy < 0.3, na.rm = TRUE),
       parameters = s,
       effects = effects[c("effect", "rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_ratio")])
}

#' Select sampled parameters without deterministic outputs
#' @noRd
.causal_sem_sampled_roots <- function(model) {
  # The Stan sources declare these roots in their parameters blocks.
  common <- c("lambda_free", "threshold_1", "threshold_2", "threshold_3",
              "item_sd_treated", "latent_raw")
  c(common, if (model == "latent_covariate") {
    c("mu_v_treated", "sigma_v", "alpha", "beta", "sigma_y")
  } else if (model == "latent_outcome") {
    c("mu_v", "sigma_v", "alpha_treated", "beta", "sigma_y")
  } else if (model == "latent_mediator") {
    c("mu_q", "sigma_q", "L_q", "alpha_m_treated", "B", "alpha_y", "D", "d", "sigma_m", "sigma_y")
  } else {
    c("mu_q1_treated", "mu_q2", "sigma_q", "L_q", "alpha_m", "B", "alpha_y", "D", "d", "sigma_m", "sigma_y")
  })
}

#' Draw replicated observations from a causal SEM
#'
#' @param object A `dcvar_causal_sem_fit` object.
#' @param newdata Not supported. Replicates use the fit's manifest predictors.
#' @param type `"outcome"`, `"indicators"`, `"latent"`, or `"mediator"`.
#' @param summary Return means and credible intervals if `TRUE`.
#' @param level Probability covered by the central interval.
#' @param ... Additional arguments, currently unused.
#' @details
#' Replicates draw new person factors. The latent covariate model draws a new
#' covariate and outcome. The other models condition on the supplied manifest
#' predictors. The baseline mediator model conditions on the observed outcome
#' baseline and generates a new mediator baseline and mediator. Indicator draws
#' contain integer category codes; use the stored levels to recover labels.
#' Outcome and manifest mediator draws use the input units. Latent draws use
#' the identified factor scale. Intervals for category codes summarize those
#' codes; raw draws support category probabilities and predictive checks.
#' @return A `posterior::draws_array`, or a data frame when `summary = TRUE`.
#' @export
predict.dcvar_causal_sem_fit <- function(object, newdata = NULL,
                                         type = c("outcome", "indicators", "latent", "mediator"),
                                         summary = FALSE, level = 0.95, ...) {
  if (!is.null(newdata)) cli_abort("New data predictions are not available for causal SEM fits.")
  type <- match.arg(type)
  .validate_interval_level(level, "level")
  if (!is.logical(summary) || length(summary) != 1L || is.na(summary)) {
    cli_abort("{.arg summary} must be TRUE or FALSE.")
  }
  if (type == "mediator" && !object$model %in% c("latent_mediator", "latent_mediator_baseline")) {
    cli_abort("Mediator predictions require a mediation model.")
  }
  variable <- switch(type, outcome = "y_rep", indicators = "u_rep",
                     latent = "latent_rep", mediator = "m_rep")
  d <- draws(object, variable = variable)
  role <- switch(type, outcome = "outcome", mediator = "mediator", NULL)
  if (!is.null(role)) {
    tr <- object$meta$scales[[role]]
    d[] <- d * tr$scale + tr$center
  }
  if (summary) .causal_sem_summarise(d, level, "variable") else d
}
