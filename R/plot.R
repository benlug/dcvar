# ============================================================================
# Plotting Functions
# ============================================================================

#' Plot the rho trajectory with credible intervals
#'
#' @param object A `dcvar_fit` or `dcvar_hmm_fit` object.
#' @param show_ci Logical; show credible interval ribbons (default: `TRUE`).
#' @param ci_level Credible interval level for the outer ribbon (default: 0.95).
#' @param inner_level Credible interval level for the inner ribbon (default:
#'   0.80). Set to `NULL` to disable the inner ribbon.
#' @param true_rho Optional numeric vector of true rho values for overlay
#'   (useful for simulation studies).
#' @param title Plot title.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
plot_rho <- function(object, show_ci = TRUE, ci_level = 0.95, inner_level = 0.80,
                     true_rho = NULL, title = NULL, ...) {
  .validate_interval_level(ci_level, arg_name = "ci_level")
  if (!is.null(inner_level)) {
    .validate_interval_level(inner_level, arg_name = "inner_level")
    if (inner_level >= ci_level) {
      cli_abort(c(
        "{.arg inner_level} must be smaller than {.arg ci_level}.",
        "i" = "The inner ribbon must be narrower than the outer ribbon."
      ))
    }
  }

  # Compute the quantile probs needed for the requested CI levels
  alpha_outer <- (1 - ci_level) / 2
  probs <- c(alpha_outer, 1 - alpha_outer, 0.5)
  if (!is.null(inner_level)) {
    alpha_inner <- (1 - inner_level) / 2
    probs <- c(alpha_outer, alpha_inner, 0.5, 1 - alpha_inner, 1 - alpha_outer)
  }
  probs <- sort(unique(probs))

  rho_df <- rho_trajectory(object, probs = probs)

  if (is.null(title)) {
    model_name <- switch(object$model,
      dcvar = "DC-VAR",
      hmm = "HMM Copula",
      constant = "Constant Copula",
      object$model
    )
    title <- paste0("Time-Varying Copula Dependence (", model_name, ")")
  }

  q_outer_low <- paste0("q", alpha_outer * 100)
  q_outer_high <- paste0("q", (1 - alpha_outer) * 100)

  p <- ggplot2::ggplot(rho_df, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")

  if (show_ci && q_outer_low %in% names(rho_df)) {
    # Outer CI
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[q_outer_low]], ymax = .data[[q_outer_high]]),
      fill = "steelblue", alpha = 0.2
    )
    # Inner CI
    if (!is.null(inner_level)) {
      q_inner_low <- paste0("q", alpha_inner * 100)
      q_inner_high <- paste0("q", (1 - alpha_inner) * 100)
      if (q_inner_low %in% names(rho_df) && q_inner_high %in% names(rho_df)) {
        p <- p + ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data[[q_inner_low]], ymax = .data[[q_inner_high]]),
          fill = "steelblue", alpha = 0.3
        )
      }
    }
  }

  p <- p +
    ggplot2::geom_line(ggplot2::aes(y = .data$mean), color = "steelblue", linewidth = 1)

  if ("q50" %in% names(rho_df)) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = .data$q50), color = "darkblue", linetype = "dashed")
  }

  if (!is.null(true_rho)) {
    time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)
    if (!is.numeric(true_rho)) {
      cli_abort("{.arg true_rho} must be a numeric vector.")
    }
    if (length(true_rho) != length(time_values)) {
      cli_abort(c(
        "{.arg true_rho} must have the same length as the observed rho time axis.",
        "i" = "Use the fit's stored time values for simulation overlays."
      ))
    }
    true_df <- data.frame(time = time_values, rho = true_rho)
    p <- p + ggplot2::geom_line(
      data = true_df,
      ggplot2::aes(x = .data$time, y = .data$rho),
      color = "red", linewidth = 1, linetype = "solid"
    )
  }

  # Build subtitle dynamically from the actual CI levels
  subtitle <- if (show_ci) {
    if (!is.null(inner_level)) {
      sprintf("%d%% and %d%% credible intervals",
              round(ci_level * 100), round(inner_level * 100))
    } else {
      sprintf("%d%% credible interval", round(ci_level * 100))
    }
  } else {
    NULL
  }

  p <- p +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::labs(
      x = "Time",
      y = expression(rho[t]),
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )

  p
}


#' Plot VAR(1) coefficient matrix as a heatmap
#'
#' @param object A fitted model object.
#' @param var_names Character vector of variable names for axis labels.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
plot_phi <- function(object, var_names = NULL, ...) {
  summ <- .fit_summary(
    object$fit, variables = NULL, backend = object$backend,
    required = .stan_output_group_pattern("Phi"),
    required_type = "pattern",
    context = "plot_phi()",
    output_type = "parameter group",
    mean,
    ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  )

  phi_rows <- grep("^Phi\\[", summ$variable)
  phi_df <- data.frame(
    variable = summ$variable[phi_rows],
    mean = summ$mean[phi_rows],
    q2.5 = summ$q2.5[phi_rows],
    q97.5 = summ$q97.5[phi_rows]
  )

  idx <- regmatches(phi_df$variable, regexec("Phi\\[(\\d+),(\\d+)\\]", phi_df$variable))
  if (any(lengths(idx) != 3)) {
    cli_abort("Unexpected Phi variable name format in Stan output.")
  }
  phi_df$row <- as.integer(vapply(idx, `[`, character(1), 2))
  phi_df$col <- as.integer(vapply(idx, `[`, character(1), 3))

  if (is.null(var_names)) {
    var_names <- if (!is.null(object$vars)) object$vars else paste0("Y", 1:max(phi_df$row))
  }

  phi_df$from <- var_names[phi_df$col]
  phi_df$to <- var_names[phi_df$row]
  phi_df$label <- sprintf("%.2f\n[%.2f, %.2f]", phi_df$mean, phi_df$q2.5, phi_df$q97.5)

  # Symmetric limits around zero based on actual coefficient range
  max_abs <- max(abs(phi_df$mean), na.rm = TRUE)
  fill_limit <- max(max_abs * 1.1, 0.1)  # at least 0.1 to avoid degenerate scale

  ggplot2::ggplot(phi_df, ggplot2::aes(x = .data$from, y = .data$to)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$mean), color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3) +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-fill_limit, fill_limit)
    ) +
    ggplot2::labs(
      x = "From (t-1)", y = "To (t)",
      fill = expression(Phi),
      title = "VAR(1) Coefficient Matrix"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::coord_fixed()
}


#' Plot MCMC diagnostics
#'
#' Creates a combined panel with trace plots, Rhat, and ESS diagnostics.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments (unused).
#'
#' @return A combined ggplot object (via patchwork).
#' @export
plot_diagnostics <- function(object, ...) {
  margins <- object$margins %||% "normal"

  trace_pars <- if (object$model == "multilevel") {
    c("phi_bar[1]", "phi_bar[2]", "phi_bar[3]", "phi_bar[4]", "rho")
  } else if (object$model == "sem") {
    if (margins == "exponential") {
      c("mu[1]", "mu[2]", "phi11", "phi22", "sigma_exp[1]", "sigma_exp[2]", "rho")
    } else {
      c("mu[1]", "mu[2]", "phi11", "phi22", "sigma[1]", "sigma[2]", "rho")
    }
  } else if (margins == "normal") {
    c("mu[1]", "mu[2]", "sigma_eps[1]", "sigma_eps[2]")
  } else if (margins == "exponential") {
    c("mu[1]", "mu[2]", "sigma_exp[1]", "sigma_exp[2]")
  } else if (margins == "skew_normal") {
    c("mu[1]", "mu[2]", "omega[1]", "omega[2]", "delta[1]", "delta[2]")
  } else if (margins == "gamma") {
    c("mu[1]", "mu[2]", "sigma_gam[1]", "sigma_gam[2]", "shape_gam")
  } else {
    c("mu[1]", "mu[2]")
  }
  if (object$model == "dcvar") trace_pars <- c(trace_pars, "sigma_omega")
  if (object$model %in% c("dcvar_covariate", "dcvar_covariate_nodrift")) {
    trace_pars <- c(trace_pars, "beta_0", paste0("beta[", seq_len(object$stan_data$P), "]"))
    if (object$model == "dcvar_covariate") trace_pars <- c(trace_pars, "sigma_omega")
  }
  if (object$model == "hmm") {
    K <- object$K
    trace_pars <- c(trace_pars, paste0("rho_state[", 1:K, "]"))
  }
  if (object$model == "constant") trace_pars <- c(trace_pars, "rho")

  required_patterns <- unique(vapply(
    trace_pars,
    .stan_output_group_pattern,
    character(1)
  ))

  draws_array <- .fit_draws(
    object$fit,
    format = "draws_array",
    backend = object$backend,
    required = required_patterns,
    required_type = "pattern",
    context = "plot_diagnostics()",
    output_type = "parameter group"
  )

  p1 <- bayesplot::mcmc_trace(draws_array, pars = trace_pars) +
    ggplot2::ggtitle("Trace Plots")

  summ <- .fit_summary(
    object$fit,
    backend = object$backend,
    required = required_patterns,
    required_type = "pattern",
    context = "plot_diagnostics()",
    output_type = "parameter group"
  )
  rhats <- summ$rhat[is.finite(summ$rhat)]
  p2 <- bayesplot::mcmc_rhat(rhats) +
    ggplot2::ggtitle("R-hat Diagnostics")

  n_eff_ratio <- summ$ess_bulk[is.finite(summ$ess_bulk)]
  total_draws <- object$meta$iter_sampling * object$meta$chains
  p3 <- bayesplot::mcmc_neff(n_eff_ratio / total_draws) +
    ggplot2::ggtitle("Effective Sample Size Ratio")

  bottom_row <- patchwork::wrap_plots(p2, p3, ncol = 2)
  (p1 / bottom_row) +
    patchwork::plot_annotation(
      title = "MCMC Diagnostics",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14)
      )
    )
}


#' Posterior predictive check for residual correlations
#'
#' @param object A fitted model object.
#' @param n_sample Number of posterior draws to use (default: 100).
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#'
#' @details Posterior predictive checks are currently available for normal and
#'   exponential margins. Gamma and skew-normal fits store copula-level
#'   replicated z-scores in `eps_rep`, so their replicated draws are not on the
#'   same residual scale as `eps`.
#' @export
plot_ppc <- function(object, n_sample = 100, ...) {
  if (inherits(object, "dcvar_multilevel_fit") || inherits(object, "dcvar_sem_fit")) {
    cli_abort("Posterior predictive checks are not yet supported for {.cls {class(object)[1]}} models.")
  }

  margins <- object$margins %||% "normal"
  if (margins %in% c("gamma", "skew_normal")) {
    cli_abort(
      c(
        "Posterior predictive checks are not supported for {.val {margins}} margins.",
        "i" = "The stored {.field eps_rep} draws are copula-level z-scores rather than replicated residuals on the observed margin scale."
      )
    )
  }

  if (!margins %in% c("normal", "exponential")) {
    cli_abort("Posterior predictive checks are not implemented for {.val {margins}} margins.")
  }

  eps_rep_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "eps_rep", backend = object$backend,
    required = .stan_output_group_pattern("eps_rep"),
    required_type = "pattern",
    context = "plot_ppc()",
    output_type = "generated quantity"
  ))
  eps_obs_draws <- posterior::as_draws_matrix(.fit_draws(
    object$fit, "eps", backend = object$backend,
    required = .stan_output_group_pattern("eps"),
    required_type = "pattern",
    context = "plot_ppc()",
    output_type = "transformed parameter"
  ))
  eps_rep_cols_1 <- grep("eps_rep\\[.*,1\\]", colnames(eps_rep_draws))
  eps_rep_cols_2 <- grep("eps_rep\\[.*,2\\]", colnames(eps_rep_draws))
  eps_obs_cols_1 <- grep("eps\\[.*,1\\]", colnames(eps_obs_draws))
  eps_obs_cols_2 <- grep("eps\\[.*,2\\]", colnames(eps_obs_draws))

  # Observed residuals (posterior mean)
  eps_obs_1 <- colMeans(eps_obs_draws[, eps_obs_cols_1, drop = FALSE])
  eps_obs_2 <- colMeans(eps_obs_draws[, eps_obs_cols_2, drop = FALSE])
  cor_obs <- cor(eps_obs_1, eps_obs_2, use = "complete.obs")

  # Replicated correlations
  n_draws <- nrow(eps_rep_draws)
  n_sample <- min(n_sample, n_draws)
  sample_idx <- sample(seq_len(n_draws), n_sample)

  cor_rep <- numeric(n_sample)
  for (i in seq_along(sample_idx)) {
    idx <- sample_idx[i]
    rep_1 <- as.numeric(eps_rep_draws[idx, eps_rep_cols_1, drop = TRUE])
    rep_2 <- as.numeric(eps_rep_draws[idx, eps_rep_cols_2, drop = TRUE])
    cor_rep[i] <- cor(rep_1, rep_2, use = "complete.obs")
  }

  ggplot2::ggplot(data.frame(cor = cor_rep), ggplot2::aes(x = .data$cor)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = cor_obs, color = "red", linewidth = 1.5) +
    ggplot2::labs(
      x = "Residual Correlation",
      y = "Count",
      title = "Posterior Predictive Check",
      subtitle = "Red line: observed; Histogram: replicated data"
    ) +
    ggplot2::theme_minimal()
}


#' Plot HMM state posteriors
#'
#' @param object A `dcvar_hmm_fit` object.
#' @param show_viterbi Logical; overlay the Viterbi (MAP) state sequence
#'   (default: `TRUE`).
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
plot_hmm_states <- function(object, show_viterbi = TRUE, ...) {
  states <- hmm_states(object)
  K <- object$K
  T_eff <- nrow(states$gamma)
  time_values <- .observed_time_values(object$stan_data, drop_first = TRUE)

  # Build long-format data for gamma
  gamma_df <- data.frame(
    time = rep(time_values, K),
    state = rep(paste0("State ", 1:K), each = T_eff),
    probability = as.vector(states$gamma)
  )

  p <- ggplot2::ggplot(gamma_df, ggplot2::aes(x = .data$time, y = .data$probability,
                                                fill = .data$state)) +
    ggplot2::geom_area(alpha = 0.7, position = "stack") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      x = "Time", y = "State Probability",
      title = "HMM State Posteriors",
      fill = "State"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (show_viterbi) {
    viterbi_df <- data.frame(
      time = time_values,
      state = paste0("State ", states$viterbi),
      y = 0.5
    )
    p <- p + ggplot2::geom_point(
      data = viterbi_df,
      ggplot2::aes(x = .data$time, y = .data$y, color = .data$state),
      size = 0.5, alpha = 0.5, inherit.aes = FALSE
    ) +
      ggplot2::scale_color_brewer(palette = "Set2", guide = "none")
  }

  p
}


#' Plot random effects (caterpillar plot)
#'
#' Displays unit-specific VAR coefficients with credible intervals.
#'
#' @param object A `dcvar_multilevel_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
plot_random_effects <- function(object, ...) {
  re <- random_effects(object)

  re$unit <- factor(re$unit)

  ggplot2::ggplot(re, ggplot2::aes(
    x = .data$unit, y = .data$mean,
    ymin = .data$q2.5, ymax = .data$q97.5
  )) +
    ggplot2::geom_pointrange(size = 0.3, color = "steelblue") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::labs(
      x = "Unit",
      y = "Posterior Mean [95% CI]",
      title = "Unit-Specific VAR Coefficients"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(size = 6)
    )
}


#' Plot estimated latent states with credible intervals
#'
#' @param object A `dcvar_sem_fit` object.
#' @param true_states Optional T x 2 matrix of true latent states for overlay.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object.
#' @export
plot_latent_states <- function(object, true_states = NULL, ...) {
  states_df <- latent_states(object)

  p <- ggplot2::ggplot(states_df, ggplot2::aes(x = .data$time)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$q2.5, ymax = .data$q97.5),
      fill = "steelblue", alpha = 0.2
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$mean),
      color = "steelblue", linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::labs(
      x = "Time", y = "Latent State",
      title = "Estimated Latent States",
      subtitle = "95% credible intervals"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
  )

  if (!is.null(true_states)) {
    time_values <- .observed_time_values(object$stan_data)
    if (!is.matrix(true_states) || ncol(true_states) != 2) {
      cli_abort("{.arg true_states} must be a matrix with 2 columns.")
    }
    if (nrow(true_states) != length(time_values)) {
      cli_abort(c(
        "{.arg true_states} must have one row per observed time point.",
        "i" = "Use the fit's stored time values for simulation overlays."
      ))
    }
    true_df <- data.frame(
      time = rep(time_values, 2),
      variable = rep(object$vars, each = nrow(true_states)),
      value = c(true_states[, 1], true_states[, 2])
    )
    p <- p + ggplot2::geom_line(
      data = true_df,
      ggplot2::aes(x = .data$time, y = .data$value),
      color = "red", linewidth = 0.5, linetype = "dashed"
    )
  }

  p
}


#' Internal: plot transition matrix as a heatmap
#' @noRd
.plot_transition_matrix <- function(A, K) {
  A_df <- expand.grid(from = 1:K, to = 1:K)
  A_df$prob <- as.vector(A)
  A_df$label <- sprintf("%.3f", A_df$prob)
  A_df$from <- paste0("State ", A_df$from)
  A_df$to <- paste0("State ", A_df$to)

  ggplot2::ggplot(A_df, ggplot2::aes(x = .data$from, y = .data$to)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$prob), color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 4) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
    ggplot2::labs(
      x = "From State", y = "To State",
      fill = "P(transition)",
      title = "Transition Matrix"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::coord_fixed()
}
