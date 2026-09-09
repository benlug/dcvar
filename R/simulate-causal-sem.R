#' Simulate data from a causal SEM model
#'
#' Generate independent persons with one latent factor and three ordinal
#' indicators. All parameters use the identified model scale. The first
#' loading is one. The control factor mean is zero. Each control item SD
#' is one. The returned fit arguments turn off standardization.
#'
#' @param n Total number of persons, or two group sizes. A total uses groups
#'   of equal size when possible. Each group needs at least one person.
#' @param model A causal SEM model name. See [prepare_causal_sem_data()].
#' @param parameters Named list of parameter values. For the first two models,
#'   use length-two vectors `mu_v`, `sigma_v`, `alpha`, `beta`, and `sigma_y`.
#'   For mediation, use 2 by 2 matrices `mu_q`, `sigma_q`, `B`, and `D`.
#'   Also use length-two vectors `rho_q`, `alpha_m`, `alpha_y`, `d`,
#'   `sigma_m`, and `sigma_y`. Matrix rows identify control and treatment.
#'   Matrix columns identify the mediator baseline and outcome baseline.
#'   All models accept `lambda` of length three and `item_sd` of size 2 by 3.
#'   An explicit value must agree with each identification restriction.
#' @param thresholds A list of three threshold vectors, or one vector for
#'   all items. Each vector has two to four strictly increasing values.
#'   Alternatively, set `threshold_1`, `threshold_2`, and `threshold_3` in
#'   `parameters`. The default is `c(-0.7, 0.7)` for each item.
#' @param seed Random seed.
#' @param target_weights Target population weights in control, treatment order.
#'   A sum within `1e-8` of one is accepted. The weights are then divided by
#'   their sum to remove roundoff.
#' @param parameter_prior Logical. Draw parameters from the default model prior.
#'   This option uses three categories for each item. To use other category
#'   counts, supply a full parameter draw through `parameters`.
#' @param priors Prior settings for parameter draws. See [prepare_causal_sem_data()].
#' @return A list with `data`, `parameters`, `effects`, `args`, `roles`, and
#'   `latent`. The last element contains the true person factors. Use
#'   `do.call(dcvar_causal_sem, sim$args)` to fit the returned data.
#'   Prior draws can produce a constant indicator. Data preparation reports
#'   this case as an error. The simulator does not replace such a draw.
#' @export
simulate_dcvar_causal_sem <- function(n = 400, model = "latent_covariate",
                                      parameters = list(), thresholds = NULL,
                                      seed = NULL, target_weights = c(0.5, 0.5),
                                      parameter_prior = FALSE, priors = list()) {
  model <- .causal_sem_validate_model(model)
  .causal_sem_validate_flag(parameter_prior, "parameter_prior")
  target_weights <- .causal_sem_validate_weights(target_weights)
  priors <- .causal_sem_priors(priors)
  if (!is.numeric(n) || !is.null(dim(n)) || !length(n) %in% c(1L, 2L) ||
      any(!is.finite(n)) ||
      any(n < 1) || any(n != floor(n)) || any(n > .Machine$integer.max)) {
    cli_abort("{.arg n} must be one total size or two positive integer group sizes.")
  }
  sizes <- if (length(n) == 1L) c(floor(n / 2), n - floor(n / 2)) else n
  if (any(sizes < 1) || sum(sizes) > .Machine$integer.max) {
    cli_abort("The total size must allow at least one person in each group.")
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || !is.null(dim(seed)) ||
        length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
      cli_abort("{.arg seed} must be NULL or a non-negative integer.")
    }
    set.seed(seed)
  }
  if (parameter_prior) {
    if (length(parameters) || !is.null(thresholds)) {
      cli_abort("Use {.arg parameter_prior} without explicit parameters or thresholds.")
    }
    parameters <- .causal_sem_draw_prior(model, priors)
  }
  p <- .causal_sem_sim_parameters(model, parameters, thresholds)
  group <- rep(seq_len(2L), times = sizes)
  N <- length(group)
  roles <- list(outcome = "y", covariate = NULL, mediator = NULL,
                mediator_baseline = NULL, outcome_baseline = NULL)
  data <- data.frame(x = group - 1L)
  if (model %in% c("latent_covariate", "latent_outcome")) {
    v <- stats::rnorm(N, p$mu_v[group], p$sigma_v[group])
    y <- stats::rnorm(N, p$alpha[group] + p$beta[group] * v, p$sigma_y[group])
    if (model == "latent_covariate") {
      latent <- v
      latent_name <- "xi"
      roles$covariate <- latent_name
      data$y <- y
    } else {
      latent <- y
      latent_name <- "eta"
      roles$outcome <- latent_name
      roles$covariate <- "z"
      data$z <- v
    }
  } else {
    noise <- matrix(stats::rnorm(N * 2L), N, 2L)
    q1 <- p$mu_q[group, 1] + p$sigma_q[group, 1] * noise[, 1]
    q2 <- p$mu_q[group, 2] + p$sigma_q[group, 2] * (
      p$rho_q[group] * noise[, 1] + sqrt(1 - p$rho_q[group]^2) * noise[, 2]
    )
    m <- stats::rnorm(N,
      p$alpha_m[group] + p$B[group, 1] * q1 + p$B[group, 2] * q2,
      p$sigma_m[group]
    )
    data$y <- stats::rnorm(N,
      p$alpha_y[group] + p$D[group, 1] * q1 + p$D[group, 2] * q2 + p$d[group] * m,
      p$sigma_y[group]
    )
    data$ypre <- q2
    roles$mediator <- "m"
    roles$mediator_baseline <- "mpre"
    roles$outcome_baseline <- "ypre"
    if (model == "latent_mediator") {
      latent <- m
      latent_name <- "eta_m"
      roles$mediator <- latent_name
      data$mpre <- q1
    } else {
      latent <- q1
      latent_name <- "xi_mpre"
      roles$mediator_baseline <- latent_name
      data$m <- m
    }
  }
  item_names <- paste0("u", seq_len(3L))
  for (j in seq_len(3L)) {
    threshold <- p[[paste0("threshold_", j)]]
    response <- stats::rnorm(N, p$lambda[j] * latent, p$item_sd[group, j])
    data[[item_names[j]]] <- cut(
      response, breaks = c(-Inf, threshold, Inf),
      labels = as.character(seq_len(length(threshold) + 1L)), ordered_result = TRUE
    )
  }
  indicators <- stats::setNames(list(item_names), latent_name)
  args <- c(list(data = data, model = model, treatment = "x", control = 0, treated = 1),
    roles, list(indicators = indicators, target_weights = target_weights,
                standardize = FALSE, priors = priors)
  )
  list(data = data, parameters = p,
       effects = .causal_sem_sim_effects(model, p, target_weights),
       args = args, roles = roles, latent = as.numeric(latent))
}

#' Internal: complete parameters on the identified scale
#' @noRd
.causal_sem_sim_parameters <- function(model, parameters, thresholds = NULL) {
  if (!is.list(parameters) || (length(parameters) &&
      (is.null(names(parameters)) || anyNA(names(parameters)) ||
       any(!nzchar(names(parameters))) || anyDuplicated(names(parameters))))) {
    cli_abort("{.arg parameters} must be a list with unique non-empty names.")
  }
  p <- list(lambda = c(1, 0.9, 1.1), item_sd = rbind(rep(1, 3), c(1.1, 0.9, 1.2)),
            threshold_1 = c(-0.7, 0.7), threshold_2 = c(-0.7, 0.7),
            threshold_3 = c(-0.7, 0.7))
  regression <- model %in% c("latent_covariate", "latent_outcome")
  if (regression) {
    p <- c(p, list(mu_v = c(0, 0.4), sigma_v = c(1, 1.1),
                   alpha = c(0, 0.5), beta = c(0.4, 0.7), sigma_y = c(0.8, 1)))
  } else {
    p <- c(p, list(mu_q = matrix(c(0, 0, 0.3, -0.2), 2, 2, byrow = TRUE),
      sigma_q = matrix(c(1, 0.9, 1.1, 1.2), 2, 2, byrow = TRUE),
      rho_q = c(0.25, -0.15), alpha_m = c(0, 0.4),
      B = matrix(c(0.45, 0.1, 0.6, 0.2), 2, 2, byrow = TRUE),
      alpha_y = c(0.1, 0.6), D = matrix(c(0.1, 0.45, 0.2, 0.3), 2, 2, byrow = TRUE),
      d = c(0.5, 0.75), sigma_m = c(0.8, 0.9), sigma_y = c(0.9, 1)))
  }
  unknown <- setdiff(names(parameters), names(p))
  if (length(unknown)) cli_abort("Unknown parameter name{?s}: {.val {unknown}}.")
  if (!is.null(thresholds)) {
    if (any(paste0("threshold_", seq_len(3L)) %in% names(parameters))) {
      cli_abort("Supply thresholds in {.arg thresholds} or {.arg parameters}, not both.")
    }
    if (is.numeric(thresholds)) thresholds <- rep(list(thresholds), 3L)
    if (!is.list(thresholds) || length(thresholds) != 3L) {
      cli_abort("{.arg thresholds} must be one vector or a list of three vectors.")
    }
    for (j in seq_len(3L)) p[paste0("threshold_", j)] <- thresholds[j]
  }
  for (name in names(parameters)) p[name] <- parameters[name]
  matrix_dims <- list(item_sd = c(2L, 3L))
  if (!regression) {
    matrix_dims <- c(matrix_dims, list(mu_q = c(2L, 2L), sigma_q = c(2L, 2L),
                                      B = c(2L, 2L), D = c(2L, 2L)))
  }
  for (name in names(p)) {
    value <- p[[name]]
    if (!is.numeric(value) || any(!is.finite(value))) {
      cli_abort("Parameter {.val {name}} must contain finite numeric values.")
    }
    if (name %in% names(matrix_dims)) {
      if (!is.matrix(value) || !identical(dim(value), matrix_dims[[name]])) {
        cli_abort("Parameter {.val {name}} has an invalid matrix size.")
      }
    } else if (startsWith(name, "threshold_")) {
      if (!is.null(dim(value)) || !length(value) %in% 2:4 || any(diff(value) <= 0)) {
        cli_abort("Parameter {.val {name}} needs two to four strictly increasing thresholds.")
      }
    } else {
      expected <- if (name == "lambda") 3L else 2L
      if (!is.null(dim(value)) || length(value) != expected) {
        cli_abort("Parameter {.val {name}} must have length {expected}.")
      }
    }
    if ((startsWith(name, "sigma_") || name == "item_sd") && any(value <= 0)) {
      cli_abort("Parameter {.val {name}} must contain positive SDs.")
    }
    if (name == "rho_q" && any(abs(value) >= 1)) {
      cli_abort("Parameter {.val {name}} must lie strictly between -1 and 1.")
    }
    p[[name]] <- unname(value)
  }
  if (p$lambda[1] != 1) cli_abort("The first loading must equal one.")
  if (any(p$item_sd[1, ] != 1)) cli_abort("Each control item SD must equal one.")
  check_anchor <- function(name, current, required) {
    if (name %in% names(parameters) && !isTRUE(all.equal(current, required, tolerance = 1e-12))) {
      cli_abort("Parameter {.val {name}} conflicts with the control factor mean of zero.")
    }
  }
  if (model == "latent_covariate") {
    check_anchor("mu_v", p$mu_v[1], 0)
    p$mu_v[1] <- 0
  } else if (model == "latent_outcome") {
    anchor <- -p$beta[1] * p$mu_v[1]
    check_anchor("alpha", p$alpha[1], anchor)
    p$alpha[1] <- anchor
  } else if (model == "latent_mediator") {
    anchor <- -sum(p$B[1, ] * p$mu_q[1, ])
    check_anchor("alpha_m", p$alpha_m[1], anchor)
    p$alpha_m[1] <- anchor
  } else {
    check_anchor("mu_q", p$mu_q[1, 1], 0)
    p$mu_q[1, 1] <- 0
  }
  p
}

#' Internal: calculate true effects without using the posterior extractor
#' @noRd
.causal_sem_sim_effects <- function(model, parameters, weights) {
  p <- parameters
  if (model %in% c("latent_covariate", "latent_outcome")) {
    covariate_mean <- sum(weights * p$mu_v)
    responses <- p$alpha + p$beta * covariate_mean
    return(c(ATE = responses[2] - responses[1], interaction = p$beta[2] - p$beta[1]))
  }
  target_q <- as.vector(crossprod(weights, p$mu_q))
  actual_m <- p$alpha_m + rowSums(p$B * p$mu_q)
  target_m <- sum(weights * actual_m)
  total <- direct <- numeric(2L)
  for (g in seq_len(2L)) {
    model_m <- p$alpha_m[g] + sum(p$B[g, ] * target_q)
    baseline_y <- p$alpha_y[g] + sum(p$D[g, ] * target_q)
    total[g] <- baseline_y + p$d[g] * model_m
    direct[g] <- baseline_y + p$d[g] * target_m
  }
  ate <- total[2] - total[1]
  ade <- direct[2] - direct[1]
  c(ATE = ate, ADE = ade, AIE = ate - ade)
}
