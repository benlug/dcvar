#' Prepare data for a causal SEM model
#'
#' Validate roles and ordinal indicators. Each row describes one person.
#' This function keeps all rows. It uses both groups to standardize each
#' manifest continuous variable when requested.
#'
#' @param data A data frame with complete observations for all model columns.
#' @param model One of `"latent_covariate"`, `"latent_outcome"`,
#'   `"latent_mediator"`, or `"latent_mediator_baseline"`.
#' @param treatment Name of the column that identifies the two groups.
#' @param control,treated Distinct values that identify the control and
#'   treatment groups.
#' @param outcome Name of the outcome role.
#' @param covariate Name of the covariate role for the first two models.
#' @param mediator Name of the mediator role for the mediation models.
#' @param mediator_baseline,outcome_baseline Names of the baseline roles
#'   for the mediation models.
#' @param indicators A named list with one element. Its name is the latent
#'   role name. Its value names three ordered factor columns. Each column
#'   must declare three to five levels.
#' @param target_weights Numeric vector with control and treatment weights.
#'   The values must be non-negative and sum to one.
#'   A sum within `1e-8` of one is accepted. The weights are then divided by
#'   their sum to remove roundoff.
#' @param standardize Logical. Standardize manifest continuous variables
#'   with their pooled sample mean and sample SD.
#' @param priors Named list of prior settings. Supported names are `beta_sd`,
#'   `intercept_sd`, `scale_sd`, `loading_mean`, `loading_sd`, `threshold_sd`,
#'   `item_log_sd`, and `lkj_shape`. All values except `loading_mean` must
#'   be positive. The defaults are 1, 2.5, 1.5, 1, 1, 2.5, 0.5, and 2.
#' @param prior_only Logical. Exclude the observed data likelihood when true.
#' @return A list with `stan_data` and `meta`. The metadata stores roles,
#'   levels, groups, weights, priors, and each role's `center` and `scale`.
#' @export
prepare_causal_sem_data <- function(data, model, treatment, control = 0,
                                   treated = 1, outcome, covariate = NULL,
                                   mediator = NULL, mediator_baseline = NULL,
                                   outcome_baseline = NULL, indicators,
                                   target_weights = c(0.5, 0.5),
                                   standardize = TRUE, priors = list(),
                                   prior_only = FALSE) {
  model <- .causal_sem_validate_model(model)
  priors <- .causal_sem_priors(priors)
  if (!is.data.frame(data) || nrow(data) < 2L) {
    cli_abort("{.arg data} must be a data frame with at least two rows.")
  }
  if (anyNA(names(data)) || any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    cli_abort("Data columns must have unique non-empty names.")
  }
  data <- as.data.frame(data)
  .causal_sem_validate_flag(standardize, "standardize")
  .causal_sem_validate_flag(prior_only, "prior_only")
  target_weights <- .causal_sem_validate_weights(target_weights)
  .causal_sem_validate_name(treatment, "treatment")
  roles <- list(outcome = outcome, covariate = covariate, mediator = mediator,
                mediator_baseline = mediator_baseline, outcome_baseline = outcome_baseline)
  regression <- model %in% c("latent_covariate", "latent_outcome")
  required <- if (regression) c("outcome", "covariate") else {
    c("outcome", "mediator", "mediator_baseline", "outcome_baseline")
  }
  for (role in names(roles)) {
    if (role %in% required) {
      .causal_sem_validate_name(roles[[role]], role)
    } else if (!is.null(roles[[role]])) {
      cli_abort("Role {.val {role}} is not used by model {.val {model}}.")
    }
  }
  role_names <- unlist(roles[required], use.names = FALSE)
  if (anyDuplicated(c(treatment, role_names))) {
    cli_abort("Treatment and model roles must use distinct names.")
  }
  latent_role <- switch(model,
    latent_covariate = "covariate", latent_outcome = "outcome",
    latent_mediator = "mediator", latent_mediator_baseline = "mediator_baseline"
  )
  latent_name <- roles[[latent_role]]
  if (!is.list(indicators) || length(indicators) != 1L ||
      !identical(names(indicators), latent_name)) {
    cli_abort("{.arg indicators} must be a list named {.val {latent_name}}.")
  }
  item_names <- indicators[[1]]
  if (!is.character(item_names) || length(item_names) != 3L || anyNA(item_names) ||
      any(!nzchar(item_names)) || anyDuplicated(item_names)) {
    cli_abort("{.arg indicators} must name three distinct columns.")
  }
  if (any(item_names %in% c(treatment, role_names))) {
    cli_abort("Indicator columns cannot also identify treatment or model roles.")
  }
  if (latent_name %in% names(data)) {
    cli_abort("Latent role {.val {latent_name}} must not also be a data column.")
  }
  manifest_roles <- setdiff(required, latent_role)
  model_columns <- c(treatment, unlist(roles[manifest_roles], use.names = FALSE), item_names)
  missing_columns <- setdiff(model_columns, names(data))
  if (length(missing_columns)) {
    cli_abort("Model column{?s} not found: {.val {missing_columns}}.")
  }
  for (label in list(control, treated)) {
    if (!is.atomic(label) || length(label) != 1L || is.na(label) ||
        (is.numeric(label) && !is.finite(label))) {
      cli_abort("{.arg control} and {.arg treated} must be finite single values.")
    }
  }
  group_values <- data[[treatment]]
  if (!is.atomic(group_values) || !is.null(dim(group_values)) || anyNA(group_values)) {
    cli_abort("The treatment column must be a complete vector.")
  }
  control_key <- as.character(control)
  treated_key <- as.character(treated)
  if (identical(control_key, treated_key)) {
    cli_abort("{.arg control} and {.arg treated} must differ.")
  }
  group <- match(as.character(group_values), c(control_key, treated_key))
  if (anyNA(group)) cli_abort("The treatment column contains an unknown group value.")
  n_group <- stats::setNames(tabulate(group, nbins = 2L), c("control", "treated"))
  if (any(n_group == 0L)) cli_abort("Both treatment groups must have observations.")
  item_levels <- stats::setNames(vector("list", 3L), item_names)
  u <- matrix(0L, nrow(data), 3L)
  K <- integer(3)
  for (j in seq_len(3L)) {
    item <- data[[item_names[j]]]
    if (!is.ordered(item)) {
      cli_abort("Indicator {.val {item_names[j]}} must be an ordered factor.")
    }
    if (anyNA(item)) cli_abort("Indicator {.val {item_names[j]}} contains missing values.")
    K[j] <- nlevels(item)
    if (K[j] < 3L || K[j] > 5L) {
      cli_abort("Indicator {.val {item_names[j]}} must declare three to five levels.")
    }
    if (length(unique(item)) < 2L) {
      cli_abort("Indicator {.val {item_names[j]}} is constant across all persons.")
    }
    item_levels[[j]] <- levels(item)
    u[, j] <- as.integer(item)
  }
  stan_data <- list(N = nrow(data), group = as.integer(group), J = 3L, K = K, u = u)
  role_fields <- c(outcome = "y", covariate = "z", mediator = "m",
                   mediator_baseline = "q1", outcome_baseline = "q2")
  for (field in unname(role_fields)) stan_data[[field]] <- numeric(nrow(data))
  scales <- stats::setNames(vector("list", length(required)), required)
  for (role in required) {
    scales[[role]] <- list(center = 0, scale = 1)
    if (role == latent_role) next
    column <- roles[[role]]
    value <- data[[column]]
    if (!is.numeric(value) || !is.null(dim(value)) || any(!is.finite(value))) {
      cli_abort("Role {.val {role}} must use a complete finite numeric column.")
    }
    if (standardize) {
      center <- mean(value)
      scale <- stats::sd(value)
      if (!is.finite(center) || !is.finite(scale) || scale <= 0) {
        cli_abort("Role {.val {role}} needs a positive finite SD for standardization.")
      }
      scales[[role]] <- list(center = center, scale = scale)
      value <- (value - center) / scale
    }
    stan_data[[role_fields[[role]]]] <- as.numeric(value)
  }
  stan_data <- c(stan_data, .causal_sem_stan_priors(priors), list(
    prior_only = as.integer(prior_only), target_weights = target_weights
  ))
  list(stan_data = stan_data, meta = list(
    model = model, roles = roles, indicators = indicators,
    group_labels = list(control = control, treated = treated),
    target_weights = target_weights, scales = scales, levels = item_levels,
    n_group = n_group, standardized = standardize, priors = priors
  ))
}

#' Internal: validate a causal SEM role name
#' @noRd
.causal_sem_validate_name <- function(value, role) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    cli_abort("Role {.val {role}} must be one non-empty column or latent name.")
  }
}

#' Internal: validate a causal SEM flag
#' @noRd
.causal_sem_validate_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    cli_abort("{.arg {name}} must be TRUE or FALSE.")
  }
}

#' Internal: validate target population weights
#' @noRd
.causal_sem_validate_weights <- function(target_weights) {
  if (!is.numeric(target_weights) || !is.null(dim(target_weights)) ||
      length(target_weights) != 2L ||
      any(!is.finite(target_weights)) || any(target_weights < 0) ||
      abs(sum(target_weights) - 1) > 1e-8) {
    cli_abort("{.arg target_weights} must contain two non-negative weights that sum to one.")
  }
  unname(as.numeric(target_weights / sum(target_weights)))
}
