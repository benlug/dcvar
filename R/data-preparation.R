# ============================================================================
# Data Preparation Functions
# ============================================================================

#' Internal: coerce time-like vectors to a numeric ordering scale
#' @noRd
.time_to_numeric <- function(x) {
  if (inherits(x, "POSIXt") || inherits(x, "Date") || inherits(x, "difftime")) {
    return(as.numeric(x))
  }
  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  xtfrm(x)
}

#' Internal: validate a positive finite scalar
#' @noRd
.prep_validate_positive_scalar <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    cli_abort("{.arg {arg_name}} must be a single positive finite numeric value.")
  }
}

#' Internal: validate a scalar logical flag
#' @noRd
.prep_validate_scalar_logical <- function(x, arg_name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    cli_abort("{.arg {arg_name}} must be a single logical value.")
  }
}

#' Internal: validate a finite numeric vector
#' @noRd
.prep_validate_numeric_vector <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x))) {
    cli_abort("{.arg {arg_name}} must be a non-empty finite numeric vector.")
  }
}

#' Internal: validate that a variable vector has no duplicates
#' @noRd
.prep_validate_unique_vars <- function(vars, context = "vars") {
  if (anyDuplicated(vars)) {
    cli_abort("{.arg {context}} must contain distinct variable names.")
  }
}

#' Internal: validate ordered time values for VAR-style models
#' @noRd
.validate_time_values <- function(time_values, allow_gaps = FALSE, context = "data") {
  if (anyNA(time_values)) {
    cli_abort("{.arg time_var} contains missing values in {context}.")
  }
  if (length(time_values) < 2) {
    return(invisible(time_values))
  }
  if (anyDuplicated(time_values)) {
    cli_abort(c(
      "Duplicate {.arg time_var} values found in {context}.",
      "i" = "VAR-style models require a unique observation time for each row."
    ))
  }

  steps <- diff(.time_to_numeric(time_values))
  if (any(!is.finite(steps)) || any(steps <= 0)) {
    cli_abort(c(
      "{.arg time_var} in {context} could not be ordered into a strictly increasing sequence.",
      "i" = "Use a numeric, Date, POSIXct, or otherwise orderable time index."
    ))
  }

  regular_steps <- isTRUE(all.equal(
    steps,
    rep(steps[[1]], length(steps)),
    tolerance = sqrt(.Machine$double.eps)
  ))
  if (!regular_steps) {
    if (!allow_gaps) {
      cli_abort(c(
        "{.arg time_var} in {context} is not evenly spaced.",
        "!" = "The VAR(1) model assumes consecutive, equally spaced observations.",
        "i" = "Reindex the series to consecutive time steps or set {.arg allow_gaps = TRUE} to proceed with a warning."
      ))
    }
    cli_warn(c(
      "{.arg time_var} in {context} is not evenly spaced.",
      "!" = "The VAR(1) model still treats adjacent rows as consecutive observations.",
      "i" = "Interpret lagged quantities with caution."
    ))
  }

  invisible(time_values)
}

#' Internal: derive the common multilevel time grid
#' @noRd
.validate_panel_time_grid <- function(data, ids, id_var, time_var) {
  time_grid <- NULL

  for (id in ids) {
    unit_times <- data[data[[id_var]] == id, time_var]
    .validate_time_values(unit_times, allow_gaps = FALSE, context = paste0("unit ", id))

    if (is.null(time_grid)) {
      time_grid <- unit_times
    } else {
      same_grid <- isTRUE(all.equal(time_grid, unit_times, check.attributes = TRUE))
      if (!same_grid) {
        cli_abort(c(
          "All units must share the same {.arg time_var} grid.",
          "!" = "The multilevel model assumes aligned measurement occasions across units."
        ))
      }
    }
  }

  time_grid
}

#' Internal: shared data preparation logic for all copula VAR models
#'
#' @param data A data frame.
#' @param vars Character vector of exactly two variable names.
#' @param time_var Name of the time column.
#' @param standardize Logical scalar; whether to z-score the variables.
#' @param allow_gaps Logical scalar; if `FALSE` (default), interior missing values
#'   cause an error. If `TRUE`, they produce a warning and are removed.
#'
#' @return A list with elements `Y`, `T`, `D`, and standardization metadata.
#' @noRd
.prepare_var_data <- function(data, vars, time_var = "time",
                              standardize = TRUE, allow_gaps = FALSE) {
  # Validate inputs
  if (!is.data.frame(data)) {
    cli_abort("{.arg data} must be a data frame.")
  }
  if (nrow(data) == 0) {
    cli_abort("{.arg data} must not be empty.")
  }
  if (length(vars) != 2) {
    cli_abort("Exactly 2 variables required (bivariate model). Got {.val {length(vars)}}.")
  }
  .prep_validate_unique_vars(vars)
  .prep_validate_scalar_logical(standardize, "standardize")
  .prep_validate_scalar_logical(allow_gaps, "allow_gaps")
  missing_vars <- setdiff(c(vars, time_var), names(data))
  if (length(missing_vars) > 0) {
    cli_abort("Column{?s} not found in data: {.val {missing_vars}}")
  }
  for (v in vars) {
    if (!is.numeric(data[[v]])) {
      cli_abort("Column {.val {v}} must be numeric, got {.cls {class(data[[v]])[1]}}.")
    }
  }

  # Sort by time
  data <- data[order(data[[time_var]]), , drop = FALSE]
  time_values <- data[[time_var]]
  .validate_time_values(time_values, allow_gaps = allow_gaps)

  # Extract variables
  Y <- as.matrix(data[, vars])

  # Check for non-finite values
  if (any(is.nan(Y) | is.infinite(Y))) {
    cli_abort("{.arg data} contains NaN or Inf values. Remove or impute them first.")
  }

  # Handle missing values
  cc <- complete.cases(Y)
  if (sum(!cc) > 0) {
    n_missing <- sum(!cc)
    missing_idx <- which(!cc)
    interior <- missing_idx[missing_idx > 1 & missing_idx < nrow(Y)]

    if (length(interior) > 0) {
      if (!allow_gaps) {
        cli_abort(c(
          "Found {n_missing} interior missing value{?s} that would break time series adjacency.",
          "!" = "The VAR(1) model assumes consecutive observations.",
          "i" = "Consider imputation or explicit missing-data handling.",
          "i" = "Set {.arg allow_gaps = TRUE} to remove these rows with a warning instead."
        ))
      }
      cli_warn(c(
        "Removing {n_missing} row{?s} with missing values.",
        "!" = "This breaks time series adjacency: the VAR(1) model assumes consecutive observations.",
        "i" = "Consider imputation or explicit missing-data handling for real data."
      ))
    } else {
      cli_warn("Removing {n_missing} row{?s} with missing values (edges only).")
    }
    Y <- Y[cc, ]
    time_values <- time_values[cc]
  }

  # Require at least 3 observations for a VAR(1) model
  if (nrow(Y) < 3) {
    cli_abort("At least 3 observations required for a VAR(1) model, got {.val {nrow(Y)}}.")
  }

  # Standardize if requested
  Y_means <- NULL
  Y_sds <- NULL
  if (standardize) {
    Y_means <- colMeans(Y)
    Y_sds <- apply(Y, 2, sd)
    if (any(Y_sds < 1e-10)) {
      cli_abort("Cannot standardize: one or more variables have zero (or near-zero) variance.")
    }
    Y <- scale(Y)
  }

  list(
    Y = Y,
    T_obs = nrow(Y),
    D = ncol(Y),
    standardized = standardize,
    time_values = time_values,
    Y_means = Y_means,
    Y_sds = Y_sds,
    vars = vars
  )
}


#' Prepare data for the DC-VAR model
#'
#' Transforms a data frame into a list suitable for the DC-VAR Stan model.
#' Handles sorting, missing values, and optional standardization.
#'
#' @section Prior naming conventions:
#' Parameters ending in `_sd` specify normal prior standard deviations
#' (location parameters). Parameters ending in `_rate` specify exponential
#' prior means (scale parameters), where the exponential rate is
#' `1/prior_*_rate`. The constant and HMM models use `prior_z_rho_sd`
#' (normal prior on the Fisher-z scale), while the DC-VAR model uses
#' `prior_sigma_omega_rate` (exponential prior on the random-walk SD)
#' because the two quantities have fundamentally different roles.
#'
#' @param data A data frame with time series observations.
#' @param vars Character vector of two variable names to model.
#' @param time_var Name of the time column (default: `"time"`).
#' @param standardize Logical; whether to z-score variables (default: `TRUE`).
#' @param margins Character string specifying the marginal distribution.
#'   One of `"normal"` (default), `"exponential"`, `"skew_normal"`, or `"gamma"`.
#' @param skew_direction Integer vector of length D indicating skew direction
#'   for asymmetric margins. Each element must be `1` (right-skewed) or `-1`
#'   (left-skewed). Required for `"exponential"` and `"gamma"` margins.
#' @param prior_mu_sd Prior SD for intercepts: `mu ~ normal(0, prior_mu_sd)`.
#' @param prior_phi_sd Prior SD for VAR coefficients: `Phi ~ normal(0, prior_phi_sd)`.
#' @param prior_sigma_eps_rate Prior mean for innovation SDs:
#'   `sigma_eps ~ exponential(1/prior_sigma_eps_rate)`. Default `1` gives
#'   `exponential(1)` with prior mean 1.
#' @param prior_sigma_omega_rate Prior mean for rho process SD:
#'   `sigma_omega ~ exponential(1/prior_sigma_omega_rate)`. Default `0.1` gives
#'   `exponential(10)` with prior mean 0.1.
#' @param prior_rho_init_sd Prior SD for initial rho on Fisher-z scale.
#' @param allow_gaps Logical; if `FALSE` (default), interior missing values
#'   cause an error. If `TRUE`, they produce a warning and are removed.
#'
#' @return A named list suitable as Stan data input.
#' @export
prepare_dcvar_data <- function(data, vars, time_var = "time",
                               standardize = TRUE,
                               margins = "normal",
                               skew_direction = NULL,
                               prior_mu_sd = 2,
                               prior_phi_sd = 0.5,
                               prior_sigma_eps_rate = 1,
                               prior_sigma_omega_rate = 0.1,
                               prior_rho_init_sd = 1,
                               allow_gaps = FALSE) {
  .validate_margins(margins, skew_direction)
  .prep_validate_positive_scalar(prior_mu_sd, "prior_mu_sd")
  .prep_validate_positive_scalar(prior_phi_sd, "prior_phi_sd")
  .prep_validate_positive_scalar(prior_sigma_eps_rate, "prior_sigma_eps_rate")
  .prep_validate_positive_scalar(prior_sigma_omega_rate, "prior_sigma_omega_rate")
  .prep_validate_positive_scalar(prior_rho_init_sd, "prior_rho_init_sd")
  prep <- .prepare_var_data(data, vars, time_var, standardize, allow_gaps)

  stan_data <- list(
    T = prep$T_obs,
    D = prep$D,
    Y = prep$Y,
    sigma_mu_prior = prior_mu_sd,
    sigma_phi_prior = prior_phi_sd,
    sigma_omega_prior = prior_sigma_omega_rate,
    rho_init_prior_sd = prior_rho_init_sd
  )

  # Normal margins need sigma_eps_prior; non-Gaussian do not
  if (margins == "normal") {
    stan_data$sigma_eps_prior <- prior_sigma_eps_rate
  }

  # Add margin-specific Stan data
  if (margins %in% c("exponential", "gamma")) {
    if (is.null(skew_direction)) {
      cli_abort(c(
        "Margin {.val {margins}} requires {.arg skew_direction}.",
        "i" = "{.arg skew_direction} must be a numeric vector of length {.val {prep$D}} with values 1 or -1."
      ))
    }
    stan_data$skew_direction <- as.numeric(skew_direction)
  }

  attr(stan_data, "vars") <- prep$vars
  attr(stan_data, "standardized") <- prep$standardized
  attr(stan_data, "margins") <- margins
  attr(stan_data, "time_values") <- prep$time_values
  if (!is.null(skew_direction)) attr(stan_data, "skew_direction") <- skew_direction
  if (prep$standardized) {
    attr(stan_data, "Y_means") <- prep$Y_means
    attr(stan_data, "Y_sds") <- prep$Y_sds
  }

  stan_data
}


#' Prepare data for the HMM copula model
#'
#' Transforms a data frame into a list suitable for the HMM copula Stan model.
#' Includes HMM-specific prior hyperparameters.
#'
#' @inheritParams prepare_dcvar_data
#' @param K Number of hidden states (default: 2).
#' @param prior_kappa Sticky Dirichlet self-transition concentration (default: 10).
#' @param prior_alpha_off Sticky Dirichlet off-diagonal concentration (default: 1).
#' @param prior_z_rho_sd Prior SD for state-specific z_rho values (default: 1.0).
#'
#' @return A named list suitable as Stan data input.
#' @export
prepare_hmm_data <- function(data, vars, K = 2, time_var = "time",
                             standardize = TRUE,
                             margins = "normal",
                             skew_direction = NULL,
                             prior_mu_sd = 2,
                             prior_phi_sd = 0.5,
                             prior_sigma_eps_rate = 1,
                             prior_kappa = 10,
                             prior_alpha_off = 1,
                             prior_z_rho_sd = 1.0,
                             allow_gaps = FALSE) {
  .validate_margins(margins, skew_direction)
  .prep_validate_positive_scalar(prior_mu_sd, "prior_mu_sd")
  .prep_validate_positive_scalar(prior_phi_sd, "prior_phi_sd")
  .prep_validate_positive_scalar(prior_sigma_eps_rate, "prior_sigma_eps_rate")
  .prep_validate_positive_scalar(prior_kappa, "prior_kappa")
  .prep_validate_positive_scalar(prior_alpha_off, "prior_alpha_off")
  .prep_validate_positive_scalar(prior_z_rho_sd, "prior_z_rho_sd")
  if (!is.numeric(K) || length(K) != 1 || K < 2 || K != as.integer(K)) {
    cli_abort("{.arg K} must be an integer >= 2, got {.val {K}}.")
  }
  prep <- .prepare_var_data(data, vars, time_var, standardize, allow_gaps)

  stan_data <- list(
    T = prep$T_obs,
    D = prep$D,
    Y = prep$Y,
    K = as.integer(K),
    sigma_mu_prior = prior_mu_sd,
    sigma_phi_prior = prior_phi_sd,
    kappa = prior_kappa,
    alpha_off = prior_alpha_off,
    z_rho_prior_sd = prior_z_rho_sd
  )

  if (margins == "normal") {
    stan_data$sigma_eps_prior <- prior_sigma_eps_rate
  }

  if (margins %in% c("exponential", "gamma")) {
    if (is.null(skew_direction)) {
      cli_abort(c(
        "Margin {.val {margins}} requires {.arg skew_direction}.",
        "i" = "{.arg skew_direction} must be a numeric vector of length {.val {prep$D}} with values 1 or -1."
      ))
    }
    stan_data$skew_direction <- as.numeric(skew_direction)
  }

  attr(stan_data, "vars") <- prep$vars
  attr(stan_data, "standardized") <- prep$standardized
  attr(stan_data, "margins") <- margins
  attr(stan_data, "time_values") <- prep$time_values
  if (!is.null(skew_direction)) attr(stan_data, "skew_direction") <- skew_direction
  if (prep$standardized) {
    attr(stan_data, "Y_means") <- prep$Y_means
    attr(stan_data, "Y_sds") <- prep$Y_sds
  }

  stan_data
}


#' Prepare data for the multilevel copula VAR model
#'
#' @param data A data frame in long (panel) format.
#' @param vars Character vector of two variable names.
#' @param id_var Name of the unit/person ID column.
#' @param time_var Name of the time column.
#' @param center Logical scalar; person-mean center (default: `TRUE`).
#' @param prior_phi_bar_sd Prior SD for phi_bar.
#' @param prior_tau_phi_scale Prior scale for tau_phi.
#' @param prior_sigma_sd Prior SD for sigma.
#' @param prior_rho_sd Prior SD for rho.
#'
#' @return A named list suitable as Stan data input.
#' @export
prepare_multilevel_data <- function(data, vars, id_var = "id",
                                    time_var = "time", center = TRUE,
                                    prior_phi_bar_sd = 0.5,
                                    prior_tau_phi_scale = 0.2,
                                    prior_sigma_sd = 1,
                                    prior_rho_sd = 0.5) {
  if (!is.data.frame(data)) cli_abort("{.arg data} must be a data frame.")
  .prep_validate_scalar_logical(center, "center")
  if (length(vars) != 2) {
    cli_abort("Exactly 2 variables required for multilevel models. Got {.val {length(vars)}}.")
  }
  .prep_validate_unique_vars(vars)
  .prep_validate_positive_scalar(prior_phi_bar_sd, "prior_phi_bar_sd")
  .prep_validate_positive_scalar(prior_tau_phi_scale, "prior_tau_phi_scale")
  .prep_validate_positive_scalar(prior_sigma_sd, "prior_sigma_sd")
  .prep_validate_positive_scalar(prior_rho_sd, "prior_rho_sd")
  missing_cols <- setdiff(c(vars, id_var, time_var), names(data))
  if (length(missing_cols) > 0) {
    cli_abort("Column{?s} not found: {.val {missing_cols}}")
  }
  for (v in vars) {
    if (!is.numeric(data[[v]])) {
      cli_abort("Column {.val {v}} must be numeric, got {.cls {class(data[[v]])[1]}}.")
    }
  }

  id_values <- data[[id_var]]
  if (anyNA(id_values)) {
    cli_abort(c(
      "{.arg id_var} contains missing values.",
      "i" = "Remove or impute missing unit IDs before fitting the multilevel model."
    ))
  }
  if (is.factor(id_values)) {
    id_values <- droplevels(id_values)
    data[[id_var]] <- id_values
  }

  # Split by unit
  ids <- unique(id_values)
  N <- length(ids)

  # Sort within each unit
  data <- data[order(data[[id_var]], data[[time_var]]), , drop = FALSE]

  # Check balanced panel
  counts <- table(id_values)
  if (length(unique(counts)) != 1) {
    cli_abort("Unbalanced panels not yet supported. All units must have the same number of observations.")
  }
  T_obs <- as.integer(unique(counts))
  time_grid <- .validate_panel_time_grid(data, ids, id_var, time_var)

  # Build 3D array
  y_array <- array(NA_real_, dim = c(N, T_obs, 2))
  person_means <- matrix(NA_real_, N, 2)

  for (i in seq_len(N)) {
    unit_data <- data[data[[id_var]] == ids[i], ]
    y_i <- as.matrix(unit_data[, vars])

    # Check for non-finite values
    if (any(is.nan(y_i) | is.infinite(y_i))) {
      cli_abort(c(
        "Unit {.val {ids[i]}} contains NaN or Inf values.",
        "i" = "Multilevel models require complete cases. Remove or impute non-finite values first."
      ))
    }

    # Check for NAs
    if (any(is.na(y_i))) {
      cli_abort(c(
        "Unit {.val {ids[i]}} contains missing values (NA).",
        "i" = "Multilevel models require complete cases. Consider imputation before fitting."
      ))
    }

    if (center) {
      person_means[i, ] <- colMeans(y_i)
      y_i <- sweep(y_i, 2, person_means[i, ])
    }

    y_array[i, , ] <- y_i
  }

  # Convert to list of matrices (cmdstanr array format)
  y_list <- lapply(seq_len(N), function(i) y_array[i, , ])

  stan_data <- list(
    N = N,
    T = T_obs,
    y = y_list,
    prior_phi_bar_sd = prior_phi_bar_sd,
    prior_tau_phi_scale = prior_tau_phi_scale,
    prior_sigma_sd = prior_sigma_sd,
    prior_rho_sd = prior_rho_sd
  )

  attr(stan_data, "vars") <- vars
  attr(stan_data, "person_means") <- person_means
  attr(stan_data, "ids") <- ids
  attr(stan_data, "time_values") <- time_grid

  stan_data
}


#' Prepare data for the SEM copula VAR model
#'
#' Transforms a data frame of indicator variables into a list suitable for
#' the SEM copula Stan model. The measurement model parameters (lambda,
#' sigma_e) are fixed and passed through to Stan.
#'
#' @param data A data frame with time series of indicator variables.
#' @param indicators A list of two character vectors, each naming J indicator
#'   columns per latent variable.
#' @param J Number of indicators per latent variable.
#' @param lambda Numeric vector of length J with fixed factor loadings.
#' @param sigma_e Fixed measurement error SD (scalar).
#' @param time_var Name of the time column (default: `"time"`).
#' @param prior_mu_sd Prior SD for intercepts.
#' @param prior_phi_sd Prior SD for VAR coefficients.
#' @param prior_sigma_sd Prior SD for lognormal on innovation SDs.
#' @param prior_rho_sd Prior SD for rho_raw.
#'
#' @return A named list suitable as Stan data input.
#' @export
prepare_sem_data <- function(data, indicators, J, lambda, sigma_e,
                              time_var = "time",
                              prior_mu_sd = 0.25,
                              prior_phi_sd = 0.5,
                              prior_sigma_sd = 0.5,
                              prior_rho_sd = 0.75) {
  if (!is.data.frame(data)) {
    cli_abort("{.arg data} must be a data frame.")
  }
  if (!is.list(indicators) || length(indicators) != 2) {
    cli_abort("{.arg indicators} must be a list of two character vectors.")
  }
  if (!is.numeric(J) || length(J) != 1L || !is.finite(J) || J < 1 || J != as.integer(J)) {
    cli_abort("{.arg J} must be a positive integer.")
  }
  if (length(indicators[[1]]) != J || length(indicators[[2]]) != J) {
    cli_abort("Each element of {.arg indicators} must have {.val {J}} indicator names.")
  }
  .prep_validate_numeric_vector(lambda, "lambda")
  if (length(lambda) != J) {
    cli_abort("{.arg lambda} must have length {.val {J}}.")
  }
  .prep_validate_positive_scalar(sigma_e, "sigma_e")
  .prep_validate_positive_scalar(prior_mu_sd, "prior_mu_sd")
  .prep_validate_positive_scalar(prior_phi_sd, "prior_phi_sd")
  .prep_validate_positive_scalar(prior_sigma_sd, "prior_sigma_sd")
  .prep_validate_positive_scalar(prior_rho_sd, "prior_rho_sd")

  all_ind <- c(indicators[[1]], indicators[[2]])
  if (anyDuplicated(all_ind)) {
    cli_abort(c(
      "Indicator names must be unique across both latent variables.",
      "i" = "No indicator may appear in both measurement blocks."
    ))
  }
  missing_cols <- setdiff(c(all_ind, time_var), names(data))
  if (length(missing_cols) > 0) {
    cli_abort("Column{?s} not found in data: {.val {missing_cols}}")
  }
  for (col in all_ind) {
    if (!is.numeric(data[[col]])) {
      cli_abort("Indicator column {.val {col}} must be numeric, got {.cls {class(data[[col]])[1]}}.")
    }
  }

  # Sort by time
  data <- data[order(data[[time_var]]), , drop = FALSE]
  time_values <- data[[time_var]]
  .validate_time_values(time_values, allow_gaps = FALSE, context = "SEM data")

  # Build T x 2J matrix (latent 1 indicators first, then latent 2)
  y <- as.matrix(data[, all_ind])
  T_obs <- nrow(y)
  if (any(is.nan(y) | is.infinite(y))) {
    cli_abort("{.arg data} contains NaN or Inf values in the indicator columns.")
  }
  if (any(is.na(y))) {
    cli_abort("{.arg data} contains missing values in the indicator columns.")
  }

  if (T_obs < 3) {
    cli_abort("At least 3 observations required, got {.val {T_obs}}.")
  }

  stan_data <- list(
    T = T_obs,
    J = as.integer(J),
    y = y,
    lambda = lambda,
    sigma_e = sigma_e,
    prior_mu_sd = prior_mu_sd,
    prior_phi_sd = prior_phi_sd,
    prior_sigma_sd = prior_sigma_sd,
    prior_rho_sd = prior_rho_sd
  )

  # Use list names as latent variable names, or default
  latent_names <- names(indicators)
  if (is.null(latent_names) || any(latent_names == "")) {
    latent_names <- c("latent1", "latent2")
  }

  attr(stan_data, "vars") <- latent_names
  attr(stan_data, "indicators") <- indicators
  attr(stan_data, "time_values") <- time_values

  stan_data
}


#' Prepare data for the constant copula model
#'
#' Transforms a data frame into a list suitable for the constant copula Stan model.
#'
#' @inheritParams prepare_dcvar_data
#' @param prior_z_rho_sd Prior SD for rho on Fisher-z scale (default: 1.0).
#'
#' @return A named list suitable as Stan data input.
#' @export
prepare_constant_data <- function(data, vars, time_var = "time",
                                  standardize = TRUE,
                                  margins = "normal",
                                  skew_direction = NULL,
                                  prior_mu_sd = 2,
                                  prior_phi_sd = 0.5,
                                  prior_sigma_eps_rate = 1,
                                  prior_z_rho_sd = 1.0,
                                  allow_gaps = FALSE) {
  .prep_validate_unique_vars(vars)
  .prep_validate_positive_scalar(prior_mu_sd, "prior_mu_sd")
  .prep_validate_positive_scalar(prior_phi_sd, "prior_phi_sd")
  .prep_validate_positive_scalar(prior_sigma_eps_rate, "prior_sigma_eps_rate")
  .prep_validate_positive_scalar(prior_z_rho_sd, "prior_z_rho_sd")
  .validate_margins(margins, skew_direction)
  prep <- .prepare_var_data(data, vars, time_var, standardize, allow_gaps)

  stan_data <- list(
    T = prep$T_obs,
    D = prep$D,
    Y = prep$Y,
    sigma_mu_prior = prior_mu_sd,
    sigma_phi_prior = prior_phi_sd,
    z_rho_prior_sd = prior_z_rho_sd
  )

  if (margins == "normal") {
    stan_data$sigma_eps_prior <- prior_sigma_eps_rate
  }

  if (margins %in% c("exponential", "gamma")) {
    if (is.null(skew_direction)) {
      cli_abort(c(
        "Margin {.val {margins}} requires {.arg skew_direction}.",
        "i" = "{.arg skew_direction} must be a numeric vector of length {.val {prep$D}} with values 1 or -1."
      ))
    }
    stan_data$skew_direction <- as.numeric(skew_direction)
  }

  attr(stan_data, "vars") <- prep$vars
  attr(stan_data, "standardized") <- prep$standardized
  attr(stan_data, "margins") <- margins
  attr(stan_data, "time_values") <- prep$time_values
  if (!is.null(skew_direction)) attr(stan_data, "skew_direction") <- skew_direction
  if (prep$standardized) {
    attr(stan_data, "Y_means") <- prep$Y_means
    attr(stan_data, "Y_sds") <- prep$Y_sds
  }

  stan_data
}
