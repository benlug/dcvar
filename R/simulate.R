# ============================================================================
# Data Simulation
# ============================================================================

#' Internal: validate a positive finite scalar for simulation inputs
#' @noRd
.simulate_validate_positive_scalar <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    cli_abort("{.arg {arg_name}} must be a single positive finite numeric value.")
  }
}

#' Internal: validate a finite numeric vector for simulation inputs
#' @noRd
.simulate_validate_numeric_vector <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x))) {
    cli_abort("{.arg {arg_name}} must be a non-empty finite numeric vector.")
  }
}

#' Simulate data from a copula VAR(1) model
#'
#' Generates bivariate time series data with correlated innovations
#' driven by a specified rho trajectory.
#'
#' @param T Number of time points.
#' @param rho_trajectory Numeric vector of length T-1 specifying the
#'   correlation at each time step. Use [rho_constant()], [rho_decreasing()],
#'   etc.
#' @param mu Intercept vector of length 2 (default: `c(0, 0)`).
#' @param Phi VAR(1) coefficient matrix, 2x2 (default:
#'   `matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2)`).
#' @param sigma_eps Innovation standard deviations, length 2 (default:
#'   `c(1, 1)`). Used for normal margins.
#' @param margins Character: `"normal"` (default), `"exponential"`,
#'   `"skew_normal"`, or `"gamma"`.
#' @param skew_direction Length-2 integer vector of +1/-1. Required for
#'   `"exponential"` and `"gamma"` margins.
#' @param skew_params Named list of margin-specific parameters. For
#'   `"skew_normal"`: `alpha` (length-2 vector of skew-normal shape params).
#'   For `"gamma"`: `shape` (scalar gamma shape parameter).
#' @param seed Random seed for reproducibility.
#'
#' @return A named list with:
#'   - `Y`: T x 2 observation matrix
#'   - `Y_df`: data frame with columns `time`, `y1`, `y2` (ready for
#'     [dcvar()])
#'   - `true_params`: list of true parameter values
#' @export
#'
#' @examples
#' sim <- simulate_dcvar(T = 100, rho_trajectory = rho_decreasing(100))
#' head(sim$Y_df)
#' plot(sim$true_params$rho, type = "l")
simulate_dcvar <- function(T,
                           rho_trajectory,
                           mu = c(0, 0),
                           Phi = matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2),
                           sigma_eps = c(1, 1),
                           margins = "normal",
                           skew_direction = NULL,
                           skew_params = NULL,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.numeric(T) || length(T) != 1 || T != as.integer(T) || T < 2) {
    cli_abort("{.arg T} must be an integer >= 2, got {.val {T}}.")
  }

  if (length(rho_trajectory) != T - 1) {
    cli_abort(
      "{.arg rho_trajectory} must have length T-1 = {T - 1}, got {length(rho_trajectory)}."
    )
  }
  if (any(abs(rho_trajectory) > 1)) {
    cli_abort("{.arg rho_trajectory} values must be in [-1, 1].")
  }

  D <- length(mu)
  if (D != 2) {
    cli_abort(
      "{.fun simulate_dcvar} currently supports bivariate (D = 2) models only, got D = {D}."
    )
  }
  if (!is.matrix(Phi) || !all(dim(Phi) == c(D, D))) {
    cli_abort("{.arg Phi} must be a {D}x{D} matrix.")
  }

  .validate_margins(margins, skew_direction)
  if (margins == "normal" && length(sigma_eps) != D) {
    cli_abort("{.arg sigma_eps} must have length {.val {D}}, got {.val {length(sigma_eps)}}.")
  }
  if (margins == "normal") {
    .simulate_validate_numeric_vector(sigma_eps, "sigma_eps")
    if (any(sigma_eps <= 0)) {
      cli_abort("{.arg sigma_eps} values must be positive.")
    }
  }
  if (margins == "gamma") {
    if (!is.null(skew_params) && !is.list(skew_params)) {
      cli_abort("{.arg skew_params} must be a list when using gamma margins.")
    }
    gamma_shape <- if (is.null(skew_params) || is.null(skew_params$shape)) 1 else skew_params$shape
    .simulate_validate_positive_scalar(gamma_shape, "skew_params$shape")
  }
  if (margins == "skew_normal") {
    if (!is.null(skew_params) && !is.list(skew_params)) {
      cli_abort("{.arg skew_params} must be a list when using skew_normal margins.")
    }
    alpha <- if (is.null(skew_params) || is.null(skew_params$alpha)) c(0, 0) else skew_params$alpha
    .simulate_validate_numeric_vector(alpha, "skew_params$alpha")
    if (length(alpha) != D) {
      cli_abort("{.arg skew_params$alpha} must have length {.val {D}}, got {.val {length(alpha)}}.")
    }
    skew_params <- list(alpha = alpha)
  }
  if (margins == "gamma") {
    skew_params <- list(shape = gamma_shape)
  }

  Y <- matrix(0, T, D)
  Y[1, ] <- mu

  for (t in 2:T) {
    rho_t <- rho_trajectory[t - 1]

    # Generate copula uniforms via Gaussian copula
    L <- matrix(c(1, rho_t, 0, sqrt(1 - rho_t^2)), 2, 2)
    z <- rnorm(D)
    w <- L %*% z  # correlated standard normals

    # Transform through marginal quantiles
    eps <- .sim_marginal_quantile(w, margins, sigma_eps, skew_direction, skew_params)

    # VAR(1) update
    Y[t, ] <- mu + Phi %*% (Y[t - 1, ] - mu) + eps
  }

  Y_df <- data.frame(
    time = 1:T,
    y1 = Y[, 1],
    y2 = Y[, 2]
  )

  true_params <- list(
    rho = rho_trajectory,
    Phi = Phi,
    mu = mu,
    margins = margins
  )

  # Add margin-specific true params
  if (margins == "normal") {
    true_params$sigma_eps <- sigma_eps
  } else if (margins == "exponential") {
    true_params$skew_direction <- skew_direction
  } else if (margins == "skew_normal") {
    true_params$skew_params <- skew_params
  } else if (margins == "gamma") {
    true_params$skew_direction <- skew_direction
    true_params$skew_params <- skew_params
  }

  list(
    Y = Y,
    Y_df = Y_df,
    true_params = true_params
  )
}


#' Internal: transform correlated standard normals to marginal quantiles
#' @noRd
.sim_marginal_quantile <- function(w, margins, sigma_eps, skew_direction, skew_params) {
  D <- length(w)

  switch(margins,
    normal = {
      # w are already standard normal, just scale
      w * sigma_eps
    },
    exponential = {
      # Convert to uniforms, then to standardized exponential
      u <- stats::pnorm(w)
      eps <- numeric(D)
      for (i in seq_len(D)) {
        x_raw <- stats::qexp(u[i], rate = 1)
        x_std <- x_raw - 1  # standardize: mean=1, sd=1 for Exp(1)
        eps[i] <- if (skew_direction[i] < 0) -x_std else x_std
      }
      eps
    },
    skew_normal = {
      if (!requireNamespace("sn", quietly = TRUE)) {
        cli_abort("Package {.pkg sn} is required for skew-normal simulation.")
      }
      alpha <- skew_params$alpha %||% c(0, 0)
      u <- stats::pnorm(w)
      eps <- numeric(D)
      for (i in seq_len(D)) {
        delta <- alpha[i] / sqrt(1 + alpha[i]^2)
        omega_i <- sqrt(1 / (1 - 2 * delta^2 / pi))
        xi_i <- -omega_i * delta * sqrt(2 / pi)
        eps[i] <- sn::qsn(u[i], xi = xi_i, omega = omega_i, alpha = alpha[i])
      }
      eps
    },
    gamma = {
      shape <- skew_params$shape %||% 1
      u <- stats::pnorm(w)
      eps <- numeric(D)
      for (i in seq_len(D)) {
        x_raw <- stats::qgamma(u[i], shape = shape, rate = sqrt(shape))
        x_std <- x_raw - sqrt(shape)  # standardize: mean=sqrt(shape)*1/sqrt(shape)=1
        eps[i] <- if (skew_direction[i] < 0) -x_std else x_std
      }
      eps
    },
    cli_abort("Unknown margin type: {.val {margins}}")
  )
}
