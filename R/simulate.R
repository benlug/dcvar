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
#' @param n_time Number of time points.
#' @param rho_trajectory Numeric vector of length `n_time - 1` specifying the
#'   correlation at each time step. Use [rho_constant()], [rho_decreasing()],
#'   etc.
#' @param mu Intercept vector of length 2 (default: `c(0, 0)`).
#' @param Phi VAR(1) coefficient matrix, 2x2 (default:
#'   `matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2)`).
#' @param sigma_eps Innovation standard deviations, length 2 (default:
#'   `c(1, 1)`). Used for normal margins.
#' @param margins Marginal family. Either a single string applied to both
#'   variables, or a length-2 character vector for per-variable (mixed) margins,
#'   e.g. `c("normal", "exponential")`. Each entry is one of `"normal"`
#'   (default), `"exponential"`, `"skew_normal"`, or `"gamma"`.
#' @param skew_direction Length-2 integer vector of +1/-1. Required whenever any
#'   dimension uses an `"exponential"` or `"gamma"` margin; only those
#'   dimensions consult it.
#' @param skew_params Named list of margin-specific parameters. `alpha`
#'   (length-2 vector of skew-normal shape params) is used by skew-normal
#'   dimensions; `shape` (scalar gamma shape parameter) is used by gamma
#'   dimensions. Both may be supplied together for mixed margins.
#' @param seed Random seed for reproducibility.
#'
#' @return A named list with:
#'   - `Y`: `n_time x 2` observation matrix
#'   - `Y_df`: data frame with columns `time`, `y1`, `y2` (ready for
#'     [dcvar()])
#'   - `true_params`: list of true parameter values
#' @export
#'
#' @examples
#' sim <- simulate_dcvar(n_time = 100, rho_trajectory = rho_decreasing(100))
#' head(sim$Y_df)
#' plot(sim$true_params$rho, type = "l")
simulate_dcvar <- function(n_time,
                           rho_trajectory,
                           mu = c(0, 0),
                           Phi = matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2),
                           sigma_eps = c(1, 1),
                           margins = "normal",
                           skew_direction = NULL,
                           skew_params = NULL,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.numeric(n_time) || length(n_time) != 1 || n_time != as.integer(n_time) || n_time < 2) {
    cli_abort("{.arg n_time} must be an integer >= 2, got {.val {n_time}}.")
  }

  if (length(rho_trajectory) != n_time - 1L) {
    cli_abort(
      "{.arg rho_trajectory} must have length n_time - 1 = {n_time - 1}, got {length(rho_trajectory)}."
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

  margins <- .normalize_margins_spec(margins)
  .validate_margins(margins, skew_direction)
  margins_vec <- if (length(margins) == 1L) rep(margins, D) else margins

  if (!is.null(skew_params) && !is.list(skew_params)) {
    cli_abort("{.arg skew_params} must be a list.")
  }
  skew_params <- if (is.list(skew_params)) skew_params else list()

  # Normal dimensions use sigma_eps
  if (any(margins_vec == "normal")) {
    .simulate_validate_numeric_vector(sigma_eps, "sigma_eps")
    if (length(sigma_eps) != D) {
      cli_abort("{.arg sigma_eps} must have length {.val {D}}, got {.val {length(sigma_eps)}}.")
    }
    if (any(sigma_eps <= 0)) {
      cli_abort("{.arg sigma_eps} values must be positive.")
    }
  }
  # Skew-normal dimensions use a per-dimension alpha
  if (any(margins_vec == "skew_normal")) {
    alpha <- skew_params$alpha %||% rep(0, D)
    .simulate_validate_numeric_vector(alpha, "skew_params$alpha")
    if (length(alpha) != D) {
      cli_abort("{.arg skew_params$alpha} must have length {.val {D}}, got {.val {length(alpha)}}.")
    }
    skew_params$alpha <- alpha
  }
  # Gamma dimensions use a shared shape
  if (any(margins_vec == "gamma")) {
    gamma_shape <- skew_params$shape %||% 1
    .simulate_validate_positive_scalar(gamma_shape, "skew_params$shape")
    skew_params$shape <- gamma_shape
  }

  Y <- matrix(0, n_time, D)
  Y[1, ] <- mu

  for (time_index in 2:n_time) {
    rho_t <- rho_trajectory[time_index - 1L]

    # Generate copula uniforms via Gaussian copula
    L <- matrix(c(1, rho_t, 0, sqrt(1 - rho_t^2)), 2, 2)
    z <- rnorm(D)
    w <- L %*% z  # correlated standard normals

    # Transform through marginal quantiles
    eps <- .sim_marginal_quantile(w, margins, sigma_eps, skew_direction, skew_params)

    # VAR(1) update
    Y[time_index, ] <- mu + Phi %*% (Y[time_index - 1L, ] - mu) + eps
  }

  Y_df <- data.frame(
    time = seq_len(n_time),
    y1 = Y[, 1],
    y2 = Y[, 2]
  )

  true_params <- list(
    rho = rho_trajectory,
    Phi = Phi,
    mu = mu,
    margins = margins
  )

  # Add margin-specific true params for every family present (independent
  # checks so mixed margins record each family's parameters).
  if (any(margins_vec == "normal")) {
    true_params$sigma_eps <- sigma_eps
  }
  if (any(margins_vec %in% c("exponential", "gamma"))) {
    true_params$skew_direction <- skew_direction
  }
  if (any(margins_vec %in% c("skew_normal", "gamma"))) {
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
  margins_vec <- if (length(margins) == 1L) rep(margins, D) else margins
  eps <- numeric(D)

  for (i in seq_len(D)) {
    fam <- margins_vec[[i]]
    if (identical(fam, "normal")) {
      # w[i] is already standard normal, just scale
      eps[i] <- w[i] * sigma_eps[i]
    } else if (identical(fam, "exponential")) {
      # Convert to uniform, then to standardized exponential (Exp(1), mean/sd 1).
      # The Stan likelihoods use u = 1 - F(x_shifted) for left-skewed dimensions,
      # so the uniform must be flipped before the quantile to keep the simulated
      # eps comonotone with the latent copula normal w.
      u <- stats::pnorm(w[i])
      if (skew_direction[i] < 0) u <- 1 - u
      x_std <- stats::qexp(u, rate = 1) - 1
      eps[i] <- skew_direction[i] * x_std
    } else if (identical(fam, "skew_normal")) {
      if (!requireNamespace("sn", quietly = TRUE)) {
        cli_abort("Package {.pkg sn} is required for skew-normal simulation.")
      }
      alpha_i <- (skew_params$alpha %||% rep(0, D))[i]
      delta <- alpha_i / sqrt(1 + alpha_i^2)
      omega_i <- sqrt(1 / (1 - 2 * delta^2 / pi))
      xi_i <- -omega_i * delta * sqrt(2 / pi)
      eps[i] <- sn::qsn(stats::pnorm(w[i]), xi = xi_i, omega = omega_i, alpha = alpha_i)
    } else if (identical(fam, "gamma")) {
      shape <- skew_params$shape %||% 1
      # Same uniform flip as the exponential branch (see comment there).
      u <- stats::pnorm(w[i])
      if (skew_direction[i] < 0) u <- 1 - u
      x_std <- stats::qgamma(u, shape = shape, rate = sqrt(shape)) - sqrt(shape)
      eps[i] <- skew_direction[i] * x_std
    } else {
      cli_abort("Unknown margin type: {.val {fam}}")
    }
  }

  eps
}
