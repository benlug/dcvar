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

#' Internal: normalise a (possibly time-varying) trajectory specification
#'
#' Accepts an `(n_time - 1) x k` matrix, a list of `k` length-`(n_time - 1)`
#' vectors (so the `rho_*` generators can be reused per column), or a constant
#' specification (`k`-vector, or `2 x 2` matrix when `k = 4`, expanded over
#' time). Returns an `(n_time - 1) x k` matrix. For `k = 4` the column order
#' is row-major: phi11, phi12, phi21, phi22.
#' @noRd
.tv_resolve_trajectory <- function(x, n_time, k, arg_name) {
  n_eff <- n_time - 1L

  if (is.list(x)) {
    if (length(x) != k) {
      cli_abort("{.arg {arg_name}} as a list must have {.val {k}} elements, got {.val {length(x)}}.")
    }
    lens <- lengths(x)
    if (!all(lens == n_eff)) {
      cli_abort("Each element of {.arg {arg_name}} must have length {.val {n_eff}} (n_time - 1).")
    }
    x <- do.call(cbind, x)
  } else if (is.matrix(x)) {
    if (k == 4L && all(dim(x) == c(2L, 2L))) {
      # Constant 2x2 Phi, expanded row-major over time
      x <- matrix(rep(c(t(x)), each = n_eff), n_eff, k)
    } else if (!all(dim(x) == c(n_eff, k))) {
      cli_abort(c(
        "{.arg {arg_name}} must be a {.val {n_eff}} x {.val {k}} matrix.",
        "i" = "Rows follow the transitions 2..n_time; for Phi the columns are phi11, phi12, phi21, phi22."
      ))
    }
  } else if (is.numeric(x) && length(x) == k) {
    x <- matrix(rep(x, each = n_eff), n_eff, k)
  } else {
    cli_abort(c(
      "{.arg {arg_name}} must be a {.val {n_eff}} x {.val {k}} matrix, a list of {.val {k}} length-{.val {n_eff}} vectors, or a constant length-{.val {k}} specification.",
      "i" = "For Phi a constant {.val 2 x 2} matrix is also accepted."
    ))
  }

  x <- unname(as.matrix(x))
  if (!is.numeric(x) || any(!is.finite(x))) {
    cli_abort("{.arg {arg_name}} must contain only finite numeric values.")
  }
  x
}

#' Internal: per-t spectral radius of a row-major (phi11, phi12, phi21, phi22) path
#' @noRd
.tv_spectral_radius <- function(phi_mat) {
  apply(phi_mat, 1, function(p) {
    tr <- p[1] + p[4]
    det_phi <- p[1] * p[4] - p[2] * p[3]
    disc <- tr^2 - 4 * det_phi
    if (disc >= 0) {
      max(abs(0.5 * (tr + sqrt(disc))), abs(0.5 * (tr - sqrt(disc))))
    } else {
      sqrt(abs(det_phi))
    }
  })
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
#' @param phi_trajectory Optional time-varying VAR coefficient paths for data
#'   matching `dcvar(tv_phi = TRUE)`: an `(n_time - 1) x 4` matrix with
#'   columns in row-major order (phi11, phi12, phi21, phi22), a list of 4
#'   length-`(n_time - 1)` vectors (the `rho_*` trajectory generators can be
#'   reused per coefficient), or a constant `2 x 2` matrix. Mutually exclusive
#'   with `Phi`. `NULL` (default) keeps the constant `Phi`.
#' @param sigma_trajectory Optional time-varying residual scale paths for data
#'   matching `dcvar(tv_sigma = TRUE)`: an `(n_time - 1) x 2` positive matrix,
#'   a list of 2 length-`(n_time - 1)` vectors, or a constant length-2 vector.
#'   The scale is each family's natural scale (innovation SD for normal,
#'   residual SD for skew-normal, `sigma_exp` / `sigma_gam` for
#'   exponential / gamma). Exponential and gamma dimensions must have a
#'   constant path (matching the fitted model's restriction). Mutually
#'   exclusive with `sigma_eps`. `NULL` (default) keeps constant scales.
#' @param seed Random seed for reproducibility.
#'
#' @return A named list with:
#'   - `Y`: `n_time x 2` observation matrix
#'   - `Y_df`: data frame with columns `time`, `y1`, `y2` (ready for
#'     [dcvar()])
#'   - `true_params`: list of true parameter values. These are on the raw
#'     data scale; the fitting functions standardize by default, so round-trip
#'     comparisons of `mu`, `Phi`, and `sigma_eps` require fitting with
#'     `standardize = FALSE` (the `rho` trajectory is scale-invariant).
#'     Without `sigma_trajectory`, exponential and gamma dimensions are
#'     simulated with unit-SD standardized innovations, so their implied true
#'     scale (`sigma_exp`/`sigma_gam`) is 1. With trajectories supplied,
#'     `Phi` is the `(n_time - 1) x 4` coefficient path and `sigma` the
#'     `(n_time - 1) x 2` scale path.
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
                           phi_trajectory = NULL,
                           sigma_trajectory = NULL,
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
  if (!is.null(phi_trajectory) && !missing(Phi)) {
    cli_abort("Supply either {.arg Phi} (constant) or {.arg phi_trajectory} (time-varying), not both.")
  }
  if (!is.null(sigma_trajectory) && !missing(sigma_eps)) {
    cli_abort("Supply either {.arg sigma_eps} (constant) or {.arg sigma_trajectory} (time-varying), not both.")
  }
  if (!is.matrix(Phi) || !all(dim(Phi) == c(D, D))) {
    cli_abort("{.arg Phi} must be a {D}x{D} matrix.")
  }

  margins <- .normalize_margins_spec(margins)
  .validate_margins(margins, skew_direction)
  margins_vec <- if (length(margins) == 1L) rep(margins, D) else margins

  # Resolve the (possibly constant) per-transition parameter paths. Row t
  # governs the transition Y[t] -> Y[t+1], matching the Stan models' eps[t].
  phi_mat <- if (is.null(phi_trajectory)) {
    .tv_resolve_trajectory(Phi, n_time, 4L, "Phi")
  } else {
    .tv_resolve_trajectory(phi_trajectory, n_time, 4L, "phi_trajectory")
  }
  sr <- .tv_spectral_radius(phi_mat)
  if (any(sr >= 1)) {
    cli_warn(c(
      "{sum(sr >= 1)} time point{?s} of the VAR coefficient path {?is/are} nonstationary (spectral radius >= 1).",
      "i" = "Brief excursions are fine; persistent ones make the series explode."
    ))
  }

  scale_mat <- if (is.null(sigma_trajectory)) {
    NULL
  } else {
    m <- .tv_resolve_trajectory(sigma_trajectory, n_time, D, "sigma_trajectory")
    if (any(m <= 0)) {
      cli_abort("{.arg sigma_trajectory} values must be positive.")
    }
    for (d in seq_len(D)) {
      if (margins_vec[d] %in% c("exponential", "gamma") && length(unique(m[, d])) > 1L) {
        cli_abort(c(
          "{.arg sigma_trajectory} must be constant on dimension {d} ({.val {margins_vec[d]}} margin).",
          "i" = "The fitted model holds shifted exponential/gamma scales constant; see {.fun dcvar}."
        ))
      }
    }
    m
  }

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

  # Constant per-dimension scales used when no sigma_trajectory is supplied:
  # the innovation SD for normal dims, unit scale otherwise.
  base_scales <- rep(1, D)
  is_normal <- margins_vec == "normal"
  if (any(is_normal)) base_scales[is_normal] <- sigma_eps[is_normal]

  Y <- matrix(0, n_time, D)
  Y[1, ] <- mu

  for (time_index in 2:n_time) {
    t_eff <- time_index - 1L
    rho_t <- rho_trajectory[t_eff]

    # Generate copula uniforms via Gaussian copula
    L <- matrix(c(1, rho_t, 0, sqrt(1 - rho_t^2)), 2, 2)
    z <- rnorm(D)
    w <- L %*% z  # correlated standard normals

    # Transform through marginal quantiles with the per-transition scales
    scales_t <- if (is.null(scale_mat)) base_scales else scale_mat[t_eff, ]
    eps <- .sim_marginal_quantile_scaled(w, margins_vec, scales_t, skew_direction, skew_params)

    # VAR(1) update with the per-transition coefficients (byrow: the path
    # columns are row-major phi11, phi12, phi21, phi22)
    Phi_t <- matrix(phi_mat[t_eff, ], 2, 2, byrow = TRUE)
    Y[time_index, ] <- mu + Phi_t %*% (Y[time_index - 1L, ] - mu) + eps
  }

  Y_df <- data.frame(
    time = seq_len(n_time),
    y1 = Y[, 1],
    y2 = Y[, 2]
  )

  true_params <- list(
    rho = rho_trajectory,
    Phi = if (is.null(phi_trajectory)) {
      Phi
    } else {
      colnames(phi_mat) <- c("phi11", "phi12", "phi21", "phi22")
      phi_mat
    },
    mu = mu,
    margins = margins
  )

  # Add margin-specific true params for every family present (independent
  # checks so mixed margins record each family's parameters).
  if (!is.null(scale_mat)) {
    colnames(scale_mat) <- c("y1", "y2")
    true_params$sigma <- scale_mat
  } else if (any(margins_vec == "normal")) {
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
  scales <- rep(1, D)
  is_normal <- margins_vec == "normal"
  if (any(is_normal)) scales[is_normal] <- sigma_eps[is_normal]
  .sim_marginal_quantile_scaled(w, margins_vec, scales, skew_direction, skew_params)
}


#' Internal: scale-aware marginal quantile transform
#'
#' Core shared by the constant-parameter and time-varying simulators. `scales`
#' is the per-dimension scale on each family's natural reporting scale: the
#' innovation SD for normal dims, `sigma_exp` for exponential dims (mean and
#' SD of the shifted exponential), the SD (`sigma_gam`) for gamma dims, and
#' the residual SD for skew-normal dims. All residuals have mean zero by
#' construction.
#' @noRd
.sim_marginal_quantile_scaled <- function(w, margins_vec, scales, skew_direction, skew_params) {
  D <- length(w)
  eps <- numeric(D)

  for (i in seq_len(D)) {
    fam <- margins_vec[[i]]
    s <- scales[i]
    if (identical(fam, "normal")) {
      # w[i] is already standard normal, just scale
      eps[i] <- w[i] * s
    } else if (identical(fam, "exponential")) {
      # Convert to uniform, then to a shifted exponential with mean/sd s.
      # The Stan likelihoods use u = 1 - F(x_shifted) for left-skewed dimensions,
      # so the uniform must be flipped before the quantile to keep the simulated
      # eps comonotone with the latent copula normal w.
      u <- stats::pnorm(w[i])
      if (skew_direction[i] < 0) u <- 1 - u
      x_std <- stats::qexp(u, rate = 1 / s) - s
      eps[i] <- skew_direction[i] * x_std
    } else if (identical(fam, "skew_normal")) {
      if (!requireNamespace("sn", quietly = TRUE)) {
        cli_abort("Package {.pkg sn} is required for skew-normal simulation.")
      }
      alpha_i <- (skew_params$alpha %||% rep(0, D))[i]
      delta <- alpha_i / sqrt(1 + alpha_i^2)
      omega_i <- s * sqrt(1 / (1 - 2 * delta^2 / pi))
      xi_i <- -omega_i * delta * sqrt(2 / pi)
      eps[i] <- sn::qsn(stats::pnorm(w[i]), xi = xi_i, omega = omega_i, alpha = alpha_i)
    } else if (identical(fam, "gamma")) {
      shape <- skew_params$shape %||% 1
      # Same uniform flip as the exponential branch (see comment there).
      # rate = sqrt(shape)/s gives SD s; the mean shift keeps E[eps] = 0.
      u <- stats::pnorm(w[i])
      if (skew_direction[i] < 0) u <- 1 - u
      x_std <- stats::qgamma(u, shape = shape, rate = sqrt(shape) / s) - sqrt(shape) * s
      eps[i] <- skew_direction[i] * x_std
    } else {
      cli_abort("Unknown margin type: {.val {fam}}")
    }
  }

  eps
}
