# ============================================================================
# SEM Simulation
# ============================================================================

#' Internal: validate a positive finite scalar for SEM simulation inputs
#' @noRd
.simulate_sem_validate_positive_scalar <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    cli_abort("{.arg {arg_name}} must be a single positive finite numeric value.")
  }
}

#' Internal: validate a finite numeric vector for SEM simulation inputs
#' @noRd
.simulate_sem_validate_numeric_vector <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x))) {
    cli_abort("{.arg {arg_name}} must be a non-empty finite numeric vector.")
  }
}

#' Simulate data from a SEM copula VAR(1) model
#'
#' Generates indicator-level time series data from a latent VAR(1) process
#' with Gaussian copula dependence and a fixed measurement model.
#'
#' @param n_time Number of time points.
#' @param J Number of indicators per latent variable.
#' @param lambda Numeric vector of length J with factor loadings.
#' @param sigma_e Measurement error SD (scalar).
#' @param Phi 2x2 VAR coefficient matrix.
#' @param mu Length-2 intercept vector.
#' @param margins Latent-innovation marginal family. Either a single string
#'   applied to both latent variables, or a length-2 character vector for
#'   per-variable (mixed) margins, e.g. `c("normal", "gamma")`, where each entry
#'   is one of `"normal"`, `"exponential"`, `"skew_normal"`, or `"gamma"`.
#' @param sigma Length-2 latent innovation scale vector used by normal,
#'   skew-normal, gamma, and mixed-margin specifications. Values are on each
#'   family's natural residual scale: innovation SD for normal/skew-normal,
#'   and `sigma_gam` for gamma. Homogeneous exponential margins use
#'   `sigma_exp` instead.
#' @param sigma_exp Length-2 shifted-exponential scale vector for the
#'   single-family exponential path.
#' @param skew_direction Integer vector of length 2 of `1`/`-1`. Required
#'   whenever any dimension uses an `"exponential"` or `"gamma"` margin.
#' @param skew_params Named list of margin-specific parameters for mixed
#'   margins: `alpha` (length-2 skew-normal shape) and/or `shape` (scalar gamma
#'   shape).
#' @param rho Copula correlation.
#' @param rho_trajectory Optional numeric vector of length `n_time - 1` for a
#'   time-varying copula correlation. Mutually exclusive with `rho`.
#' @param phi_trajectory Optional time-varying latent VAR coefficient paths,
#'   accepted in the same forms as [simulate_dcvar()].
#' @param sigma_trajectory Optional time-varying latent innovation scale paths,
#'   accepted in the same forms as [simulate_dcvar()]. The supplied value is
#'   each family's natural scale: innovation SD (normal), residual SD
#'   (skew-normal), `sigma_exp` (exponential), `sigma_gam` (gamma).
#' @param tv_sigma_k Soft-barrier sharpness for time-varying exponential/gamma
#'   scales.
#' @param burnin Retained for backward compatibility but ignored. Default `0`
#'   keeps the default simulation path aligned with the fitted SEM model,
#'   which conditions on `x_0 = 0` and treats the first returned state as
#'   observed rather than drawn after a burn-in period.
#' @param seed Random seed for reproducibility.
#'
#' @return A named list with:
#'   - `data`: data frame with columns `time`, `y1_1`, ..., `y1_J`,
#'     `y2_1`, ..., `y2_J`
#'   - `true_params`: list of true parameter values
#'   - `latent_states`: `n_time x 2` matrix of true latent states
#'   - `innovations`: `n_time x 2` matrix of true innovations
#' @export
simulate_dcvar_sem <- function(n_time = 200, J = 3,
                                lambda = rep(sqrt(0.8), 3),
                                sigma_e = sqrt(0.2),
                                Phi = matrix(c(0.5, 0.15, 0.15, 0.3), 2, 2),
                                mu = c(0, 0),
                                margins = "normal",
                                sigma = c(1, 1),
                                sigma_exp = c(1, 1),
                                skew_direction = NULL,
                                skew_params = NULL,
                                rho = 0.3,
                                rho_trajectory = NULL,
                                phi_trajectory = NULL,
                                sigma_trajectory = NULL,
                                tv_sigma_k = NULL,
                                burnin = 0,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  margins <- .normalize_margins_spec(margins)
  .validate_sem_margins(margins, skew_direction)
  margins_vec <- if (length(margins) == 1L) rep(margins, 2L) else margins
  mixed <- .is_mixed_margins(margins)
  mixed_engine <- .uses_sem_multilevel_mixed_engine("sem", margins)
  engine_margins <- .sem_multilevel_engine_margins("sem", margins, 2L)
  skew_params <- if (is.list(skew_params)) skew_params else list()
  if (any(margins_vec == "skew_normal")) skew_params$alpha <- skew_params$alpha %||% rep(0, 2L)
  if (any(margins_vec == "gamma")) skew_params$shape <- skew_params$shape %||% 1
  if (!is.numeric(n_time) || length(n_time) != 1L || n_time != as.integer(n_time) || n_time < 1) {
    cli_abort("{.arg n_time} must be an integer >= 1, got {.val {n_time}}.")
  }
  if (!is.numeric(J) || length(J) != 1L || J != as.integer(J) || J < 1) {
    cli_abort("{.arg J} must be an integer >= 1, got {.val {J}}.")
  }

  if (length(lambda) != J) {
    cli_abort("{.arg lambda} must have length {.val {J}}, got {.val {length(lambda)}}.")
  }
  .simulate_sem_validate_numeric_vector(lambda, "lambda")
  .simulate_sem_validate_positive_scalar(sigma_e, "sigma_e")
  if (!is.matrix(Phi) || !all(dim(Phi) == c(2L, 2L)) || any(!is.finite(Phi))) {
    cli_abort("{.arg Phi} must be a finite 2x2 matrix.")
  }
  if (!is.null(phi_trajectory) && !missing(Phi)) {
    cli_abort("Supply either {.arg Phi} (constant) or {.arg phi_trajectory} (time-varying), not both.")
  }
  .simulate_sem_validate_numeric_vector(mu, "mu")
  if (length(mu) != 2L) {
    cli_abort("{.arg mu} must have length 2, got {.val {length(mu)}}.")
  }
  # sigma is the constant scale for every non-homogeneous-exponential path.
  if (!(!mixed && identical(margins, "exponential"))) {
    .simulate_sem_validate_numeric_vector(sigma, "sigma")
    if (length(sigma) != 2L) {
      cli_abort("{.arg sigma} must have length 2, got {.val {length(sigma)}}.")
    }
    if (any(sigma <= 0)) {
      cli_abort("{.arg sigma} values must be positive.")
    }
  }
  if (!mixed && identical(margins, "exponential")) {
    .simulate_sem_validate_numeric_vector(sigma_exp, "sigma_exp")
    if (length(sigma_exp) != 2L) {
      cli_abort("{.arg sigma_exp} must have length 2, got {.val {length(sigma_exp)}}.")
    }
    if (any(sigma_exp <= 0)) {
      cli_abort("{.arg sigma_exp} values must be positive.")
    }
  }
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < -1 || rho > 1) {
    cli_abort("{.arg rho} must be a single finite numeric value in [-1, 1].")
  }
  if (!is.null(rho_trajectory) && !missing(rho)) {
    cli_abort("Supply either {.arg rho} (constant) or {.arg rho_trajectory} (time-varying), not both.")
  }
  if (!is.null(sigma_trajectory) && !missing(sigma)) {
    cli_abort("Supply either {.arg sigma} (constant) or {.arg sigma_trajectory} (time-varying), not both.")
  }
  if (!is.null(sigma_trajectory) && !missing(sigma_exp)) {
    cli_abort("Supply either {.arg sigma_exp} (constant) or {.arg sigma_trajectory} (time-varying), not both.")
  }
  if (!is.null(tv_sigma_k)) {
    .simulate_validate_positive_scalar(tv_sigma_k, "tv_sigma_k")
  }

  if (!identical(burnin, 0L) && !identical(burnin, 0)) {
    cli_warn(
      c(
        "{.arg burnin} is ignored for SEM simulation.",
        "i" = "The fitted SEM model treats the first returned state as the observed draw implied by `x_0 = 0`."
      )
    )
  }

  rho_path <- if (is.null(rho_trajectory)) {
    rep(rho, max(n_time - 1L, 1L))
  } else {
    if (length(rho_trajectory) != n_time - 1L) {
      cli_abort("{.arg rho_trajectory} must have length {.val {n_time - 1L}}.")
    }
    if (any(!is.finite(rho_trajectory)) || any(abs(rho_trajectory) > 1)) {
      cli_abort("{.arg rho_trajectory} values must be finite and in [-1, 1].")
    }
    rho_trajectory
  }
  phi_mat <- if (is.null(phi_trajectory)) {
    .tv_resolve_trajectory(Phi, n_time, 4L, "Phi")
  } else {
    .tv_resolve_trajectory(phi_trajectory, n_time, 4L, "phi_trajectory")
  }
  scale_mat <- if (is.null(sigma_trajectory)) {
    NULL
  } else {
    m <- .tv_resolve_trajectory(sigma_trajectory, n_time, 2L, "sigma_trajectory")
    if (any(m <= 0)) {
      cli_abort("{.arg sigma_trajectory} values must be positive.")
    }
    if (is.null(tv_sigma_k)) {
      for (d in seq_len(2L)) {
        if (margins_vec[d] %in% c("exponential", "gamma") && length(unique(m[, d])) > 1L) {
          cli_abort(c(
            "{.arg sigma_trajectory} is not constant on latent dimension {d} ({.val {margins_vec[d]}} margin).",
            "i" = "Set {.arg tv_sigma_k} to simulate from the soft-barrier model, or keep the scale constant."
          ))
        }
      }
    }
    m
  }

  base_scales <- if (!mixed && identical(margins, "exponential")) {
    sigma_exp
  } else {
    sigma
  }

  # Generate correlated innovations via Gaussian copula
  zeta <- matrix(NA_real_, n_time, 2)
  for (time_index in seq_len(n_time)) {
    tt <- if (time_index == 1L) 1L else time_index - 1L
    rho_t <- rho_path[tt]
    L <- matrix(c(1, rho_t, 0, sqrt(1 - rho_t^2)), 2, 2)
    z <- rnorm(2)
    w <- drop(L %*% z)
    scales_t <- if (is.null(scale_mat)) base_scales else scale_mat[tt, ]
    sim_margins <- if (mixed_engine) engine_margins else margins_vec
    zeta[time_index, ] <- .sim_marginal_quantile_scaled(
      w, sim_margins, scales_t, skew_direction, skew_params,
      barrier_k = tv_sigma_k
    )
  }

  # Latent VAR(1) recursion matching the SEM Stan model, which conditions on x_0 = 0.
  state <- matrix(0, n_time, 2)
  state[1, ] <- mu + zeta[1, ]
  for (time_index in seq_len(n_time - 1L) + 1L) {
    Phi_t <- matrix(phi_mat[time_index - 1L, ], 2, 2, byrow = TRUE)
    state[time_index, ] <- mu + as.vector(Phi_t %*% state[time_index - 1L, ]) + zeta[time_index, ]
  }

  # Measurement model: y_{ij,t} = lambda_j * state_{i,t} + e_{ij,t}
  y <- matrix(NA_real_, n_time, 2 * J)
  for (j in seq_len(J)) {
    y[, j]     <- lambda[j] * state[, 1] + rnorm(n_time, 0, sigma_e)
    y[, J + j] <- lambda[j] * state[, 2] + rnorm(n_time, 0, sigma_e)
  }
  colnames(y) <- c(paste0("y1_", seq_len(J)), paste0("y2_", seq_len(J)))

  data <- data.frame(time = seq_len(n_time), y, check.names = FALSE)

  true_params <- list(
    Phi = if (is.null(phi_trajectory)) {
      Phi
    } else {
      colnames(phi_mat) <- c("phi11", "phi12", "phi21", "phi22")
      phi_mat
    },
    mu = mu,
    margins = margins,
    rho = if (is.null(rho_trajectory)) rho else rho_path,
    lambda = lambda,
    sigma_e = sigma_e,
    J = J
  )
  if (mixed_engine) {
    true_params$sigma <- if (is.null(scale_mat)) sigma else scale_mat
    if (any(engine_margins %in% c("exponential", "gamma"))) {
      true_params$skew_direction <- skew_direction
    }
    if (any(engine_margins %in% c("skew_normal", "gamma"))) {
      true_params$skew_params <- skew_params
    }
  } else if (identical(margins, "normal")) {
    true_params$sigma <- if (is.null(scale_mat)) sigma else scale_mat
  } else {
    true_params$sigma_exp <- if (is.null(scale_mat)) sigma_exp else scale_mat
    true_params$skew_direction <- skew_direction
  }

  list(
    data = data,
    true_params = true_params,
    latent_states = state,
    innovations = zeta
  )
}
