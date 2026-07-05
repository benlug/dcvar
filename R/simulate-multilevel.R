# ============================================================================
# Multilevel Simulation
# ============================================================================

#' Simulate data from a multilevel copula VAR(1) model
#'
#' Generates panel data with unit-specific VAR coefficients drawn from
#' a population distribution and a global copula correlation. The simulator
#' matches the fitted multilevel model support by leaving unit-level VAR
#' matrices unconstrained; nonstationary draws are possible.
#'
#' @param N Number of units.
#' @param n_time Number of time points per unit.
#' @param phi_bar Population mean for VAR coefficients (length-4 vector:
#'   phi11, phi12, phi21, phi22).
#' @param tau_phi Population SD for each VAR coefficient (length-4 vector).
#' @param sigma Innovation scale vector (length 2). Values are on each
#'   family's natural residual scale: innovation SD for normal/skew-normal,
#'   `sigma_exp` for exponential, and `sigma_gam` for gamma.
#' @param rho Global copula correlation.
#' @param margins Marginal family. Either a single string applied to both
#'   variables, or a length-2 character vector for per-variable (mixed) margins,
#'   e.g. `c("normal", "exponential")`. Each entry is one of `"normal"`
#'   (default), `"exponential"`, `"skew_normal"`, or `"gamma"`.
#' @param skew_direction Length-2 `1`/`-1` vector. Required whenever any
#'   dimension uses an `"exponential"` or `"gamma"` margin.
#' @param skew_params Named list of margin-specific parameters: `alpha`
#'   (length-2 skew-normal shape) and/or `shape` (scalar gamma shape).
#' @param phi_trajectory Optional population-level time-varying VAR coefficient
#'   path. Unit-specific coefficients follow this shared drift around their
#'   sampled baselines.
#' @param sigma_trajectory Optional shared time-varying innovation scale path.
#'   The supplied value is each family's natural scale (innovation SD for
#'   normal, residual SD for skew-normal, `sigma_exp` / `sigma_gam` for
#'   exponential / gamma).
#' @param tv_sigma_k Soft-barrier sharpness for time-varying exponential/gamma
#'   scales.
#' @param burnin Number of burn-in observations to discard (default: 30).
#' @param center Logical; person-mean center the data (default: `TRUE`).
#' @param seed Random seed for reproducibility.
#'
#' @return A named list with:
#'   - `data`: panel data frame with columns `id`, `time`, `y1`, `y2`
#'   - `true_params`: list of true parameter values, including `phi_bar`,
#'     `tau_phi`, `sigma`, `rho`, `margins`, `skew_direction`, `skew_params`,
#'     the per-unit VAR coefficients `Phi_mat` (an `N x 4` matrix) and
#'     `Phi_list` (a length-`N` list of `2 x 2` matrices), `Phi_population`
#'     (the shared population VAR path: `phi_bar` when constant, or an
#'     `(n_time - 1) x 4` matrix when `phi_trajectory` is supplied), and
#'     `Phi_unit_paths` (a length-`N` list of per-unit effective coefficient
#'     paths, each `(n_time - 1) x 4`)
#'   - `person_means`: N x 2 matrix of person means (before centering)
#' @export
simulate_dcvar_multilevel <- function(N = 40, n_time = 100,
                                      phi_bar = c(0.3, 0.1, 0.1, 0.3),
                                      tau_phi = c(0.1, 0.05, 0.05, 0.1),
                                      sigma = c(1, 1),
                                      rho = 0.3,
                                      margins = "normal",
                                      skew_direction = NULL,
                                      skew_params = NULL,
                                      phi_trajectory = NULL,
                                      sigma_trajectory = NULL,
                                      tv_sigma_k = NULL,
                                      burnin = 30,
                                      center = TRUE,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  margins <- .normalize_margins_spec(margins)
  .validate_margins(margins, skew_direction)
  margins_vec <- if (length(margins) == 1L) rep(margins, 2L) else margins
  skew_params <- if (is.list(skew_params)) skew_params else list()
  if (any(margins_vec == "skew_normal")) {
    skew_params$alpha <- skew_params$alpha %||% rep(0, 2L)
  }
  if (any(margins_vec == "gamma")) {
    skew_params$shape <- skew_params$shape %||% 1
  }
  if (!is.numeric(N) || length(N) != 1L || N != as.integer(N) || N < 1) {
    cli_abort("{.arg N} must be an integer >= 1, got {.val {N}}.")
  }
  if (!is.numeric(n_time) || length(n_time) != 1L || n_time != as.integer(n_time) || n_time < 2) {
    cli_abort("{.arg n_time} must be an integer >= 2, got {.val {n_time}}.")
  }
  .simulate_validate_numeric_vector(phi_bar, "phi_bar")
  if (length(phi_bar) != 4L) {
    cli_abort("{.arg phi_bar} must have length 4, got {.val {length(phi_bar)}}.")
  }
  .simulate_validate_numeric_vector(tau_phi, "tau_phi")
  if (length(tau_phi) != 4L) {
    cli_abort("{.arg tau_phi} must have length 4, got {.val {length(tau_phi)}}.")
  }
  if (any(tau_phi < 0)) {
    cli_abort("{.arg tau_phi} values must be non-negative.")
  }
  .simulate_validate_numeric_vector(sigma, "sigma")
  if (length(sigma) != 2L) {
    cli_abort("{.arg sigma} must have length 2, got {.val {length(sigma)}}.")
  }
  if (any(sigma <= 0)) {
    cli_abort("{.arg sigma} values must be positive.")
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    cli_abort("{.arg center} must be TRUE or FALSE.")
  }
  if (!is.numeric(burnin) || length(burnin) != 1L ||
      burnin != as.integer(burnin) || burnin < 0) {
    cli_abort("{.arg burnin} must be a non-negative integer, got {.val {burnin}}.")
  }
  if (!is.numeric(rho) || length(rho) != 1 || !is.finite(rho) || rho < -1 || rho > 1) {
    cli_abort("{.arg rho} must be a single finite numeric value in [-1, 1], got {.val {rho}}.")
  }
  if (!is.null(sigma_trajectory) && !missing(sigma)) {
    cli_abort("Supply either {.arg sigma} (constant) or {.arg sigma_trajectory} (time-varying), not both.")
  }
  if (!is.null(tv_sigma_k)) {
    .simulate_validate_positive_scalar(tv_sigma_k, "tv_sigma_k")
  }

  n_time_sim <- n_time + burnin
  phi_pop <- if (is.null(phi_trajectory)) {
    matrix(rep(phi_bar, each = n_time - 1L), n_time - 1L, 4L)
  } else {
    .tv_resolve_trajectory(phi_trajectory, n_time, 4L, "phi_trajectory")
  }
  phi_dev_ret <- sweep(phi_pop, 2, phi_bar, "-")
  phi_dev_sim <- if (burnin > 0L) {
    rbind(matrix(rep(phi_dev_ret[1, ], each = burnin), burnin, 4L), phi_dev_ret)
  } else {
    phi_dev_ret
  }

  scale_ret <- if (is.null(sigma_trajectory)) {
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
            "{.arg sigma_trajectory} is not constant on dimension {d} ({.val {margins_vec[d]}} margin).",
            "i" = "Set {.arg tv_sigma_k} to simulate from the soft-barrier model, or keep the scale constant."
          ))
        }
      }
    }
    m
  }
  scale_sim <- if (is.null(scale_ret)) {
    NULL
  } else if (burnin > 0L) {
    rbind(matrix(rep(scale_ret[1, ], each = burnin), burnin, 2L), scale_ret)
  } else {
    scale_ret
  }
  base_scales <- sigma

  Phi_mat <- matrix(NA_real_, N, 4)
  Phi_list <- vector("list", N)
  Phi_unit_paths <- vector("list", N)

  # Generate unit-specific VAR matrices
  for (i in seq_len(N)) {
    Phi_mat[i, ] <- rnorm(4, phi_bar, tau_phi)
    Phi_i <- matrix(Phi_mat[i, ], 2, 2, byrow = TRUE)
    Phi_list[[i]] <- Phi_i
    Phi_unit_paths[[i]] <- sweep(phi_pop, 2, phi_bar - Phi_mat[i, ], "-")
  }

  # Simulate data
  y_raw_list <- vector("list", N)
  for (i in seq_len(N)) {
    Y <- matrix(0, n_time_sim, 2)
    for (time_index in 2:n_time_sim) {
      # Correlated standard normals via Cholesky, then per-dimension margins.
      L <- matrix(c(1, rho, 0, sqrt(1 - rho^2)), 2, 2)
      z <- rnorm(2)
      w <- as.numeric(L %*% z)
      scales_t <- if (is.null(scale_sim)) base_scales else scale_sim[time_index - 1L, ]
      eps <- .sim_marginal_quantile_scaled(
        w, margins_vec, scales_t, skew_direction, skew_params,
        barrier_k = tv_sigma_k
      )

      Phi_t <- matrix(Phi_mat[i, ] + phi_dev_sim[time_index - 1L, ], 2, 2, byrow = TRUE)
      Y[time_index, ] <- Phi_t %*% Y[time_index - 1L, ] + eps
    }
    y_raw_list[[i]] <- Y[(burnin + 1L):n_time_sim, , drop = FALSE]
  }

  # Person-mean center
  person_means <- matrix(NA_real_, N, 2)
  y_list <- vector("list", N)
  for (i in seq_len(N)) {
    person_means[i, ] <- colMeans(y_raw_list[[i]])
    y_list[[i]] <- if (center) {
      sweep(y_raw_list[[i]], 2, person_means[i, ])
    } else {
      y_raw_list[[i]]
    }
  }

  # Build panel data frame
  all_rows <- vector("list", N)
  for (i in seq_len(N)) {
    all_rows[[i]] <- data.frame(
      id = i,
      time = seq_len(n_time),
      y1 = y_list[[i]][, 1],
      y2 = y_list[[i]][, 2]
    )
  }
  panel <- do.call(rbind, all_rows)

  list(
    data = panel,
    true_params = list(
      phi_bar = phi_bar,
      tau_phi = tau_phi,
      sigma = if (is.null(scale_ret)) sigma else scale_ret,
      rho = rho,
      margins = margins,
      skew_direction = skew_direction,
      skew_params = skew_params,
      Phi_mat = Phi_mat,
      Phi_list = Phi_list,
      Phi_population = if (is.null(phi_trajectory)) phi_bar else phi_pop,
      Phi_unit_paths = Phi_unit_paths
    ),
    person_means = person_means
  )
}

