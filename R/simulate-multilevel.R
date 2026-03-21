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
#' @param T Number of time points per unit.
#' @param phi_bar Population mean for VAR coefficients (length-4 vector:
#'   phi11, phi12, phi21, phi22).
#' @param tau_phi Population SD for each VAR coefficient (length-4 vector).
#' @param sigma Innovation SDs (length-2 vector).
#' @param rho Global copula correlation.
#' @param burnin Number of burn-in observations to discard (default: 30).
#' @param center Logical; person-mean center the data (default: `TRUE`).
#' @param seed Random seed for reproducibility.
#'
#' @return A named list with:
#'   - `data`: panel data frame with columns `id`, `time`, `y1`, `y2`
#'   - `true_params`: list of true parameter values
#'   - `person_means`: N x 2 matrix of person means (before centering)
#' @export
simulate_dcvar_multilevel <- function(N = 40, T = 100,
                                      phi_bar = c(0.3, 0.1, 0.1, 0.3),
                                      tau_phi = c(0.1, 0.05, 0.05, 0.1),
                                      sigma = c(1, 1),
                                      rho = 0.3,
                                      burnin = 30,
                                      center = TRUE,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.numeric(N) || length(N) != 1L || N != as.integer(N) || N < 1) {
    cli_abort("{.arg N} must be an integer >= 1, got {.val {N}}.")
  }
  if (!is.numeric(T) || length(T) != 1L || T != as.integer(T) || T < 2) {
    cli_abort("{.arg T} must be an integer >= 2, got {.val {T}}.")
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
  if (!is.numeric(rho) || length(rho) != 1 || rho < -1 || rho > 1) {
    cli_abort("{.arg rho} must be a single numeric value in [-1, 1], got {.val {rho}}.")
  }

  T_sim <- T + burnin
  Phi_mat <- matrix(NA_real_, N, 4)
  Phi_list <- vector("list", N)

  # Generate unit-specific VAR matrices
  for (i in seq_len(N)) {
    Phi_mat[i, ] <- rnorm(4, phi_bar, tau_phi)
    Phi_i <- matrix(Phi_mat[i, ], 2, 2, byrow = TRUE)
    Phi_list[[i]] <- Phi_i
  }

  # Simulate data
  y_raw_list <- vector("list", N)
  for (i in seq_len(N)) {
    Y <- matrix(0, T_sim, 2)
    for (t in 2:T_sim) {
      # Correlated innovations via Cholesky
      L <- matrix(c(1, rho, 0, sqrt(1 - rho^2)), 2, 2)
      z <- rnorm(2)
      eps <- L %*% z * sigma

      Y[t, ] <- Phi_list[[i]] %*% Y[t - 1, ] + eps
    }
    y_raw_list[[i]] <- Y[(burnin + 1):T_sim, , drop = FALSE]
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
      time = seq_len(T),
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
      sigma = sigma,
      rho = rho,
      Phi_mat = Phi_mat,
      Phi_list = Phi_list
    ),
    person_means = person_means
  )
}


#' Compute the spectral radius of a 2x2 VAR coefficient matrix
#' @noRd
.spectral_radius <- function(Phi) {
  max(Mod(eigen(Phi, only.values = TRUE)$values))
}


#' Project VAR matrix to ensure stationarity
#' @noRd
.project_if_needed <- function(Phi, alpha = 0.995) {
  ev <- eigen(Phi, only.values = TRUE)$values
  r <- max(Mod(ev))
  if (is.na(r) || r < 1) return(Phi)
  Phi * (alpha / r)
}
