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
#' @param T Number of time points.
#' @param J Number of indicators per latent variable.
#' @param lambda Numeric vector of length J with factor loadings.
#' @param sigma_e Measurement error SD (scalar).
#' @param Phi 2x2 VAR coefficient matrix.
#' @param mu Length-2 intercept vector.
#' @param sigma Length-2 innovation SD vector.
#' @param rho Copula correlation.
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
#'   - `latent_states`: T x 2 matrix of true latent states
#'   - `innovations`: T x 2 matrix of true innovations
#' @export
simulate_dcvar_sem <- function(T = 200, J = 3,
                                lambda = rep(sqrt(0.8), 3),
                                sigma_e = sqrt(0.2),
                                Phi = matrix(c(0.5, 0.15, 0.15, 0.3), 2, 2),
                                mu = c(0, 0),
                                sigma = c(1, 1),
                                rho = 0.3,
                                burnin = 0,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.numeric(T) || length(T) != 1L || T != as.integer(T) || T < 1) {
    cli_abort("{.arg T} must be an integer >= 1, got {.val {T}}.")
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
  .simulate_sem_validate_numeric_vector(mu, "mu")
  if (length(mu) != 2L) {
    cli_abort("{.arg mu} must have length 2, got {.val {length(mu)}}.")
  }
  .simulate_sem_validate_numeric_vector(sigma, "sigma")
  if (length(sigma) != 2L) {
    cli_abort("{.arg sigma} must have length 2, got {.val {length(sigma)}}.")
  }
  if (any(sigma <= 0)) {
    cli_abort("{.arg sigma} values must be positive.")
  }
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho < -1 || rho > 1) {
    cli_abort("{.arg rho} must be a single finite numeric value in [-1, 1].")
  }

  if (!identical(burnin, 0L) && !identical(burnin, 0)) {
    cli_warn(
      c(
        "{.arg burnin} is ignored for SEM simulation.",
        "i" = "The fitted SEM model treats the first returned state as the observed draw implied by `x_0 = 0`."
      )
    )
  }

  # Generate correlated innovations via Gaussian copula
  Sigma_copula <- matrix(c(1, rho, rho, 1), 2, 2)
  L <- t(chol(Sigma_copula))
  zeta <- matrix(NA_real_, T, 2)
  for (t in seq_len(T)) {
    z <- rnorm(2)
    zeta[t, ] <- (L %*% z) * sigma
  }

  # Latent VAR(1) recursion matching the SEM Stan model, which conditions on x_0 = 0.
  state <- matrix(0, T, 2)
  state[1, ] <- mu + zeta[1, ]
  for (t in 2:T) {
    state[t, ] <- mu + as.vector(Phi %*% state[t - 1, ]) + zeta[t, ]
  }

  # Measurement model: y_{ij,t} = lambda_j * state_{i,t} + e_{ij,t}
  y <- matrix(NA_real_, T, 2 * J)
  for (j in seq_len(J)) {
    y[, j]     <- lambda[j] * state[, 1] + rnorm(T, 0, sigma_e)
    y[, J + j] <- lambda[j] * state[, 2] + rnorm(T, 0, sigma_e)
  }
  colnames(y) <- c(paste0("y1_", seq_len(J)), paste0("y2_", seq_len(J)))

  data <- data.frame(time = seq_len(T), y, check.names = FALSE)

  list(
    data = data,
    true_params = list(
      Phi = Phi,
      mu = mu,
      sigma = sigma,
      rho = rho,
      lambda = lambda,
      sigma_e = sigma_e,
      J = J
    ),
    latent_states = state,
    innovations = zeta
  )
}
