# ============================================================================
# Breakpoint Data Simulation Convenience
# ============================================================================

#' Simulate data with a breakpoint rho trajectory
#'
#' Convenience wrapper that combines [rho_step()] or [rho_double_step()] with
#' [simulate_dcVar()] for quick breakpoint simulation studies.
#'
#' @param T Number of time points.
#' @param type Character; one of `"single"` (single breakpoint) or `"double"`
#'   (double breakpoint / relapse pattern).
#' @param rho_before Rho before breakpoint (default: 0.7).
#' @param rho_after Rho after breakpoint (default: 0.3).
#' @param rho_levels Numeric vector of three rho levels for double breakpoint
#'   (default: `c(0.7, 0.3, 0.7)`). Only used when `type = "double"`.
#' @param breakpoint Breakpoint location as proportion of T-1 (default: 0.5).
#' @param breakpoints Numeric vector of two breakpoints for double type
#'   (default: `c(1/3, 2/3)`).
#' @param transition_width Number of time points for smooth transition
#'   (default: 0 = abrupt).
#' @param mu Intercept vector of length 2 (default: `c(0, 0)`).
#' @param Phi VAR(1) coefficient matrix, 2x2.
#' @param sigma_eps Innovation SDs, length 2 (default: `c(1, 1)`).
#' @param seed Random seed.
#'
#' @return A named list as returned by [simulate_dcVar()].
#' @export
#'
#' @examples
#' sim <- simulate_breakpoint_data(T = 100, type = "single", seed = 42)
#' plot(sim$true_params$rho, type = "l")
simulate_breakpoint_data <- function(T,
                                     type = c("single", "double"),
                                     rho_before = 0.7,
                                     rho_after = 0.3,
                                     rho_levels = c(0.7, 0.3, 0.7),
                                     breakpoint = 0.5,
                                     breakpoints = c(1 / 3, 2 / 3),
                                     transition_width = 0,
                                     mu = c(0, 0),
                                     Phi = matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2),
                                     sigma_eps = c(1, 1),
                                     seed = NULL) {
  type <- match.arg(type)

  rho_traj <- if (type == "single") {
    rho_step(T, rho_before = rho_before, rho_after = rho_after,
             breakpoint = breakpoint, transition_width = transition_width)
  } else {
    rho_double_step(T, rho_levels = rho_levels,
                    breakpoints = breakpoints,
                    transition_width = transition_width)
  }

  simulate_dcVar(T = T, rho_trajectory = rho_traj,
                 mu = mu, Phi = Phi, sigma_eps = sigma_eps, seed = seed)
}
