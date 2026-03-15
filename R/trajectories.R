# ============================================================================
# Trajectory Generation Functions
# ============================================================================

#' Generate a constant rho trajectory
#'
#' @param T Number of time points.
#' @param rho Constant correlation value (default: 0.5). Must be in \[-1, 1\].
#'
#' @return Numeric vector of length T-1.
#' @export
#'
#' @examples
#' rho_constant(100, rho = 0.5)
rho_constant <- function(T, rho = 0.5) {
  if (T < 2) cli_abort("{.arg T} must be >= 2, got {.val {T}}.")
  if (rho < -1 || rho > 1) {
    cli_abort("{.arg rho} must be between -1 and 1, got {rho}.")
  }
  rep(rho, T - 1)
}


#' Generate a logistically decreasing rho trajectory
#'
#' Mimics a therapy effect where coupling decreases from high to low.
#'
#' @param T Number of time points.
#' @param rho_start Starting rho value (default: 0.7).
#' @param rho_end Ending rho value (default: 0.3).
#' @param midpoint Time point of inflection (default: `T/2`).
#' @param steepness Controls transition sharpness (default: 0.05).
#'
#' @return Numeric vector of length T-1.
#' @export
#'
#' @examples
#' rho_decreasing(100)
rho_decreasing <- function(T, rho_start = 0.7, rho_end = 0.3,
                            midpoint = NULL, steepness = 0.05) {
  if (T < 2) cli_abort("{.arg T} must be >= 2, got {.val {T}}.")
  if (is.null(midpoint)) midpoint <- T / 2
  t_vec <- 1:(T - 1)
  rho_end + (rho_start - rho_end) / (1 + exp(steepness * (t_vec - midpoint)))
}


#' Generate a logistically increasing rho trajectory
#'
#' Mimics deterioration where coupling increases from low to high.
#'
#' @inheritParams rho_decreasing
#' @param rho_start Starting rho value (default: 0.3).
#' @param rho_end Ending rho value (default: 0.7).
#'
#' @return Numeric vector of length T-1.
#' @export
#'
#' @examples
#' rho_increasing(100)
rho_increasing <- function(T, rho_start = 0.3, rho_end = 0.7,
                            midpoint = NULL, steepness = 0.05) {
  if (T < 2) cli_abort("{.arg T} must be >= 2, got {.val {T}}.")
  if (is.null(midpoint)) midpoint <- T / 2
  t_vec <- 1:(T - 1)
  rho_start + (rho_end - rho_start) / (1 + exp(-steepness * (t_vec - midpoint)))
}


#' Generate a random walk rho trajectory on the Fisher-z scale
#'
#' Stochastic trajectory matching the DC-VAR data-generating process.
#'
#' @param T Number of time points.
#' @param z_init Initial value on Fisher-z scale (default: 0.5, corresponding
#'   to rho = 0.46).
#' @param sigma_omega Innovation SD for the random walk (default: 0.05).
#' @param seed Random seed for reproducibility.
#'
#' @return Numeric vector of length T-1.
#' @export
#'
#' @examples
#' rho_random_walk(100, seed = 42)
rho_random_walk <- function(T, z_init = 0.5, sigma_omega = 0.05, seed = NULL) {
  if (T < 2) cli_abort("{.arg T} must be >= 2, got {.val {T}}.")
  if (!is.null(seed)) set.seed(seed)
  if (sigma_omega < 0) {
    cli_abort("{.arg sigma_omega} must be non-negative, got {.val {sigma_omega}}.")
  }

  z <- numeric(T - 1)
  z[1] <- z_init + rnorm(1, 0, sigma_omega)
  if (T > 2) {
    for (t in 2:(T - 1)) {
      z[t] <- z[t - 1] + rnorm(1, 0, sigma_omega)
    }
  }

  # Clamp extreme z values to avoid rho numerically equal to +/- 1
  z <- pmax(pmin(z, 5), -5)

  tanh(z)
}


#' Generate a single-breakpoint (step function) rho trajectory
#'
#' Abrupt change from one rho level to another at a specified time.
#'
#' @param T Number of time points.
#' @param rho_before Rho before breakpoint (default: 0.7).
#' @param rho_after Rho after breakpoint (default: 0.3).
#' @param breakpoint Breakpoint location as a proportion of T-1 (if <= 1)
#'   or an absolute time index (default: 0.5).
#' @param transition_width Number of time points for smooth transition.
#'   0 = abrupt (default: 0).
#'
#' @return Numeric vector of length T-1.
#' @export
#'
#' @examples
#' rho_step(100, rho_before = 0.7, rho_after = 0.3)
rho_step <- function(T, rho_before = 0.7, rho_after = 0.3,
                      breakpoint = 0.5, transition_width = 0) {
  if (T < 2) cli_abort("{.arg T} must be >= 2, got {.val {T}}.")
  T_eff <- T - 1

  if (breakpoint <= 1) {
    bp_time <- round(breakpoint * T_eff)
  } else {
    bp_time <- round(breakpoint)
  }
  bp_time <- max(1, min(T_eff - 1, bp_time))

  t_vec <- 1:T_eff

  if (transition_width == 0) {
    ifelse(t_vec <= bp_time, rho_before, rho_after)
  } else {
    steep <- 10 / transition_width
    rho_after + (rho_before - rho_after) / (1 + exp(steep * (t_vec - bp_time)))
  }
}


#' Generate a double-breakpoint (relapse pattern) rho trajectory
#'
#' Three-phase trajectory: level A -> level B -> level C.
#'
#' @param T Number of time points.
#' @param rho_levels Numeric vector of three rho levels (default:
#'   `c(0.7, 0.3, 0.7)`).
#' @param breakpoints Numeric vector of two breakpoint positions (default:
#'   `c(1/3, 2/3)`). Interpreted as proportions of T-1 if <= 1.
#' @param transition_width Number of time points for smooth transitions
#'   (default: 0).
#'
#' @return Numeric vector of length T-1.
#' @export
#'
#' @examples
#' rho_double_step(100, rho_levels = c(0.7, 0.3, 0.7))
rho_double_step <- function(T, rho_levels = c(0.7, 0.3, 0.7),
                             breakpoints = c(1 / 3, 2 / 3),
                             transition_width = 0) {
  if (T < 2) cli_abort("{.arg T} must be >= 2, got {.val {T}}.")
  T_eff <- T - 1

  bp1 <- if (breakpoints[1] <= 1) round(breakpoints[1] * T_eff) else round(breakpoints[1])
  bp2 <- if (breakpoints[2] <= 1) round(breakpoints[2] * T_eff) else round(breakpoints[2])
  bp1 <- max(1, min(T_eff - 2, bp1))
  bp2 <- max(bp1 + 1, min(T_eff - 1, bp2))

  t_vec <- 1:T_eff

  if (transition_width == 0) {
    ifelse(t_vec <= bp1, rho_levels[1],
           ifelse(t_vec <= bp2, rho_levels[2], rho_levels[3]))
  } else {
    steep <- 10 / transition_width
    trans1 <- rho_levels[2] + (rho_levels[1] - rho_levels[2]) / (1 + exp(steep * (t_vec - bp1)))
    trans2 <- rho_levels[3] + (rho_levels[2] - rho_levels[3]) / (1 + exp(steep * (t_vec - bp2)))
    ifelse(t_vec <= (bp1 + bp2) / 2, trans1, trans2)
  }
}


#' Get a named trajectory scenario
#'
#' Convenience function to retrieve a standard scenario by name.
#'
#' @param scenario Character string. One of:
#'   - Smooth: `"constant"`, `"decreasing"`, `"increasing"`, `"random_walk"`
#'   - Step: `"single_middle"`, `"large_change"`, `"small_change"`,
#'     `"increase"`, `"double_relapse"`
#' @param T Number of time points.
#' @param ... Additional arguments passed to the generator.
#'
#' @return Numeric vector of length T-1.
#' @export
#'
#' @examples
#' rho_scenario("decreasing", T = 100)
#' rho_scenario("double_relapse", T = 150)
rho_scenario <- function(scenario, T, ...) {
  scenarios <- list(
    constant = list(fn = rho_constant, defaults = list(rho = 0.5)),
    decreasing = list(fn = rho_decreasing, defaults = list(rho_start = 0.7, rho_end = 0.3)),
    increasing = list(fn = rho_increasing, defaults = list(rho_start = 0.3, rho_end = 0.7)),
    random_walk = list(fn = rho_random_walk, defaults = list(z_init = 0.5, sigma_omega = 0.05)),
    single_middle = list(fn = rho_step, defaults = list(rho_before = 0.7, rho_after = 0.3, breakpoint = 0.5)),
    large_change = list(fn = rho_step, defaults = list(rho_before = 0.8, rho_after = 0.2, breakpoint = 0.5)),
    small_change = list(fn = rho_step, defaults = list(rho_before = 0.6, rho_after = 0.4, breakpoint = 0.5)),
    increase = list(fn = rho_step, defaults = list(rho_before = 0.3, rho_after = 0.7, breakpoint = 0.5)),
    double_relapse = list(fn = rho_double_step, defaults = list(rho_levels = c(0.7, 0.3, 0.7)))
  )

  if (!scenario %in% names(scenarios)) {
    cli_abort("Unknown scenario: {.val {scenario}}. Available: {.val {names(scenarios)}}")
  }

  scenario_spec <- scenarios[[scenario]]
  args <- utils::modifyList(scenario_spec$defaults, list(...))
  do.call(scenario_spec$fn, c(list(T = T), args))
}
