# ============================================================================
# Trajectory Generation Functions
# ============================================================================

#' Internal: validate a rho scalar in [-1, 1]
#' @noRd
.trajectory_validate_rho_scalar <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < -1 || x > 1) {
    cli_abort("{.arg {arg_name}} must be a single finite numeric value in [-1, 1].")
  }
}

#' Internal: validate a rho level vector in [-1, 1]
#' @noRd
.trajectory_validate_rho_levels <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 3L || any(!is.finite(x)) || any(x < -1 | x > 1)) {
    cli_abort("{.arg {arg_name}} must be a length-3 finite numeric vector in [-1, 1].")
  }
}

#' Internal: validate and resolve a breakpoint specification
#' @noRd
.trajectory_resolve_breakpoint <- function(breakpoint, T_eff, arg_name = "breakpoint") {
  if (!is.numeric(breakpoint) || length(breakpoint) != 1L || !is.finite(breakpoint)) {
    cli_abort("{.arg {arg_name}} must be a single finite numeric value.")
  }

  if (breakpoint >= 0 && breakpoint <= 1) {
    bp_time <- round(breakpoint * T_eff)
  } else {
    if (breakpoint != as.integer(breakpoint) || breakpoint < 1) {
      cli_abort(c(
        "{.arg {arg_name}} must be either a proportion in [0, 1] or a positive integer breakpoint.",
        "i" = "Use a whole number for an absolute breakpoint index."
      ))
    }
    bp_time <- as.integer(breakpoint)
  }

  max(1L, min(T_eff - 1L, as.integer(bp_time)))
}

#' Internal: validate a non-negative transition width
#' @noRd
.trajectory_validate_transition_width <- function(x, arg_name = "transition_width") {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
    cli_abort("{.arg {arg_name}} must be a single non-negative finite numeric value.")
  }
}

#' Generate a constant rho trajectory
#'
#' @param n_time Number of time points.
#' @param rho Constant correlation value (default: 0.5). Must be in \[-1, 1\].
#'
#' @return Numeric vector of length `n_time - 1`.
#' @export
#'
#' @examples
#' rho_constant(100, rho = 0.5)
rho_constant <- function(n_time, rho = 0.5) {
  if (n_time < 2) cli_abort("{.arg n_time} must be >= 2, got {.val {n_time}}.")
  .trajectory_validate_rho_scalar(rho, "rho")
  rep(rho, n_time - 1)
}


#' Generate a logistically decreasing rho trajectory
#'
#' Mimics a therapy effect where coupling decreases from high to low.
#'
#' @param n_time Number of time points.
#' @param rho_start Starting rho value (default: 0.7).
#' @param rho_end Ending rho value (default: 0.3).
#' @param midpoint Time point of inflection (default: `n_time/2`).
#' @param steepness Controls transition sharpness (default: 0.05).
#'
#' @return Numeric vector of length `n_time - 1`.
#' @export
#'
#' @examples
#' rho_decreasing(100)
rho_decreasing <- function(n_time, rho_start = 0.7, rho_end = 0.3,
                            midpoint = NULL, steepness = 0.05) {
  if (n_time < 2) cli_abort("{.arg n_time} must be >= 2, got {.val {n_time}}.")
  if (is.null(midpoint)) midpoint <- n_time / 2
  t_vec <- seq_len(n_time - 1L)
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
#' @return Numeric vector of length `n_time - 1`.
#' @export
#'
#' @examples
#' rho_increasing(100)
rho_increasing <- function(n_time, rho_start = 0.3, rho_end = 0.7,
                            midpoint = NULL, steepness = 0.05) {
  if (n_time < 2) cli_abort("{.arg n_time} must be >= 2, got {.val {n_time}}.")
  if (is.null(midpoint)) midpoint <- n_time / 2
  t_vec <- seq_len(n_time - 1L)
  rho_start + (rho_end - rho_start) / (1 + exp(-steepness * (t_vec - midpoint)))
}


#' Generate a random walk rho trajectory on the Fisher-z scale
#'
#' Stochastic trajectory matching the DC-VAR data-generating process.
#'
#' @param n_time Number of time points.
#' @param z_init Initial value on Fisher-z scale (default: 0.5, corresponding
#'   to rho = 0.46).
#' @param sigma_omega Innovation SD for the random walk (default: 0.05).
#' @param seed Random seed for reproducibility.
#'
#' @return Numeric vector of length `n_time - 1`.
#' @export
#'
#' @examples
#' rho_random_walk(100, seed = 42)
rho_random_walk <- function(n_time, z_init = 0.5, sigma_omega = 0.05, seed = NULL) {
  if (n_time < 2) cli_abort("{.arg n_time} must be >= 2, got {.val {n_time}}.")
  if (!is.null(seed)) set.seed(seed)
  if (sigma_omega < 0) {
    cli_abort("{.arg sigma_omega} must be non-negative, got {.val {sigma_omega}}.")
  }

  z <- numeric(n_time - 1L)
  z[1] <- z_init + rnorm(1, 0, sigma_omega)
  if (n_time > 2) {
    for (time_index in 2:(n_time - 1L)) {
      z[time_index] <- z[time_index - 1L] + rnorm(1, 0, sigma_omega)
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
#' @param n_time Number of time points.
#' @param rho_before Rho before breakpoint (default: 0.7).
#' @param rho_after Rho after breakpoint (default: 0.3).
#' @param breakpoint Breakpoint location as a proportion of `n_time - 1` (if <= 1)
#'   or an absolute time index (default: 0.5).
#' @param transition_width Number of time points for smooth transition.
#'   0 = abrupt (default: 0).
#'
#' @return Numeric vector of length `n_time - 1`.
#' @export
#'
#' @examples
#' rho_step(100, rho_before = 0.7, rho_after = 0.3)
rho_step <- function(n_time, rho_before = 0.7, rho_after = 0.3,
                      breakpoint = 0.5, transition_width = 0) {
  if (n_time < 2) cli_abort("{.arg n_time} must be >= 2, got {.val {n_time}}.")
  .trajectory_validate_rho_scalar(rho_before, "rho_before")
  .trajectory_validate_rho_scalar(rho_after, "rho_after")
  n_time_eff <- n_time - 1L
  .trajectory_validate_transition_width(transition_width)
  bp_time <- .trajectory_resolve_breakpoint(breakpoint, n_time_eff)

  t_vec <- seq_len(n_time_eff)

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
#' @param n_time Number of time points.
#' @param rho_levels Numeric vector of three rho levels (default:
#'   `c(0.7, 0.3, 0.7)`).
#' @param breakpoints Numeric vector of two breakpoint positions (default:
#'   `c(1/3, 2/3)`). Interpreted as proportions of `n_time - 1` if <= 1.
#' @param transition_width Number of time points for smooth transitions
#'   (default: 0).
#'
#' @return Numeric vector of length `n_time - 1`.
#' @export
#'
#' @examples
#' rho_double_step(100, rho_levels = c(0.7, 0.3, 0.7))
rho_double_step <- function(n_time, rho_levels = c(0.7, 0.3, 0.7),
                             breakpoints = c(1 / 3, 2 / 3),
                             transition_width = 0) {
  if (n_time < 2) cli_abort("{.arg n_time} must be >= 2, got {.val {n_time}}.")
  .trajectory_validate_rho_levels(rho_levels, "rho_levels")
  if (!is.numeric(breakpoints) || length(breakpoints) != 2L || any(!is.finite(breakpoints))) {
    cli_abort("{.arg breakpoints} must be a length-2 finite numeric vector.")
  }
  .trajectory_validate_transition_width(transition_width)
  n_time_eff <- n_time - 1L

  bp1 <- .trajectory_resolve_breakpoint(breakpoints[1], n_time_eff, "breakpoints[1]")
  bp2 <- .trajectory_resolve_breakpoint(breakpoints[2], n_time_eff, "breakpoints[2]")
  if (bp2 <= bp1) {
    cli_abort(c(
      "{.arg breakpoints} must resolve to two strictly increasing breakpoint positions.",
      "i" = "Choose distinct breakpoint values so the middle phase is non-empty."
    ))
  }

  t_vec <- seq_len(n_time_eff)

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
#' @param n_time Number of time points.
#' @param ... Additional arguments passed to the generator.
#'
#' @return Numeric vector of length `n_time - 1`.
#' @export
#'
#' @examples
#' rho_scenario("decreasing", n_time = 100)
#' rho_scenario("double_relapse", n_time = 150)
rho_scenario <- function(scenario, n_time, ...) {
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
  do.call(scenario_spec$fn, c(list(n_time = n_time), args))
}
