# ============================================================================
# Trajectory Comparison Plot
# ============================================================================

#' Plot and compare multiple rho trajectory shapes
#'
#' Visualises several named trajectory scenarios side by side for comparison.
#'
#' @param n_time Number of time points.
#' @param scenarios Character vector of scenario names (see [rho_scenario()]).
#'   Default: all built-in scenarios.
#' @param ... Additional arguments passed to [rho_scenario()]. Each argument is
#'   forwarded only to the scenarios whose generator accepts it; arguments not
#'   used by any requested scenario produce a warning.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' plot_trajectories(100)
#' plot_trajectories(100, scenarios = c("decreasing", "single_middle"))
plot_trajectories <- function(n_time,
                              scenarios = c("constant", "decreasing", "increasing",
                                            "random_walk", "single_middle",
                                            "large_change", "double_relapse"),
                              ...) {
  unknown <- setdiff(scenarios, names(.rho_scenarios))
  if (length(unknown) > 0) {
    cli_abort("Unknown scenario{?s}: {.val {unknown}}. Available: {.val {names(.rho_scenarios)}}")
  }

  extra <- list(...)
  accepted <- lapply(scenarios, function(s) names(formals(.rho_scenarios[[s]]$fn)))
  unused <- setdiff(names(extra), unlist(accepted))
  if (length(unused) > 0) {
    cli_warn("Ignoring argument{?s} {.arg {unused}}: not used by any requested scenario.")
  }

  rows <- lapply(seq_along(scenarios), function(i) {
    s <- scenarios[[i]]
    # Forward only the arguments this scenario's generator accepts; passing
    # the shared dots verbatim would error on any scenario-specific argument.
    args <- extra[intersect(names(extra), accepted[[i]])]
    rho <- do.call(rho_scenario, c(list(scenario = s, n_time = n_time), args))
    data.frame(
      time = 2:(length(rho) + 1),
      rho = rho,
      scenario = s,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df$scenario <- factor(df$scenario, levels = scenarios)

  ggplot2::ggplot(df, ggplot2::aes(
    x = .data$time, y = .data$rho, color = .data$scenario
  )) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::labs(
      x = "Time",
      y = expression(rho[t]),
      color = "Scenario",
      title = "Rho Trajectory Comparison"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}
