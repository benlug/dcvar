# ============================================================================
# Trajectory Comparison Plot
# ============================================================================

#' Plot and compare multiple rho trajectory shapes
#'
#' Visualises several named trajectory scenarios side by side for comparison.
#'
#' @param T Number of time points.
#' @param scenarios Character vector of scenario names (see [rho_scenario()]).
#'   Default: all built-in scenarios.
#' @param ... Additional arguments passed to [rho_scenario()].
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' plot_trajectories(100)
#' plot_trajectories(100, scenarios = c("decreasing", "single_middle"))
plot_trajectories <- function(T,
                              scenarios = c("constant", "decreasing", "increasing",
                                            "random_walk", "single_middle",
                                            "large_change", "double_relapse"),
                              ...) {
  rows <- lapply(scenarios, function(s) {
    rho <- rho_scenario(s, T = T, ...)
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
