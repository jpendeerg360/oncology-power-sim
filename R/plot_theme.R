# Shared ggplot theme and colours. Categorical hues are used in fixed slot order.

ARM_COLOURS <- c("Obs" = "#2a78d6", "Lev" = "#eb6834", "Lev+5FU" = "#1baf7a")
SERIES_COLOURS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")

theme_sim <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(colour = "grey35"),
      axis.title = ggplot2::element_text(colour = "grey25"),
      plot.title.position = "plot"
    )
}
