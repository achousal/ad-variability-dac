# theme_research.R -- ggplot2 theme matching _code/styles/STYLE_GUIDE.md + PLOT_DESIGN.md
#
# Provides theme_research() for consistent styling across all R plots.
# Anchored to the canonical style guide and plot design specs in the vault.

# -- Standard figure sizes (inches) -------------------------------------------

#' Canonical figure sizes by plot type
FIGURE_SIZES <- list(
  box_grouped     = c(width = 8,  height = 6),
  box_single      = c(width = 6,  height = 6),
  violin_grouped  = c(width = 8,  height = 6),
  violin_single   = c(width = 6,  height = 6),
  scatter_multi   = c(width = 18, height = 7),
  scatter_single  = c(width = 14, height = 8),
  scatter_bivar   = c(width = 10, height = 6),
  heatmap         = c(width = 10, height = 6),
  forest          = c(width = 10, height = 8),
  volcano         = c(width = 10, height = 8),
  roc             = c(width = 7,  height = 7),
  bar             = c(width = 8,  height = 6)
)

#' Research theme for ggplot2
#'
#' Clean publication-ready theme based on theme_classic with:
#' - Bold titles and strip text
#' - Configurable strip background color
#' - Configurable legend position
#' - 14pt base font size
#'
#' @param base_size Base font size (default 14).
#' @param strip_color Fill color for facet strip backgrounds (default "grey90").
#' @param legend_position Legend position (default "bottom").
#' @return A ggplot2 theme object.
#' @export
theme_research <- function(base_size = 14,
                           strip_color = "grey90",
                           legend_position = "bottom") {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      # Title
      plot.title = ggplot2::element_text(face = "bold"),
      # Strip (facet labels)
      strip.background = ggplot2::element_rect(fill = strip_color, colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      # Legend
      legend.position = legend_position,
      # Clean axes
      axis.line = ggplot2::element_line(colour = "black"),
      panel.grid = ggplot2::element_blank()
    )
}

#' Get canonical figure size for a plot type
#'
#' @param plot_type Character string matching a key in FIGURE_SIZES.
#' @return Named numeric vector with "width" and "height".
#' @export
get_figure_size <- function(plot_type) {
  if (!plot_type %in% names(FIGURE_SIZES)) {
    stop(
      "Unknown plot type: '", plot_type, "'. Choose from: ",
      paste(names(FIGURE_SIZES), collapse = ", ")
    )
  }
  FIGURE_SIZES[[plot_type]]
}
