# palettes.R -- Canonical color palettes for all research projects
#
# Provides universal semantic palettes (binary, direction, significance),
# lab-specific categorical accent palettes, and ggplot2 scale helpers.
# Replaces per-project color definitions with a single source of truth.
#
# Palettes are loaded from _code/styles/palettes.yaml when available, with
# hardcoded fallbacks for offline / CI use. Domain profiles add lab-specific
# and domain-specific palettes via merge.
#
# Usage:
#   source("palettes.R")
#   ggplot(df, aes(x, y, color = direction)) + scale_color_direction()
#   ggplot(df, aes(x, y, fill = group)) + scale_fill_lab("your-lab")

# -- YAML loading with hardcoded fallback --------------------------------------

.load_palettes_yaml <- function() {
  # Find _code/styles/palettes.yaml relative to this file
  # Expected layout: _code/R/palettes.R -> _code/ is 1 level up
  script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
  )
  candidates <- c(
    file.path(script_dir, "..", "styles", "palettes.yaml"),
    file.path(getwd(), "_code", "styles", "palettes.yaml")
  )
  for (path in candidates) {
    path <- normalizePath(path, mustWork = FALSE)
    if (file.exists(path)) {
      if (requireNamespace("yaml", quietly = TRUE)) {
        tryCatch({
          data <- yaml::read_yaml(path)
          return(data)
        }, error = function(e) {
          message("Warning: failed to parse ", path, ": ", conditionMessage(e))
        })
      }
    }
  }
  return(NULL)
}

.yaml_data <- .load_palettes_yaml()

# -- Semantic palettes (universal, all labs) -----------------------------------

.get_semantic <- function(name, fallback) {
  if (!is.null(.yaml_data) && !is.null(.yaml_data$semantic[[name]])) {
    unlist(.yaml_data$semantic[[name]])
  } else {
    fallback
  }
}

#' Binary outcome palette (control vs case)
BINARY_COLORS <- .get_semantic("binary", c(
  "Control" = "#4DAF4A",
  "Case"    = "#E41A1C"
))

#' Direction palette (up/down regulation)
DIRECTION_COLORS <- .get_semantic("direction", c(
  "Up"   = "#E41A1C",
  "Down" = "#377EB8",
  "NS"   = "#999999"
))

#' Significance palette
SIG_COLORS <- .get_semantic("sig", c(
  "sig"     = "#E41A1C",
  "not sig" = "#999999"
))

#' Diverging heatmap palette name (RdBu for centered-at-zero data)
DIVERGING_PALETTE <- if (!is.null(.yaml_data$diverging)) .yaml_data$diverging else "RdBu"

#' Sequential heatmap palette name (Blues for density/counts)
SEQUENTIAL_PALETTE <- if (!is.null(.yaml_data$sequential)) .yaml_data$sequential else "Blues"

# -- Lab-specific categorical palettes -----------------------------------------
# Lab palettes loaded from palettes.yaml. Domain profiles add labs via merge.

#' Lab accent palettes for categorical variables
#' Loaded from _code/styles/palettes.yaml if available, otherwise empty.
LAB_PALETTES <- if (!is.null(.yaml_data$labs)) {
  lapply(.yaml_data$labs, unlist)
} else {
  list()
}

# -- Helper functions ----------------------------------------------------------

#' Get a lab-specific categorical palette
#'
#' @param lab Lab name (case-insensitive).
#' @param n Number of colors to return (default: all 8).
#' @return Character vector of hex colors.
#' @export
lab_palette <- function(lab, n = NULL) {
  lab_lower <- tolower(lab)
  if (!lab_lower %in% names(LAB_PALETTES)) {
    stop(
      "Unknown lab: '", lab, "'. Choose from: ",
      paste(names(LAB_PALETTES), collapse = ", ")
    )
  }
  pal <- LAB_PALETTES[[lab_lower]]
  if (!is.null(n)) {
    if (n > length(pal)) {
      stop("Requested ", n, " colors but palette has only ", length(pal))
    }
    pal <- pal[seq_len(n)]
  }
  pal
}

# -- ggplot2 scale helpers -----------------------------------------------------

#' Discrete color scale for direction (up/down/NS)
#' @param ... Additional arguments passed to scale_colour_manual.
#' @export
scale_color_direction <- function(...) {
  ggplot2::scale_colour_manual(values = DIRECTION_COLORS, ...)
}

#' Discrete fill scale for direction (up/down/NS)
#' @param ... Additional arguments passed to scale_fill_manual.
#' @export
scale_fill_direction <- function(...) {
  ggplot2::scale_fill_manual(values = DIRECTION_COLORS, ...)
}

#' Discrete color scale for significance
#' @param ... Additional arguments passed to scale_colour_manual.
#' @export
scale_color_sig <- function(...) {
  ggplot2::scale_colour_manual(values = SIG_COLORS, ...)
}

#' Discrete fill scale for significance
#' @param ... Additional arguments passed to scale_fill_manual.
#' @export
scale_fill_sig <- function(...) {
  ggplot2::scale_fill_manual(values = SIG_COLORS, ...)
}

#' Discrete color scale using a lab's categorical palette
#'
#' @param lab Lab name (case-insensitive).
#' @param ... Additional arguments passed to scale_colour_manual.
#' @export
scale_color_lab <- function(lab, ...) {
  ggplot2::scale_colour_manual(values = lab_palette(lab), ...)
}

#' Discrete fill scale using a lab's categorical palette
#'
#' @param lab Lab name (case-insensitive).
#' @param ... Additional arguments passed to scale_fill_manual.
#' @export
scale_fill_lab <- function(lab, ...) {
  ggplot2::scale_fill_manual(values = lab_palette(lab), ...)
}
