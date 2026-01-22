#' BRC Theme and Styling Helpers
#'
#' Reusable theme components, labels, and scales for consistent visualization
#' across BRC core analysis plots.
#'
#' @name brc_themes
#' @author Bolívar Aponte Rolón
#' @date 2025-01-20

# Theme Components ----

#' Standard Title and Text Sizing Theme
#'
#' Provides consistent text sizing for titles, axes, and legends
#'
#' @param title_size Numeric, size for plot titles (default: 12)
#' @param legend_title_size Numeric, size for legend titles (default: 12)
#' @param axis_title_size Numeric, size for axis titles (default: 14)
#' @param axis_text_size Numeric, size for axis text (default: 12)
#' @param legend_text_size Numeric, size for legend text (default: 12)
#'
#' @return A ggplot2 theme object
#' # For scripts, not packaged
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   brc_theme_title_size()
brc_theme_title_size <- function(
  title_size = 12,
  legend_title_size = 12,
  axis_title_size = 14,
  axis_text_size = 12,
  legend_text_size = 12
) {
  ggplot2::theme(
    title = ggplot2::element_text(size = title_size),
    legend.title = ggplot2::element_text(
      size = legend_title_size,
      face = "bold"
    ),
    axis.title.x = ggplot2::element_text(size = axis_title_size, face = "bold"),
    axis.title.y = ggplot2::element_text(size = axis_title_size, face = "bold"),
    axis.text.x = ggplot2::element_text(size = axis_text_size),
    axis.text.y = ggplot2::element_text(size = axis_text_size),
    legend.text = ggplot2::element_text(size = legend_text_size)
  )
}

#' Paper-Ready Theme Package
#'
#' Comprehensive theme for publication-ready figures
#'
#' @return A list of ggplot2 theme components
#' # For scripts, not packaged
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   brc_theme_paper()
brc_theme_paper <- function() {
  list(
    ggplot2::theme_bw(base_size = 12),
    brc_theme_title_size()
  )
}

# Label Helpers ----

#' Standard BRC Crop Labels
#'
#' Returns standardized crop labels for consistent naming across plots
#'
#' @return Character vector of crop labels
#' # For scripts, not packaged
#'
#' @examples
#' brc_crop_labels()
brc_crop_labels <- function() {
  c(
    "Corn",
    "Miscanthus",
    "Poplar",
    "Restored Prairie",
    "Sorghum",
    "Sorghum + Rye",
    "Soy",
    "Switchgrass"
  )
}

#' BRC Name Labels (Abbreviated to Full Names)
#'
#' Convert abbreviated BRC names to full names
#'
#' @return Named character vector for use with as_labeller()
#' # For scripts, not packaged
#'
#' @examples
#' library(ggplot2)
#' # Use with facet_wrap
#' # facet_wrap(~brc, labeller = as_labeller(brc_name_labels()))
brc_name_labels <- function() {
  c(
    cabbi = "CABBI",
    cbi = "CBI",
    jbei = "JBEI",
    glbrc = "GLBRC"
  )
}

#' Full BRC and Crop Labels for Faceting
#'
#' Combined labeller for nested faceting with BRC and crop
#'
#' @return Named character vector for use with as_labeller()
#' # For scripts, not packaged
#'
#' @examples
#' library(ggplot2)
#' # Use with facet_nested
#' # facet_nested(~brc + crop, labeller = as_labeller(brc_full_labels()))
brc_full_labels <- function() {
  c(
    cabbi = "CABBI",
    cbi = "CBI",
    jbei = "JBEI",
    glbrc = "GLBRC",
    Sorghum = "Sorghum",
    poplar = "Poplar",
    Restored_Prairie = "Restored Prairie",
    Switchgrass = "Switchgrass",
    Corn = "Corn",
    Miscanthus = "Miscanthus",
    Soy = "Soy",
    `Sorghum + Rye` = "Sorghum + Rye"
  )
}

# Scale Helpers ----

#' Percent Scale for Y-Axis (0-100%)
#'
#' Standardized y-axis scale for percentage data
#'
#' @param limits Numeric vector of length 2 (default: c(0, 1))
#' @param breaks Numeric vector of break points (default: waiver())
#' @param suffix Character to append after percentage (default: "")
#'
#' @return A ggplot2 scale object
#' # For scripts, not packaged
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(x = 1:10, y = seq(0, 1, length.out = 10))
#' ggplot(df, aes(x, y)) +
#'   geom_line() +
#'   brc_scale_y_percent()
brc_scale_y_percent <- function(
  limits = c(0, 1),
  breaks = ggplot2::waiver(),
  suffix = ""
) {
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(suffix = suffix),
    limits = limits,
    breaks = breaks
  )
}

#' Color Scale for BRC Set2 Palette
#'
#' Standard color palette for BRC crops
#'
#' @param labels Optional labels for the color scale
#' @param ... Additional arguments passed to scale_color_brewer
#'
#' @return A ggplot2 scale object
#' # For scripts, not packaged
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt, color = factor(cyl))) +
#'   geom_point() +
#'   brc_scale_color_crops()
brc_scale_color_crops <- function(labels = brc_crop_labels(), ...) {
  ggplot2::scale_color_brewer(
    palette = "Set2",
    labels = labels,
    ...
  )
}

#' Fill Scale for Core vs Non-Core ASVs
#'
#' Two-color fill scale for distinguishing core from non-core ASVs
#'
#' @param core_color Color for core ASVs (default: "#377EB8")
#' @param noncore_color Color for non-core ASVs (default: "#E41A1C")
#' @param labels Labels for legend (default: c("Core", "Non-core"))
#'
#' @return A ggplot2 scale object
#' # For scripts, not packaged
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(x = 1:10, y = 1:10, type = rep(c("Core", "Non-core"), 5))
#' ggplot(df, aes(x, y, fill = type)) +
#'   geom_point(shape = 21, size = 3) +
#'   brc_scale_fill_core()
brc_scale_fill_core <- function(
  core_color = "#E41A1C",
  noncore_color = "#377EB8",
  labels = c("Core", "Non-core")
) {
  ggplot2::scale_fill_manual(
    labels = labels,
    values = c(core_color, noncore_color)
  )
}

# Pre-configured Theme Lists ----

#' Pre-configured Theme List for 50 Core ASVs (Crops)
#'
#' @return A list of ggplot2 components
#' # For scripts, not packaged
brc_theme_core50_crops <- function() {
  list(
    ggplot2::ggtitle("50 core ASVs in BRC crops"),
    brc_scale_color_crops(),
    brc_theme_title_size()
  )
}

#' Pre-configured Theme List for 50 Core ASVs (BRCs)
#'
#' @return A list of ggplot2 components
#' # For scripts, not packaged
brc_theme_core50_brc <- function() {
  list(
    ggplot2::ggtitle("50 core ASVs in BRCs"),
    brc_theme_title_size()
  )
}

#' Pre-configured Theme List for 60% Threshold Core
#'
#' @return A list of ggplot2 components
#' # For scripts, not packaged
brc_theme_threshold60_core <- function() {
  list(
    ggplot2::ggtitle("ASVs in >60% Samples Across Crops and BRCs"),
    ggplot2::guides(fill = ggplot2::guide_legend(title = "ASV Association")),
    brc_theme_title_size()
  )
}

#' Pre-configured Theme List for 60% Threshold (Crops)
#'
#' @return A list of ggplot2 components
#' # For scripts, not packaged
brc_theme_threshold60_crops <- function() {
  list(
    ggplot2::ggtitle("12 high prevalence ASVs in BRC crops"),
    brc_scale_color_crops(),
    brc_theme_title_size()
  )
}

#' Pre-configured Theme List for 60% Threshold (BRCs)
#'
#' @return A list of ggplot2 components
#' # For scripts, not packaged
brc_theme_threshold60_brc <- function() {
  list(
    ggplot2::ggtitle("12 high prevalence ASVs in BRCs"),
    brc_theme_title_size()
  )
}

#' Pre-configured Theme List for 100% Threshold (All ASVs)
#'
#' @return A list of ggplot2 components
#' # For scripts, not packaged
brc_theme_threshold100_palette <- function() {
  list(
    brc_scale_color_crops(),
    brc_theme_title_size()
  )
}

#' Override Scale for Core/Non-Core with Single Color
#'
#' Special scale that uses same color for both core and non-core
#' (used in specific visualizations where distinction isn't needed)
#'
#' @return A list of ggplot2 components
#' # For scripts, not packaged
brc_scale_override_single <- function() {
  list(
    brc_scale_y_percent(),
    ggplot2::scale_fill_manual(
      labels = c("Core", "Non-core"),
      values = c("#377EB8", "#377EB8")
    )
  )
}
