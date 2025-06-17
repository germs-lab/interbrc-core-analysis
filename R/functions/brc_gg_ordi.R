#' Standard Ordination Plot with ggplot2
#'
#' This function creates a standardized ordination plot using `ggplot2`.
#' It includes points, ellipses, and reference lines for visualizing NMDS or PCoA results. The function allows
#' customization of color, shape, and handling of missing values.
#'
#' @param .data A data frame containing ordination coordinates and additional metadata.
#' @param ordi A character string specifying the ordination method. Options are "NMDS" and "PcoA".
#' @param .color A column in `.data` to use for coloring points and ellipses. This should be a categorical variable.
#' @param .drop_na A column in `.data` used to filter out rows with missing values. Rows with `NA` in this column will be removed.
#'
#' @return A `ggplot` object representing the NMDS plot with points, ellipses, and reference lines.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import rlang
#' @export
#'
#' @examples
#' # Load example data
#' library(vegan)
#' data(dune, package = "vegan")
#'
#' # Run NMDS
#' nmds_result <- metaMDS(dune, k = 2)
#' nmds_scores <- scores(nmds_result, display = "sites")
#'
#' # Adding factor to data
#' dune$Management <- sample(c("A", "B"), nrow(dune), replace = TRUE)
#' nmds_data <- cbind(nmds_scores, dune)
#'
#' # Create NMDS plot
#' brc_gg_ordi(nmds_data, .color = Management, ordi = "NMDS", .drop_na = Management)
# Standard Ordination plot
brc_gg_ordi <- function(
    .data,
    .color = NULL,
    .shape = NULL,
    color_values = NULL,
    .labels = NULL,
    .size = 1,
    .alpha = 1,
    ordi,
    .drop_na = NULL
) {
    # Input validation
    if (!inherits(.data, "data.frame")) {
        cli::cli_abort("{.arg .data} must be a data frame")
    }

    # Capture quosures safely
    .color <- rlang::enquo(.color)
    .shape <- rlang::enquo(.shape)
    .drop_na <- rlang::enquo(.drop_na)

    ORDIS <- c("NMDS", "PCoA")

    ordi <- match.arg(ordi, ORDIS)

    # Base ggplot elements
    base_plot <- .data %>%
        {
            if (!rlang::quo_is_null(.drop_na)) {
                tidyr::drop_na(., !!.drop_na)
            } else {
                .
            }
        } %>%
        ggplot() +
        geom_hline(
            yintercept = 0,
            colour = "grey70",
            linewidth = 0.65
        ) +
        geom_vline(
            xintercept = 0,
            colour = "grey70",
            linewidth = 0.65
        ) +
        theme_bw(base_size = 12) +
        theme(
            legend.position = "right",
            legend.title = element_text()
        )

    # Add coordinates based on ordination type
    if (ordi == "NMDS") {
        base_plot <- base_plot +
            aes(x = NMDS1, y = NMDS2)
    } else if (ordi == "PCoA") {
        base_plot <- base_plot +
            aes(x = Dim1, y = Dim2)
    }

    # Build aesthetics conditionally
    point_aes <- aes()
    ellipse_aes <- aes()

    if (!rlang::quo_is_null(.color)) {
        point_aes$colour <- rlang::quo_get_expr(.color)
        ellipse_aes$colour <- rlang::quo_get_expr(.color)
    }

    if (!rlang::quo_is_null(.shape)) {
        point_aes$shape <- rlang::quo_get_expr(.shape)
    }

    # Add points
    final_plot <- base_plot +
        geom_point(
            mapping = point_aes,
            stroke = 1,
            alpha = .alpha,
            size = .size,
            na.rm = TRUE
        )

    # Add ellipses only if color is provided
    if (!rlang::quo_is_null(.color)) {
        final_plot <- final_plot +
            stat_ellipse(
                mapping = ellipse_aes,
                geom = "path",
                linewidth = 1.3,
                type = "t",
                level = 0.95,
                show.legend = TRUE
            )
    }

    # Add guide titles
    if (!rlang::quo_is_null(.color) && rlang::quo_is_null(.shape)) {
        # Only color, no shape
        .color2 <- rlang::as_label(.color)
        guide_params <- list(title = str_to_title(.color2))

        if (!is.null(.labels)) {
            guide_params$labels <- .labels
        }

        final_plot <- final_plot +
            guides(color = do.call(guide_legend, guide_params))

        if (!is.null(color_values)) {
            final_plot <- final_plot +
                scale_color_manual(values = color_values, labels = .labels)
        }
    }

    if (!rlang::quo_is_null(.shape) && !rlang::quo_is_null(.color)) {
        # Both color and shape - combine into one legend
        .shape2 <- rlang::as_label(.shape)
        .color2 <- rlang::as_label(.color)

        guide_params <- list(
            title = str_to_upper(.color2),
            override.aes = list(
                shape = c(16, 17), # Will use first n shapes based on data
                alpha = 1
            )
        )

        if (!is.null(.labels)) {
            guide_params$labels <- .labels
        }

        final_plot <- final_plot +
            guides(
                color = do.call(guide_legend, guide_params),
                shape = "none" # Hide separate shape legend
            )

        if (!is.null(color_values)) {
            final_plot <- final_plot +
                scale_color_manual(values = color_values, labels = .labels)
        }
    }

    if (!rlang::quo_is_null(.shape) && rlang::quo_is_null(.color)) {
        # Only shape, no color
        .shape2 <- rlang::as_label(.shape)
        guide_params <- list(title = str_to_upper(.shape2))

        if (!is.null(.labels)) {
            guide_params$labels <- .labels
        }

        final_plot <- final_plot +
            guides(shape = do.call(guide_legend, guide_params))
    }
    return(final_plot)
}
