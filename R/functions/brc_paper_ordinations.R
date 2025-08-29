brc_paper_ordinations <- function(
  nmds_data,
  pcoa_data,
  nmds_colors = NULL,
  pcoa_colors = NULL,
  nmds_labels = NULL,
  pcoa_labels = NULL,
  point_size = 2.5,
  point_alpha = 0.5,
  crop_theme = NULL,
  brc_theme = NULL
) {
  # Prepare data and ordination types
  data <- list(nmds_data, pcoa_data)
  ordi_type <- list("NMDS", "PCoA")

  # Create plots using map2
  plots <- purrr::map2(
    ordi_type,
    data,
    function(ord_type, plot_data) {
      list(
        crops = brc_gg_ordi(
          .data = plot_data,
          ordi = ord_type,
          .color = crop,
          .size = point_size,
          .alpha = point_alpha,
          .drop_na = brc
        ) +
          crop_theme,

        brc = brc_gg_ordi(
          .data = plot_data,
          ordi = ord_type,
          .color = brc,
          color_values = if (ord_type == "NMDS") nmds_colors else pcoa_colors,
          .labels = if (ord_type == "NMDS") nmds_labels else pcoa_labels,
          .size = point_size,
          .alpha = point_alpha,
          .drop_na = brc
        ) +
          brc_theme +
          guides(color = guide_legend(title = "BRC"))
      )
    }
  ) %>%
    purrr::set_names(c("nmds", "pcoa"))

  return(plots)
}
