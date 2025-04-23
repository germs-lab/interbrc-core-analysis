brc_bc_occ_curve <- function(
  core_summary_list,
  point_size = 2.5,
  point_alpha = 1,
  point_shape = 21,
  x_lab = "log10(Mean Relative Abundance)",
  y_lab = "Occupancy",
  theme_settings = NULL
) {
  # Extract occupancy-abundance data from the list
  occ_abun_data <- core_summary_list[[4]]

  # Create base plot
  p <- ggplot(
    occ_abun_data,
    aes(
      x = log10(otu_rel),
      y = otu_occ,
      fill = fill
    )
  ) +
    geom_point(
      pch = point_shape,
      alpha = point_alpha,
      size = point_size
    ) +
    labs(x = x_lab, y = y_lab) +
    theme_bw()

  # Apply custom theme settings if provided, otherwise use defaults
  if (is.null(theme_settings)) {
    p <- p +
      theme(
        axis.title.x = element_text(size = 12),
        title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.margin = unit(c(.5, 1, .5, .5), "cm")
      )
  } else {
    p <- p + theme_settings
  }

  return(p)
}
