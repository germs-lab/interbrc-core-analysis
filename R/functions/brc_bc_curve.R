brc_bc_curve <- function(
  core_summary_list,
  max_otus = 100,
  threshold = 1.02
) {
  # Extract Bray-Curtis ranked data from the list
  BC_ranked <- core_summary_list[[2]] %>%
    dplyr::mutate(rank = factor(rank, levels = rank)) %>%
    tidyr::drop_na()

  # Subset to max OTUs
  BC_ranked_max <- BC_ranked[1:max_otus, ]

  # Find last OTU meeting threshold
  last_otu <- last(as.numeric(BC_ranked$rank[
    (BC_ranked$IncreaseBC >= threshold)
  ]))

  # Create plot
  p <- BC_ranked_max %>%
    ggplot(aes(x = rank[1:max_otus], y = proportionBC)) +
    geom_point(
      pch = 21,
      col = "black",
      alpha = 0.85,
      size = 3.5
    ) +
    geom_vline(
      xintercept = last_otu,
      lty = 4,
      col = "darkred",
      cex = 0.5
    ) +
    scale_x_discrete(
      limits = BC_ranked$rank[1:max_otus],
      breaks = BC_ranked$rank[1:max_otus][seq(
        1,
        length(BC_ranked$rank[1:max_otus]),
        by = 10
      )]
    ) +
    xlab("Ranked OTUs") +
    ylab("Contribution to Bray-Curtis Similarity") +
    annotate(
      geom = "text",
      x = last_otu,
      y = 0.18,
      label = paste(
        "Last",
        (threshold - 1) * 100,
        "% increase\n(",
        last_otu,
        " OTUs)",
        sep = ""
      ),
      color = "darkred",
      size = 4
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 14),
      title = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      strip.text.x = element_text(size = 10),
      strip.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      plot.margin = unit(c(.5, 1, .5, .5), "cm")
    )

  return(p)
}
