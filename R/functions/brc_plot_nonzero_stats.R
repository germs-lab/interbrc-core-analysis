plot_nonzero_stats <- function(phyloseq_obj) {
  # Extract ASV table
  asv_mat <- as.matrix(otu_table(phyloseq_obj))

  # Calculate ASVs with NO zeros (i.e., present in all samples)
  nonzero_stats <- data.frame(
    Sample = colnames(asv_mat),
    total_asvs = nrow(asv_mat),
    asvs_always_present = sum(rowSums(asv_mat == 0) == 0) # Count ASVs with no zeros
  ) %>%
    mutate(
      percent_always_present = (asvs_always_present / total_asvs) * 100
    )

  # Create visualization
  p <- ggplot(nonzero_stats, aes(x = Sample, y = percent_always_present)) +
    geom_bar(stat = "identity", fill = "darkgreen") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "ASVs Present in ALL Samples",
      x = "Sample",
      y = "Percentage of ASVs Present in All Samples (%)"
    )

  # Add summary statistics
  summary_stats <- data.frame(
    Metric = c(
      "Total ASVs",
      "ASVs present in all samples",
      "Percentage of ASVs present in all samples"
    ),
    Value = c(
      unique(nonzero_stats$total_asvs),
      unique(nonzero_stats$asvs_always_present),
      round(unique(nonzero_stats$percent_always_present), 2)
    )
  )

  return(list(
    plot = p,
    statistics = nonzero_stats,
    summary = summary_stats
  ))
}
