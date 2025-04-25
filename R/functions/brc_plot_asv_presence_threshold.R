plot_asv_presence_threshold <- function(
  phyloseq_obj,
  thresholds = c(60, 70, 80, 90)
) {
  # Extract ASV table
  asv_mat <- as.matrix(otu_table(phyloseq_obj))
  n_samples <- ncol(asv_mat)

  # Create empty list to store plots
  plots <- list()

  # Calculate statistics for each threshold
  for (threshold in thresholds) {
    # Calculate the minimum number of samples needed for this threshold
    min_samples_needed <- ceiling((threshold / 100) * n_samples)

    # For each ASV, calculate in how many samples it has non-zero reads
    presence_counts <- rowSums(asv_mat > 0)

    # Calculate ASVs present in at least threshold% of samples
    presence_stats <- data.frame(
      Sample = colnames(asv_mat),
      total_asvs = nrow(asv_mat),
      asvs_above_threshold = sum(presence_counts >= min_samples_needed)
    ) %>%
      mutate(
        percent_asvs = (asvs_above_threshold / total_asvs) * 100
      )

    # Create plot for this threshold
    p <- ggplot(presence_stats, aes(x = Sample, y = percent_asvs)) +
      geom_bar(stat = "identity", fill = "darkgreen") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(
        title = sprintf("ASVs Present in ≥%d%% of Samples", threshold),
        x = "Sample",
        y = sprintf("Perc. of ASVs present in ≥%d%% of samples", threshold)
      )

    plots[[sprintf("threshold_%d", threshold)]] <- p
  }

  # Combine all plots using cowplot
  combined_plot <- plot_grid(
    plotlist = plots,
    ncol = 2
    #labels = sprintf("≥%d%% presence", thresholds)
  )

  # Create summary statistics for each threshold
  summary_stats <- lapply(thresholds, function(threshold) {
    min_samples_needed <- ceiling((threshold / 100) * n_samples)
    presence_counts <- rowSums(asv_mat > 0)
    asvs_above_threshold <- sum(presence_counts >= min_samples_needed)

    data.frame(
      Threshold = threshold,
      Total_ASVs = nrow(asv_mat),
      ASVs_Meeting_Threshold = asvs_above_threshold,
      Percentage = round((asvs_above_threshold / nrow(asv_mat)) * 100, 2)
    )
  }) %>%
    bind_rows()

  return(list(
    #individual_plots = plots,
    combined_plot = combined_plot,
    summary_statistics = summary_stats
  ))
}
