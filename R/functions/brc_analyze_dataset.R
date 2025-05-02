# Function to generate and display all analyses for a dataset
analyze_dataset <- function(phyloseq_obj, label) {
  cat("\n### ", label, "\n")

  # Generate read pattern analysis
  patterns <- analyze_read_patterns(phyloseq_obj) +
    ggtitle(label)
  print(patterns)

  # Generate non-zero stats
  nonzero_results <- plot_nonzero_stats(phyloseq_obj)
  print(nonzero_results$plot)

  # Generate presence threshold analysis
  presence_results <- plot_asv_presence_threshold(
    phyloseq_obj,
    thresholds = thresholds
  )
  print(presence_results$combined_plot)
}
