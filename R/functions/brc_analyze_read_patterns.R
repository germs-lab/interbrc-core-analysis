library(phyloseq)
library(ggplot2)
library(tidyverse)
library(cowplot)

analyze_read_patterns <- function(phyloseq_obj) {
  # Extract ASV table and convert to long format
  asv_mat <- as.matrix(otu_table(phyloseq_obj))
  sample_data_df <- as.data.frame(sample_data(phyloseq_obj))

  # 1. Read count distribution plot - FIXED data transformation
  read_dist <- data.frame(
    ASV = rep(rownames(asv_mat), ncol(asv_mat)),
    Sample = rep(colnames(asv_mat), each = nrow(asv_mat)),
    reads = as.vector(asv_mat) # Changed order to ensure 'reads' exists
  ) %>%
    filter(reads > 0) # Now 'reads' exists when filtering

  p1 <- ggplot(read_dist, aes(x = reads)) +
    geom_histogram(binwidth = 1) +
    scale_x_log10() +
    theme_bw() +
    labs(
      title = "Distribution of Read Counts",
      x = "Number of Reads (log10)",
      y = "Frequency"
    ) +
    theme(plot.title = element_text(hjust = 0.5))

  # 2. Scatter plot to check for diagonal patterns
  # Get pairs of consecutive samples
  read_pairs <- data.frame(
    sample1 = as.vector(asv_mat[, 1:(ncol(asv_mat) - 1)]),
    sample2 = as.vector(asv_mat[, 2:ncol(asv_mat)])
  ) %>%
    filter(sample1 > 0 | sample2 > 0)

  p2 <- ggplot(read_pairs, aes(x = sample1, y = sample2)) +
    geom_point(alpha = 0.1) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    labs(
      title = "Read Count Correlation Between Consecutive Samples",
      x = "Sample N Reads (log10)",
      y = "Sample N+1 Reads (log10)"
    ) +
    theme(plot.title = element_text(hjust = 0.5))

  # 3. Pattern analysis by sample characteristics
  # Assuming 'run_id' is in sample_data
  if ("run_id" %in% colnames(sample_data_df)) {
    reads_by_run <- merge(
      read_dist,
      sample_data_df,
      by.x = "Sample",
      by.y = "row.names"
    )

    p3 <- ggplot(reads_by_run, aes(x = reads, fill = run_id)) +
      geom_histogram(position = "dodge", bins = 30) +
      scale_x_log10() +
      theme_bw() +
      labs(
        title = "Read Distribution by Sequencing Run",
        x = "Number of Reads (log10)",
        y = "Frequency",
        fill = "Run ID"
      ) +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    p3 <- NULL
  }

  # 4. Heatmap of low read counts
  low_reads_mat <- asv_mat
  low_reads_mat[low_reads_mat > 100] <- NA
  low_reads_df <- reshape2::melt(low_reads_mat, varnames = c("ASV", "Sample"))
  colnames(low_reads_df)[3] <- "Reads" # Explicitly name the value column

  p4 <- ggplot(
    low_reads_df %>% filter(!is.na(Reads)),
    aes(x = Sample, y = ASV, fill = Reads)
  ) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank()
    ) +
    labs(title = "Heatmap of Low Read Counts (<100)", x = "Sample", y = "ASV") +
    theme(plot.title = element_text(hjust = 0.5))

  # Combine plots
  if (!is.null(p3)) {
    plot_grid(p1, p2, p3, p4, ncol = 2)
  } else {
    plot_grid(p1, p2, p4, ncol = 2)
  }
}
