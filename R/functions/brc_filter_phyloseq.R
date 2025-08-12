filter_phyloseq <- function(
  phyloseq_obj,
  min_sample_reads = 100,
  min_asv_reads = 20
) {
  # Input validation
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("phyloseq package is required for this function")
  }

  if (!inherits(phyloseq_obj, "phyloseq")) {
    stop("Input must be a phyloseq object")
  }

  # Remove samples with low number of reads
  filtered_obj <- phyloseq::prune_samples(
    phyloseq::sample_sums(phyloseq_obj) >= min_sample_reads,
    phyloseq_obj
  )

  # Filter ASVs based on minimum read count
  filtered_obj <- phyloseq::filter_taxa(
    filtered_obj,
    function(x) {
      sum(x > min_asv_reads) > (0.00 * length(x))
    },
    TRUE
  )

  return(filtered_obj)
}

filter_zero_asvs <- function(phyloseq_obj, max_zero_percent = 90) {
  otu_mat <- as.matrix(otu_table(phyloseq_obj))
  
  # Calculate percentage of zeros per ASV (row-wise)
  zero_percentages <- apply(otu_mat, 1, function(x) sum(x == 0) / length(x) * 100)
  
  # Identify ASVs to keep
  asvs_to_keep <- names(zero_percentages[zero_percentages <= max_zero_percent])
  
  # Filter phyloseq object
  filtered_obj <- prune_taxa(asvs_to_keep, phyloseq_obj)
  
  cat("Zero threshold:", max_zero_percent, "% - ASVs remaining:", 
      ntaxa(filtered_obj), "out of", ntaxa(phyloseq_obj), "\n")
  
  return(list(
    phyloseq = filtered_obj,
    zero_percentages = zero_percentages,
    asvs_removed = setdiff(taxa_names(phyloseq_obj), asvs_to_keep)
  ))
}