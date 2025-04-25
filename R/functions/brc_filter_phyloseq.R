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
