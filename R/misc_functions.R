# Miscellaneous functions for Inter-BRC Core Microbiome Analysis
# By Bolívar Aponte Rolón

# ExtracMatrix(): extract/subset a ASV/OTU matrix based on a vector of strings and a rowSums() condition.

ExtractMatrix <- function(physeq, .vec, keep_rows_sums = 0) {
  # Error handling
  if (!inherits(physeq, c("phyloseq", "matrix"))) {
    cli::cli_abort(
      "{.arg physeq} must be a 'phyloseq' or 'matrix' object.\nYou've supplied a {class(physeq)[1]} class."
    )
  }

  if (inherits(physeq, "phyloseq")) {
    cli::cli_alert_info("Detected a 'phyloseq' object. Input object is valid!")
    physeq <- otu_table(physeq)
  } else {
    cli::cli_alert_info("Detected a 'matrix' object. Proceeding with the input matrix.")
  }

  if (!inherits(physeq, "matrix")) {
    cli::cli_abort(
      "The extracted OTU/ASV table is not a matrix. Something went wrong."
    )
  }

  cli::cli_alert_success("Extracting the ASV/OTU matrix...")

  # Main function
  physeq %>%
    t() %>% # We need samples as rows and ASV as columns
    as.data.frame() %>%
    select(contains(.vec)) %>%
    .[rowSums(.) > keep_rows_sums, ] %>% # Keep only samples with a non-zero sum. Not all samples have the "core".
    as.matrix()

  # cli::cli_alert_success("Yay! You got a matrix.")
}
