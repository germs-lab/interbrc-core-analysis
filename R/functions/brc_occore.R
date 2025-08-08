brc_occore <- function(physeq, threshold = 0.6, as = "rows") {
  otu <- otu_table(physeq) %>% as("matrix")

  # Helper function to calculate occupancy and abundance
  occ_calc <- function() {
    otu_PA <- 1 * ((otu > 0) == 1) # presence-absence data
    otu_occ <- rowSums(otu_PA) / ncol(otu_PA) # occupancy calculation
    otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean) # mean relative abundance

    occ_abun <- tibble::rownames_to_column(
      as.data.frame(cbind(otu_occ, otu_rel)),
      "otu"
    ) %>%
      tibble::as_tibble()

    return(occ_abun)
  }

  # Helper function to handle empty OTU tables
  handle_empty_otu <- function() {
    warning("Empty OTU table detected. Returning empty result.")
    return(tibble::tibble(
      otu = character(0),
      otu_occ = numeric(0),
      otu_rel = numeric(0),
      fill = character(0)
    ))
  }

  # Check for zero dimensions first
  if (nrow(otu) == 0 || ncol(otu) == 0) {
    return(handle_empty_otu())
  }

  # Handle edge case: threshold = 0 (all samples)
  if (threshold == 0) {
    occ_abun <- occ_calc()
    occ_abun$fill <- "no"
    message("Threshold = 0: All ASVs included but none considered 'core'")
    return(occ_abun)
  }

  # Normal cases (threshold > 0)
  threshold_asv <- BRCore::filter_core(physeq, threshold = threshold, as = as)
  occ_abun <- occ_calc()
  
  # Get ASV IDs that passed the threshold filter
  high_occ_asvs <- rownames(threshold_asv$physeq_high_occ@otu_table)

  # Check if all ASVs passed the threshold (edge case)
  if (length(high_occ_asvs) == nrow(otu)) {
    occ_abun$fill <- "no"
    message("All ASVs passed the threshold. None considered 'core'.")
  } else {
    occ_abun$fill <- ifelse(occ_abun$otu %in% high_occ_asvs, "core", "no")
  }

  return(occ_abun)
}