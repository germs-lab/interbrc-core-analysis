abunocc_core <- function(physeq, threshold = 0.6, as = "rows") {
  threshold_asv <- BRCore::filter_core(physeq, threshold = threshold, as = as)

  otu <- otu_table(physeq) %>% as("matrix")

  # calculating occupancy and abundance
  otu_PA <-
    1 * ((otu > 0) == 1) # presence-absence data
  otu_occ <-
    rowSums(otu_PA) / ncol(otu_PA) # occupancy calculation
  otu_rel <-
    apply(decostand(otu, method = "total", MARGIN = 2), 1, mean) # mean relative abundance

  occ_abun <- tibble::rownames_to_column(
    as.data.frame(cbind(otu_occ, otu_rel)),
    "otu"
  ) %>%
    tibble::as_tibble()

  # Get ASV IDs that passed the threshold filter
  high_occ_asvs <- rownames(threshold_asv$physeq_high_occ@otu_table)

  # Add fill column based on threshold
  occ_abun$fill <- ifelse(occ_abun$otu %in% high_occ_asvs, "core", "no")

  return(occ_abun)
}
