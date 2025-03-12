#' Select High and Low Occurrence ASVs from a Phyloseq Object
#'
#' This function selects a specified number of high and low occurrence ASVs from a phyloseq object,
#' ensuring no empty samples remain. It returns pruned phyloseq objects for both high and low occurrence ASVs.
#'
#' @param physeq A phyloseq object.
#' @param high_occurrence_asvs A character vector of high occurrence ASV names.
#' @param low_occurrence_asvs A character vector of low occurrence ASV names.
#' @param as The matrix dimension that is desired for taxa. Must be
#'   "rows" for rows and "columns" or "cols" for columns.
#' @return A named list containing:
#'   \itemize{
#'     \item \code{ps_high_occ}: Phyloseq object with selected high occurrence ASVs.
#'     \item \code{ps_low_occ}: Phyloseq object with selected low occurrence ASVs.
#'   }
#' @export
select_asvs <- function(physeq, high_occurrence_asvs, low_occurrence_asvs, as = "rows") {
  # Validate inputs
  if (!inherits(physeq, "phyloseq")) {
    cli::cli_abort("{.arg physeq} must be a phyloseq object. Got {class(physeq)[1]} instead.")
  }

  if (length(high_occurrence_asvs) == 0 || length(low_occurrence_asvs) == 0) {
    cli::cli_abort("Both {.arg high_occurrence_asvs} and {.arg low_occurrence_asvs} must be non-empty.")
  }

  # Ensure num_asvs is not greater than available ASVs
  num_asvs <- min(length(high_occurrence_asvs), length(low_occurrence_asvs))

  cli::cli_alert_info("Selecting {num_asvs} ASVs from each group.")

  # Helper function to select ASVs while avoiding empty samples
  select_valid_asvs <- function(asv_list, physeq, num_asvs, as = "rows") {
    stopifnot(as %in% c("rows", "columns", "cols")) # Borrowed from https://github.com/mikemc/speedyseq/blob/0057652ff7a4244ccef2b786dca58d901ec2fc62/R/utils.R#L114
    if (identical(as, "rows")) {
      if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
      }
    } else {
      if (taxa_are_rows(physeq)) {
        physeq <- t(physeq)
      }
    }

    return(physeq)
  }
  asv_table <- as.data.frame(otu_table(physeq))

  select_valid_asvs <- function(asv_table, high_occurrence_asvs, low_occurrence_asvs, num_asvs = 12) {
    asv_table <- as.data.frame(asv_table)
    asv_list <- list(
      high_asvs = high_occurrence_asvs,
      low_asvs = low_occurrence_asvs
    )

    # Validate inputs
    if (!is.list(asv_list)) {
      cli::cli_abort("{.arg asv_list} must be a named list of ASV groups.")
    }

    if (!is.data.frame(asv_table)) {
      cli::cli_abort("{.arg asv_table} must be a data frame or matrix.")
    }

    if (num_asvs < 1) {
      cli::cli_abort("{.arg num_asvs} must be at least 1.")
    }

    # Initialize output list
    selected_asvs <- list()

    # Process each ASV group
    for (group_name in names(asv_list)) {
      asvs <- asv_list[[group_name]]

      # Ensure ASV names match column names
      valid_asvs <- intersect(asvs, colnames(asv_table))
      if (length(valid_asvs) == 0) {
        cli::cli_alert_warning("No matching ASVs found for group '{group_name}'.")
        selected_asvs[[group_name]] <- character(0)
        next
      }

      # Sort ASVs by occurrence frequency
      sorted_asvs <- valid_asvs[order(-colSums(asv_table[, valid_asvs] > 0))]

      # Select top 'num_asvs' ASVs
      selected_asvs[[group_name]] <- sorted_asvs[1:min(num_asvs, length(sorted_asvs))]

      if (length(selected_asvs[[group_name]]) < num_asvs) {
        cli::cli_alert_warning(
          "Only {length(selected_asvs[[group_name]])} ASVs available for group '{group_name}'."
        )
      }
    }

    return(selected_asvs)
  }

  # Select ASVs
  selected_high_occ_asvs <- select_valid_asvs(high_occurrence_asvs, physeq, num_asvs)
  selected_low_occ_asvs <- select_valid_asvs(low_occurrence_asvs, physeq, num_asvs)

  # Prune the phyloseq object to include only selected ASVs
  ps_high_occ <- prune_taxa(selected_high_occ_asvs, physeq)
  ps_low_occ <- prune_taxa(selected_low_occ_asvs, physeq)

  # Check for empty samples
  if (any(sample_sums(ps_high_occ) == 0)) {
    cli::cli_alert_warning("Some samples in the high occurrence phyloseq object have zero reads.")
  }
  if (any(sample_sums(ps_low_occ) == 0)) {
    cli::cli_alert_warning("Some samples in the low occurrence phyloseq object have zero reads.")
  }

  # Return results
  return(list(
    ps_high_occ = ps_high_occ,
    ps_low_occ = ps_low_occ
  ))
}
