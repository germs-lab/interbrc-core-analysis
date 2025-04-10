# Miscellaneous functions for Inter-BRC Core Microbiome Analysis
# By Bolívar Aponte Rolón

# extract_matrix(): extract/subset a ASV/OTU matrix based on a vector of strings and a rowSums() condition.

extract_matrix <- function(physeq, .vec, keep_rows_sums = 0) {
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



# Borrowed subset.fasta from https://github.com/GuillemSalazar/FastaUtils/blob/master/R/FastaUtils_functions.R

#' Select a subset of sequences from a fasta file
#'
#' This function loads a fasta file, selects a subset on sequences based on the sequence's names and writes a new fsata file.
#' @param file Fasta file
#' @param subset Vector with (exact) names of the sequences to be retrieved.
#' @param out Path to output file. If absent, the '.mod' is added to the input file's name (path is thus conserved).
#' @keywords FastaUtils
#' @return Writes a fasta file with the selected sequences.
#' @export
#' @author Guillem Salazar <salazar@@icm.csic.es>
#' @examples
#' subset.fasta(file = "http://greengenes.lbl.gov/Data/JD_Tutorial/UnAligSeq24606.txt", subset = c("24.6jsd1.Tut", "24.6jsd2.Tut ", "24.6jsd3.Tut "), out = "out.fasta")
subset_fasta <- function(file = NULL, subset = NULL, out = paste(file, ".subset", sep = "")) {
  library(Biostrings)
  sequences <- readDNAStringSet(file)
  if (all(as.character(subset) %in% names(sequences)) == FALSE) stop("There are names in 'subset' not present in the fasta file")
  pos <- match(as.character(subset), names(sequences))
  writeXStringSet(sequences[pos], filepath = out)
}


# Parallelized function for choosing dimensions for NMDS
nmds_screen_parallel <- function(x, ncores = parallel::detectCores() - 1) {
  # Function to calculate stress for a given number of dimensions
  calculate_stress <- function(k) {
    replicate(10, metaMDS(x, autotransform = FALSE, k = k, maxit = 100, trymax = 10)$stress)
  }

  # Use mclapply to parallelize the stress calculation
  stress_values <- parallel::mclapply(1:10, calculate_stress, mc.cores = ncores)

  # Plot the results
  plot(rep(1, 10), stress_values[[1]],
    xlim = c(1, 10),
    ylim = c(0, 0.30),
    xlab = "# of Dimensions",
    ylab = "Stress",
    main = "NMDS stress plot"
  )

  for (i in 1:9) {
    points(rep(i + 1, 10), stress_values[[i + 1]])
  }
}

# Standard Ordination plot
gg_ordi <- function(.data, .color, ordi, .drop_na = NULL) {
  # Input validation
  if (!inherits(.data, "data.frame")) {
    cli::cli_abort("{.arg .data} must be a data frame")
  }

  # Capture quosures safely
  .color <- rlang::enquo(.color)
  .drop_na <- rlang::enquo(.drop_na)

  ORDIS <- c("NMDS", "PCoA")

  ordi <- match.arg(ordi, ORDIS)

  # Base ggplot elements
  base_plot <- .data %>%
    {
      if (!rlang::quo_is_null(.drop_na)) {
        tidyr::drop_na(., !!.drop_na)
      } else {
        .
      }
    } %>%
    ggplot() +
    geom_hline(
      yintercept = 0,
      colour = "grey70",
      linewidth = 0.65
    ) +
    geom_vline(
      xintercept = 0,
      colour = "grey70",
      linewidth = 0.65
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text()
    )

  # Add coordinates based on ordination type
  if (ordi == "NMDS") {
    base_plot <- base_plot +
      aes(x = NMDS1, y = NMDS2)
  } else if (ordi == "PCoA") {
    base_plot <- base_plot +
      aes(x = Dim1, y = Dim2)
  }

  # Add points and ellipses
  final_plot <- base_plot +
    geom_point(
      aes(color = !!.color),
      stroke = 1,
      alpha = 0.5,
      na.rm = TRUE
    ) +
    stat_ellipse(
      aes(color = !!.color),
      geom = "path",
      linewidth = 1.3,
      type = "t",
      level = 0.95,
      show.legend = TRUE
    )

  return(final_plot)
}


#' Subset phyloseq object based on core analysis results
#'
#' Filters a phyloseq object to retain only samples and taxa identified in a core microbiome analysis
#' using `extract_core()`.
#' Returns both the subsetted phyloseq object and corresponding ASV matrix.
#'
#' @param core_obj A list object containing core microbiome analysis results.
#'        Expected to contain a dataframe at position [[4]] with OTU information.
#' @param physeq A phyloseq object to be subsetted.
#' @param .var Character string specifying the column name in `core_obj[[4]]`
#'        containing sample/group identifiers to filter by.
#' @param type Value to filter by in the 'fill' column of `core_obj[[4]]`.
#'        Determines which OTUs are retained.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \code{expected_physeq} - Subsetted phyloseq object
#'   \item \code{asv_matrix} - Corresponding ASV abundance matrix
#' }
#'
#' @details
#' This function performs the following operations:
#' \enumerate{
#'   \item Extracts OTUs of interest from the core analysis object based on the specified type
#'   \item Filters the phyloseq object to retain only these OTUs and samples where they appear
#'   \item Verifies the filtering was successful
#'   \item Returns both the subsetted phyloseq object and ASV matrix
#' }
#'
#' The function uses \code{prune_samples} and \code{prune_taxa} to safely subset the phyloseq object
#' while preserving all associated data (sample data, taxonomy table, etc.). The ASV matrix is
#' processed to remove samples with zero row sums.
#'
#' @examples
#' \dontrun{
#' # Subset for "persistent" core OTUs
#' result <- subset_physeq(
#'   core_obj = my_core_results,
#'   physeq = my_phyloseq,
#'   .var = "body_site",
#'   type = "persistent"
#' )
#'
#' # Access results
#' subset_physeq <- result$expected_physeq
#' asv_matrix <- result$asv_matrix
#' }
#'
#' @import phyloseq
#' @import dplyr
#' @import tibble
#' @import cli
#' @export
subset_physeq <- function(core_obj, physeq, .var, type) {
  asv_matrix <-
    core_obj[[4]] %>% # Get strings from all OTUs
    dplyr::filter(., .$fill == type) %>%
    tibble::column_to_rownames(., var = .var) %>%
    rownames() %>%
    extract_matrix(physeq, .vec = .) # Remove samples where row sum 0 through extract_matrix()

  sample_strings <- rownames(asv_matrix)

  ## Phyloseqs
  expected_physeq <- phyloseq::prune_samples(
    sort(phyloseq::sample_names(physeq)) %in% sort(sample_strings),
    physeq
  ) %>%
    phyloseq::prune_taxa(rownames(.@otu_table) %in% colnames(asv_matrix), .)

  ## Checking that core is filtered out
  if (any(rownames(expected_physeq@otu_table) %in% colnames(asv_matrix))) {
    cli::cli_alert_info("Successfully selected ASVs of interest")
  }

  return(
    list(
      expected_physeq = expected_physeq,
      asv_matrix = asv_matrix
    )
  )
}



brc_ggsave <- function(objects, save_path) {
  ggnames <- c(names(objects))


  plot_paths <- str_glue("{save_path}{tolower(ggnames)}.png")

  purrr::walk2(
    plot_paths, objects,
    \(path, plot) ggsave(
      filename = path,
      plot = plot,
      dpi = 300,
      width = 200,
      height = 200,
      units = "mm"
    )
  )
}
