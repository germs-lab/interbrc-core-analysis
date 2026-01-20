#' Unified Ordination Interface for ASV Communities
#'
#' This function provides a unified interface for multiple ordination methods
#' including NMDS, PCoA, dbRDA, and capscale. It standardizes parameter names
#' and return structures across different ordination types.
#'
#' @param data ASV count matrix (samples x ASVs), phyloseq object, or distance matrix
#' @param physeq Phyloseq object containing sample metadata (required unless data is phyloseq)
#' @param method Character string specifying ordination method. One of:
#'   \itemize{
#'     \item "NMDS" - Non-metric Multidimensional Scaling
#'     \item "PCoA" - Principal Coordinates Analysis
#'     \item "dbRDA" - Distance-based Redundancy Analysis
#'     \item "capscale" - Constrained Analysis of Principal Coordinates
#'   }
#' @param distance Distance metric for calculations. Default is "bray".
#'   Ignored if data is already a distance matrix.
#' @param transform Data transformation method. Options: "hellinger", "none".
#'   Default is "hellinger" for NMDS, "none" for PCoA.
#' @param formula Formula for constrained ordinations (dbRDA, capscale).
#'   Example: ~ treatment + site
#' @param k Number of dimensions/axes. Default is 2.
#' @param ncores Number of cores for parallel processing. Default uses available cores - 1.
#' @param trymax For NMDS: Maximum number of random starts. Default is 100.
#' @param maxit For NMDS: Maximum iterations. Default is 999.
#' @param previous.best For NMDS: Previous best solution to use as starting point.
#' @param ... Additional arguments passed to underlying ordination functions
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{scores_df}: Ordination scores merged with sample metadata
#'     \item \code{ordi}: The ordination object from vegan
#'     \item \code{method}: The ordination method used
#'     \item \code{metadata}: Additional metadata about the ordination
#'   }
#'
#' @details
#' This function consolidates the functionality of \code{brc_nmds()} and
#' \code{brc_pcoa()} into a single interface. It maintains backward compatibility
#' while providing a more intuitive API for ordination analysis.
#'
#' For NMDS, the function performs Hellinger transformation (by default),
#' calculates Bray-Curtis distances, and computes NMDS ordination with
#' convergence checking.
#'
#' For PCoA, the function uses weighted classical multidimensional scaling
#' and can accept pre-computed distance matrices.
#'
#' @section Backward Compatibility:
#' The original \code{brc_nmds()} and \code{brc_pcoa()} functions are maintained
#' as wrappers around this function for backward compatibility.
#'
#' # For scripts, not packaged
#'
#' @examples
#' \dontrun{
#' # NMDS ordination
#' nmds_result <- brc_ordination(
#'   data = asv_matrix,
#'   physeq = ps_obj,
#'   method = "NMDS",
#'   k = 2
#' )
#'
#' # PCoA ordination
#' pcoa_result <- brc_ordination(
#'   data = distance_matrix,
#'   physeq = ps_obj,
#'   method = "PCoA"
#' )
#'
#' # dbRDA ordination
#' dbrda_result <- brc_ordination(
#'   data = asv_matrix,
#'   physeq = ps_obj,
#'   method = "dbRDA",
#'   formula = ~ treatment + block
#' )
#' }
brc_ordination <- function(
    data,
    physeq = NULL,
    method = c("NMDS", "PCoA", "dbRDA", "capscale"),
    distance = "bray",
    transform = NULL,
    formula = NULL,
    k = 2,
    ncores = parallel::detectCores() - 1,
    trymax = 100,
    maxit = 999,
    previous.best = NULL,
    ...) {
  
  # Match method argument
  method <- match.arg(method)
  
  # Set default transform based on method
  if (is.null(transform)) {
    transform <- ifelse(method == "NMDS", "hellinger", "none")
  }
  
  # Extract phyloseq if data is phyloseq object
  if (inherits(data, "phyloseq")) {
    physeq <- data
    data <- phyloseq::otu_table(physeq)
    if (!phyloseq::taxa_are_rows(physeq)) {
      data <- t(data)
    }
  }
  
  # Validate inputs
  if (is.null(physeq)) {
    cli::cli_abort("{.arg physeq} must be provided or {.arg data} must be a phyloseq object.")
  }
  
  if (!inherits(physeq, "phyloseq")) {
    cli::cli_abort("{.arg physeq} must be a phyloseq object.")
  }
  
  # Check for constrained methods requiring formula
  if (method %in% c("dbRDA", "capscale") && is.null(formula)) {
    cli::cli_abort("{.arg formula} is required for constrained ordinations (dbRDA, capscale).")
  }
  
  # Dispatch to appropriate method
  result <- switch(
    method,
    NMDS = brc_ordination_nmds(
      data = data,
      physeq = physeq,
      distance = distance,
      transform = transform,
      k = k,
      ncores = ncores,
      trymax = trymax,
      maxit = maxit,
      previous.best = previous.best,
      ...
    ),
    PCoA = brc_ordination_pcoa(
      data = data,
      physeq = physeq,
      distance = distance,
      transform = transform,
      k = k,
      ncores = ncores,
      ...
    ),
    dbRDA = brc_ordination_dbrda(
      data = data,
      physeq = physeq,
      formula = formula,
      distance = distance,
      transform = transform,
      k = k,
      ...
    ),
    capscale = brc_ordination_capscale(
      data = data,
      physeq = physeq,
      formula = formula,
      distance = distance,
      transform = transform,
      k = k,
      ...
    )
  )
  
  # Add method to result
  result$method <- method
  
  return(result)
}

# Internal helper functions ----

#' @keywords internal
brc_ordination_nmds <- function(
    data,
    physeq,
    distance,
    transform,
    k,
    ncores,
    trymax,
    maxit,
    previous.best,
    seed = 54641,
    ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  if (!is.matrix(data) && !is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a matrix or data frame.")
  }
  
  if (any(is.na(data)) || any(is.nan(data))) {
    cli::cli_abort("{.arg data} contains NA or NaN values.")
  }
  
  # Apply transformation
  if (transform == "hellinger") {
    transformed_data <- vegan::decostand(t(data), MARGIN = 1, method = "hellinger")
  } else {
    transformed_data <- t(data)
  }
  
  # Remove columns that sum to 0
  transformed_data <- transformed_data %>%
    as.data.frame() %>%
    dplyr::select(where(~ is.numeric(.) && sum(.) > 0)) %>%
    as.matrix()
  
  # Check for empty matrix after filtering
  if (ncol(transformed_data) == 0) {
    cli::cli_abort("No valid columns remaining after filtering.")
  }
  
  # Calculate distances
  dist_matrix <- vegan::vegdist(
    t(transformed_data),
    method = distance,
    upper = FALSE,
    binary = FALSE,
    na.rm = TRUE
  )
  
  # Perform NMDS
  ordi <- vegan::metaMDS(
    as.matrix(dist_matrix),
    distance = distance,
    display = "sites",
    noshare = TRUE,
    autotransform = FALSE,
    wascores = TRUE,
    tidy = TRUE,
    k = k,
    maxit = maxit,
    trymax = trymax,
    parallel = ncores,
    previous.best = previous.best
  )
  
  # Check NMDS convergence
  if (ordi$converged == FALSE) {
    cli::cli_alert_warning(
      "NMDS did not converge. Consider increasing {.arg trymax}."
    )
  }
  
  # Add species scores
  vegan::sppscores(ordi) <- t(transformed_data)
  
  # Extract NMDS scores
  scores_raw <- as.data.frame(vegan::scores(ordi)$sites)
  scores_raw <- scores_raw %>%
    tibble::rownames_to_column(var = "unique_id")
  
  # Merge with sample metadata
  sample_df <- physeq@sam_data %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "unique_id")
  
  scores_df <- dplyr::right_join(sample_df, scores_raw, by = "unique_id")
  
  # Return results
  return(list(
    scores_df = scores_df,
    ordi = ordi,
    metadata = list(
      stress = ordi$stress,
      converged = ordi$converged,
      distance = distance,
      transform = transform
    )
  ))
}

#' @keywords internal
brc_ordination_pcoa <- function(
    data,
    physeq,
    distance,
    transform,
    k,
    ncores,
    eig = TRUE,
    add = FALSE,
    x.ret = FALSE,
    seed = 54641,
    ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Check if data is already a distance matrix
  is_dist <- "dist" %in% class(data)
  
  if (!is_dist) {
    # Validate matrix/data frame
    if (!is.matrix(data) && !is.data.frame(data)) {
      cli::cli_abort("{.arg data} must be a matrix, data frame, or distance matrix.")
    }
    
    if (any(is.na(data)) || any(is.nan(data))) {
      cli::cli_abort("{.arg data} contains NA or NaN values.")
    }
    
    # Apply transformation if needed
    if (transform == "hellinger") {
      transformed_data <- vegan::decostand(t(data), MARGIN = 1, method = "hellinger")
      data <- vegan::vegdist(t(transformed_data), method = distance)
    } else if (transform != "none") {
      cli::cli_alert_warning("Unknown transform '{transform}'. Using 'none'.")
    }
  } else {
    # Data is already a distance matrix
    if (any(is.na(data)) || any(is.nan(data))) {
      cli::cli_abort("{.arg data} (distance matrix) contains NA or NaN values.")
    }
  }
  
  # Perform PCoA
  ordi <- vegan::wcmdscale(data, k = k, eig = eig, add = add, x.ret = x.ret)
  
  # Extract PCoA scores and add metadata
  scores_df <- data.frame(
    Dim1 = ordi$points[, 1],
    Dim2 = ordi$points[, 2],
    brc = factor(physeq@sam_data$brc),
    crop = factor(physeq@sam_data$crop)
  )
  
  # Add all sample data
  sample_df <- physeq@sam_data %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "sample_id")
  
  scores_df <- scores_df %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::left_join(sample_df, by = "sample_id")
  
  # Calculate variance explained
  variance_explained <- NULL
  if (eig) {
    variance_explained <- ordi$eig / sum(ordi$eig)
  }
  
  # Return results
  return(list(
    scores_df = scores_df,
    ordi = ordi,
    metadata = list(
      variance_explained = variance_explained,
      distance = distance,
      transform = transform
    )
  ))
}

#' @keywords internal
brc_ordination_dbrda <- function(
    data,
    physeq,
    formula,
    distance,
    transform,
    k,
    ...) {
  
  cli::cli_alert_info("dbRDA ordination is currently under development.")
  cli::cli_abort("Please use vegan::dbrda() directly or wait for implementation.")
}

#' @keywords internal
brc_ordination_capscale <- function(
    data,
    physeq,
    formula,
    distance,
    transform,
    k,
    ...) {
  
  cli::cli_alert_info("capscale ordination is currently under development.")
  cli::cli_abort("Please use vegan::capscale() directly or wait for implementation.")
}

# Backward compatibility wrappers ----

#' Calculate NMDS for ASV Communities (Legacy)
#'
#' This is a backward-compatible wrapper around \code{brc_ordination()}.
#' For new code, please use \code{brc_ordination(method = "NMDS")}.
#'
#' @inheritParams brc_ordination
#' @param asv_matrix A matrix or data frame of ASV counts (samples x ASVs).
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{nmds_scores}: NMDS scores for samples.
#'     \item \code{nmds_df}: NMDS scores merged with sample metadata.
#'     \item \code{ordi}: NMDS ordination object.
#'   }
#' # For scripts, not packaged
#'
#' @examples
#' \dontrun{
#' nmds_result <- brc_nmds(asv_matrix, physeq_obj)
#' }
brc_nmds_new <- function(
    asv_matrix,
    physeq,
    ncores = parallel::detectCores(),
    k = 2,
    maxit = 999,
    trymax = 100,
    previous.best = NULL) {
  
  result <- brc_ordination(
    data = asv_matrix,
    physeq = physeq,
    method = "NMDS",
    k = k,
    ncores = ncores,
    trymax = trymax,
    maxit = maxit,
    previous.best = previous.best
  )
  
  # Reformat to match legacy output structure
  return(list(
    nmds_scores = result$scores_df %>%
      dplyr::select(unique_id, NMDS1, NMDS2),
    nmds_df = result$scores_df,
    ordi = result$ordi
  ))
}

#' Calculate PCoA for ASV Communities (Legacy)
#'
#' This is a backward-compatible wrapper around \code{brc_ordination()}.
#' For new code, please use \code{brc_ordination(method = "PCoA")}.
#'
#' @inheritParams brc_ordination
#' @param asv_matrix A matrix, data frame, or distance object
#' @param eig Logical, return eigenvalues
#' @param add Logical, add constant to make matrix Euclidean
#' @param x.ret Logical, return distance matrix
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{pcoa_df}: PCoA scores merged with sample metadata.
#'     \item \code{ordi}: PCoA ordination object.
#'     \item \code{distance_matrix}: Distance matrix used.
#'   }
#' # For scripts, not packaged
#'
#' @examples
#' \dontrun{
#' pcoa_result <- brc_pcoa(distance_matrix, physeq_obj)
#' }
brc_pcoa_new <- function(
    asv_matrix,
    physeq,
    ncores = parallel::detectCores() - 1,
    k = 2,
    eig = TRUE,
    add = FALSE,
    x.ret = FALSE) {
  
  result <- brc_ordination(
    data = asv_matrix,
    physeq = physeq,
    method = "PCoA",
    k = k,
    ncores = ncores,
    eig = eig,
    add = add,
    x.ret = x.ret
  )
  
  # Reformat to match legacy output structure
  return(list(
    pcoa_df = result$scores_df,
    ordi = result$ordi,
    distance_matrix = asv_matrix
  ))
}
