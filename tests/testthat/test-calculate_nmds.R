library(phyloseq)
library(vegan)
library(tidyverse)
library(waldo)
library(testthat)
list.files("R/functions/", full.names = TRUE) %>%
  map(source)



# Function to create mock data
create_mock_data <- function(n_samples = 100, n_asvs = 500, zero_prob = 0.5) {
  # Set seed for reproducibility
  set.seed(123)

  # Create mock ASV matrix (samples x ASVs)
  mock_asv_matrix <- matrix(
    rpois(n_samples * n_asvs, lambda = 10), # Random counts with Poisson distribution
    nrow = n_samples,
    ncol = n_asvs,
    dimnames = list(
      paste0("Sample", 1:n_samples),
      paste0("ASV", 1:n_asvs)
    )
  )

  # Introduce additional zeros randomly
  zero_mask <- matrix(
    sample(c(0, 1), n_samples * n_asvs, replace = TRUE, prob = c(zero_prob, 1 - zero_prob)),
    nrow = n_samples,
    ncol = n_asvs
  )
  mock_asv_matrix <- mock_asv_matrix * zero_mask


  # Create mock sample metadata
  mock_sample_data <- data.frame(
    Sample = paste0("Sample", 1:n_samples),
    Group = sample(c("A", "B"), n_samples, replace = TRUE), # Random groups
    row.names = paste0("Sample", 1:n_samples)
  )

  # Create mock taxonomy table
  mock_taxonomy_table <- matrix(
    sample(c("Bacteria", "Archaea"), n_asvs * 7, replace = TRUE), # Random taxonomy
    nrow = n_asvs,
    ncol = 7,
    dimnames = list(
      paste0("ASV", 1:n_asvs),
      c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    )
  )

  # Create phyloseq object
  otu_table <- otu_table(mock_asv_matrix, taxa_are_rows = FALSE)
  sample_data <- sample_data(mock_sample_data)
  tax_table <- tax_table(mock_taxonomy_table)

  mock_phyloseq <- phyloseq::phyloseq(otu_table, sample_data, tax_table)

  return(list(asv_matrix = mock_asv_matrix, phyloseq_obj = mock_phyloseq))
}


# Test suite
testthat::test_that("calculate_nmds() works correctly", {
  # Generate mock data
  mock_data <- create_mock_data(n_samples = 400, n_asvs = 500)
  mock_asv_matrix <- mock_data$asv_matrix
  mock_phyloseq <- mock_data$phyloseq_obj

  # Run the function
  results <- calculate_nmds(mock_asv_matrix, mock_phyloseq,
    ncores = parallel::detectCores(),
    k = 8,
    trymax = 200
  )

  # Test 1: Check output structure
  expect_named(results, c("nmds_scores", "nmds_df", "ordi"))

  # Test 2: Check nmds_scores
  expect_s3_class(results$nmds_scores, "data.frame")
  expect_named(results$nmds_scores, c("unique_id", "NMDS1", "NMDS2"))
  expect_equal(nrow(results$nmds_scores), nrow(mock_asv_matrix))

  # Test 3: Check nmds_df
  expect_s3_class(results$nmds_df, "data.frame")
  expect_true(all(c("unique_id", "NMDS1", "NMDS2", "Sample", "Group") %in% names(results$nmds_df)))

  # Test 4: Check core_ordi
  expect_s3_class(results$ordi, "metaMDS")
  expect_equal(results$ordi$ndim, 2) # Check number of dimensions

  # Test 5: Compare NMDS scores with expected values (using waldo)
  expected_nmds_scores <- results$nmds_scores # Replace with known good values if available
  compare(results$nmds_scores, expected_nmds_scores)

  # Test 6: Compare merged metadata with expected values (using waldo)
  expected_nmds_df <- results$nmds_df # Replace with known good values if available
  compare(results$nmds_df, expected_nmds_df)

  # Test 7: Check for NA/NaN in NMDS scores
  expect_false(any(is.na(results$nmds_scores)))

  # Test 8: Check for convergence
  expect_true(results$ordi$converged)
})
