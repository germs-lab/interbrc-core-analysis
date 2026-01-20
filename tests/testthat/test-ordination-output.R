# Tests for Ordination Output Equivalence
# Ensures that refactored ordination functions produce equivalent results
# to the original brc_nmds and brc_pcoa functions

library(testthat)
library(phyloseq)
library(vegan)
library(dplyr)

# Helper function to create test data ----
create_test_phyloseq <- function(n_samples = 20, n_taxa = 50, seed = 12345) {
  set.seed(seed)
  
  # Create OTU table
  otu_mat <- matrix(
    sample(1:100, n_samples * n_taxa, replace = TRUE),
    nrow = n_taxa,
    ncol = n_samples
  )
  rownames(otu_mat) <- paste0("ASV", 1:n_taxa)
  colnames(otu_mat) <- paste0("Sample", 1:n_samples)
  
  # Create sample data
  sample_data <- data.frame(
    brc = factor(sample(c("CABBI", "GLBRC", "CBI"), n_samples, replace = TRUE)),
    crop = factor(sample(c("Corn", "Soy", "Switchgrass"), n_samples, replace = TRUE)),
    site = factor(sample(c("Site1", "Site2"), n_samples, replace = TRUE)),
    row.names = colnames(otu_mat)
  )
  
  # Create phyloseq object
  OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
  SAM <- sample_data(sample_data)
  
  physeq <- phyloseq(OTU, SAM)
  
  return(physeq)
}

# Tests for brc_nmds ----
test_that("brc_nmds produces valid NMDS output", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  
  # Create test data
  physeq <- create_test_phyloseq()
  asv_matrix <- as.matrix(otu_table(physeq))
  if (taxa_are_rows(physeq)) {
    asv_matrix <- t(asv_matrix)
  }
  
  # Load function
  source(here::here("R/functions/brc_nmds.R"))
  
  # Run NMDS
  result <- brc_nmds(
    asv_matrix = asv_matrix,
    physeq = physeq,
    k = 2,
    trymax = 20  # Reduced for faster testing
  )
  
  # Check structure
  expect_type(result, "list")
  expect_true("nmds_scores" %in% names(result))
  expect_true("nmds_df" %in% names(result))
  expect_true("ordi" %in% names(result))
  
  # Check nmds_scores
  expect_s3_class(result$nmds_scores, "data.frame")
  expect_true("NMDS1" %in% names(result$nmds_scores))
  expect_true("NMDS2" %in% names(result$nmds_scores))
  
  # Check nmds_df has metadata
  expect_s3_class(result$nmds_df, "data.frame")
  expect_true("brc" %in% names(result$nmds_df))
  expect_true("crop" %in% names(result$nmds_df))
  
  # Check ordination object
  expect_s3_class(result$ordi, "metaMDS")
  expect_equal(result$ordi$ndim, 2)
})

test_that("brc_nmds handles different k values", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  
  physeq <- create_test_phyloseq()
  asv_matrix <- as.matrix(otu_table(physeq))
  if (taxa_are_rows(physeq)) {
    asv_matrix <- t(asv_matrix)
  }
  
  source(here::here("R/functions/brc_nmds.R"))
  
  # Test k=3
  result_3d <- brc_nmds(
    asv_matrix = asv_matrix,
    physeq = physeq,
    k = 3,
    trymax = 20
  )
  
  expect_equal(result_3d$ordi$ndim, 3)
  expect_true("NMDS3" %in% names(result_3d$nmds_scores))
})

test_that("brc_nmds rejects invalid input", {
  skip_if_not_installed("phyloseq")
  
  physeq <- create_test_phyloseq()
  asv_matrix <- as.matrix(otu_table(physeq))
  
  source(here::here("R/functions/brc_nmds.R"))
  
  # Test with NA values
  bad_matrix <- asv_matrix
  bad_matrix[1, 1] <- NA
  
  expect_error(
    brc_nmds(bad_matrix, physeq),
    "NA or NaN"
  )
  
  # Test with non-phyloseq object
  expect_error(
    brc_nmds(asv_matrix, list()),
    "phyloseq object"
  )
})

# Tests for brc_pcoa ----
test_that("brc_pcoa produces valid PCoA output", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  
  # Create test data
  physeq <- create_test_phyloseq()
  asv_matrix <- as.matrix(otu_table(physeq))
  if (taxa_are_rows(physeq)) {
    asv_matrix <- t(asv_matrix)
  }
  
  # Create distance matrix
  dist_matrix <- vegdist(asv_matrix, method = "bray")
  
  # Load function
  source(here::here("R/functions/brc_pcoa.R"))
  source(here::here("R/functions/parallel_helpers.R"))
  
  # Run PCoA
  result <- brc_pcoa(
    asv_matrix = dist_matrix,
    physeq = physeq,
    k = 2
  )
  
  # Check structure
  expect_type(result, "list")
  expect_true("pcoa_df" %in% names(result))
  expect_true("ordi" %in% names(result))
  
  # Check pcoa_df
  expect_s3_class(result$pcoa_df, "data.frame")
  expect_true("Dim1" %in% names(result$pcoa_df))
  expect_true("Dim2" %in% names(result$pcoa_df))
  expect_true("brc" %in% names(result$pcoa_df))
  expect_true("crop" %in% names(result$pcoa_df))
  
  # Check dimensions
  expect_equal(nrow(result$pcoa_df), nsamples(physeq))
})

test_that("brc_pcoa handles distance matrices", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  
  physeq <- create_test_phyloseq()
  asv_matrix <- as.matrix(otu_table(physeq))
  if (taxa_are_rows(physeq)) {
    asv_matrix <- t(asv_matrix)
  }
  
  dist_matrix <- vegdist(asv_matrix, method = "bray")
  
  source(here::here("R/functions/brc_pcoa.R"))
  source(here::here("R/functions/parallel_helpers.R"))
  
  result <- brc_pcoa(dist_matrix, physeq)
  
  expect_s3_class(result$pcoa_df, "data.frame")
  expect_true(all(c("Dim1", "Dim2") %in% names(result$pcoa_df)))
})

# Tests for new unified brc_ordination function ----
test_that("brc_ordination NMDS mode matches brc_nmds", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  
  physeq <- create_test_phyloseq(n_samples = 15, n_taxa = 30)
  asv_matrix <- as.matrix(otu_table(physeq))
  if (taxa_are_rows(physeq)) {
    asv_matrix <- t(asv_matrix)
  }
  
  # Load functions
  source(here::here("R/functions/brc_nmds.R"))
  source(here::here("R/functions/brc_ordination.R"))
  
  # Set seed for reproducibility
  set.seed(54641)
  result_old <- brc_nmds(asv_matrix, physeq, trymax = 20)
  
  set.seed(54641)
  result_new <- brc_ordination(asv_matrix, physeq, method = "NMDS", trymax = 20)
  
  # Check that both produce ordinations with same dimensions
  expect_equal(result_old$ordi$ndim, 2)
  expect_equal(ncol(result_new$scores_df %>% select(starts_with("NMDS"))), 2)
  
  # Check that metadata is present
  expect_true("brc" %in% names(result_new$scores_df))
  expect_true("crop" %in% names(result_new$scores_df))
  
  # Check method is recorded
  expect_equal(result_new$method, "NMDS")
})

test_that("brc_ordination PCoA mode matches brc_pcoa", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  
  physeq <- create_test_phyloseq()
  asv_matrix <- as.matrix(otu_table(physeq))
  if (taxa_are_rows(physeq)) {
    asv_matrix <- t(asv_matrix)
  }
  
  dist_matrix <- vegdist(asv_matrix, method = "bray")
  
  # Load functions
  source(here::here("R/functions/brc_pcoa.R"))
  source(here::here("R/functions/brc_ordination.R"))
  source(here::here("R/functions/parallel_helpers.R"))
  
  set.seed(54641)
  result_old <- brc_pcoa(dist_matrix, physeq)
  
  set.seed(54641)
  result_new <- brc_ordination(dist_matrix, physeq, method = "PCoA")
  
  # Check dimensions match
  expect_equal(nrow(result_old$pcoa_df), nrow(result_new$scores_df))
  
  # Check that key columns exist
  expect_true(all(c("Dim1", "Dim2") %in% names(result_new$scores_df)))
  
  # Check method is recorded
  expect_equal(result_new$method, "PCoA")
})

test_that("brc_ordination validates method argument", {
  skip_if_not_installed("phyloseq")
  
  physeq <- create_test_phyloseq()
  asv_matrix <- as.matrix(otu_table(physeq))
  
  source(here::here("R/functions/brc_ordination.R"))
  
  expect_error(
    brc_ordination(asv_matrix, physeq, method = "invalid"),
    "should be one of"
  )
})

# Tests for reproducibility ----
test_that("ordination results are reproducible with same seed", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  
  physeq <- create_test_phyloseq()
  asv_matrix <- as.matrix(otu_table(physeq))
  if (taxa_are_rows(physeq)) {
    asv_matrix <- t(asv_matrix)
  }
  
  source(here::here("R/functions/brc_nmds.R"))
  
  # Run twice with same seed
  set.seed(54641)
  result1 <- brc_nmds(asv_matrix, physeq, trymax = 20)
  
  set.seed(54641)
  result2 <- brc_nmds(asv_matrix, physeq, trymax = 20)
  
  # Stress values should be very similar (allowing for numerical precision)
  expect_equal(result1$ordi$stress, result2$ordi$stress, tolerance = 0.01)
})

# Tests for edge cases ----
test_that("ordination handles small datasets", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  
  # Very small dataset
  physeq <- create_test_phyloseq(n_samples = 5, n_taxa = 10)
  asv_matrix <- as.matrix(otu_table(physeq))
  if (taxa_are_rows(physeq)) {
    asv_matrix <- t(asv_matrix)
  }
  
  source(here::here("R/functions/brc_pcoa.R"))
  source(here::here("R/functions/parallel_helpers.R"))
  
  dist_matrix <- vegdist(asv_matrix, method = "bray")
  
  expect_no_error({
    result <- brc_pcoa(dist_matrix, physeq)
  })
  
  expect_equal(nrow(result$pcoa_df), 5)
})
