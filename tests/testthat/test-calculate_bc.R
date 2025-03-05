library(testthat)

# Helper function for Bray-Curtis calculation (as defined in ExtracCore.R)
calculate_bc <- function(matrix, nReads) {
  if (nrow(matrix) == 0) {
    cli::cli_alert_warning("{.arg matrix} is empty. Enter a non-empty matrix.")
    return(list(values = numeric(0), names = character(0)))
  }

  bc_values <- apply(combn(ncol(matrix), 2), 2, function(cols) {
    sum(abs(matrix[, cols[1]] - matrix[, cols[2]])) / (2 * nReads)
  })

  x_names <- apply(combn(ncol(matrix), 2), 2, function(cols) {
    paste(colnames(matrix)[cols], collapse = "-")
  })

  list(values = bc_values, names = x_names)
}



# Unit Tests
test_that("calculate_bc() correctly computes Bray-Curtis between samples", {
  # Test 1: 1 OTU x 2 samples
  test_matrix <- matrix(
    c(10, 20),
    nrow = 1, # 1 OTU (row)
    ncol = 2, # 2 samples (columns)
    dimnames = list("OTU1", c("SampleA", "SampleB"))
  )
  nReads <- 30 # Total reads per sample after rarefaction
  result <- calculate_bc(test_matrix, nReads)

  # Expected BC: |10-20|/(2*30) = 10/60 = 0.1666667
  expect_equal(result$values, 10 / 60)
  expect_equal(result$names, "SampleA-SampleB")

  # Test 2: 2 OTUs x 3 samples
  test_matrix <- matrix(
    c(10, 5, 20, 15, 30, 25),
    nrow = 2,
    ncol = 3,
    dimnames = list(
      c("OTU1", "OTU2"),
      c("SampleA", "SampleB", "SampleC")
    )
  )
  result <- calculate_bc(test_matrix, nReads = 30)

  # Expected pairwise comparisons (3 pairs)
  expect_equal(length(result$values), 3) # combn(3 samples, 2) = 3 pairs

  # Verify BC for SampleA vs SampleB:
  # OTU1: |10-20| = 10, OTU2: |5-15| = 10 → sum = 20 → 20/(2*30) = 0.333...
  expect_equal(result$values[1], 20 / 60, tolerance = 0.001)
  expect_equal(result$names, c("SampleA-SampleB", "SampleA-SampleC", "SampleB-SampleC"))

  # Test 3: Empty matrix (0 OTUs x 3 samples)
  empty_matrix <- matrix(nrow = 0, ncol = 3)
  result <- calculate_bc(empty_matrix, nReads = 30)
  expect_equal(length(result$values), 0) # No OTUs → BC = 0 for all pairs
})
