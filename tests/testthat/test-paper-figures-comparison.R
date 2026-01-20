# Tests for Paper Figures Refactoring
# Ensures that the refactored 009_paper_figures_refactored.R script
# produces equivalent outputs to the original 009_paper_figures.R

library(testthat)
library(ggplot2)
library(dplyr)

# Source brc_themes.R once at the top for all tests
source(here::here("R/functions/brc_themes.R"))

# Tests for theme helpers ----
test_that("brc_theme_title_size returns valid theme", {
  skip_if_not_installed("ggplot2")
  
  theme_obj <- brc_theme_title_size()
  
  expect_s3_class(theme_obj, "theme")
  expect_type(theme_obj, "list")
})

test_that("brc_crop_labels returns correct labels", {
  
  labels <- brc_crop_labels()
  
  expect_type(labels, "character")
  expect_length(labels, 8)
  expect_true("Corn" %in% labels)
  expect_true("Switchgrass" %in% labels)
  expect_true("Restored Prairie" %in% labels)
})

test_that("brc_name_labels returns correct BRC labels", {
  
  labels <- brc_name_labels()
  
  expect_type(labels, "character")
  expect_true("CABBI" %in% labels)
  expect_true("GLBRC" %in% labels)
  expect_equal(labels["cabbi"], c(cabbi = "CABBI"))
})

test_that("brc_full_labels includes both BRC and crop labels", {
  
  labels <- brc_full_labels()
  
  expect_type(labels, "character")
  # Should have BRC labels
  expect_true("CABBI" %in% labels)
  # Should have crop labels
  expect_true("Corn" %in% labels)
  expect_true("Switchgrass" %in% labels)
})

test_that("brc_scale_y_percent creates valid scale", {
  skip_if_not_installed("ggplot2")
  
  
  scale_obj <- brc_scale_y_percent()
  
  expect_s3_class(scale_obj, "ScaleContinuous")
})

test_that("brc_scale_color_crops creates valid color scale", {
  skip_if_not_installed("ggplot2")
  
  
  scale_obj <- brc_scale_color_crops()
  
  expect_s3_class(scale_obj, "ScaleDiscrete")
})

test_that("brc_scale_fill_core creates valid fill scale", {
  skip_if_not_installed("ggplot2")
  
  
  scale_obj <- brc_scale_fill_core()
  
  expect_s3_class(scale_obj, "ScaleDiscrete")
})

test_that("pre-configured theme lists are valid", {
  skip_if_not_installed("ggplot2")
  
  
  # Test each pre-configured theme
  theme_core50_crops <- brc_theme_core50_crops()
  expect_type(theme_core50_crops, "list")
  expect_true(length(theme_core50_crops) > 0)
  
  theme_core50_brc <- brc_theme_core50_brc()
  expect_type(theme_core50_brc, "list")
  
  theme_t60_core <- brc_theme_threshold60_core()
  expect_type(theme_t60_core, "list")
  
  theme_t60_crops <- brc_theme_threshold60_crops()
  expect_type(theme_t60_crops, "list")
  
  theme_t60_brc <- brc_theme_threshold60_brc()
  expect_type(theme_t60_brc, "list")
  
  theme_t100 <- brc_theme_threshold100_palette()
  expect_type(theme_t100, "list")
})

# Tests for theme application ----
test_that("themes can be applied to ggplot objects", {
  skip_if_not_installed("ggplot2")
  
  
  # Create a simple plot
  p <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
  
  # Apply theme
  expect_no_error({
    p + brc_theme_title_size()
  })
  
  # Apply scale
  expect_no_error({
    p + brc_scale_y_percent(limits = c(0, 1))
  })
})

test_that("theme lists can be added to plots", {
  skip_if_not_installed("ggplot2")
  
  
  p <- ggplot(mtcars, aes(mpg, wt, color = factor(cyl))) + geom_point()
  
  theme_list <- brc_theme_core50_crops()
  
  expect_no_error({
    p + theme_list
  })
})

# Tests for consistency between original and refactored ----
test_that("refactored script uses same theme values as original", {
  
  # Get theme from helper
  theme_obj <- brc_theme_title_size()
  
  # Original values from 009_paper_figures.R
  expect_equal(theme_obj$title$size, 12)
  expect_equal(theme_obj$legend.title$size, 12)
  expect_equal(theme_obj$axis.title.x$size, 14)
  expect_equal(theme_obj$axis.title.y$size, 14)
  expect_equal(theme_obj$axis.text.x$size, 12)
  expect_equal(theme_obj$axis.text.y$size, 12)
})

test_that("crop labels match original values", {
  
  labels <- brc_crop_labels()
  
  # Original labels from 009_paper_figures.R lines 74-83
  expected_labels <- c(
    "Corn",
    "Miscanthus",
    "Poplar",
    "Restored Prairie",
    "Sorghum",
    "Sorghum + Rye",
    "Soy",
    "Switchgrass"
  )
  
  expect_equal(labels, expected_labels)
})

# Tests for plot generation helpers ----
test_that("brc_gg_ordi can use theme helpers", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  
  # Create simple test data
  test_df <- data.frame(
    NMDS1 = rnorm(20),
    NMDS2 = rnorm(20),
    brc = factor(sample(c("CABBI", "GLBRC"), 20, replace = TRUE)),
    crop = factor(sample(c("Corn", "Soy"), 20, replace = TRUE))
  )
  
  source(here::here("R/functions/brc_gg_ordi.R"))
  
  # Create plot with theme
  expect_no_error({
    p <- brc_gg_ordi(
      .data = test_df,
      ordi = "NMDS",
      .color = brc,
      .size = 2
    ) + brc_theme_title_size()
  })
})

# Test for script structure ----
test_that("refactored script has correct structure", {
  refactored_path <- here::here("R/analysis/009_paper_figures_refactored.R")
  
  expect_true(file.exists(refactored_path))
  
  # Read the script
  script_lines <- readLines(refactored_path)
  
  # Check that it sources brc_themes.R
  expect_true(any(grepl("brc_themes", script_lines)))
  
  # Check that syntax errors are fixed (no invalid ---- patterns)
  # Original had "----raw" and "----box_2" which are invalid
  invalid_patterns <- grepl("^-----$", script_lines) | 
                      grepl("----raw", script_lines) |
                      grepl("----box_2", script_lines)
  expect_false(any(invalid_patterns))
  
  # Check that key sections exist
  expect_true(any(grepl("THEMES AND LABELS", script_lines)))
  expect_true(any(grepl("PRE-PROCESSING", script_lines)))
  expect_true(any(grepl("ARRANGEMENT", script_lines)))
})

test_that("refactored script uses helper functions", {
  refactored_path <- here::here("R/analysis/009_paper_figures_refactored.R")
  script_lines <- readLines(refactored_path)
  script_text <- paste(script_lines, collapse = "\n")
  
  # Check usage of theme helpers
  expect_true(grepl("brc_crop_labels\\(\\)", script_text))
  expect_true(grepl("brc_theme_core50_crops\\(\\)", script_text))
  expect_true(grepl("brc_theme_title_size\\(\\)", script_text))
  expect_true(grepl("brc_full_labels\\(\\)", script_text))
})

# Integration test for complete workflow ----
test_that("theme helpers integrate with ordination workflow", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("vegan")
  
  # Create test ordination data
  test_nmds <- data.frame(
    NMDS1 = rnorm(30),
    NMDS2 = rnorm(30),
    brc = factor(rep(c("CABBI", "GLBRC", "CBI"), 10)),
    crop = factor(rep(c("Corn", "Soy", "Switchgrass"), 10))
  )
  
  test_pcoa <- data.frame(
    Dim1 = rnorm(30),
    Dim2 = rnorm(30),
    brc = factor(rep(c("CABBI", "GLBRC", "CBI"), 10)),
    crop = factor(rep(c("Corn", "Soy", "Switchgrass"), 10))
  )
  
  source(here::here("R/functions/brc_gg_ordi.R"))
  source(here::here("R/functions/brc_paper_ordinations.R"))
  
  # Test that workflow runs without error
  expect_no_error({
    plots <- brc_paper_ordinations(
      nmds_data = test_nmds,
      pcoa_data = test_pcoa,
      nmds_colors = c("#E7B800", "#FC4E07"),
      pcoa_colors = c("#E7B800", "#383961", "#FC4E07"),
      nmds_labels = c("CABBI", "GLBRC"),
      pcoa_labels = c("CABBI", "CBI", "GLBRC"),
      nmds_crop_labels = brc_crop_labels(),
      pcoa_crop_labels = brc_crop_labels(),
      point_size = 1.25,
      crop_theme = brc_theme_core50_crops(),
      brc_theme = brc_theme_core50_brc()
    )
  })
  
  # Check output structure
  expect_type(plots, "list")
  expect_true("nmds" %in% names(plots))
  expect_true("pcoa" %in% names(plots))
})

# Test for backward compatibility ----
test_that("refactored code maintains backward compatibility", {
  
  # The old inline definition should match the new function
  old_title_size <- list(theme(
    title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  ))
  
  new_title_size <- brc_theme_title_size()
  
  # Check that key properties match
  expect_equal(old_title_size[[1]]$title$size, new_title_size$title$size)
  expect_equal(old_title_size[[1]]$axis.title.x$size, new_title_size$axis.title.x$size)
})
