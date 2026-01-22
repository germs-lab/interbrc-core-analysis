# Quick Start Guide: Using Refactored BRC Functions

**Last Updated:** 2025-01-20

This guide shows you how to use the newly refactored functions for BRC core analysis.

---

## 🎨 Using Theme Helpers

### Basic Usage

```r
# Load the theme helpers
source("R/functions/brc_themes.R")

library(ggplot2)

# Create a plot with standard BRC styling
my_plot <- ggplot(data, aes(x = abundance, y = occupancy)) +
  geom_point() +
  brc_theme_title_size() +
  ggtitle("My Analysis")
```

### Using Pre-configured Themes

```r
# For 50 core ASVs analysis (crops)
core_plot <- ggplot(core_data, aes(x, y, color = crop)) +
  geom_point() +
  brc_theme_core50_crops()  # Includes title, colors, and sizing

# For threshold-based analysis (60%)
threshold_plot <- ggplot(threshold_data, aes(x, y, fill = type)) +
  geom_bar(stat = "identity") +
  brc_theme_threshold60_core()  # Includes percent scale and title
```

### Custom Combinations

```r
# Mix and match components
custom_plot <- ggplot(data, aes(log10(abundance), occupancy, fill = core_type)) +
  geom_point(shape = 21, size = 3) +
  brc_scale_fill_core() +           # Core vs non-core colors
  brc_scale_y_percent() +            # 0-100% y-axis
  brc_theme_title_size() +           # Standard text sizes
  ggtitle("Custom Analysis") +
  labs(x = "log10(Abundance)")
```

### Using Label Helpers

```r
# Get standardized crop labels
crop_labels <- brc_crop_labels()
# Returns: c("Corn", "Miscanthus", "Poplar", ...)

# Use in plots
ggplot(data, aes(x = crop, y = value)) +
  geom_boxplot() +
  scale_x_discrete(labels = crop_labels)

# For faceting with nested labels
library(ggh4x)

ggplot(data, aes(x, y)) +
  geom_point() +
  facet_nested(~ brc + crop, labeller = as_labeller(brc_full_labels()))
```

---

## 🔬 Using Unified Ordination Interface

### NMDS Analysis

```r
# Load the function
source("R/functions/brc_ordination.R")

library(phyloseq)
library(vegan)

# Prepare data
asv_matrix <- as.matrix(otu_table(my_phyloseq))
if (taxa_are_rows(my_phyloseq)) {
  asv_matrix <- t(asv_matrix)
}

# Run NMDS
nmds_result <- brc_ordination(
  data = asv_matrix,
  physeq = my_phyloseq,
  method = "NMDS",
  k = 2,
  trymax = 100
)

# Access results
plot_data <- nmds_result$scores_df  # Data with metadata
stress <- nmds_result$metadata$stress  # Stress value
converged <- nmds_result$metadata$converged  # Convergence status

# Check convergence
if (nmds_result$metadata$converged) {
  cat("NMDS converged! Stress:", stress, "\n")
} else {
  cat("Warning: NMDS did not converge\n")
}
```

### PCoA Analysis

```r
# Calculate distance matrix
dist_matrix <- vegdist(asv_matrix, method = "bray")

# Run PCoA
pcoa_result <- brc_ordination(
  data = dist_matrix,
  physeq = my_phyloseq,
  method = "PCoA",
  k = 2
)

# Access results
plot_data <- pcoa_result$scores_df
variance_explained <- pcoa_result$metadata$variance_explained

# Check variance explained
cat("Axis 1:", round(variance_explained[1] * 100, 1), "%\n")
cat("Axis 2:", round(variance_explained[2] * 100, 1), "%\n")
```

### Using Legacy Functions (Backward Compatibility)

```r
# Old way still works!
source("R/functions/brc_nmds.R")
source("R/functions/brc_pcoa.R")

# These work exactly as before
nmds_old <- brc_nmds(asv_matrix, my_phyloseq)
pcoa_old <- brc_pcoa(dist_matrix, my_phyloseq)
```

---

## 📊 Creating Paper-Ready Figures

### Complete Workflow Example

```r
# Load all necessary functions
source("R/functions/brc_themes.R")
source("R/functions/brc_nmds.R")
source("R/functions/brc_pcoa.R")
source("R/functions/brc_gg_ordi.R")
source("R/functions/brc_paper_ordinations.R")

# 1. Calculate ordinations
nmds_result <- brc_nmds(asv_matrix, my_phyloseq)
pcoa_result <- brc_pcoa(dist_matrix, my_phyloseq)

# 2. Create matching plot sets
ordination_plots <- brc_paper_ordinations(
  nmds_data = nmds_result$nmds_df,
  pcoa_data = pcoa_result$pcoa_df,
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

# 3. Access individual plots
ordination_plots$nmds$crops   # NMDS colored by crop
ordination_plots$nmds$brc     # NMDS colored by BRC
ordination_plots$pcoa$crops   # PCoA colored by crop
ordination_plots$pcoa$brc     # PCoA colored by BRC

# 4. Save plots
library(ggpubr)

combined <- ggarrange(
  ordination_plots$pcoa$crops,
  ordination_plots$pcoa$brc,
  ncol = 2,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "bottom"
)

ggsave("my_ordination_figure.png", combined, width = 300, height = 150, units = "mm", dpi = 300)
```

### Single Ordination Plot

```r
# Create a standalone ordination plot
source("R/functions/brc_gg_ordi.R")
source("R/functions/brc_themes.R")

single_plot <- brc_gg_ordi(
  .data = nmds_result$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .size = 2,
  .alpha = 0.7,
  color_values = c("#E7B800", "#383961", "#FC4E07"),
  .labels = c("CABBI", "CBI", "GLBRC")
) +
  brc_theme_title_size() +
  ggtitle("NMDS of ASV Communities")
```

---

## 🧪 Running Tests

### Run All Tests

```r
# Install testthat if needed
if (!require("testthat")) install.packages("testthat")

# Run all tests
testthat::test_dir("tests/testthat")
```

### Run Specific Test Files

```r
# Test ordination functions
testthat::test_file("tests/testthat/test-ordination-output.R")

# Test theme helpers
testthat::test_file("tests/testthat/test-paper-figures-comparison.R")
```

### Expected Output

```
✓ | F W S  OK | Context
✓ |        15 | test-ordination-output
✓ |        17 | test-paper-figures-comparison

══ Results ══════════════════════════════════════════
Duration: 15.3 s

[ FAIL 0 | WARN 0 | SKIP 0 | PASS 32 ]
```

---

## 📁 Complete Analysis Script Template

```r
#!/usr/bin/env Rscript
#########################################################
# BRC Core Microbiome Analysis Template
# Using refactored functions
#########################################################

# Load packages
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggpubr)

# Source BRC functions
source("R/functions/brc_themes.R")
source("R/functions/brc_ordination.R")
source("R/functions/brc_gg_ordi.R")
source("R/functions/brc_paper_ordinations.R")

# Load data
my_phyloseq <- readRDS("data/my_phyloseq.rds")

# Filter phyloseq (example: >500 reads/sample, >20 reads/ASV)
ps_filtered <- my_phyloseq %>%
  prune_samples(sample_sums(.) >= 500, .) %>%
  filter_taxa(function(x) sum(x) >= 20, prune = TRUE)

# Prepare ASV matrix
asv_matrix <- as.matrix(otu_table(ps_filtered))
if (taxa_are_rows(ps_filtered)) {
  asv_matrix <- t(asv_matrix)
}

# Calculate distance matrix
dist_matrix <- vegdist(asv_matrix, method = "bray")

# Run ordinations
nmds_result <- brc_ordination(
  data = asv_matrix,
  physeq = ps_filtered,
  method = "NMDS",
  k = 2,
  trymax = 100
)

pcoa_result <- brc_ordination(
  data = dist_matrix,
  physeq = ps_filtered,
  method = "PCoA",
  k = 2
)

# Create plots
plots <- brc_paper_ordinations(
  nmds_data = nmds_result$scores_df,
  pcoa_data = pcoa_result$scores_df,
  nmds_colors = c("#E7B800", "#FC4E07"),
  pcoa_colors = c("#E7B800", "#FC4E07"),
  nmds_labels = c("Site A", "Site B"),
  pcoa_labels = c("Site A", "Site B"),
  crop_theme = list(
    ggtitle("Ordination by Crop"),
    brc_theme_title_size()
  ),
  brc_theme = list(
    ggtitle("Ordination by Site"),
    brc_theme_title_size()
  )
)

# Arrange and save
final_figure <- ggarrange(
  plots$pcoa$crops,
  plots$pcoa$brc,
  ncol = 2,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  "output/ordination_figure.png",
  final_figure,
  width = 300,
  height = 150,
  units = "mm",
  dpi = 300
)

# Print summary
cat("\n=== Analysis Complete ===\n")
cat("NMDS Stress:", nmds_result$metadata$stress, "\n")
cat("NMDS Converged:", nmds_result$metadata$converged, "\n")
cat("PCoA Variance Explained (Axis 1):", 
    round(pcoa_result$metadata$variance_explained[1] * 100, 1), "%\n")
cat("PCoA Variance Explained (Axis 2):", 
    round(pcoa_result$metadata$variance_explained[2] * 100, 1), "%\n")
cat("Plots saved to: output/ordination_figure.png\n")
```

---

## 🔍 Common Patterns

### Pattern 1: Standard PCoA Workflow

```r
# Distance → PCoA → Plot
dist <- vegdist(asv_matrix, "bray")
pcoa <- brc_ordination(dist, physeq, method = "PCoA")
plot <- brc_gg_ordi(pcoa$scores_df, ordi = "PCoA", .color = treatment)
```

### Pattern 2: Standard NMDS Workflow

```r
# Matrix → NMDS → Plot
nmds <- brc_ordination(asv_matrix, physeq, method = "NMDS")
plot <- brc_gg_ordi(nmds$scores_df, ordi = "NMDS", .color = site)
```

### Pattern 3: Paper Figure Set

```r
# Calculate → Create pairs → Arrange → Save
nmds <- brc_nmds(asv_matrix, physeq)
pcoa <- brc_pcoa(dist_matrix, physeq)
plots <- brc_paper_ordinations(nmds$nmds_df, pcoa$pcoa_df, ...)
final <- ggarrange(plots$pcoa$crops, plots$pcoa$brc, ...)
ggsave("figure.png", final)
```

---

## 💡 Tips and Best Practices

### 1. Always Check Convergence (NMDS)

```r
nmds <- brc_ordination(data, physeq, method = "NMDS", trymax = 100)
if (!nmds$metadata$converged) {
  # Try again with more random starts
  nmds <- brc_ordination(data, physeq, method = "NMDS", trymax = 500)
}
```

### 2. Use Consistent Themes

```r
# Define once, reuse everywhere
my_theme <- list(
  brc_theme_title_size(),
  theme(legend.position = "bottom")
)

plot1 + my_theme
plot2 + my_theme
plot3 + my_theme
```

### 3. Leverage Pre-configured Themes

```r
# Instead of manually building themes every time
plot + brc_theme_core50_crops()  # Much cleaner!
```

### 4. Save Ordination Results

```r
# Save for reuse
saveRDS(nmds_result, "output/nmds_result.rds")
saveRDS(pcoa_result, "output/pcoa_result.rds")

# Load later
nmds <- readRDS("output/nmds_result.rds")
```

---

## 📚 Additional Resources

- **Work Plan:** `docs/REFACTORING_WORK_PLAN.md`
- **Summary:** `docs/REFACTORING_SUMMARY.md`
- **Validation:** `docs/REFACTORING_CHECKLIST.md`
- **Original Script:** `R/analysis/009_paper_figures.R`
- **Refactored Script:** `R/analysis/009_paper_figures_refactored.R`

---

## ❓ Troubleshooting

### Error: "NMDS did not converge"
**Solution:** Increase `trymax` parameter
```r
nmds <- brc_ordination(data, physeq, method = "NMDS", trymax = 500)
```

### Error: "data contains NA or NaN values"
**Solution:** Remove NAs before analysis
```r
asv_matrix[is.na(asv_matrix)] <- 0
```

### Error: "No valid columns remaining after filtering"
**Solution:** Check for all-zero samples/taxa
```r
# Remove zero-sum samples
ps <- prune_samples(sample_sums(ps) > 0, ps)
# Remove zero-sum taxa
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
```

---

**Questions?** Refer to the full documentation in `docs/` or examine the test files in `tests/testthat/` for more examples.
