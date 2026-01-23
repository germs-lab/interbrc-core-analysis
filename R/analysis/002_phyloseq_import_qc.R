#########################################################
# PHYLOSEQ OBJECT CREATION AND FILTERING
# Import and organize sequence data into phyloseq format
#
# Project:  Inter-BRC-Core-Microbiome
# Author: B. Kristy
# Modified by: Bolívar Aponte Rolón
# Date: 2025-02-20
# Last modified: 2026-01-16
#########################################################

# DESCRIPTION:
# This script creates phyloseq objects from taxonomy, metadata, and feature tables.
# It performs initial filtering steps including singleton removal and plant contaminant filtering.

# SETUP AND DEPENDENCIES ----

source("R/utils/000_setup.R")


# DATA IMPORT ----

# Import taxonomy, metadata, and feature tables
taxonomy <- read.delim("data/input/taxonomy.tsv", comment.char = "#")
metadata <- read.delim("data/input/metadata.tsv", row.names = 1)
feature.table.modified <- read.delim("data/input/feature-table-modified.tsv")

# Format feature table for phyloseq
feature.table.modified <- feature.table.modified %>%
  remove_rownames() %>%
  column_to_rownames(var = "OTU.ID")


# TAXONOMY PROCESSING ----

# Define taxonomy ranks
ranks <- c(
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species"
)

# Clean and format taxonomy data
taxonomy_modified <- taxonomy %>%
  mutate_at("Taxon", str_replace_all, "[a-z]__", "") %>%
  separate(Taxon, sep = ";", into = ranks, remove = TRUE) %>%
  column_to_rownames(var = "Feature.ID") %>%
  as.matrix()


# PHYLOSEQ OBJECT CREATION ----

# Convert processed data to phyloseq components
TAX <- tax_table(taxonomy_modified)
feature.table.modified <- as.matrix(feature.table.modified)
OTU <- otu_table(feature.table.modified, taxa_are_rows = T)
metadata <- sample_data(metadata)

# Combine into single phyloseq object
unfiltered_phyloseq <- phyloseq(OTU, TAX, metadata)

# Clean up column names
colnames(sample_data(unfiltered_phyloseq)) <- unfiltered_phyloseq@sam_data |>
  janitor::clean_names() |>
  colnames()


# PHYLOSEQ FILTERING ----

# Remove singletons
phyloseq_removed_singletons <- prune_taxa(
  taxa_sums(unfiltered_phyloseq) > 1,
  phyloseq
)
tax.remove <- ntaxa(unfiltered_phyloseq) - ntaxa(phyloseq_removed_singletons) # 1,675 ASVs

# Remove plant contaminants (Eukaryota and unassigned)
phyloseq_removed_singletons_plant <- subset_taxa(
  phyloseq_removed_singletons,
  kingdom != "Eukaryota" & kingdom != "Unassigned"
)

n.filtered <- ntaxa(phyloseq_removed_singletons) -
  ntaxa(phyloseq_removed_singletons_plant) # 1,989 ASVs removed


# SEQUENCE COVERAGE ANALYSIS ----

# Calculate basic sequence coverage statistics
sorted <- sort(sample_sums(phyloseq_removed_singletons_plant))
min <- min(sample_sums(phyloseq_removed_singletons_plant)) # 0
max <- max(sample_sums(phyloseq_removed_singletons_plant)) # 237,153
mean <- mean(sample_sums(phyloseq_removed_singletons_plant)) # 19,051.85
median <- median(sample_sums(phyloseq_removed_singletons_plant)) # 16934

# Remove low-coverage samples (<100 reads)
phyloseq_removed_singletons_plant <- prune_samples(
  sample_sums(phyloseq_removed_singletons_plant) >= 100,
  phyloseq_removed_singletons_plant
)

# Filter by read abundance (minimum 20 reads per ASV)
filtered_phyloseq <- filter_taxa(
  phyloseq_removed_singletons_plant,
  function(x) {
    sum(x > 20) > (0.00 * length(x))
  },
  TRUE
)


# EXPORT PHYLOSEQ OBJECTS ----

# Save filtered and unfiltered phyloseq objects for further analysis
save(
  filtered_phyloseq,
  file = "data/output/phyloseq_objects/filtered_phyloseq_002.rda"
)
save(
  unfiltered_phyloseq,
  file = "data/output/phyloseq_objects/unfiltered_phyloseq_002.rda"
)
