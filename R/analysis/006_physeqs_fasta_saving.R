#########################################################
# PHYLOSEQ AND FASTA FILE GENERATION
# Export phyloseq objects and corresponding FASTA files
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-05-08
#########################################################

# DESCRIPTION:
# This script extracts core and non-core community components from phyloseq
# objects and exports them as separate phyloseq objects and FASTA files for
# downstream analyses.

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
source("R/utils/000_setup.R")
if (exists("phyloseq")) remove(phyloseq)

#--------------------------------------------------------
# EXTRACT CORE AND NON-CORE TAXA
#--------------------------------------------------------
# Get ASV identifiers for core and non-core taxa
core_asv_strings <- core_summary_lists[[4]] %>%
  dplyr::filter(., .$fill == "core") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

non_core_asv_strings <- core_summary_lists[[4]] %>%
  dplyr::filter(., .$fill == "no") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

# Extract ASV matrices for core and non-core communities
core_asv_matrix <- extract_matrix(filtered_phyloseq, .vec = core_asv_strings)
non_core_asv_matrix <- extract_matrix(
  filtered_phyloseq,
  .vec = non_core_asv_strings
)

# Get sample names for core and non-core communities
core_sample_strings <- rownames(core_asv_matrix)
non_core_sample_strings <- rownames(non_core_asv_matrix)

#--------------------------------------------------------
# COMMUNITY COMPOSITION SUMMARY
#--------------------------------------------------------
# Note on community composition:
# - Core dimensions: 50 ASVs out of a total of 23473 non-core ASVs in 1813 samples
# - When non-zero sums samples are gathered, 50 core ASVs are present in 1733 samples
# - Core ASVs are present in 1733/1813 samples (96%)

#--------------------------------------------------------
# CREATE PHYLOSEQ OBJECTS
#--------------------------------------------------------
# Create core and non-core phyloseq objects
core_brc_phyloseq <- prune_samples(
  sort(sample_names(filtered_phyloseq)) %in% sort(core_sample_strings),
  filtered_phyloseq
) %>%
  prune_taxa(rownames(.@otu_table) %in% core_asv_strings, .)

non_core_brc_phyloseq <- prune_samples(
  sort(sample_names(filtered_phyloseq)) %in% sort(non_core_sample_strings),
  filtered_phyloseq
) %>%
  prune_taxa(rownames(.@otu_table) %in% non_core_asv_strings, .)

# Verify proper filtering (should return all FALSE)
rownames(core_brc_phyloseq@otu_table) %in% non_core_asv_strings

#--------------------------------------------------------
# SAVE PHYLOSEQ OBJECTS
#--------------------------------------------------------
# Save core and non-core phyloseq objects
save(
  core_brc_phyloseq,
  file = "data/output/phyloseq_objects/core_brc_phyloseq.rda"
)
save(
  non_core_brc_phyloseq,
  file = "data/output/phyloseq_objects/non_core_brc_phyloseq.rda"
)

#--------------------------------------------------------
# GENERATE FASTA FILES
#--------------------------------------------------------
# Create FASTA files for core community
subset_fasta(
  file = "data/output/fasta_files/rep_asv_seqs.fasta",
  subset = core_asv_strings,
  out = "data/output/fasta_files/core_asv_seqs.fasta"
)

# Create FASTA files for non-core community
subset_fasta(
  file = "data/output/fasta_files/rep_asv_seqs.fasta",
  subset = non_core_asv_strings,
  out = "data/output/fasta_files/non_core_asv_seqs.fasta"
)

#--------------------------------------------------------
# THRESHOLD-BASED CORE SELECTION
#--------------------------------------------------------
# Extract high and low occupancy ASVs from threshold-based approach
physeq_high_occ <- core_asvs_threshold$physeq_high_occ
physeq_low_occ <- core_asvs_threshold$physeq_low_occ

# Create matrices for high and low occupancy ASVs
high_occ_matrix <- physeq_high_occ@otu_table %>%
  t() %>% # Samples as rows
  as.data.frame() %>%
  .[rowSums(.) > 0, ] %>% # Keep only samples with a non-zero sum
  as.matrix()

low_occ_matrix <- physeq_low_occ@otu_table %>%
  t() %>% # Samples as rows
  as.data.frame() %>%
  .[rowSums(.) > 0, ] %>% # Keep only samples with a non-zero sum
  as.matrix()

# Generate FASTA files for high occupancy ASVs
subset_fasta(
  file = "data/output/fasta_files/rep_asv_seqs.fasta",
  subset = colnames(high_occ_matrix),
  out = "data/output/fasta_files/high_occ_asv_seqs.fasta"
)

# Generate FASTA files for low occupancy ASVs
subset_fasta(
  file = "data/output/fasta_files/rep_asv_seqs.fasta",
  subset = colnames(low_occ_matrix),
  out = "data/output/fasta_files/low_occ_asv_seqs.fasta"
)

#--------------------------------------------------------
# SAVE MATRICES FOR DOWNSTREAM ANALYSES
#--------------------------------------------------------
# Save all matrices as a single list object
asv_matrices <- list(
  core = core_asv_matrix,
  non_core = non_core_asv_matrix,
  high_occ = high_occ_matrix,
  low_occ = low_occ_matrix
)

save(asv_matrices, file = "data/output/asv_matrices.rda")
