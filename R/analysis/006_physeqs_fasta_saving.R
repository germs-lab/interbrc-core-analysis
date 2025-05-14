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
bc_core_asv_ids <- core_summary_lists[[4]] %>%
  dplyr::filter(., .$fill == "core") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

bc_noncore_asv_ids <- core_summary_lists[[4]] %>%
  dplyr::filter(., .$fill == "no") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

# Extract ASV matrices for core and non-core communities
bc_core_matrix <- extract_matrix(
  filtered_phyloseq,
  .vec = bc_core_ids
)
bc_noncore_matrix <- extract_matrix(
  filtered_phyloseq,
  .vec = bc_noncore_ids
)

# Get sample names for core and non-core communities
bc_core_sample_ids <- rownames(bc_core_matrix)
bc_noncore_sample_ids <- rownames(bc_noncore_matrix)

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
braycurt_core <- prune_samples(
  sort(sample_names(filtered_phyloseq)) %in% sort(bc_core_sample_ids),
  filtered_phyloseq
) %>%
  prune_taxa(rownames(otu_table(.)) %in% bc_core_ids, .)

braycurt_noncore <- prune_samples(
  sort(sample_names(filtered_phyloseq)) %in%
    sort(bc_noncore_sample_ids),
  filtered_phyloseq
) %>%
  prune_taxa(rownames(otu_table(.)) %in% bc_noncore_asv_ids, .)

# Verify proper filtering (should return all FALSE)
rownames(braycurt_core@otu_table) %in% bc_noncore_sample_ids

#--------------------------------------------------------
# SAVE PHYLOSEQ OBJECTS
#--------------------------------------------------------
# Save core and non-core phyloseq objects
save(
  braycurt_core,
  file = "data/output/phyloseq_objects/braycurt_core.rda"
)
save(
  braycurt_noncore,
  file = "data/output/phyloseq_objects/braycurt_noncore.rda"
)

#--------------------------------------------------------
# GENERATE FASTA FILES
#--------------------------------------------------------
# Create FASTA files for core community
subset_fasta(
  file = "data/input/rep_asv_seqs.fasta",
  subset = bc_core_ids,
  out = "data/output/fasta_files/braycurt_core_seqs.fasta"
)

# Create FASTA files for non-core community
subset_fasta(
  file = "data/input/rep_asv_seqs.fasta",
  subset = bc_noncore_ids,
  out = "data/output/fasta_files/braycurt_noncore_seqs.fasta"
)

#--------------------------------------------------------
# THRESHOLD-BASED CORE SELECTION
#--------------------------------------------------------
# Extract high and low occupancy ASVs from threshold-based approach
physeq_high_occ <- prevalence_core$physeq_high_occ
physeq_low_occ <- prevalence_core$physeq_low_occ

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
  file = "data/input/rep_asv_seqs.fasta",
  subset = colnames(high_occ_matrix),
  out = "data/output/fasta_files/high_occ_seqs.fasta"
)

# Generate FASTA files for low occupancy ASVs
subset_fasta(
  file = "data/input/rep_asv_seqs.fasta",
  subset = colnames(low_occ_matrix),
  out = "data/output/fasta_files/low_occ_seqs.fasta"
)

#--------------------------------------------------------
# SAVE MATRICES FOR DOWNSTREAM ANALYSES
#--------------------------------------------------------
# Save all matrices as a single list object
asv_matrices <- list(
  bc_core = bc_core_matrix,
  bc_noncore = bc_noncore_matrix,
  high_occ = high_occ_matrix,
  low_occ = low_occ_matrix
)

save(asv_matrices, file = "data/output/asv_matrices.rda")
