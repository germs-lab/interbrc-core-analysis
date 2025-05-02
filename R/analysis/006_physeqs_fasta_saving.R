## Inter-BRC Core Analysis
## Core & non_core Phyloseqs and FASTA files
## By Bolívar Aponte Rolón

# Setup
source("R/utils/000_setup.R")
if (exists("phyloseq")) remove(phyloseq)

#################################################
### Core and Non_Core objects and fasta files ###
#################################################
### extract_core() ###
#######################

# Subset by "core" and "non-core" names using "core_summary_lists" object

## Name strings to subset
core_asv_strings <- core_summary_lists[[4]] %>%
  dplyr::filter(., .$fill == "core") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

non_core_asv_strings <- core_summary_lists[[4]] %>%
  dplyr::filter(., .$fill == "no") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

# "core and "non_core" communities as matrices and phyloseq objects
## ASV (OTU) Matrices

## Remove samples where row sum 0 through extract_matrix()
core_asv_matrix <- extract_matrix(filtered_phyloseq, .vec = core_asv_strings) # Useful for manual ordinations
non_core_asv_matrix <- extract_matrix(
  filtered_phyloseq,
  .vec = non_core_asv_strings
)

## Sample names and data
core_sample_strings <- rownames(core_asv_matrix)
non_core_sample_strings <- rownames(non_core_asv_matrix)


#################################################
## Core dimensions is 50 ASVs out of a total of 23473 non-core ASVs in 1813 samples
## When non-zero sums samples are gathered we en up with 50 ASVs present in 1733 samples
## Core ASVs are present in 1733 / 1813 samples (96%)
#################################################

## Phyloseqs
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

## Checking that core is filtered out
rownames(core_brc_phyloseq@otu_table) %in% non_core_asv_strings # Should return FALSE for all 50 ASVs

# Save phyloseqs
save(
  core_brc_phyloseq,
  file = "data/output/phyloseq_objects/core_brc_phyloseq.rda"
)
save(
  non_core_brc_phyloseq,
  file = "data/output/phyloseq_objects/non_core_brc_phyloseq.rda"
)

# Creating Core and non_core Fasta files
## Subsetting fasta based on `subset.fasta` function from
## https://github.com/GuillemSalazar/FastaUtils/blob/master/R/FastaUtils_functions.R

## Core
subset_fasta(
  file = "data/output/fasta_files/rep_asv_seqs.fasta",
  subset = core_asv_strings,
  out = "data/output/fasta_files/core_asv_seqs.fasta"
)

## Non-Core
subset_fasta(
  file = "data/output/fasta_files/rep_asv_seqs.fasta",
  subset = non_core_asv_strings,
  out = "data/output/fasta_files/non_core_asv_seqs.fasta"
)


##############################################
### Microbiome Core Selection by Threshold ###
##############################################
physeq_high_occ <- core_asvs_threshold$physeq_high_occ
physeq_low_occ <- core_asvs_threshold$physeq_low_occ

high_occ_matrix <- physeq_high_occ@otu_table %>%
  t() %>% # Samples as rows
  as.data.frame() %>%
  .[rowSums(.) > 0, ] %>% # Keep only samples with a non-zero sum. Not all samples have the "core".
  as.matrix()


low_occ_matrix <- physeq_low_occ@otu_table %>%
  t() %>% # Samples as rows
  as.data.frame() %>%
  .[rowSums(.) > 0, ] %>% # Keep only samples with a non-zero sum. Not all samples have the "core".
  as.matrix()

## High occupancy
subset_fasta(
  file = "data/output/fasta_files/rep_asv_seqs.fasta",
  subset = colnames(high_occ_matrix),
  out = "data/output/fasta_files/high_occ_asv_seqs.fasta"
)

## Low occupancy
subset_fasta(
  file = "data/output/fasta_files/rep_asv_seqs.fasta",
  subset = colnames(low_occ_matrix),
  out = "data/output/fasta_files/low_occ_asv_seqs.fasta"
)

################################
### Save as list of matrices ###
################################
asv_matrices <- list(
  core = core_asv_matrix,
  non_core = non_core_asv_matrix,
  high_occ = high_occ_matrix,
  low_occ = low_occ_matrix
)

save(asv_matrices, file = "data/output/asv_matrices.rda")
