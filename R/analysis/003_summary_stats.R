#########################################################
# SUMMARY STATISTICS
# Generate summary statistics for phyloseq objects
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-02-18
# Last modified: 2026-01-16
#########################################################

# DESCRIPTION:
# This script generates summary statistics for filtered phyloseq objects
# and exports the data as a CSV file for further analysis.

# SETUP AND DEPENDENCIES ----

source("R/utils/000_setup.R")
if (exists("unfiltered_phyloseq")) {
  remove(unfiltered_phyloseq)
}

# GENERATE SUMMARY STATISTICS ----

# Basic phyloseq summary
metagMisc::phyloseq_summary(filtered_phyloseq, more_stats = F, long = F)

# Calculate phylum distribution across samples
percent_phyla_samples <- metagMisc::phyloseq_ntaxa_by_tax(
  filtered_phyloseq,
  TaxRank = "phylum",
  relative = F,
  add_meta_data = F
) |>
  as.data.frame() |>
  mutate(sum = sum(N.OTU)) |>
  group_by(phylum) |>
  summarise(occurance_in_samples = n())


# DATA EXPORT ----

# Convert phyloseq object to dataframe and export as CSV
metagMisc::phyloseq_to_df(
  filtered_phyloseq,
  addtax = T,
  addtot = F,
  addmaxrank = F,
  sorting = "abundance"
) %>%
  # rename(ASV = OTU) %>%
  write.csv(
    .,
    file.path("data/output/phyloseq_objects/filtered_phyloseq_df_003.csv")
  )
