#########################################################
# CORE SELECTION PER BRC
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-05-08
#########################################################

# DESCRIPTION:
# This script identifies the core microbiome across samples for individual
# BRCs using multiple approaches:
# 1. extract_core() method based on Bray-Curtis dissimilarity
# 2. Threshold-based approach

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
source("R/utils/000_setup.R")
if (exists("phyloseq")) remove(phyloseq)

#--------------------------------------------------------
# BRC SELECTION AND PARAMETERS
#--------------------------------------------------------
# Define BRC of interest and output naming parameters
BRC <- "jbei"
CORE <- "core"
NOCORE <- "nocore"

# Filter phyloseq object for the selected BRC
physeq <- subset_samples(filtered_phyloseq, brc == BRC)

#--------------------------------------------------------
# CORE MICROBIOME EXTRACTION
#--------------------------------------------------------
# Extract core microbiome using Bray-Curtis dissimilarity
braycore_summary <- extract_core(
  physeq,
  Var = "site",
  method = "increase",
  increase_value = 2
)

# Subset by "core" and "non-core" taxa
core <- subset_physeq(braycore_summary, physeq, .var = "otu", type = "core")
nocore <- subset_physeq(braycore_summary, physeq, .var = "otu", type = "no")

#--------------------------------------------------------
# SAVE BRC-SPECIFIC PHYLOSEQ OBJECTS
#--------------------------------------------------------
# Helper function to create and name phyloseq objects
name_object <- function(brc, physeq, type) {
  obj_name <- stringr::str_glue("{brc}_{type}_phyloseq")
  assign(obj_name, value = physeq, envir = .GlobalEnv)
}

# Save core and non-core phyloseq objects
name_object(brc = BRC, physeq = physeq, type = "core") %>%
  saveRDS(
    .,
    file = stringr::str_glue(
      "data/output/phyloseq_objects/{BRC}_{CORE}_phyloseq.rda"
    )
  )

name_object(brc = BRC, physeq = physeq, type = "nocore") %>%
  saveRDS(
    .,
    file = stringr::str_glue(
      "data/output/phyloseq_objects/{BRC}_{NOCORE}_phyloseq.rda"
    )
  )
