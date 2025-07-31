#########################################################
# COMMUNITY ORDINATION ANALYSIS
# NMDS and PCoA ordinations and dbRDA analysis and ordinations
# of core and non-core communities from extract_core()
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-07-10
#########################################################

# This is a version of the 007_ordinations.R focuse on  performing an NMDS on the full dataset
# using a Docker/Singularity container in Nova HPC at Iowa State


# Print the current library paths to verify
print(.libPaths())

# List of packages to load
packages <- c(
  'conflicted',
  'phyloseq',
  'vegan',
  'tidyverse',
  'minpack.lm',
  'Hmisc',
  'stats4',
  'BRCore',
  'furrr',
  'parallelly',
  'doParallel',
  'future',
  'here'
)


# Load packages using lapply
invisible(lapply(packages, library, character.only = TRUE))

cat("Session info")
print(sessionInfo())

cat("Loading data")
load(here::here("data/output/phyloseq_objects/filtered_phyloseq.rda"))
load(here::here("data/output/asv_matrices.rda"))
source(here::here("R/functions/brc_pcoa.R"))
cat("Data ready")



#--------------------------------------------------------
# ORDINATIONS ANALYSIS: FULL COMMUNITY
#--------------------------------------------------------
# workers <- parallelly::availableWorkers()
# cl <- parallelly::makeClusterPSOCK(length(workers) - 1L, autoStop = TRUE)
# doParallel::registerDoParallel(cl)
# plan(multicore, workers = parallel::detectCores() - 1)

#--------------------------------------------------------
# DATA TRANSFORMATION
#--------------------------------------------------------
# Transform community matrices using Hellinger transformation
cat("Start data transformations")
set.seed(54645)

hell_matrices <- purrr::map(
  asv_matrices,
  ~ {
    decostand(t(.x), method = "hellinger", MARGIN = 1)
  }
)

# Calculate distance matrices in parallel
plan(multicore, workers = parallel::detectCores() - 1)

distance_matrices <- hell_matrices %>%
  future_map(
    ~ vegdist(
      t(.x),
      method = "bray",
      upper = FALSE,
      binary = FALSE,
      na.rm = TRUE
    ),
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  ) %>%
  purrr::set_names(names(hell_matrices))

plan(sequential) # Close parallel processing

save(distance_matrices, file = here::here("data/output/distance_matrices.rda"))

#----------------------
# PCoA
#----------------------
cat("Starting PCoA calculations")

all_brc_pcoa <- brc_pcoa(
  distance_matrices$full_asv_matrix,
  filtered_phyloseq
)

cat("Finished PCoA")
save(all_brc_pcoa, file = here::here("data/output/all_brc_pcoa.rda"))
cat("Saved PCoA")

#----------------------
# NMDS
#----------------------

cat("Starting NMDS calculations")

all_brc_nmds <- BRCore::brc_nmds(
  asv_matrix = asv_matrices$full_asv_matrix,
  physeq = filtered_phyloseq,
  ncores = parallel::detectCores() - 1,
  k = 2,
  trymax = 100
)

cat("Finished NMDS")
save(all_brc_nmds, file = here::here("data/output/all_brc_nmds.rda"))
cat("Process complete")
print