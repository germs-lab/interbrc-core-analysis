#########################################################
# COMMUNITY ORDINATION ANALYSIS
# NMDS and PCoA ordinations and dbRDA analysis and ordinations
# of core and non-core communities from extract_core()
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-07-10
#########################################################

# THis is a version of the 007_ordinations.R focuse on  performing an NMDS on the full dataset
# using a Docker/Singularity container in Nova HPC at Iowa State

# First, install pak if not already installed
#if (!requireNamespace(c("renv"), quietly = TRUE)) {
#  install.packages(c("renv"))
#}

#options(renv.config.pak.enabled = TRUE)
#renv::install("pak@0.9.0")
#renv::restore(rebuild = TRUE)

# Clear existing library paths and set only the renvcache path
.libPaths("/opt/renvcache/v5/linux-ubuntu-jammy/R-4.4/x86_64-pc-linux-gnu")

# Print the current library paths to verify
print(.libPaths())

# List of packages to load
packages <- c(
  #'conflicted',
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

#purrr::walk(packages, library, character.only = TRUE)

## List files and source each
# list.files(here::here("R/functions"), pattern = "brc_", full.names = TRUE) %>%
#   purrr::map(source)

cat("Loading data")
load(here::here("data/output/phyloseq_objects/filtered_phyloseq.rda"))
load(here::here("data/output/asv_matrices.rda"))
cat("Data ready")

cat("Session info")
print(sessionInfo())

#--------------------------------------------------------
# NMDS ANALYSIS: FULL COMMUNITY
#--------------------------------------------------------
# workers <- parallelly::availableWorkers()
# cl <- parallelly::makeClusterPSOCK(length(workers) - 1L, autoStop = TRUE)
# doParallel::registerDoParallel(cl)
# plan(multicore, workers = parallel::detectCores() - 1)
cat("Starting NMDS calculations")
all_brc_nmds <- BRCore::brc_nmds(
  asv_matrix = asv_matrices$full_asv_matrix,
  physeq = filtered_phyloseq,
  ncores = parallel::detectCores() - 1,
  k = 2,
  trymax = 100
)

# plan(sequential)
cat("Finished NMDS")
#save(all_brc_nmds, file = "data/output/all_brc_nmds.rda")
cat("Process complete")
print