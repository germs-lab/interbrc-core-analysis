#----------------------------------------------------
# Package installation for R environment
# Necessary for execution of 004_core_selection_HPC.R
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-02-29
#----------------------------------------------------

# Package and Environment setup
# First, install pak if not already installed
if (!requireNamespace(c("renv"), quietly = TRUE)) {
  install.packages(c("renv"))
}

# Use pak to install and load all required packages
# pak::pkg_install(c(
#   "conflicted",
#   "styler",
#   "phyloseq",
#   "vegan",
#   "tidyverse",
#   "minpack.lm",
#   "Hmisc",
#   "stats4",
#   "vmikk/metagMisc",
#   "germs-lab/BRCore@b391575",
#   "furrr",
#   "parallelly",
#   "doParallel",
#   "future"
#   # "future.batchtools",
#   # "batchtools" # Installed by future.batchtools. Being explicit here.
# ))

#options(renv.config.ppm.enabled = FALSE)
options(renv.config.pak.enabled = TRUE)
renv::install("pak@0.9.0")
renv::restore(rebuild = TRUE)
