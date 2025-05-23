# ==============================================================================
# SLURM Job Submission Script for Microbiome Core Extraction
# ==============================================================================
#
# Description:
#   This SLURM script configures and executes microbiome core extraction analysis
#   on a high-performance computing cluster. It allocates, what I thought, were
#   appropriate compute resources for parallel processing. It sets up the necessary environment,
#   and runs the core selection analysis (004_core_selection_HPC.R) using the
#   'interbrc_env' environment.
#   See main SLURM configuration in /scripts/core_selection_HPC.sh'
#
# Resources:
#   - 1 compute node with 32 processor cores
#   - 256GB RAM
#   - 2 days walltime
#   - Output directed to the hpc directory with job ID in filenames
#
# Usage:
#   Submit interactively from R or VS Code session or via SLURM batch job.
#   The R and VS COde sessions need installation of packages prior to running the job.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-05-16
# ==============================================================================

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
# Package and Environment setup
if (!require("pacman")) {
    install.packages("pacman")
}

pacman::p_load(
    conflicted,
    styler,
    phyloseq,
    vegan,
    tidyverse,
    minpack.lm,
    Hmisc,
    stats4,
    metagMisc,
    BRCore,
    furrr,
    future,
    future.batchtools,
    batchtools # Installed by future.batchtools. Being explicit here.
)

source(
    "/work/adina/bolivar/interbrc-core-analysis/R/functions/extract_core_parallel.R"
)
load(
    "/work/adina/bolivar/interbrc-core-analysis/data/output/phyloseq_objects/filtered_phyloseq.rda"
)


# Solve known conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("survival", "cluster")


#--------------------------------------------------------
# PARALLEL PROCESSING SETUP FOR HPC
#--------------------------------------------------------

nCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # https://research.it.iastate.edu/guides/pronto/r/#using-the-parallel-library
print(paste("nCores", nCores))
myCluster <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(myCluster)

# plan(multisession, workers = 16)

#--------------------------------------------------------
# CORE EXTRACTION USING EXTRACT_CORE()
#--------------------------------------------------------

# Extract core microbiome across all sites (with minimum 2% increase in Bray-Curtis)
# Set .parallel = TRUE to use the future framework
braycore_summary <- extract_core_parallel(
    test_phyloseq,
    Var = "site",
    method = "increase",
    increase_value = 2,
    .parallel = TRUE,
    ncores = 2
)

# plan(sequential) # Close parallel processing

# Save results to avoid recomputation
save(
    braycore_summary,
    file = "/work/adina/bolivar/interbrc-core-analysis/data/output/braycore_summary_batch.rda"
)
