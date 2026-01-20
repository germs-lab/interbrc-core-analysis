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
# Last modified: 2026-01-16
# ==============================================================================

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
# Load the packages
# Installed via /scripts/package_install.R
invisible(lapply(
    c(
        "conflicted",
        "styler",
        "phyloseq",
        "vegan",
        "tidyverse",
        "minpack.lm",
        "Hmisc",
        "stats4",
        "metagMisc",
        "BRCore",
        "furrr",
        "parallelly",
        "doParallel",
        "future"
    ),
    library,
    character.only = TRUE
))


load(
    here::here("data/output/phyloseq_objects/filtered_phyloseq.rda")
)

# Solve known conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("survival", "cluster")


#--------------------------------------------------------
# CORE EXTRACTION USING IDENTIFY_CORE()
#--------------------------------------------------------
# Ensure minimum sample quality
filtered_phyloseq <- prune_samples(
    sample_sums(filtered_phyloseq) >= 100,
    filtered_phyloseq
)

# Perform multiple rarefaction
interbrc_rarefied <- multi_rarefy(
    physeq = filtered_phyloseq,
    depth_level = 7000,
    num_iter = 50,
    threads = 2,
    set_seed = 7642
)

# Update phyloseq object with rarefied data
new_interbrc_rarefied <- update_otu_table(
    physeq = filtered_phyloseq,
    otu_rare = interbrc_rarefied
)


braycore_summary <- BRCore::identify_core(
    new_interbrc_rarefied,
    priority_var = "site",
    increase_value = 0.02,
    abundance_weight = 0,
    seed = 7895
)


# Save results to avoid recomputation
save(
    braycore_summary,
    file = here::here("data/output/braycore_summary_batch.rda")
)
