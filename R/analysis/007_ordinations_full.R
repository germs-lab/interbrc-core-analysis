#########################################################
# COMMUNITY ORDINATION ANALYSIS WITH FILTERING
# NMDS and PCoA ordinations with zero-filtering optimization
#
# This pipeline is focused on calculating PCoA and NMDS for
# the full ASV community. The distance matrices calculated are
# for thresholds (subsets) of the whole community, not core and non-core
# communities from extract_core().
#
# This is a version of the 007_ordinations.R focused on performing an NMDS on the full dataset
# using a Docker/Singularity container in Nova HPC at Iowa State
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-08-01
# Last modified: 2026-01-16
#########################################################

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

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

cat("Session info\n")
print(sessionInfo())

cat("Loading data\n")
load(here::here("data/output/phyloseq_objects/filtered_phyloseq.rda"))
#load(here::here("data/output/asv_matrices.rda"))
source(here::here("R/functions/parallel_helpers.R"))
source(here::here("R/functions/brc_filter_phyloseq.R"))
source(here::here("R/functions/brc_pcoa.R"))
cat("Data ready\n")


# FILTERING PARAMETER OPTIMIZATION ----

# Define thresholds for different filtering parameters
min_sample_reads_thresholds <- c(500, 1000, 5000)
min_asv_reads_thresholds <- c(20, 50, 100)

filtered_results <- list()

cat("Starting filtering optimization\n")

# Loop through sample read and ASV thresholds
for (sample_reads in min_sample_reads_thresholds) {
    for (asv_reads in min_asv_reads_thresholds) {
        # Create unique identifier for this combination
        result_name <- paste0(
            "sample_reads_",
            sample_reads,
            "_asv_reads_",
            asv_reads
        )

        # Apply phyloseq filtering
        tryCatch(
            {
                highly_filtered <- filter_phyloseq(
                    filtered_phyloseq,
                    min_sample_reads = sample_reads,
                    min_asv_reads = asv_reads
                )

                # Store results with dimensions info
                filtered_results[[result_name]] <- list(
                    phyloseq = highly_filtered,
                    dimensions = c(
                        nsamples(highly_filtered),
                        ntaxa(highly_filtered)
                    ),
                    sample_threshold = sample_reads,
                    asv_threshold = asv_reads
                )

                cat(
                    "  ",
                    result_name,
                    ":",
                    nsamples(highly_filtered),
                    "samples,",
                    ntaxa(highly_filtered),
                    "ASVs\n"
                )
            },
            error = function(e) {
                cat("  Error with", result_name, ":", e$message, "\n")
            }
        )
    }
}


# DATA TRANSFORMATION FOR FILTERED MATRICES ----

cat("Start data transformations for filtered matrices\n")
set.seed(54645)

# Create ASV matrices from filtered phyloseq objects
filtered_asv_matrices <- list()
for (name in names(filtered_results)) {
    phyloseq_obj <- filtered_results[[name]]$phyloseq

    # Skip if no data
    if (nsamples(phyloseq_obj) == 0 || ntaxa(phyloseq_obj) == 0) {
        next
    }

    # Extract ASV matrix
    filtered_asv_matrices[[name]] <- as.matrix(as.data.frame(t(otu_table(
        phyloseq_obj,
        taxa_are_rows = TRUE
    ))))
}

save(
    filtered_asv_matrices,
    file = here::here("data/output/filtered_asv_matrices.rda")
)

# Transform community matrices using Hellinger transformation
hell_matrices_filtered <- purrr::map(
    filtered_asv_matrices,
    ~ {
        decostand(t(.x), method = "hellinger", MARGIN = 1)
    }
)

save(
    hell_matrices_filtered,
    file = here::here("data/output/hell_matrices_filtered.rda")
)

# Clean up
remove("highly_filtered", "phyloseq_obj")

# Calculate distance matrices in parallel
plan(multicore, workers = parallel::detectCores() - 1)

distance_matrices_filtered <- hell_matrices_filtered %>%
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
    purrr::set_names(names(hell_matrices_filtered))

plan(sequential) # Close parallel processing

save(
    distance_matrices_filtered,
    file = here::here("data/output/distance_matrices_filtered.rda")
)

# Clean up
remove("hell_matrices_filtered", "filtered_phyloseq")


# PCOA ANALYSIS FOR FILTERED MATRICES ----

cat("Starting PCoA calculations for filtered matrices\n")
pcoa_results_filtered <- list()

for (name in names(distance_matrices_filtered)) {
    cat("Processing PCoA for:", name, "\n")

    # Calculate memory estimate
    matrix_size <- attr(distance_matrices_filtered[[name]], "Size")
    memory_estimate <- (matrix_size^2 * 8) / (1024^3)

    cat(
        "  Matrix size:",
        matrix_size,
        "samples, Memory estimate:",
        round(memory_estimate, 3),
        "GB\n"
    )

    # Only proceed if memory estimate is reasonable
    if (memory_estimate < 20) {
        tryCatch(
            {
                pcoa_result <- brc_pcoa(
                    distance_matrices_filtered[[name]],
                    filtered_results[[name]]$phyloseq
                )

                pcoa_results_filtered[[name]] <- pcoa_result
                cat("  PCoA completed successfully!\n")
            },
            error = function(e) {
                cat("  Error with PCoA:", e$message, "\n")
                pcoa_results_filtered[[name]] <- list(error = e$message)
            }
        )
    } else {
        cat("  Skipping - memory estimate too large\n")
        pcoa_results_filtered[[name]] <- list(
            error = "Memory estimate too large"
        )
    }
}

save(
    pcoa_results_filtered,
    file = here::here("data/output/pcoa_results_filtered.rda")
)

cat("PCoA calculations complete")

# NMDS ANALYSIS FOR OPTIMAL FILTERED MATRIX ----

# Select the best filtered matrix (balance between data retention and computational feasibility)
best_filter <- NULL
best_score <- 0

for (name in names(filtered_results)) {
    if (
        name %in%
            names(pcoa_results_filtered) &&
            !("error" %in% names(pcoa_results_filtered[[name]]))
    ) {
        # Score based on number of samples and ASVs
        score <- filtered_results[[name]]$dimensions[1] *
            filtered_results[[name]]$dimensions[2]
        if (score > best_score) {
            best_score <- score
            best_filter <- name
        }
    }
}

if (!is.null(best_filter)) {
    cat("Starting NMDS calculations for best filter:", best_filter, "\n")

    nmds_result_filtered <- BRCore::brc_nmds(
        asv_matrix = filtered_asv_matrices[[best_filter]],
        physeq = filtered_results[[best_filter]]$phyloseq,
        ncores = get_available_cores() - 1,
        k = 2,
        trymax = 100
    )

    save(
        nmds_result_filtered,
        file = here::here("data/output/nmds_result_filtered.rda")
    )
    cat("Finished NMDS for filtered data\n")
} else {
    cat("No suitable filtered matrix found for NMDS\n")
}


# ORIGINAL FULL DATASET ANALYSIS (if computationally feasible) ----

# Try original analysis if not too large
# Load original distance matrices if not already loaded
if (!exists("distance_matrices")) {
    load(here::here("data/output/distance_matrices.rda")) # loads 'distance_matrices'
}
if ("full_asv_matrix" %in% names(distance_matrices)) {
    matrix_size_full <- attr(distance_matrices$full_asv_matrix, "Size")
    memory_estimate_full <- (matrix_size_full^2 * 8) / (1024^3)

    cat(
        "Full matrix size:",
        matrix_size_full,
        "samples, Memory estimate:",
        round(memory_estimate_full, 2),
        "GB\n"
    )

    if (memory_estimate_full < 20) {
        # Higher threshold for full analysis
        cat("Starting PCoA calculations for full dataset\n")

        tryCatch(
            {
                all_brc_pcoa <- brc_pcoa(
                    distance_matrices$full_asv_matrix,
                    filtered_phyloseq
                )

                save(
                    all_brc_pcoa,
                    file = here::here("data/output/all_brc_pcoa.rda")
                )
                cat("Finished PCoA for full dataset\n")
            },
            error = function(e) {
                cat("Error with full PCoA:", e$message, "\n")
            }
        )

        cat("Starting NMDS calculations for full dataset\n")

        tryCatch(
            {
                all_brc_nmds <- BRCore::brc_nmds(
                    asv_matrix = asv_matrices$full_asv_matrix,
                    physeq = filtered_phyloseq,
                    ncores = get_available_cores() - 1,
                    k = 2,
                    trymax = 100
                )

                save(
                    all_brc_nmds,
                    file = here::here("data/output/all_brc_nmds.rda")
                )
                cat("Finished NMDS for full dataset\n")
            },
            error = function(e) {
                cat("Error with full NMDS:", e$message, "\n")
            }
        )
    } else {
        cat(
            "Full dataset too large for analysis - using filtered results only\n"
        )
    }
}

cat("Process complete\n")
