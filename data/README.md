# Joint BRC Data Directory

This directory contains the joint Bioenergy Research Center (BRC) data used for collaborative analyses across all centers. It is organized to separate raw input files from processed outputs, ensuring clarity and reproducibility in the analysis pipeline.

## Directory Structure

```
data/
├── input/
│   ├── feature-table-modified.tsv       # Modified feature table for analysis
│   ├── feature-table-modified.tsv.zip  # Compressed version of the feature table
│   ├── metadata.tsv                     # Sample metadata file
│   └── taxonomy.tsv                     # Taxonomic annotations for ASVs
└── output/
    ├── asv_matrices.rda                # ASV abundance matrices
    ├── core_ext_nmds.rda               # NMDS results for core extended community
    ├── core_summary_lists.rda          # Summary lists for core microbial community
    ├── distance_matrices.rda           # Distance matrices for beta-diversity analyses
    ├── fasta_files/                    # FASTA files for sequence data
    ├── high_occ_nmds.rda               # NMDS results for high-occurrence ASVs
    ├── low_occ_nmds.rda                # NMDS results for low-occurrence ASVs
    ├── non_core_ext_nmds.rda           # NMDS results for non-core extended community
    ├── phyloseq_objects/               # Saved phyloseq objects for downstream analysis
    ├── plots/                          # Visualizations and plots generated during analysis
    └── seq_processing_results/         # Results from sequence processing
```

## Purpose

This directory serves as the central repository for:
- **Input Data**: Raw and pre-processed data files required for joint BRC analyses.
- **Output Data**: Processed results, intermediate files, and final outputs from the analysis pipeline.

## Contents

### `input/`

This folder contains the foundational files for the analysis:
- **`feature-table-modified.tsv`**: The main feature table containing ASV abundances.
- **`metadata.tsv`**: Metadata describing the samples in the dataset.
- **`taxonomy.tsv`**: Taxonomic assignments for each ASV.
- **`feature-table-modified.tsv.zip`**: A compressed version of the feature table for easier sharing.

> **Note**: These files should remain unchanged to maintain dataset integrity. Any modifications should be applied to copies saved in the `output/` directory.

### `output/`

This folder houses the results and intermediate files from the analysis:
- **RDA Files**: Saved R objects for efficient reloading of analysis results.
- **Subdirectories**:
  - **`fasta_files/`**: Contains sequence data in FASTA format.
  - **`phyloseq_objects/`**: Stores phyloseq objects for microbiome analysis.
  - **`plots/`**: Contains figures and visualizations.
  - **`seq_processing_results/`**: Includes outputs from sequence processing pipelines.

> **Note**: Unless explicitly labeled with a specific BRC name (e.g., `brc_jbei_*`), all files are inclusive of data from all BRCs and are intended for joint analyses (e.g., `core_*` or `brc_*`).

## Usage

1. **Input Files**:
   - Use the files in the `input/` directory as the starting point for analysis.
   - Ensure that these files are not overwritten or modified directly.

2. **Output Files**:
   - Save all processed and intermediate files in the `output/` directory.
   - Use the subdirectories to organize outputs based on file type or analysis stage.

## Notes

- The `input/` directory should only contain raw or pre-processed data that is consistent across all BRCs.
- The `output/` directory is dynamic and will grow as analyses progress; ensure proper organization and naming conventions.
- Shared plots, figures, and processed results can be utilized in the broader Inter-BRC analysis and reporting.
