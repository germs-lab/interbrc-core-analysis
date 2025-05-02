# BRC-Specific Data Subdirectories
This repository contains code and data for analysis specific to the Bioenergy Research Centers (BRCs). The directory is organized to separate BRC-specific resources from common code and data applicable across all BRCs.

## Directory Structure
```
├── brcs/
│   ├── cabbi/              # R scripts and data files specific to CABBI analysis
│   │   ├── R/              
│   │   └── data/           
│   ├── cbi/                # R scripts and data files specific to CBI analysis
│   │   ├── R/              
│   │   └── data/           
│   ├── glbrc/              # R scripts and data files specific to GLBRC analysis
│   │   ├── R/              
│   │   └── data/           
│   └── jbei/               # R scripts and data files specific to JBEI analysis
│       ├── R/             
│       └── data/           
├── data/                   # Common data and phyloseq object used for joint BRC analysis
│   ├── input/             
└── ├── output/          
        ├── fasta_files
        ├── phyloseq_objects
        ├── plots
        └── seq_processing_results
```
### `brcs/`
This directory contains subdirectories for each of the Bioenergy Research Centers (BRCs):

- **`cabbi/`**: Contains analysis code and data specific to the Center for Advanced Bioenergy and Bioproducts Innovation (CABBI).
- **`cbi/`**: Contains analysis code and data specific to the Consortium for Bioenergy Research (CBI).
- **`glbrc/`**: Contains analysis code and data specific to the Great Lakes Bioenergy Research Center (GLBRC).
- **`jbei/`**: Contains analysis code and data specific to the Joint BioEnergy Institute (JBEI).

Each subdirectory includes:
- **`R/`**: R scripts for preprocessing, analysis, and visualization specific to the respective BRC.
- **`data/`**: Data files required for the respective BRC's analysis.

## Purpose

This organization aims to:
1. Separate BRC-specific code and data from shared resources for better modularity and clarity.
2. Facilitate reproducibility by maintaining clear boundaries between individual and shared resources.

## Usage

1. Navigate to the relevant subdirectory in `brcs/` for BRC-specific analyses.
2. Use scripts and data from the *main* `data/` directory for any shared or common resources.
3. Follow the analysis pipeline documented in the project workflow files.

## Notes

- Always ensure that common resources in the `data/` directory are updated and version-controlled.
- BRC-specific resources should remain within their respective subdirectories to maintain clarity.
