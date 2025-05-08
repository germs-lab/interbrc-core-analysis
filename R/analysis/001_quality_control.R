#########################################################
# QUALITY CONTROL ANALYSIS
# Initial sequence quality assessment
# Project:  Inter-BRC-Core-Microbiome
#
# Author: Bolívar Aponte Rolón
# Date: 2025-02-12
#########################################################

# DESCRIPTION:
# This script assesses the quality of sequence reads using FastQC and MultiQC CLI tools
# through R's system2() function. This approach enables quality control processing
# directly from RStudio/Positron environments.

#--------------------------------------------------------
# SETUP AND CONFIGURATION
#--------------------------------------------------------
# Define paths to executables
fastqc <- "/usr/bin/fastqc"
multiqc <- "/usr/bin/multiqc"

# Verify tool versions for documentation
system2(fastqc, args = "--version")
system2(multiqc, args = "--version")

#--------------------------------------------------------
# INPUT FILE LOCATIONS
#--------------------------------------------------------
# Define paths to sequence files and output directories
seqs_1 <- "../CABBI-Miscanthus-Fertilization/'Raw sequences'/*.fastq.gz" # In CLI 'Raw sequences' = ~/Raw\ sequences/..
seqs_2 <- "../CABBI-SABR2023_Spatial/'Raw sequences'/*.fastq.gz"
out_dir1 <- "data/output/qc_reports/cabbi_mxg"
out_dir2 <- "data/output/qc_reports/cabbi_sabr"
adapters <- "data/input/adapters.txt" # File with PCR primers or adapters

#--------------------------------------------------------
# EXECUTE FASTQC ANALYSIS
#--------------------------------------------------------
# Process sequence files with FastQC
system2(fastqc, args = c("--adapters", adapters, seqs_1, "--outdir", out_dir1))
system2(fastqc, args = c("--adapters", adapters, seqs_2, "--outdir", out_dir2))

#--------------------------------------------------------
# GENERATE SUMMARY REPORTS WITH MULTIQC
#--------------------------------------------------------
# Aggregate FastQC results with MultiQC
# Note: Requires TOWER_ACCESS_TOKEN environment variable for AI summary
system2(multiqc, args = c(out_dir1, "--ai-summary-full", "--outdir", out_dir1))
system2(multiqc, args = c(out_dir2, "--ai-summary-full", "--outdir", out_dir2))

# Go and read your reports!
