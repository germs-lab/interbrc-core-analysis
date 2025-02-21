# Quality control assessment using FastQC and MultiQC
# Bolívar Aponte Rolón
# 2025-02-12

# This script is intent to assess quality of sequence reads using `fastqc` and `multiqc`
# CLI tools in an R workflow using `base::system2()` to call system tools.
# This could very well be a Shell script. This is just convenient to run from Rstudio/Positron.

### FastQC ###

# Location of executable and files
# Change accordingly

# Executables
fastqc <- "/usr/bin/fastqc"
multiqc <- "/usr/bin/multiqc"

system2(fastqc, args = "--version") # Make sure of version for posterity.
system2(multiqc, args = "--version")

# Directories and files
seqs_1 <- "../CABBI-Miscanthus-Fertilization/'Raw sequences'/*.fastq.gz" # In CLI 'Raw sequences' = ~/Raw\ sequences/..
seqs_2 <- "../CABBI-SABR2023_Spatial/'Raw sequences'/*.fastq.gz"
out_dir1 <- "data/output/qc_reports/cabbi_mxg"
out_dir2 <- "data/output/qc_reports/cabbi_sabr"
adapters <- "data/input/adapters.txt" # File with PCR primers or adapters

### QC reports
# FastQC
system2(fastqc, args = c("--adapters", adapters, seqs_1, "--outdir", out_dir1))
system2(fastqc, args = c("--adapters", adapters, seqs_2, "--outdir", out_dir2))

# MultiQC
system2(multiqc, args = c(out_dir1, "--ai-summary-full", "--outdir", out_dir1)) # Needs TOWER_ACCESS_TOKEN environment variable or change config.ai_provider
# If you run MultiQC without the appropriate key you will get a warning printed to the console, but report generation will otherwise proceed without the summary.
# MultiQC will not return an error exit code.
system2(multiqc, args = c(out_dir2, "--ai-summary-full", "--outdir", out_dir2))


# Go and read your reports!
