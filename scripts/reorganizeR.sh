#!/bin/bash

# Create new subdirectories in R/
mkdir -p R/functions R/analysis R/references R/utils

# Move function scripts
mv R/identify_core.R R/functions/
mv R/misc_functions.R R/functions/
mv R/model_fit.R R/functions/

# Move analysis scripts
mv R/001_quality_control.R R/analysis/
mv R/002_phyloseq_import.R R/analysis/
mv R/003_summary_stats.R R/analysis/
mv R/004_core_analysis.R R/analysis/
mv R/005_core_physeqs_fasta.R R/analysis/
mv R/006_core_ordinations.R R/analysis/

# Move reference scripts
mv R/REF_corn_dada2-phyloseq.R R/references/
mv R/REF_example_analysis.R R/references/
mv R/REF_selected_ASVs_high-low.R R/references/

# Move utility scripts
mv R/000_setup.R R/utils/

# Move DADA2 template (if it's a reference)
mv R/dada2_template.R R/references/

# Print confirmation
echo "R/ directory reorganization complete!"
echo "New structure:"
tree R/
