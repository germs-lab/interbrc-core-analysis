#!/bin/bash

# Function to generate BRC report
brc_report() {
    local BRC="$1"
    
    # Render the Quarto document
quarto render brc_report.qmd \
        --execute-param BRC=$BRC \
        --to html \
        --output "${BRC}.html"
    
}

# Usage example
# brc_report "my_brc_analysis"