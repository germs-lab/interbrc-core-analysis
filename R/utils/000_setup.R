# Package and Environment setup
# Bolívar Aponte Rolón

# Packages
library(phyloseq) # Bioconductor v1.50.0
library(vegan) # CRAN v2.6-10
library(tidyverse) # Posit RPSM v2.0.0
library(ggsci) # [standard::NA/NA] v3.2.0
library(ggpubr) # [standard::NA/NA] v0.6.0
# library(hillR) # [standard::NA/NA] v0.5.2
library(microViz) # [::NA/NA] v0.12.6
# library(magrittr) # Posit RPSM v2.0.3
# library(lazyeval) # Posit RPSM v0.2.2 # Deprecated. Change to rlang
# library(minpack.lm) # [standard::NA/NA] v1.2-4
# library(Hmisc) # [standard::NA/NA] v5.2-2
# library(stats4) # local install v4.4.1
library(metagMisc) # [github::vmikk/metagMisc] v0.5.0

# Scripts
## List files and source each
list.files("R/functions/", full.names = TRUE) %>% map(source)

# Objects
load("data/output/core_summary_lists.rda")
load("data/output/phyloseq_objects/filtered_phyloseq.rda")
load("data/output/phyloseq_objects/core_brc_phyloseq.rda")
load("data/output/phyloseq_objects/non_core_brc_phyloseq.rda")
load("data/output/core_asv_matrix.rda")
load("data/output/non_core_asv_matrix.rda")
load("data/output/phyloseq_objects/core_asvs_threshold.rda")
load("data/output/high_occ_matrix.rda")
load("data/output/low_occ_matrix.rda")
