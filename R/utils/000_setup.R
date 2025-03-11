# Package and Environment setup
# Bolívar Aponte Rolón

library(phyloseq) # Bioconductor v1.50.0
library(vegan) # CRAN v2.6-10
library(tidyverse) # Posit RPSM v2.0.0
library(ggsci) # [standard::NA/NA] v3.2.0
library(ggpubr) # [standard::NA/NA] v0.6.0
library(hillR) # [standard::NA/NA] v0.5.2
library(microViz) # [::NA/NA] v0.12.6
library(magrittr) # Posit RPSM v2.0.3
library(lazyeval) # Posit RPSM v0.2.2 # Deprecated. Change to rlang
library(minpack.lm) # [standard::NA/NA] v1.2-4
library(Hmisc) # [standard::NA/NA] v5.2-2
library(stats4) # local install v4.4.1
library(metagMisc) # [github::vmikk/metagMisc] v0.5.0


source("R/functions/extract_core.R")
source("R/functions/model_fit.R")
source("R/functions/misc_functions.R")
