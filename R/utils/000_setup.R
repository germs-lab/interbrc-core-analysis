# Package and Environment setup

<<<<<<< Updated upstream
# Packages
library(phyloseq) # Bioconductor v1.50.0
library(vegan) # CRAN v2.6-10
library(tidyverse) # Posit RPSM v2.0.0
library(ggsci) # [standard::NA/NA] v3.2.0
library(ggpubr) # [standard::NA/NA] v0.6.0
# library(hillR) # [standard::NA/NA] v0.5.2
# library(microViz) # [::NA/NA] v0.12.6
# library(magrittr) # Posit RPSM v2.0.3
# library(lazyeval) # Posit RPSM v0.2.2 # Deprecated. Change to rlang
# library(minpack.lm) # [standard::NA/NA] v1.2-4
# library(Hmisc) # [standard::NA/NA] v5.2-2
# library(stats4) # local install v4.4.1
library(metagMisc) # [github::vmikk/metagMisc] v0.5.0

# Scripts
## List files and source each
list.files("R/functions/", full.names = TRUE) %>%
  purrr::map(source)

# Objects
list.files("data/output",
  full.names = TRUE,
  recursive = TRUE,
  pattern = "\\.rda$"
) %>%
  purrr::walk(., load, envir = .GlobalEnv)
=======
if (!require("pacman"))
    install.packages("pacman")

pacman::p_load(
    conflicted,
    styler,
    phyloseq,
    vegan,
    tidyverse,
    ggsci,
    ggpubr,
    hillR,
    microViz,
    magrittr,
    lazyeval,
    minpack.lm,
    Hmisc,
    stats4,
    metagMisc
)

source("R/functions/extract_core.R")
source("R/functions/model_fit.R")
source("R/functions/misc_functions.R")
source("R/functions/multi_rarefy.R")

# Solve known conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("survival", "cluster")

>>>>>>> Stashed changes
