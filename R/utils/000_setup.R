# Package and Environment setup
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
