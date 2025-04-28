# Package and Environment setup
if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(
  conflicted,
  styler,
  phyloseq,
  vegan,
  tidyverse,
  janitor,
  ggsci,
  # ggpubr,
  # hillR,
  # microViz,
  minpack.lm,
  Hmisc,
  stats4,
  metagMisc,
  BRCore,
  furrr
)

## List files and source each
list.files(here::here("R/functions"), pattern = "brc_", full.names = TRUE) %>%
  purrr::map(source)

# Objects
list.files(
  here::here("data/output"),
  full.names = TRUE,
  recursive = TRUE,
  pattern = "\\.rda$"
) %>%
  purrr::walk(., load, envir = .GlobalEnv)

# Solve known conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("survival", "cluster")
