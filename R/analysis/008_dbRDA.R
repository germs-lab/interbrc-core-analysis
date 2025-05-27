#########################################################
# DISTANCE-BASED REDUNDANCY ANALYSIS
# dbRDA for JBEI core and non-core communities
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-05-08
#########################################################

# DESCRIPTION:
# This script performs distance-based redundancy analysis (dbRDA) on core and non-core
# microbial communities specifically for the JBEI dataset. It analyzes the influence
# of environmental variables on community structure.

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
source("R/utils/000_setup.R")
if (exists("phyloseq")) remove(phyloseq)

#--------------------------------------------------------
# Subset by "core" and "non-core" taxa
#--------------------------------------------------------

core <- subset_physeq(braycore_summary, physeq, .var = "otu", type = "core")
nocore <- subset_physeq(braycore_summary, physeq, .var = "otu", type = "no")

#--------------------------------------------------------
# CORE COMMUNITY dbRDA ANALYSIS
#--------------------------------------------------------
# Transform data using Hellinger transformation
hell_matrix <- decostand(t(core$asv_matrix), method = "hellinger", MARGIN = 1)

# Prepare environmental data
dbrda_traits <- braycore_summary$sample_metadata %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  select(., c(drought, treatment, block, harvest, sample_id)) %>%
  dplyr::filter(rownames(.) %in% colnames(hell_matrix))

# Build dbRDA models with different constraints
dbrda_00_core <- dbrda(
  t(hell_matrix) ~ 1,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_01_core <- dbrda(
  t(hell_matrix) ~ .,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_02_core <- dbrda(
  t(hell_matrix) ~ drought + block + harvest,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

#--------------------------------------------------------
# NON-CORE COMMUNITY dbRDA ANALYSIS
#--------------------------------------------------------
# Transform non-core data using Hellinger transformation
hell_matrix <- decostand(t(nocore$asv_matrix), method = "hellinger", MARGIN = 1)

# Prepare environmental data for non-core analysis
dbrda_traits <- braycore_summary$sample_metadata %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  select(., c(drought, treatment, block, harvest, sample_id)) %>%
  dplyr::filter(rownames(.) %in% colnames(hell_matrix))

# Build dbRDA models for non-core community
dbrda_00_nocore <- dbrda(
  t(hell_matrix) ~ 1,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_01_nocore <- dbrda(
  t(hell_matrix) ~ .,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_02_nocore <- dbrda(
  t(hell_matrix) ~ drought + block + harvest,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

#--------------------------------------------------------
# STATISTICAL TESTING
#--------------------------------------------------------
# Perform permutation tests on models
set.seed(123)
anova.cca(dbrda_02_nocore, by = "margin", permutations = 999, parallel = 8)
anova(dbrda_02_core, by = "axis")
anova(dbrda_02_core, by = "axis", perm.max = 500)

#--------------------------------------------------------
# VISUALIZATION
#--------------------------------------------------------
# Generate ordination plots for core and non-core communities
core_dbrda_gg <- brc_flex_ordi(
  dbrda_02_core,
  dbrda_traits,
  color_var = "treatment",
  sample_id_col = "sample_id"
) %>%
  +ggtitle(str_glue("{BRC}: Bray-Curtis Core"))

nocore_dbrda_gg <- brc_flex_ordi(
  dbrda_02_nocore,
  dbrda_traits,
  color_var = "treatment",
  sample_id_col = "sample_id"
) %>%
  +ggtitle(str_glue("{BRC}: Bray-Curtis Non-Core "))

# Save plots to output directory
save_gg <- list(core_dbrda_gg, nocore_dbrda_gg) %>%
  purrr::set_names(
    stringr::str_glue("{BRC}_{c(CORE, NOCORE)}_dbrda")
  )

brc_ggsave(save_gg, "data/output/plots/")

#--------------------------------------------------------
# THRESHOLD-BASED ANALYSIS
#--------------------------------------------------------
# Extract high and low occupancy communities from threshold-based core analysis
high_occ <- subset_samples(core_asvs_threshold$physeq_high_occ, brc == BRC)
low_occ <- subset_samples(core_asvs_threshold$physeq_low_occ, brc == BRC)
