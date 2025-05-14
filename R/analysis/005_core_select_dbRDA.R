#########################################################
# CORE SELECTION USING dbRDA
# Distance-based redundancy analysis for BRC-specific core microbiome
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-05-08
#########################################################

# DESCRIPTION:
# This script performs distance-based redundancy analysis (dbRDA) to analyze
# the relationship between environmental variables and core/non-core microbial
# communities for a specific BRC.

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
source("R/utils/000_setup.R")
if (exists("phyloseq")) remove(phyloseq)

#--------------------------------------------------------
# BRC SELECTION AND PARAMETERS
#--------------------------------------------------------
# Define BRC of interest and output naming parameters
BRC <- "jbei"
CORE <- "core"
NOCORE <- "nocore"

# Filter phyloseq object for the selected BRC
physeq <- subset_samples(filtered_phyloseq, brc == BRC)

#--------------------------------------------------------
# CORE MICROBIOME EXTRACTION
#--------------------------------------------------------
# Extract core microbiome using Bray-Curtis dissimilarity
braycore_summary <- extract_core(
  physeq,
  Var = "site",
  method = "increase",
  increase_value = 2
)

# Subset by "core" and "non-core" taxa
core <- subset_physeq(braycore_summary, physeq, .var = "otu", type = "core")
nocore <- subset_physeq(braycore_summary, physeq, .var = "otu", type = "no")

#--------------------------------------------------------
# SAVE BRC-SPECIFIC PHYLOSEQ OBJECTS
#--------------------------------------------------------
# Helper function to create and name phyloseq objects
name_object <- function(brc, physeq, type) {
  obj_name <- stringr::str_glue("{brc}_{type}_phyloseq")
  assign(obj_name, value = physeq, envir = .GlobalEnv)
}

# Save core and non-core phyloseq objects
name_object(brc = BRC, physeq = physeq, type = "core") %>%
  saveRDS(
    .,
    file = stringr::str_glue(
      "data/output/phyloseq_objects/{BRC}_{CORE}_phyloseq.rda"
    )
  )

name_object(brc = BRC, physeq = physeq, type = "nocore") %>%
  saveRDS(
    .,
    file = stringr::str_glue(
      "data/output/phyloseq_objects/{BRC}_{NOCORE}_phyloseq.rda"
    )
  )

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
high_occ <- subset_samples(prevalence_core$physeq_high_occ, brc == BRC)
low_occ <- subset_samples(prevalence_core$physeq_low_occ, brc == BRC)
