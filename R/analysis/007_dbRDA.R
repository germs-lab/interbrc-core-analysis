## Inter-BRC Core Analysis
## Core & non_core dbRDA for JBEI
## By Bolívar Aponte Rolón

# Setup
source("R/utils/000_setup.R")
remove(phyloseq)


#####################################
### Filtering for BRC of interest ###
#####################################
###############
# Parameters
BRC <- "jbei"
CORE <- "core"
NOCORE <- "nocore"


physeq <- subset_samples(filtered_phyloseq, brc == BRC)
###############

# Bray-Curtis core
braycore_summary <- extract_core(physeq,
  Var = "site",
  method = "increase",
  increase_value = 2
)

# Subset by "core" and "non-core"

core <- subset_physeq(braycore_summary,
  physeq,
  .var = "otu",
  type = "core"
)

nocore <- subset_physeq(braycore_summary,
  physeq,
  .var = "otu",
  type = "no"
)

# Helper function to save phyloseqs of interest
name_object <- function(brc, physeq, type) {
  obj_name <- stringr::str_glue("{brc}_{type}_phyloseq")
  assign(obj_name,
    value = physeq,
    envir = .GlobalEnv
  )
}

# Save phyloseqs
name_object(brc = BRC, physeq = physeq, type = "core") %>%
  saveRDS(., file = stringr::str_glue("data/output/phyloseq_objects/{BRC}_{CORE}_phyloseq.rda")) # For some reason saving like this result sin magic number X

name_object(brc = BRC, physeq = physeq, type = "nocore") %>%
  saveRDS(., file = stringr::str_glue("data/output/phyloseq_objects/{BRC}_{NOCORE}_phyloseq.rda"))


# Transformed matrices
### BC core
hell_matrix <- decostand(t(core$asv_matrix), method = "hellinger", MARGIN = 1)

dbrda_traits <- braycore_summary$sample_metadata %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  select(., c(drought, treatment, block, harvest, sample_id)) %>%
  dplyr::filter(rownames(.) %in% colnames(hell_matrix)) # Match hell_matrix samples


# dbRDA Modelling
dbrda_00_core <- dbrda(t(hell_matrix) ~ 1,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_01_core <- dbrda(t(hell_matrix) ~ .,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
) # Model with all explanatory variables.

dbrda_02_core <- dbrda(t(hell_matrix) ~ drought + block + harvest,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8, # Passes parallelization to metaMDS function
  na.action = na.omit
)

# Anovas (PERMANOVAs)
set.seed(123)
anova.cca(dbrda_02_nocore, by = "margin", permutations = 999, parallel = 8) # .cca adds Distance based CCA functionality
anova(dbrda_02_core, by = "axis")
anova(dbrda_02_core, by = "axis", perm.max = 500)


### BC Non-Core
hell_matrix <- decostand(t(nocore$asv_matrix), method = "hellinger", MARGIN = 1)

dbrda_traits <- braycore_summary$sample_metadata %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  select(., c(drought, treatment, block, harvest, sample_id)) %>%
  dplyr::filter(rownames(.) %in% colnames(hell_matrix)) # Match hell_matrix samples


# dbRDA Modelling
dbrda_00_nocore <- dbrda(t(hell_matrix) ~ 1,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_01_nocore <- dbrda(t(hell_matrix) ~ .,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
) # Model with all explanatory variables.

dbrda_02_nocore <- dbrda(t(hell_matrix) ~ drought + block + harvest,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8, # Passes parallelization to metaMDS function
  na.action = na.omit
)

######################
### Visualizations ###
######################

core_dbrda_gg <- brc_flex_ordi(dbrda_02_core, dbrda_traits,
  color_var = "treatment",
  sample_id_col = "sample_id"
) %>%
  +ggtitle(str_glue("{BRC}: Bray-Curtis Core"))

nocore_dbrda_gg <- brc_flex_ordi(dbrda_02_nocore, dbrda_traits,
  color_var = "treatment",
  sample_id_col = "sample_id"
) %>%
  +ggtitle(str_glue("{BRC}: Bray-Curtis Non-Core "))

# Save ggplots

save_gg <- list(core_dbrda_gg, nocore_dbrda_gg) %>%
  purrr::set_names(
    stringr::str_glue("{BRC}_{c(CORE, NOCORE)}_dbrda")
  )

brc_ggsave(save_gg, "data/output/plots/")

###############################################
### Core and Non-Core selected by threshold ###
###############################################

high_occ <- subset_samples(core_asvs_threshold$physeq_high_occ, brc == BRC)
low_occ <- subset_samples(core_asvs_threshold$physeq_low_occ, brc == BRC)
