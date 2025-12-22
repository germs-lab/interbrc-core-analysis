#########################################################
# COMMUNITY ORDINATION ANALYSIS
# NMDS and PCoA ordinations and dbRDA analysis and ordinations
# of core and non-core communities from identify_core()
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-05-08
#########################################################

# DESCRIPTION:
# This script performs ordination analyses (NMDS and PCoA) to visualize and compare
# community structure of core and non-core microbiomes across different BRCs. As well as,
# performs distance-based redundancy analysis (dbRDA) on core and non-core
# microbial communities specifically for the JBEI dataset. It analyzes the influence
# of environmental variables on community structure.

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
source("R/utils/000_setup.R")
if (exists("phyloseq")) {
  remove(phyloseq)
}

#--------------------------------------------------------
# DATA TRANSFORMATION
#--------------------------------------------------------
# Transform community matrices using Hellinger transformation
hell_matrices <- purrr::map(
  asv_matrices,
  ~ {
    decostand(t(.x), method = "hellinger", MARGIN = 1)
  }
)

# Calculate distance matrices in parallel
plan(multicore, workers = parallel::detectCores() - 1)

distance_matrices <- hell_matrices %>%
  future_map(
    ~ vegdist(
      t(.x),
      method = "bray",
      upper = FALSE,
      binary = FALSE,
      na.rm = TRUE
    ),
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  ) %>%
  purrr::set_names(names(hell_matrices))

plan(sequential) # Close parallel processing

save(distance_matrices, file = here::here("data/output/distance_matrices.rda"))

# Note: Two dimensions keeps stress below 0.20 (from previous analysis)

#--------------------------------------------------------
# NMDS ANALYSIS: BC_CORE COMMUNITY
#--------------------------------------------------------
# Perform NMDS on BC_CORE COMMUNITY
bc_core_nmds <- brc_nmds(
  asv_matrix = asv_matrices$bc_core,
  physeq = braycurt_core,
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
)

save(bc_core_nmds, file = here::here("data/output/bc_core_nmds.rda"))


# Generate crop-based and BRC-based visualizations
bc_core_nmds_crops <- brc_gg_ordi(
  .data = bc_core_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle(
    "50 core ASVs in BRC crops",
    subtitle = "Core that contributes 2% to Bray-Curtis"
  )

bc_core_nmds_brc <- brc_gg_ordi(
  .data = bc_core_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle(
    "50 core ASVs in BRC crops",
    subtitle = "Core that contributes 2% to Bray-Curtis"
  )

#--------------------------------------------------------
# NMDS ANALYSIS: BC_NONCORE COMMUNITY
#--------------------------------------------------------
# Perform NMDS on BC_NONCORE COMMUNITY
bc_noncore_nmds <- brc_nmds(
  asv_matrix = asv_matrices$bc_noncore,
  physeq = braycurt_noncore,
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
)


save(bc_noncore_nmds, file = here::here("data/output/bc_noncore_nmds.rda"))

# Generate crop-based and BRC-based visualizations
bc_noncore_nmds_crops <- brc_gg_ordi(
  .data = bc_noncore_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("Non-core ASVs in BRC crops (100% samples)")

bc_noncore_nmds_brc <- brc_gg_ordi(
  .data = bc_noncore_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("Non-core ASVs in BRC (100% samples)")


#--------------------------------------------------------
# NMDS ANALYSIS: THRESHOLD-BASED CORE
#--------------------------------------------------------
# Perform NMDS on high-occupancy (threshold-based core) community
high_occ_nmds <- brc_nmds(
  asv_matrix = asv_matrices$high_occ,
  physeq = prevalence_core$physeq_high_occ,
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
)

# Generate crop-based and BRC-based visualizations
high_occ_nmds_crops <- brc_gg_ordi(
  .data = high_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("12 high prevalence ASVs in BRC crops (>60% samples)")

high_occ_nmds_brc <- brc_gg_ordi(
  .data = high_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("12 high prevalence ASVs in BRCs (>60% samples)")

#--------------------------------------------------------
# NMDS ANALYSIS: THRESHOLD-BASED NON-CORE
#--------------------------------------------------------
# Perform NMDS on low-occupancy (threshold-based non-core) community
low_occ_nmds <- brc_nmds(
  asv_matrix = asv_matrices$low_occ,
  physeq = subset_samples(prevalence_core$physeq_low_occ, brc != "jbei"),
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
)

save(low_occ_nmds, file = here::here("data/output/low_occ_nmds.rda"))

# Generate crop-based and BRC-based visualizations
low_occ_nmds_crops <- brc_gg_ordi(
  .data = low_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("Low prevalence ASVs in BRC crops (<60% samples)")

low_occ_nmds_brc <- brc_gg_ordi(
  .data = low_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("Low prevalence ASVs in BRCs (<60% samples)")

#--------------------------------------------------------
# SAVE NMDS VISUALIZATIONS
#--------------------------------------------------------
# Combine all NMDS plots and save to output directory
nmds_plots <- list(
  bc_core_nmds_brc,
  bc_core_nmds_crops,
  bc_noncore_nmds_brc,
  bc_noncore_nmds_crops,
  high_occ_nmds_crops,
  high_occ_nmds_brc,
  low_occ_nmds_crops,
  low_occ_nmds_brc
)
plot_names <- c(
  "bc_core_nmds_brc",
  "bc_core_nmds_crops",
  "bc_noncore_nmds_brc",
  "bc_noncore_nmds_crops",
  "high_occ_nmds_crops",
  "high_occ_nmds_brc",
  "low_occ_nmds_crops",
  "low_occ_nmds_brc"
)

plot_paths <- str_glue(here::here("data/output/plots/{tolower(plot_names)}.png"))

purrr::walk2(
  plot_paths,
  nmds_plots,
  \(path, plot) {
    ggsave(
      filename = path,
      plot = plot,
      dpi = 300,
      width = 200,
      height = 200,
      units = "mm"
    )
  }
)

#--------------------------------------------------------
# PCOA ANALYSIS: BC_CORE COMMUNITY
#--------------------------------------------------------
# Perform PCoA on BC_CORE COMMUNITY
bc_core_asv_pcoa <- brc_pcoa(
  distance_matrices$bc_core,
  braycurt_core
)

# Generate crop-based and BRC-based visualizations
bc_core_pcoa_brc <- brc_gg_ordi(
  .data = bc_core_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = brc,
  .drop_na = brc
) +
  ggtitle(
    "PCoA: 50 core ASVs in BRC",
    subtitle = "Core that contributes 2% to Bray-Curtis"
  )

bc_core_pcoa_crops <- brc_gg_ordi(
  .data = bc_core_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = crop,
  .drop_na = brc
) +
  ggtitle(
    "PCoA: 50 core ASVs in BRC crops",
    subtitle = "Core that contributes 2% to Bray-Curtis"
  )

#--------------------------------------------------------
# PCOA ANALYSIS: BC_NONCORE COMMUNITY
#--------------------------------------------------------
# Perform PCoA on BC_NONCORE COMMUNITY
bc_noncore_asv_pcoa <- brc_pcoa(
  distance_matrices$bc_noncore,
  braycurt_noncore
)

# Generate crop-based and BRC-based visualizations
bc_noncore_pcoa_brc <- brc_gg_ordi(
  .data = bc_noncore_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("PCoA: Non-Core ASVs in BRCs")

bc_noncore_pcoa_crops <- brc_gg_ordi(
  .data = bc_noncore_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("PCoA: Non-Core ASVs in BRC crops")

#--------------------------------------------------------
# PCOA ANALYSIS: THRESHOLD-BASED COMMUNITIES
#--------------------------------------------------------
# Perform PCoA on high-occupancy (threshold-based core) community
high_asv_pcoa <- brc_pcoa(
  distance_matrices$high_occ,
  prevalence_core$physeq_high_occ %>%
    prune_samples(sample_sums(.) > 0, .)
)

# Generate crop-based and BRC-based visualizations
high_occ_pcoa_brc <- brc_gg_ordi(
  .data = high_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("PCoA: High prevalence ASVs in BRCs")

high_occ_pcoa_crops <- brc_gg_ordi(
  .data = high_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("PCoA: High prevalence ASVs in BRC crops")

# Perform PCoA on low-occupancy (threshold-based non-core) community
low_asv_pcoa <- brc_pcoa(
  distance_matrices$low_occ,
  prevalence_core$physeq_low_occ %>%
    prune_samples(sample_sums(.) > 0, .)
)

# Generate crop-based and BRC-based visualizations
low_occ_pcoa_brc <- brc_gg_ordi(
  .data = low_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("PCoA: Low prevalence ASVs in BRCs")

low_occ_pcoa_crops <- brc_gg_ordi(
  .data = low_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("PCoA: Low prevalence ASVs in BRC crops")

#--------------------------------------------------------
# SAVE PCOA VISUALIZATIONS
#--------------------------------------------------------
# Combine all PCoA plots and save to output directory
pcoa_plots <- list(
  bc_core_pcoa_brc,
  bc_core_pcoa_crops,
  bc_noncore_pcoa_brc,
  bc_noncore_pcoa_crops,
  high_occ_pcoa_crops,
  high_occ_pcoa_brc,
  low_occ_pcoa_crops,
  low_occ_pcoa_brc
)
plot_names <- c(
  "bc_core_pcoa_brc",
  "bc_core_pcoa_crops",
  "bc_noncore_pcoa_brc",
  "bc_noncore_pcoa_crops",
  "high_occ_pcoa_crops",
  "high_occ_pcoa_brc",
  "low_occ_pcoa_crops",
  "low_occ_pcoa_brc"
)

plot_paths <- str_glue(here::here("data/output/plots/{tolower(plot_names)}.png"))

purrr::walk2(
  plot_paths,
  pcoa_plots,
  \(path, plot) {
    ggsave(
      filename = path,
      plot = plot,
      dpi = 300,
      width = 200,
      height = 200,
      units = "mm"
    )
  }
)

#--------------------------------------------------------
# DISTANCE-BASED REDUNDANCY ANALYSIS
#--------------------------------------------------------

#########################################################
# Ignore this section is performing for ALL BRCs
# Script 005_core_selection_per_BRCs.R performs the same
#--------------------------------------------------------
# BRC SELECTION AND PARAMETERS
#--------------------------------------------------------
# Define BRC of interest and output naming parameters
BRC <- "jbei"
CORE <- "bc_core"
NOCORE <- "bc_noncore"

# Filter phyloseq object for the selected BRC
physeq <- subset_samples(filtered_phyloseq, brc == BRC)


# Extract core microbiome using Bray-Curtis dissimilarity
braycore_summary <- identify_core(
  physeq,
  Var = "site",
  method = "increase",
  increase_value = 2
)

##########################################################

#--------------------------------------------------------
# Subset by "core" and "non-core" taxa
#--------------------------------------------------------

bc_core <- subset_physeq(braycore_summary, physeq, .var = "otu", type = "core")
bc_noncore <- subset_physeq(braycore_summary, physeq, .var = "otu", type = "no")


#--------------------------------------------------------
# dbRDA BC_CORE COMMUNITY ANALYSIS
#--------------------------------------------------------
# Transform data using Hellinger transformation
hell_matrix <- decostand(
  t(bc_core$asv_matrix),
  method = "hellinger",
  MARGIN = 1
)

# Prepare environmental data
dbrda_traits <- braycore_summary$sample_metadata %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  select(., c(drought, treatment, block, harvest, sample_id)) %>%
  dplyr::filter(rownames(.) %in% colnames(hell_matrix))

# Build dbRDA models with different constraints
dbrda_00_bc_core <- dbrda(
  t(hell_matrix) ~ 1,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_01_bc_core <- dbrda(
  t(hell_matrix) ~ .,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_02_bc_core <- dbrda(
  t(hell_matrix) ~ drought + block + harvest,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

#--------------------------------------------------------
# dbRDA BC_NONCORE COMMUNITY ANALYSIS
#--------------------------------------------------------
# Transform non-core data using Hellinger transformation
hell_matrix <- decostand(
  t(bc_noncore$asv_matrix),
  method = "hellinger",
  MARGIN = 1
)

# Prepare environmental data for non-core analysis
dbrda_traits <- braycore_summary$sample_metadata %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  select(., c(drought, treatment, block, harvest, sample_id)) %>%
  dplyr::filter(rownames(.) %in% colnames(hell_matrix))

# Build dbRDA models for non-core community
dbrda_00_bc_noncore <- dbrda(
  t(hell_matrix) ~ 1,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_01_bc_nocore <- dbrda(
  t(hell_matrix) ~ .,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_02_bc_noncore <- dbrda(
  t(hell_matrix) ~ drought + block + harvest,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

#--------------------------------------------------------
# dbRDA STATISTICAL TESTING
#--------------------------------------------------------
# Perform permutation tests on models
set.seed(123)
anova.cca(dbrda_02_bc_noncore, by = "margin", permutations = 999, parallel = 8)
anova(dbrda_02_bc_core, by = "axis")
anova(dbrda_02_bc_core, by = "axis", perm.max = 500)


#--------------------------------------------------------
# dbRDA VISUALIZATION
#--------------------------------------------------------
# Generate ordination plots for core and non-core communities
bc_core_dbrda_gg <- brc_flex_ordi(
  dbrda_02_bc_core,
  dbrda_traits,
  color_var = "treatment",
  sample_id_col = "sample_id"
) %>%
  +ggtitle(str_glue("All {BRC}: Bray-Curtis Core"))

bc_noncore_dbrda_gg <- brc_flex_ordi(
  dbrda_02_bc_noncore,
  dbrda_traits,
  color_var = "treatment",
  sample_id_col = "sample_id"
) %>%
  +ggtitle(str_glue("All {BRC}s: Bray-Curtis Non-Core "))

# Save plots to output directory
save_gg <- list(bc_core_dbrda_gg, bc_noncore_dbrda_gg) %>%
  purrr::set_names(
    stringr::str_glue("brc_{c(CORE, NOCORE)}_dbrda")
  )

brc_ggsave(save_gg, here::here("data/output/plots/"))

#--------------------------------------------------------
# THRESHOLD-BASED ANALYSIS
#--------------------------------------------------------
# Extract high and low occupancy communities from threshold-based core analysis
high_occ <- subset_samples(prevalence_core$physeq_high_occ, brc == BRC)
low_occ <- subset_samples(prevalence_core$physeq_low_occ, brc == BRC)
