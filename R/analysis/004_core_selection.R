##################################################################
# CORE MICROBIOME SELECTION
# Extract and analyze core microbiome across samples from all BRCs
#
# Project:  Inter-BRC-Core-Microbiome
# Original Author: Brandon Kristy
# Modified by: Bolívar Aponte Rolón
# Date: 2025-02-20
##################################################################

# DESCRIPTION:
# This script identifies the core microbiome across samples using multiple approaches:
# 1. extract_core() method based on Bray-Curtis dissimilarity
# 2. Threshold-based approach
# 3. Neutral model fitting for abundance-occupancy patterns

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
source("R/utils/000_setup.R")
if (exists("phyloseq")) {
  remove(phyloseq)
}

#--------------------------------------------------------
# CORE EXTRACTION USING EXTRACT_CORE()
#--------------------------------------------------------
# Ensure minimum sample quality
filtered_phyloseq <- prune_samples(
  sample_sums(filtered_phyloseq) >= 100,
  filtered_phyloseq
)

# Add rarefaction metrics
rare_metrics <- add_rarefaction_metrics(filtered_phyloseq)
rare_metrics_plot <- plot_rarefaction_metrics(rare_metrics)

rare_metrics_plot
sum(sample_sums(filtered_phyloseq) > 7000) # How many samples above read count threshold?
# Minimum seq depth was ~7,000 reads selected for 1646 samples.

# Perform multiple rarefaction
interbrc_rarefied <- multi_rarefy(
  physeq = filtered_phyloseq,
  depth_level = 7000,
  num_iter = 100,
  threads = 4,
  set_seed = 7642
)

# Update phyloseq object with rarefied data
interbrc_rarefied <- update_otu_table(
  physeq = filtered_phyloseq,
  otu_rare = interbrc_rarefied
)


braycore_summary <- BRCore::identify_core(
  filtered_phyloseq,
  priority_var = "site",
  increase_value = 0.02,
  abundance_weight = 0,
  seed = 7895
)


# Save results to avoid recomputation
save(braycore_summary, file = here::here("data/output/braycore_summary.rda"))


#--------------------------------------------------------
# VISUALIZATION OF BRAY-CURTIS AND OCCUPANCY PATTERNS
#--------------------------------------------------------
# Generate Bray-Curtis dissimilarity curve
bray_curtis_curve <- brc_bc_curve(core_summary_list = braycore_summary)


# Generate abundance-occupancy plot
occ_abun_plot <- brc_occ_curve(core_summary_list = braycore_summary)
occ_abun_plot


# Save abundance-occupancy plot
ggsave(
  filename = "bray_curtis_abundance_occupancy.png",
  occ_abun_plot,
  path = here::here("data/output/plots/"),
  dpi = 300,
  width = 6,
  height = 4
)

#--------------------------------------------------------
# NEUTRAL MODEL FITTING FOR ABUNDANCE-OCCUPANCY
#--------------------------------------------------------

braycore_summary_neutral_fit <- fit_neutral_model(
  otu_table = braycore_summary$otu_table,
  core_set = braycore_summary$increase_core,
  abundance_occupancy = braycore_summary$abundance_occupancy
)

#--------------------------------------------------------
# NEUTRAL MODEL VISUALIZATION
#--------------------------------------------------------
plot_neutral_fit <- plot_neutral_model(
  braycore_summary_neutral_fit
)
plot_neutral_fit


#--------------------------------------------------------
# THRESHOLD-BASED CORE SELECTION
#--------------------------------------------------------
# Extract core ASVs based on presence threshold (Jae's method)
prevalence_core <- filter_core(
  filtered_phyloseq,
  threshold = 0.6,
  as = "rows"
)

# Save threshold-based core ASVs
save(
  prevalence_core,
  file = here::here("data/output/phyloseq_objects/prevalence_core.rda")
)

#--------------------------------------------------------
# JBEI-SPECIFIC CORE ANALYSIS
#--------------------------------------------------------
# Clean up JBEI metadata
new_metadata <- drought_jbei %>%
  sample_data() %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(
    across(brc, ~ str_to_lower(.)),
    across(everything(.), ~ as.character(.)),
    new_row = x_sample_id
  ) %>%
  column_to_rownames(., var = "new_row") %>%
  sample_data()

# Update JBEI phyloseq object with cleaned metadata
sample_data(drought_jbei) <- new_metadata

# Filter JBEI samples for quality
drought_jbei <- prune_samples(
  sample_sums(drought_jbei) >= 100,
  drought_jbei
)

drought_jbei <- filter_taxa(
  drought_jbei,
  function(x) {
    sum(x > 100) > (0.00 * length(x)) # Results depend on this cut-off.
  },
  TRUE
)

# Extract JBEI-specific core
jbei_braycore_summary <- extract_core(
  drought_jbei,
  Var = "treatment",
  method = "increase",
  increase_value = 2
)

# Generate JBEI-specific visualizations
bray_curtis_curve <- brc_bc_curve(
  core_summary_list = jbei_braycore_summary,
  max_otus = 100,
  threshold = 1.02
)
bray_curtis_curve

occ_abun_plot <- brc_bc_occ_curve(core_summary_list = jbei_braycore_summary)
occ_abun_plot
