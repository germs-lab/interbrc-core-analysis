# Exploring JBEI dataset for microbial core selection

# Setup
source("R/utils/000_setup.R")
if (exists("phyloseq")) remove(phyloseq)


###################################################
### Microbiome Core Selection via extract_core() ###
###################################################

# Data set clean up

new_metadata <- drought_jbei %>%
  sample_data() %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(
    across(brc, ~ str_to_lower(.)),
    across(everything(.), ~ as.character(.)),
    new_row = x_sample_id
  ) %>%
  column_to_rownames(., var = "new_row") %>% # Workaround to inserting "sa1" type rownames
  sample_data()

# Update phyloseq object
sample_data(drought_jbei) <- new_metadata

# save(drought_jbei, file = "data/output/phyloseq_objects/jbei/drought_jbei.rda")

# Check OTU table
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
# Obeser the data and see the patterns. Samples with less than 20 reads creat a weird pattern.
# Extract core
jbei_core_summary_lists <- extract_core(
  drought_jbei,
  Var = "treatment",
  method = "increase",
  increase_value = 2
)

# Plot Bray-Curtis Dissimilarity Curve:
bray_curtis_curve <- brc_bc_curve(
  core_summary_list = jbei_core_summary_lists,
  max_otus = 100,
  threshold = 1.02
)
bray_curtis_curve

occ_abun_plot <- brc_bc_occ_curve(core_summary_list = jbei_core_summary_lists)
occ_abun_plot
