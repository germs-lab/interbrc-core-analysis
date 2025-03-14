## Inter-BRC Core Analysis
## Core & non_core Ordinations
## By Bolívar Aponte Rolón

# Setup
source("R/utils/000_setup.R")
load("data/output/core_summary_lists.rda")
load("data/output/phyloseq_objects/filtered_phyloseq.rda")
load("data/output/phyloseq_objects/core_brc_phyloseq.rda")
load("data/output/phyloseq_objects/non_core_brc_phyloseq.rda")
load("data/output/core_asv_matrix.rda")
load("data/output/non_core_asv_matrix.rda")


#################################################
###### Ordinations and Community Analysis #######
#################################################

############
### Core ###
############
# Hellinger transformation of matrices
# To be used for NMDS, dbRDA and adonis2
core_hell_matrix <- decostand(t(core_asv_matrix),
  MARGIN = 1,
  method = "hellinger"
) # Now we need samples as columns and ASV are rows

###########

# Ordinations
core_asv_dist <- vegdist(t(core_hell_matrix),
  method = "bray",
  upper = FALSE,
  binary = FALSE,
  na.rm = TRUE
)


## Choosing the number of dimensions
# set.seed(484035)
# nmds_screen_parallel(core_asv_dist, ncores = 32) # Results: Two dimensions keeps stress below 0.20


core_ordi <- metaMDS(as.matrix(core_asv_dist),
  distance = "bray",
  display = c("sites"),
  noshare = TRUE,
  autotransform = FALSE,
  wascores = TRUE,
  tidy = TRUE,
  k = 2,
  trymax = 100,
  parallel = parallel::detectCores()
)

stressplot(core_ordi)

## Scores and sample data for NMDS
vegan::sppscores(core_ordi) <- t(core_hell_matrix)

core_nmds_scores <- as.data.frame(vegan::scores(core_ordi)$sites)

rownames(core_nmds_scores) <- rownames(vegan::scores(core_ordi)$sites)

core_nmds_scores <- core_nmds_scores %>%
  rownames_to_column(., var = "unique_id")

core_brc_sample_df <- core_brc_phyloseq@sam_data %>%
  data.frame() %>%
  rownames_to_column(., var = "unique_id")

core_nmds_df <- right_join(core_brc_sample_df, core_nmds_scores, by = "unique_id")


# Core NMDS Aesthetics ####
core_nmds_crops <- core_nmds(
  .data = core_nmds_df,
  .color = crop,
  .drop_na = brc
) + ggtitle("50 core ASVs in BRC crops (100% samples)")

core_nmds_brc <- core_nmds(
  .data = core_nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("50 core ASVs in BRCs (100% samples)")


ggsave(
  filename = "core_nmds_crops.png",
  plot = core_nmds_crops,
  path = "data/output/plots/",
  dpi = 300,
  width = 200,
  height = 200,
  units = "mm"
)

ggsave(
  filename = "core_nmds_brc.png",
  plot = core_nmds_brc,
  path = "data/output/plots/",
  dpi = 300,
  width = 200,
  height = 200,
  units = "mm"
)

# Threshold analysis

# Calculate occurrence of each ASV across samples
asv_sample_counts <- colSums(core_asv_matrix > 0)

# Get total number of samples
total_samples <- nrow(core_asv_matrix)

# Define occurrence threshold (60% of samples)
threshold_60 <- total_samples * 0.6

# Select ASVs based on occurrence criteria
high_occurrence_asvs <- names(asv_sample_counts[asv_sample_counts >= threshold_60])
low_occurrence_asvs <- names(asv_sample_counts[asv_sample_counts < threshold_60])

# Print ASV counts
cat("ASVs found in ≥60% samples:", length(high_occurrence_asvs), "\n")
cat("ASVs found in <60% samples:", length(low_occurrence_asvs), "\n")





################
### Non-Core ###
################

# Hellinger transformation of matrices
non_core_hell_matrix <- decostand(t(non_core_asv_matrix),
  MARGIN = 1,
  method = "hellinger"
) # Now we need samples as columns and ASV are rows

###########

# Ordinations
non_core_asv_dist <- vegdist(t(non_core_hell_matrix),
  method = "bray",
  upper = FALSE,
  binary = FALSE,
  na.rm = TRUE
)


non_core_ordi <- metaMDS(as.matrix(non_core_asv_dist),
  distance = "bray",
  display = c("sites"),
  noshare = TRUE,
  autotransform = FALSE,
  wascores = TRUE,
  zerodist = "ignore",
  tidy = TRUE,
  k = 2,
  trymax = 100,
  parallel = parallel::detectCores()
)

stressplot(non_core_ordi)

## Scores and sample data for NMDS
vegan::sppscores(non_core_ordi) <- t(non_core_hell_matrix)

non_nmds_scores <- as.data.frame(vegan::scores(non_core_ordi)$sites)

rownames(non_nmds_scores) <- rownames(vegan::scores(non_core_ordi)$sites)

non_nmds_scores <- non_nmds_scores %>%
  rownames_to_column(., var = "unique_id")

non_core_brc_sample_df <- non_core_brc_phyloseq@sam_data %>%
  data.frame() %>%
  rownames_to_column(., var = "unique_id")

non_core_nmds_df <- right_join(non_core_brc_sample_df, non_nmds_scores, by = "unique_id")

# Non-Core NMDS Aesthetics ####
non_core_nmds_crops <- core_nmds(
  .data = non_core_nmds_df,
  .color = crop,
  .drop_na = brc
) + ggtitle("Non-core ASVs in BRC crops (100% samples)")


non_core_nmds_brc <- core_nmds(
  .data = non_core_nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("Non-core ASVs in BRC (100% samples)")



ggsave(
  filename = "non_core_nmds_crops.png",
  plot = non_core_nmds_crops,
  path = "data/output/plots/",
  dpi = 300,
  width = 200,
  height = 200,
  units = "mm"
)

ggsave(
  filename = "non_core_nmds_brc.png",
  plot = non_core_nmds_brc,
  path = "data/output/plots/",
  dpi = 300,
  width = 200,
  height = 200,
  units = "mm"
)
