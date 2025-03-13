## Inter-BRC Core Analysis
## Core & non_core Ordinations
## By Bolívar Aponte Rolón

# Setup
source("R/utils/000_setup.R")

#################################################
###### Ordinations and Community Analysis #######
#################################################

##############################
### Core by ExctractCore() ###
##############################
# Hellinger transformation of matrices
# To be used for NMDS, dbRDA and adonis2
core_hell_matrix <- decostand(t(core_asv_matrix),
  MARGIN = 1,
  method = "hellinger"
) %>% # Now we need samples as columns and ASV are rows
  as.data.frame() %>%
  select(where(~ is.numeric(.) && sum(.) > 0)) %>% # Removal of columns that sum 0
  as.matrix()

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
core_nmds_crops <- gg_nmds(
  .data = core_nmds_df,
  .color = crop,
  .drop_na = brc
) + ggtitle("50 core ASVs in BRC crops (100% samples)")

core_nmds_brc <- gg_nmds(
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


##################################
### Non-Core by ExctractCore() ###
##################################

# Hellinger transformation of matrices
non_core_hell_matrix <- decostand(t(non_core_asv_matrix),
  MARGIN = 1,
  method = "hellinger"
) %>%
  as.data.frame() %>%
  select(where(~ is.numeric(.) && sum(.) > 0)) %>% # Removal of columns that sum 0
  as.matrix()

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


######################################
### Core and Non-Core by Threshold ###
######################################

# Perform NMDS on selected ASVs

# Hellinger transformation of matrices
# High Prevalence/Occupancy

high_occ_hell_matrix <- decostand(t(physeq_high_occ_matrix),
  MARGIN = 1,
  method = "hellinger"
) %>%
  as.data.frame() %>%
  select(where(~ is.numeric(.) && sum(.) > 0)) %>% # Removal of columns that sum 0
  as.matrix()

###########

# Ordinations
high_occ_asv_dist <- vegdist(t(high_occ_hell_matrix),
  method = "bray",
  upper = FALSE,
  binary = FALSE,
  na.rm = TRUE
)

## Choosing the number of dimensions
# set.seed(484035)
# nmds_screen_parallel(core_asv_dist, ncores = 32) # Results: Two dimensions keeps stress below 0.20

high_occ_ordi <- metaMDS(as.matrix(high_occ_asv_dist),
  distance = "bray",
  display = c("sites"),
  noshare = TRUE,
  autotransform = FALSE,
  wascores = TRUE,
  # zerodist = "add",
  tidy = TRUE,
  k = 2,
  trymax = 100,
  parallel = parallel::detectCores()
)

stressplot(high_occ_ordi)

## Scores and sample data for NMDS
vegan::sppscores(high_occ_ordi) <- t(high_occ_hell_matrix)

high_occ_nmds_scores <- as.data.frame(vegan::scores(high_occ_ordi)$sites)

rownames(high_occ_nmds_scores) <- rownames(vegan::scores(high_occ_ordi)$sites)

high_occ_nmds_scores <- high_occ_nmds_scores %>%
  rownames_to_column(., var = "unique_id")

high_occ_sample_df <- core_asvs_threshold$physeq_high_occ@sam_data %>%
  data.frame() %>%
  rownames_to_column(., var = "unique_id")

high_occ_nmds_df <- right_join(high_occ_sample_df, high_occ_nmds_scores, by = "unique_id")


# Core NMDS Aesthetics ####
test_nmds_crops <- gg_nmds(
  .data = high_occ_nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("12 core ASVs in BRC crops (60% samples)")

core_nmds_brc <- core_nmds(
  .data = core_nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("50 core ASVs in BRCs (100% samples)")



# nmds_high_occ_selected <- ordinate(physeq_high_occ, method = "NMDS", distance = "bray", trymax = 100, zerodist <- "add")
# nmds_low_occ_selected <- ordinate(physeq_low_occ, method = "NMDS", distance = "bray", trymax = 100)

# # Plot NMDS results
# nmds_plot_high_occ_selected <- plot_ordination(ps_high_occ_selected, nmds_high_occ_selected, color = "Plant", shape = "Location") +
#   geom_point(size = 3, alpha = 0.8) +
#   # facet_wrap(~Date, ncol = 2) +
#   labs(title = "NMDS of Selected Highly Prevalent ASVs", x = "NMDS1", y = "NMDS2") +
#   theme_minimal()

# nmds_plot_low_occ_selected <- plot_ordination(ps_low_occ_selected, nmds_low_occ_selected, color = "Plant", shape = "Location") +
#   geom_point(size = 3, alpha = 0.8) +
#   # facet_wrap(~Date, ncol = 2) +
#   labs(title = "NMDS of Selected Low Prevalence ASVs", x = "NMDS1", y = "NMDS2") +
#   theme_minimal()

# # Print NMDS plots
# print(nmds_plot_high_occ_selected)
# print(nmds_plot_low_occ_selected)
