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
) # Now we need samples as columns and ASV are rows

###########

# Ordinations
core_asv_dist <- vegdist(t(core_hell_matrix),
  method = "bray",
  upper = FALSE,
  binary = FALSE,
  na.rm = TRUE
)

## Choosing the number of dimensions for NMDS
set.seed(484035)
nmds_screen_parallel(core_asv_dist, ncores = 8) # Results: Two dimensions keeps stress below 0.20

core_ext_core <- calculate_nmds(
  asv_matrix = core_asv_matrix,
  physeq = core_brc_phyloseq,
  ncores = parallel::detectCores(),
  k = 2,
  trymax = 100
)

# Core NMDS Aesthetics ####
core_nmds_crops <- gg_nmds(
  .data = core_ext_core$nmds_df,
  .color = crop,
  .drop_na = brc
) + ggtitle("50 core ASVs in BRC crops",
  subtitle = "Core that contributes 2% to Bray-Curtis"
)

core_nmds_brc <- gg_nmds(
  .data = core_ext_core$nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("50 core ASVs in BRC crops",
  subtitle = "Core that contributes 2% to Bray-Curtis"
)


##################################
### Non-Core by ExctractCore() ###
##################################

non_core_ext_core <- calculate_nmds(
  asv_matrix = non_core_asv_matrix,
  physeq = non_core_brc_phyloseq,
  ncores = parallel::detectCores(),
  k = 2,
  trymax = 100
)

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

test <- calculate_nmds(physeq_high_occ_matrix, physeq_high_occ)
# Core NMDS Aesthetics ####
test_nmds_crops <- gg_nmds(
  .data = high_occ_nmds_df,
  .color = crop,
  .drop_na = brc
) + ggtitle("12 core ASVs in BRC crops (60% samples)")

core_nmds_brc <- core_nmds(
  .data = core_nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("50 core ASVs in BRCs (100% samples)")



## Saving plots
core_plots <- list(core_nmds_brc, core_nmds_crops)
plot_names <- c(
  "core_nmds_brc",
  "core_nmds_crops",
  "non_core_nmds_brc",
  "non_core_nmds_crops"
)


plot_paths <- str_glue("data/output/plots/{tolower(plot_names)}.png")

purrr::walk2(
  plot_paths, core_plots,
  \(path, plot) ggsave(
    filename = path,
    plot = plot,
    dpi = 300,
    width = 200,
    height = 200,
    units = "mm"
  )
)
