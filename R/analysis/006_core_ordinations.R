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
# Choosing the number of dimensions for NMDS

# Hellinger transformation of matrices
# To be used for NMDS, dbRDA and adonis2
core_hell_matrix <- decostand(t(core_asv_matrix),
  MARGIN = 1,
  method = "hellinger"
) # Now we need samples as columns and ASV are rows


# Ordinations
core_asv_dist <- vegdist(t(core_hell_matrix),
  method = "bray",
  upper = FALSE,
  binary = FALSE,
  na.rm = TRUE
)

# set.seed(484035)
# nmds_screen_parallel(core_asv_dist, ncores = 2) # Results: Two dimensions keeps stress below 0.20
#################

core_ext_core <- calculate_nmds(
  asv_matrix = core_asv_matrix,
  physeq = core_brc_phyloseq,
  ncores = parallel::detectCores(),
  k = 3,
  trymax = 9999
) # Best solution repeated 1 times (k = 3)


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
  ncores = 1,
  k = 3,
  trymax = 9999
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

high_occ_core <- calculate_nmds(
    asv_matrix = physeq_high_occ_matrix,
    physeq = physeq_high_occ,
    ncores = parallel::detectCores(),
    k = 3,
    trymax = 9999
)


# Core NMDS Aesthetics ####
high_occ_nmds_crops <- gg_nmds(
  .data = high_occ_nmds_df,
  .color = crop,
  .drop_na = brc
) + ggtitle("12 core ASVs in BRC crops (60% samples)")

high_occ_nmds_brc <- core_nmds(
  .data = core_nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("50 core ASVs in BRCs (100% samples)")



## Saving plots
core_plots <- list(core_nmds_brc, core_nmds_crops)
plot_names <- c(
  "core_nmds_brc",
  "core_nmds_crops"
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
