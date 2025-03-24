## Inter-BRC Core Analysis
## Core & non_core Ordinations
## By Bolívar Aponte Rolón

# Setup
source("R/utils/000_setup.R")
remove(phyloseq)

#################################################
###### Ordinations and Community Analysis #######
#################################################
# List of transformed matrices
# # Transformed matrices
# hell_matrices <- purrr::map(asv_matrices, ~ {
#   decostand(t(.x), method = "hellinger", MARGIN = 1)
# })

# # Distance matrices
# plan(multisession) # Set up parallel backend

# distance_matrices <- hell_matrices %>%
#   future_map(~ vegdist(t(.x), method = "bray", na.rm = TRUE),
#     .options = furrr_options(seed = TRUE)
#   ) %>%
#   set_names(names(hell_matrices))


# distance_matrices <- purrr::map(
#   .x = hell_matrices, ~
#     vegdist(t(.x),
#       method = "bray",
#       upper = FALSE,
#       binary = FALSE,
#       na.rm = TRUE
#     )
# )


##############################
### Core by ExctractCore() ###
##############################
# Choosing the number of dimensions for NMDS

# set.seed(484035)
# nmds_screen_parallel(distance_matrices$core, ncores = 2) # Results: Two dimensions keeps stress below 0.20
#################

core_ext_nmds <- calculate_nmds(
  asv_matrix = asv_matrices$core,
  physeq = core_brc_phyloseq,
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
) # Best solution repeated 1 times (k = 3)

save(core_ext_nmds, file = "data/output/core_ext_nmds.rda")

# Core NMDS Aesthetics ####
core_nmds_crops <- gg_nmds(
  .data = core_ext_nmds$nmds_df,
  .color = crop,
  .drop_na = brc
) + ggtitle("50 core ASVs in BRC crops",
  subtitle = "Core that contributes 2% to Bray-Curtis"
)

core_nmds_brc <- gg_nmds(
  .data = core_ext_nmds$nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("50 core ASVs in BRC crops",
  subtitle = "Core that contributes 2% to Bray-Curtis"
)

##################################
### Non-Core by ExctractCore() ###
##################################

non_core_ext_nmds <- calculate_nmds(
  asv_matrix = asv_matrices$non_core,
  physeq = non_core_brc_phyloseq,
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
) # Best solution was not repeated (k = 3)

# save(non_core_ext_nmds, file = "data/output/non_core_ext_nmds.rda")

# Non-Core NMDS Aesthetics ####
non_core_nmds_crops <- gg_nmds(
  .data = non_core_ext_nmds$nmds_df,
  .color = crop,
  .drop_na = brc
) + ggtitle("Non-core ASVs in BRC crops (100% samples)")


non_core_nmds_brc <- gg_nmds(
  .data = non_core_ext_nmds$nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("Non-core ASVs in BRC (100% samples)")



######################################
### Core and Non-Core by Threshold ###
######################################

# High occupancy
high_occ_nmds <- calculate_nmds(
  asv_matrix = asv_matrices$high_occ,
  physeq = core_asvs_threshold$physeq_high_occ,
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
)

# save(high_occ_nmds, file = "data/output/high_occ_nmds.rda")


high_occ_nmds_crops <- gg_nmds(
  .data = high_occ_nmds$nmds_df,
  .color = crop,
  .drop_na = brc
) + ggtitle("12 core ASVs in BRC crops (60% samples)")

high_occ_nmds_brc <- gg_nmds(
  .data = high_occ_nmds$nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("Non-Core ASVs in BRCs (60% samples)")


# Low occupancy
low_occ_nmds <- calculate_nmds(
  asv_matrix = asv_matrices$low_occ,
  physeq = subset_samples(core_asvs_threshold$physeq_low_occ, brc != "jbei"), # Remove "jbei"
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
)

save(low_occ_nmds, file = "data/output/low_occ_nmds.rda")


low_occ_nmds_crops <- gg_nmds(
  .data = low_occ_nmds$nmds_df,
  .color = crop,
  .drop_na = brc
) + ggtitle("12 core ASVs in BRC crops (60% samples)")

low_occ_nmds_brc <- gg_nmds(
  .data = low_occ_nmds$nmds_df,
  .color = brc,
  .drop_na = brc
) + ggtitle("Non-Core ASVs in BRCs (<60% samples)")


## Saving plots

core_plots <- list(
  core_nmds_brc, core_nmds_crops,
  non_core_nmds_brc, non_core_nmds_crops,
  high_occ_nmds_crops, high_occ_nmds_brc,
  low_occ_nmds_crops, low_occ_nmds_brc
)
plot_names <- c(
  "core_nmds_brc", "core_nmds_crops",
  "non_core_nmds_brc", "non_core_nmds_crops",
  "high_occ_nmds_crops", "high_occ_nmds_brc",
  "low_occ_nmds_crops", "low_occ_nmds_brc"
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



############
### PCoA ###
############
# List of transformed matrices

hell_matrices <- purrr::map(asv_matrices, ~ {
  decostand(t(.x), method = "hellinger", MARGIN = 1)
})

# Distances
distance_matrices <- purrr::map(hell_matrices ~
  vegdist(t(.x),
    method = "bray",
    upper = FALSE,
    binary = FALSE,
    na.rm = TRUE
  ))

distance_matrices <- hell_matrices %>%
  purrr::map(~ vegdist(t(.x), method = "bray", na.rm = TRUE),
    .progress = TRUE
  ) %>%
  set_names(names(hell_matrices))



core_asv_pcoa <- wcmdscale(core_asv_dist,
  k = 2,
  eig = TRUE,
  add = FALSE,
  x.ret = FALSE
)

# Extract PCoA scores and add metadata
core_pcoa_df <- data.frame(
  Dim1 = core_asv_pcoa$points[, 1], # First PCoA axis
  Dim2 = core_asv_pcoa$points[, 2], # Second PCoA axis
  brc = factor(core_brc_phyloseq@sam_data$brc), # Convert 'brc' to a factor
  crop = factor(core_brc_phyloseq@sam_data$crop)
)


core_pcoa_gg <- ggplot(data = core_pcoa_df, aes(x = Dim1, y = Dim2, color = brc)) +
  geom_point() +
  theme_bw() +
  labs(color = "BRC") +
  ggtitle("PCoA: 50 core ASVs in BRC crops",
    subtitle = "Core that contributes 2% to Bray-Curtis"
  )
