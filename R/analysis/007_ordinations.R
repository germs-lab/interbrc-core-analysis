## Inter-BRC Core Analysis
## Core & non_core Ordinations
## By Bolívar Aponte Rolón

# Setup
source("R/utils/000_setup.R")
if (exists("phyloseq")) remove(phyloseq)

#################################################
###### Ordinations and Community Analysis #######
#################################################

#####################################
### Number of dimensions for NMDS ###
#####################################
# Transformed matrices
hell_matrices <- purrr::map(
  asv_matrices,
  ~ {
    decostand(t(.x), method = "hellinger", MARGIN = 1)
  }
)

# Distance matrices
plan(multicore, workers = parallel::detectCores() - 1) # Set up parallel backend

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

plan(sequential) # Close background work by changing plans


# set.seed(484035)
# nmds_screen_parallel(distance_matrices$core, ncores = 2)
# Results: Two dimensions keeps stress below 0.20
###############################################################

####################################
### NMDS: Core by extract_core() ###
####################################
core_ext_nmds <- brc_nmds(
  asv_matrix = asv_matrices$core, # brc_nmds() computes Hellinger transformations
  physeq = core_brc_phyloseq,
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
) # Best solution repeated 1 times (k = 3)

save(core_ext_nmds, file = "data/output/core_ext_nmds.rda")

# Core NMDS Aesthetics ####
core_nmds_crops <- brc_gg_ordi(
  .data = core_ext_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle(
    "50 core ASVs in BRC crops",
    subtitle = "Core that contributes 2% to Bray-Curtis"
  )

core_nmds_brc <- brc_gg_ordi(
  .data = core_ext_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle(
    "50 core ASVs in BRC crops",
    subtitle = "Core that contributes 2% to Bray-Curtis"
  )

########################################
### NMDS: Non-Core by extract_core() ###
########################################
non_core_ext_nmds <- brc_nmds(
  asv_matrix = asv_matrices$non_core,
  physeq = non_core_brc_phyloseq,
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
) # Best solution was not repeated (k = 3)

# save(non_core_ext_nmds, file = "data/output/non_core_ext_nmds.rda")

# Non-Core NMDS Aesthetics ####
non_core_nmds_crops <- brc_gg_ordi(
  .data = non_core_ext_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("Non-core ASVs in BRC crops (100% samples)")


non_core_nmds_brc <- brc_gg_ordi(
  .data = non_core_ext_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("Non-core ASVs in BRC (100% samples)")


############################################
### NMDS: Core and Non-Core by Threshold ###
############################################
# High occupancy
high_occ_nmds <- brc_nmds(
  asv_matrix = asv_matrices$high_occ,
  physeq = core_asvs_threshold$physeq_high_occ,
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
)

# save(high_occ_nmds, file = "data/output/high_occ_nmds.rda")

high_nmds_crops <- brc_gg_ordi(
  .data = high_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("12 high prevanlence ASVs in BRC crops (60% samples)")

high_nmds_brc <- brc_gg_ordi(
  .data = high_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("12 high prevalence ASVs in BRCs (60% samples)")


# Low occupancy
low_occ_nmds <- brc_nmds(
  asv_matrix = asv_matrices$low_occ,
  physeq = subset_samples(core_asvs_threshold$physeq_low_occ, brc != "jbei"), # Remove "jbei"
  ncores = parallel::detectCores() - 1,
  k = 3,
  trymax = 9999
)

save(low_occ_nmds, file = "data/output/low_occ_nmds.rda")


low_nmds_crops <- brc_gg_ordi(
  .data = low_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("Low prevalence ASVs in BRC crops (<60% samples)")

low_nmds_brc <- brc_gg_ordi(
  .data = low_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("Low prevalence ASVs in BRCs (<60% samples)")

########################
## Saving NMDS plots ###
########################

nmds_plots <- list(
  core_nmds_brc,
  core_nmds_crops,
  non_core_nmds_brc,
  non_core_nmds_crops,
  high_nmds_crops,
  high_nmds_brc,
  low_nmds_crops,
  low_nmds_brc
)
plot_names <- c(
  "core_nmds_brc",
  "core_nmds_crops",
  "non_core_nmds_brc",
  "non_core_nmds_crops",
  "high_occ_nmds_crops",
  "high_occ_nmds_brc",
  "low_nmds_crops",
  "low_nmds_brc"
)


plot_paths <- str_glue("data/output/plots/{tolower(plot_names)}.png")

purrr::walk2(
  plot_paths,
  nmds_plots,
  \(path, plot)
    ggsave(
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
####################################
### PCoA: Core by extract_core() ###
####################################
core_asv_pcoa <- brc_pcoa(
  distance_matrices$core,
  core_brc_phyloseq
) # Need a distance matric from vegdist()

core_pcoa_brc <- brc_gg_ordi(
  .data = core_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = brc,
  .drop_na = brc
) +
  ggtitle(
    "PCoA: 50 core ASVs in BRC",
    subtitle = "Core that contributes 2% to Bray-Curtis"
  )

core_pcoa_crops <- brc_gg_ordi(
  .data = core_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = crop,
  .drop_na = brc
) +
  ggtitle(
    "PCoA: 50 core ASVs in BRC crops",
    subtitle = "Core that contributes 2% to Bray-Curtis"
  )


########################################
### PCoA: Non-Core by extract_core() ###
########################################
non_core_asv_pcoa <- brc_pcoa(
  distance_matrices$non_core,
  non_core_brc_phyloseq
)

non_core_pcoa_brc <- brc_gg_ordi(
  .data = non_core_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("PCoA: Non-Core ASVs in BRCs")

non_core_pcoa_crops <- brc_gg_ordi(
  .data = non_core_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("PCoA: Non-Core ASVs in BRC crops")


############################################
### PCoA: Core and Non-Core by threshold ###
############################################

## High occupancy
high_asv_pcoa <- brc_pcoa(
  distance_matrices$high_occ,
  core_asvs_threshold$physeq_high_occ %>%
    prune_samples(sample_sums(.) > 0, .) # Making the matrix and phyloseq the same: 1732 samples
)

high_pcoa_brc <- brc_gg_ordi(
  .data = high_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("PCoA: High prevalence ASVs in BRCs")

high_pcoa_crops <- brc_gg_ordi(
  .data = high_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("PCoA: High prevalence ASVs in BRC crops")

## Low prevalence
low_asv_pcoa <- brc_pcoa(
  distance_matrices$low_occ,
  core_asvs_threshold$physeq_low_occ %>%
    prune_samples(sample_sums(.) > 0, .) # Making the matrix and phyloseq the same: 1732 samples
)

low_pcoa_brc <- brc_gg_ordi(
  .data = low_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("PCoA: Low prevalence ASVs in BRCs")

low_pcoa_crops <- brc_gg_ordi(
  .data = low_asv_pcoa$pcoa_df,
  ordi = "PCoA",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("PCoA: Low prevalence ASVs in BRC crops")


########################
## Saving PCoA plots ###
########################

pcoa_plots <- list(
  core_pcoa_brc,
  core_pcoa_crops,
  non_core_pcoa_brc,
  non_core_pcoa_crops,
  high_pcoa_crops,
  high_pcoa_brc,
  low_pcoa_crops,
  low_pcoa_brc
)
plot_names <- c(
  "core_pcoa_brc",
  "core_pcoa_crops",
  "non_core_pcoa_brc",
  "non_core_pcoa_crops",
  "high_pcoa_crops",
  "high_pcoa_brc",
  "low_pcoa_crops",
  "low_pcoa_brc"
)


plot_paths <- str_glue("data/output/plots/{tolower(plot_names)}.png")

purrr::walk2(
  plot_paths,
  pcoa_plots,
  \(path, plot)
    ggsave(
      filename = path,
      plot = plot,
      dpi = 300,
      width = 200,
      height = 200,
      units = "mm"
    )
)
