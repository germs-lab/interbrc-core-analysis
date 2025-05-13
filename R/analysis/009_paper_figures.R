#########################################################
# CORE MICROBIOME SELECTION
# Figures publish ready
#
# Project:  Inter-BRC-Core-Microbiome
# Modified by: Bolívar Aponte Rolón
# Date: 2025-05-13
#########################################################

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
source("R/utils/000_setup.R")
if (exists("phyloseq")) remove(phyloseq)


#--------------------------------------------------------
# All BRC Bray-Curtis core and associated
#--------------------------------------------------------

occ_abun_plot <- brc_occ_curve(core_summary_list = core_summary_lists)
occ_abun_plot

# Generate crop-based and BRC-based visualizations
core_nmds_crops <- brc_gg_ordi(
  .data = core_ext_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle(
    "50 core ASVs in BRC crops"
    #subtitle = "Core that contributes 2% to Bray-Curtis"
  )

core_nmds_brc <- brc_gg_ordi(
  .data = core_ext_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle(
    "50 core ASVs in BRCs"
    #subtitle = "Core that contributes 2% to Bray-Curtis"
  )


bc_associated <- ggpubr::ggarrange(
  # Top row
  occ_abun_plot +
    ggtitle("ASVs that contribute at least 2% to Bray-Curtis Dissimilarity") +
    guides(fill = guide_legend(title = "ASV Association")),

  # Bottom row
  ggpubr::ggarrange(
    core_nmds_crops,
    core_nmds_brc,
    ncol = 2,
    labels = c("B", "C")
  ),
  # Arrangement parameters
  nrow = 2,
  labels = c("A", ""),
  heights = c(1.5, 1),
  common.legend = FALSE,
  legend = "right"
)

#--------------------------------------------------------
# All BRCs Threshold Cores
#--------------------------------------------------------

thres_core_60 <- filtered_phyloseq %>%
  abunocc_core(.) %>%
  brc_occ_curve(core_summary_list = .)

# Generate crop-based and BRC-based visualizations
high_nmds_crops <- brc_gg_ordi(
  .data = high_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("12 high prevalence ASVs in BRC crops")

high_nmds_brc <- brc_gg_ordi(
  .data = high_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("12 high prevalence ASVs in BRCs")


ggpubr::ggarrange(
  # Top row
  thres_core_60 +
    geom_hline(yintercept = 0.6, linetype = 2, linewidth = 1, color = "red") +
    ggtitle("ASVs in >60% Samples Across Crops and BRCs") +
    guides(fill = guide_legend(title = "ASV Association")),
  # Second row
  ggpubr::ggarrange(
    high_nmds_crops,
    high_nmds_brc,
    ncol = 2,
    labels = c("B", "C")
  ),

  # Overall arrangement settings
  nrow = 2,
  heights = c(1, 1),
  labels = c("A", ""),
  common.legend = FALSE
)

# Add a common title if needed
final_plot <- ggpubr::annotate_figure(
  arranged_plots,
  top = ggpubr::text_grob(
    "Core ASVs Across Crops and BRCs",
    face = "bold",
    size = 14
  )
)

# Display the final arranged plot
final_plot
