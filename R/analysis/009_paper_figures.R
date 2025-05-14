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
# BRC BRAY-CURTIS CORE AND ASSOCIATED
#--------------------------------------------------------

occ_abun_plot <- brc_occ_curve(core_summary_list = core_summary_lists)
occ_abun_plot

# Generate crop-based and BRC-based visualizations
bc_core_nmds_crops <- brc_gg_ordi(
  .data = bc_core_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle(
    "50 core ASVs in BRC crops"
    # subtitle = "Core that contributes 2% to Bray-Curtis"
  )

bc_core_nmds_brc <- brc_gg_ordi(
  .data = bc_core_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle(
    "50 core ASVs in BRC crops"
    # subtitle = "Core that contributes 2% to Bray-Curtis"
  )


bc_plots <- ggpubr::ggarrange(
  # Top row
  occ_abun_plot +
    ggtitle("ASVs contributing >2% to Bray-Curtis Dissimilarity") +
    guides(fill = guide_legend(title = "ASV Association")),

  # Bottom row
  ggpubr::ggarrange(
    bc_core_nmds_crops,
    bc_core_nmds_brc,
    ncol = 2,
    labels = c("B", "C"),
    common.legend = TRUE,
    legend = "bottom"
  ),
  # Arrangement parameters
  nrow = 2,
  labels = c("A", ""),
  heights = c(1.5, 1),
  common.legend = TRUE,
  legend = "bottom"
)


#--------------------------------------------------------
# BRC CORES SELECTED BY THRESHOLD
#--------------------------------------------------------

thres_core_60 <- filtered_phyloseq %>%
  brc_occore(.) %>%
  brc_occ_curve(core_summary_list = .)

# Generate crop-based and BRC-based visualizations
high_occ_nmds_crops <- brc_gg_ordi(
  .data = high_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("12 high prevalence ASVs in BRC crops (60% samples)")

high_occ_nmds_brc <- brc_gg_ordi(
  .data = high_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = brc,
  .drop_na = brc
) +
  ggtitle("12 high prevalence ASVs in BRCs (60% samples)")


threshold_plots <- ggpubr::ggarrange(
  # Top row
  thres_core_60 +
    geom_hline(yintercept = 0.6, linetype = 2, linewidth = 1, color = "red") +
    ggtitle("ASVs in >60% Samples Across Crops and BRCs") +
    guides(fill = guide_legend(title = "ASV Association")),
  # Second row
  ggpubr::ggarrange(
    high_occ_nmds_crops,
    high_occ_nmds_brc,
    ncol = 2,
    labels = c("E", "F"),
    common.legend = TRUE,
    legend = "bottom"
  ),

  # Overall arrangement settings
  nrow = 2,
  labels = c("D", ""),
  heights = c(1.5, 1),
  common.legend = TRUE,
  legend = "bottom"
)

# Add a common title if needed
final_plot <- ggpubr::ggarrange(
  bc_plots,
  threshold_plots,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"

  # top = ggpubr::text_grob(
  #   "ASVs in >60% Samples Across Crops and BRCs",
  #   face = "bold",
  #   size = 14
  # )
)

# Display the final arranged plot
final_plot

ggsave(
  "data/output/plots/combined_multicores_nmds.png",
  plot = final_plot,
  dpi = 300,
  width = 300,
  height = 250,
  units = "mm"
)
