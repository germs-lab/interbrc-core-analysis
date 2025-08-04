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
if (exists("phyloseq")) {
  remove(phyloseq)
}

library(microeco)
library(file2meco)
library(patchwork)
library(ggh4x)
library(ggplot2)
#pak::pkg_install("gmteunisse/ggnested")
library(ggnested)

#--------------------------------------------------------
# BRC BRAY-CURTIS CORE AND ASSOCIATED
#--------------------------------------------------------

occ_abun_plot <- brc_occ_curve(
  core_summary_list = core_summary_lists,
  color_values = RColorBrewer::brewer.pal(3, name = "Set1")
) +
  ggtitle("ASVs contributing >2% to Bray-Curtis Dissimilarity") +
  guides(fill = guide_legend(title = "ASV Association"))

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
  ) +
  scale_color_brewer(palette = "Set2") +
  theme(title = element_text(size = 8))


bc_core_nmds_brc <- brc_gg_ordi(
  .data = bc_core_nmds$nmds_df,
  ordi = "NMDS",
  .shape = brc,
  .color = brc,
  color_values = c("#FC4E07", "#E7B800"),
  .labels = c("GLBRC", "CABBI"),
  .size = 2.5,
  .alpha = 0.5,
  .drop_na = brc
) +
  ggtitle(
    "50 core ASVs in BRCs"
    # subtitle = "Core that contributes 2% to Bray-Curtis"
  ) +
  theme(title = element_text(size = 8))

# bc_plots <- ggpubr::ggarrange(
#   # Top row
#   occ_abun_plot,

#   # Bottom row
#   ggpubr::ggarrange(
#     bc_core_nmds_crops,
#     bc_core_nmds_brc,
#     ncol = 2,
#     labels = c("B", "C"),
#     common.legend = TRUE,
#     legend = "bottom"
#   ),
#   # Arrangement parameters
#   nrow = 2,
#   labels = c("A", ""),
#   heights = c(1.5, 1),
#   common.legend = TRUE,
#   legend = "bottom"
# )

#--------------------------------------------------------
# BRC CORES SELECTED BY THRESHOLD
#--------------------------------------------------------

thres_core_60 <- filtered_phyloseq %>%
  brc_occore(.) %>%
  brc_occ_curve(
    core_summary_list = .,
    color_values = RColorBrewer::brewer.pal(3, name = "Set1")
  ) +
  geom_hline(yintercept = 0.6, linetype = 2, linewidth = 1, color = "red") +
  ggtitle("ASVs in >60% Samples Across Crops and BRCs") +
  guides(fill = guide_legend(title = "ASV Association"))

# Generate crop-based and BRC-based visualizations
high_occ_nmds_crops <- brc_gg_ordi(
  .data = high_occ_nmds$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) +
  ggtitle("12 high prevalence ASVs in BRC crops") +
  scale_color_brewer(palette = "Set2") +
  theme(title = element_text(size = 8))


high_occ_nmds_brc <- brc_gg_ordi(
  .data = high_occ_nmds$nmds_df,
  ordi = "NMDS",
  .shape = brc,
  .color = brc,
  color_values = c("#FC4E07", "#E7B800"),
  .labels = c("GLBRC", "CABBI"),
  .size = 2.5,
  .alpha = 0.5,
  .drop_na = brc
) +
  ggtitle("12 high prevalence ASVs in BRCs") +
  theme(title = element_text(size = 8))


# threshold_plots <- ggpubr::ggarrange(
#   # Top row
#   thres_core_60,
#   # Second row
#   ggpubr::ggarrange(
#     high_occ_nmds_crops,
#     high_occ_nmds_brc,
#     ncol = 2,
#     labels = c("E", "F"),
#     common.legend = FALSE,
#     legend = "bottom"
#   ),

#   # Overall arrangement settings
#   nrow = 2,
#   labels = c("D", ""),
#   heights = c(1.5, 1),
#   common.legend = TRUE,
#   legend = "bottom"
# )

final_plot <- ggpubr::ggarrange(
  # First row - 2 plots
  ggpubr::ggarrange(
    occ_abun_plot,
    thres_core_60,
    ncol = 2,
    labels = c("A", "D"),
    common.legend = TRUE,
    legend = "bottom"
  ),
  " ",
  # Second row - 4 plots
  ggpubr::ggarrange(
    bc_core_nmds_crops,
    bc_core_nmds_brc,
    high_occ_nmds_crops,
    high_occ_nmds_brc,
    ncol = 4,
    labels = c("B", "C", "E", "F"),
    common.legend = TRUE,
    legend = "bottom"
  ),

  nrow = 3,
  heights = c(1.5, 0.1, 1),
  common.legend = TRUE,
  legend = "bottom"
)

# Display the final arranged plot
final_plot

## Need final legend arrangement in inkscape

ggsave(
  "data/output/plots/combined_multicores_nmds.svg",
  plot = final_plot,
  dpi = 300,
  width = 300,
  height = 250,
  units = "mm"
)


#---------------------------------------------------
#

#---------------------------------------------------
# Relative abundance visualization per BRC and Crop
#---------------------------------------------------
# Due to the size of the dataset, the following code needs to be executed in
# machine with >32GB RAM.

# First convert phyloseq to microtable object
micro_obj <- file2meco::phyloseq2meco(filtered_phyloseq)


# Relative Abundance
abund_obj <- trans_abund$new(
  micro_obj,
  taxrank = "phylum",
  #high_level = "phylum",
  #high_level_fix_nsub = 3,
  delete_taxonomy_prefix = TRUE
)
rel_abund <- abund_obj$plot_bar(
  others_color = "grey70",
  facet = c("brc", "crop"),
  #ggnested = TRUE,
  #high_level_add_other = TRUE,
  barwidth = 1,
  xtext_keep = FALSE,
  legend_text_italic = FALSE,
)

# Fixing labels
new_labels <- as_labeller(c(
  `cabbi` = "CABBI",
  `cbi` = "CBI",
  `jbei` = "JBEI",
  `glbrc` = "GLBRC",
  `Sorghum` = "Sorghum",
  `poplar` = "Poplar",
  `Restored_Prairie` = "Restored Prairie",
  `Switchgrass` = "Switchgrass",
  `Corn` = "Corn",
  `Miscanthus` = "Miscanthus",
  `Soy` = "Soy",
  `Sorghum + Rye` = "Sorghum + Rye"
))

# Vertical faceted nested plot
rel_abund_vertical <- rel_abund +
  guides(fill = guide_legend(title = "Phylum")) +
  facet_nested(
    ~ brc + crop,
    labeller = new_labels,
    scales = "free_x",
    strip = strip_nested(
      background_x = elem_list_rect(fill = "#F0F7E6"),
      text_x = elem_list_text(face = "bold")
    )
  )

# Horizontal faceted nested rapped plot
rel_abund_horizon <- rel_abund +
  guides(fill = guide_legend(title = "Phylum")) +
  facet_nested_wrap(
    ~ brc + crop,
    labeller = new_labels,
    scales = "free_x",
    strip = strip_nested(
      background_x = elem_list_rect(fill = "#F0F7E6"),
      text_x = elem_list_text(face = "bold")
    )
  )

ggsave(
  "data/output/plots/brc_crop_relabund_vertical.png",
  plot = rel_abund_vertical,
  dpi = 300,
  width = 350,
  height = 250,
  units = "mm"
)

ggsave(
  "data/output/plots/brc_crop_relabund_horizon.png",
  plot = rel_abund_horizon,
  dpi = 300,
  width = 300,
  height = 300,
  units = "mm"
)


#----------------------------------------------
# Full data set visualizations
#----------------------------------------------
# Based on "007_ordinations_full.R" results
# sample_reads_500_asv_reads_20
# 1809 sample 59951 ASVs

test <- brc_gg_ordi(
  .data = nmds_result_filtered$nmds_df,
  ordi = "NMDS",
  .color = crop,
  .drop_na = brc
) 
test <- brc_gg_ordi(
  .data = pcoa_results_filtered$sample_reads_500_asv_reads_20$pcoa_df,
  ordi = "PCoA",
  .color = crop,
  .drop_na = brc
) 
#-------------
# Maybe useful
#-------------

# Create the trans_venn object for visualization
# venn_obj <- trans_venn$new(micro_obj)
#
# # Generate Venn diagrams for phyla
# venn_obj$plot_venn(use_taxonomy = TRUE, taxonomy_rank = "phylum",
#                    fill_alpha = 0.5,
#                    plot_text_size = 3)
# phyla_venn <- venn_obj$result_plot

# # Alpha diversity
# alpha_obj <- trans_alpha$new(micro_obj, group = "crop")
#
# # Plot alpha diversity for phyla
# alpha_obj$plot_alpha(measure = "Chao1",
#                      add = "jitter",
#                      add.params = list(alpha = 0.3),
#                      y_start = 0.1,
#                      y_increase = 0.1,
#                      add_stat = TRUE)
#
# # Create a trans_upset object for UpSet plots
# upset_obj <- trans_upset$new(micro_obj, group = "crop")
