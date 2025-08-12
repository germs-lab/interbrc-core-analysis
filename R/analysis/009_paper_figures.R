#########################################################
# CORE MICROBIOME SELECTION
# Figures publish ready
# For the most part it is ordinations already explored in
# 007_ordinations.R but with different aesthetics and arrangement.
#
# We focus on PCoA for the publication
#
# Project:  Inter-BRC-Core-Microbiome
# Modified by: Bolívar Aponte Rolón
# Date: 2025-05-13
#########################################################

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(
  conflicted,
  styler,
  phyloseq,
  vegan,
  tidyverse,
  janitor,
  ggsci,
  minpack.lm,
  Hmisc,
  stats4,
  metagMisc,
  BRCore,
  furrr,
  microeco,
  file2meco,
  patchwork,
  ggh4x,
  ggnested
)

## List files and source each
list.files(here::here("R/functions"), pattern = "brc_", full.names = TRUE) %>%
  purrr::map(source)

# Load data
# We only need these files
load(here::here("data/output/bc_core_nmds.rda"))
load(here::here("data/output/bc_noncore_nmds.rda"))
load(here::here("data/output/high_occ_nmds.rda"))
load(here::here("data/output/low_occ_nmds.rda"))
load(here::here("data/output/pcoa_results_filtered.rda"))
load(here::here("data/output/nmds_result_filtered.rda"))
load(here::here("data/output/core_summary_lists_old.rda"))
load(here::here("data/output/distance_matrices.rda"))
load(here::here("data/output/phyloseq_objects/braycurt_core.rda"))
load(here::here("data/output/phyloseq_objects/braycurt_noncore.rda"))
load(here::here("data/output/phyloseq_objects/prevalence_core.rda"))
load(here::here("data/output/phyloseq_objects/filtered_phyloseq.rda"))

#--------------------------------------------------------
# GGTITLES & THEMES
#--------------------------------------------------------
title_size <- list(theme(title = element_text(size = 10), 
                         legend.title = element_text(size = 10, face = "bold")))
crop_labels <- c(
  "Corn",
  "Miscanthus",
  "Poplar",
  "Restored Prairie",
  "Sorghum",
  "Sorghum + Rye",
  "Soy",
  "Switchgrass"
)


scale_y_limits <- list(
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))
)

core_50_crops <- list(
  ggtitle("50 core ASVs in BRC crops"),
  scale_color_brewer(palette = "Set2"),
  title_size
)

core_50_brc <- list(
  ggtitle("50 core ASVs in BRCs"),
  title_size
)

threshold_60_core <- list(
  scale_y_limits,
  ggtitle("ASVs in >60% Samples Across Crops and BRCs"),
  guides(fill = guide_legend(title = "ASV Association")),
  title_size
)

threshold_60_crops <- list(
  scale_y_limits,
  ggtitle("12 high prevalence ASVs in BRC crops"),
  scale_color_brewer(palette = "Set2"),
  title_size
)

threshold_60_brc <- list(
  scale_y_limits,
  ggtitle("12 high prevalence ASVs in BRCs"),
  title_size
)

threshold_100_palette <- list(
  scale_y_limits,
  scale_color_brewer(
    palette = "Set2",
    labels = crop_labels
  ),
  title_size
)

override_scale <- list(
  scale_y_limits,
  scale_fill_manual(
    labels = c("Core", "Non-core"),
    values = c("#377EB8", "#377EB8") #"#E41A1C")
  )
)


#--------------------------------------------------------
# Pre-processing
#--------------------------------------------------------

# Perform PCoA on BC_CORE COMMUNITY
bc_core_asv_pcoa <- brc_pcoa(
  distance_matrices$bc_core,
  braycurt_core
)

# Perform PCoA on BC_NONCORE COMMUNITY
bc_noncore_asv_pcoa <- brc_pcoa(
  distance_matrices$bc_noncore,
  braycurt_noncore
)

# Perform PCoA on high-occupancy (threshold-based core) community
high_asv_pcoa <- brc_pcoa(
  distance_matrices$high_occ,
  prevalence_core$physeq_high_occ %>%
    prune_samples(sample_sums(.) > 0, .)
)


#--------------------------------------------------------
# NMDS & PCoA: BRC BRAY-CURTIS CORE AND ASSOCIATED
#--------------------------------------------------------

# Abundance/Occupancy curves
core_abun_plot <- brc_occ_curve(
  core_summary_list = core_summary_lists,
  color_values = NULL
) +
  override_scale +
  #geom_hline(yintercept = 0.1688, linetype = 2, linewidth = 1, color = "red") +
  ggtitle("ASVs contributing >2% to Bray-Curtis Dissimilarity") +
  title_size +
  guides(fill = "none") #guide_legend(title = "ASV Association"))

#--------------------
# NMDS & PCoA plots
#--------------------

core_nmds_pcoa_plots <- brc_paper_ordinations(
  nmds_data = bc_core_nmds$nmds_df,
  pcoa_data = bc_core_asv_pcoa$pcoa_df,
  nmds_colors = c("#E7B800", "#FC4E07"),
  pcoa_colors = c("#E7B800", "#383961", "#FC4E07"),
  nmds_labels = c("CABBI", "GLBRC"),
  pcoa_labels = c("CABBI", "CBI", "GLBRC"),
  nmds_crop_labels = crop_labels,
  pcoa_crop_labels = crop_labels,
  point_size = 1.25,
  crop_theme = core_50_crops,
  brc_theme = core_50_brc
)

# core_nmds_pcoa_plots$nmds$crops
# core_nmds_pcoa_plots$nmds$brc
# core_nmds_pcoa_plots$pcoa$crops
# core_nmds_pcoa_plots$pcoa$brc

#--------------------------------------------------------
# NMDS & PCoA: BRC CORES SELECTED BY THRESHOLD
#--------------------------------------------------------
thres_core_60_curve <- filtered_phyloseq %>%
  brc_occore(.) %>%
  brc_occ_curve(
    core_summary_list = .,
    color_values = NULL
  ) +
  override_scale +
  geom_hline(yintercept = 0.6, linetype = 2, linewidth = 1, color = "red") +
  threshold_60_core


high_nmds_pcoa_plots <- brc_paper_ordinations(
  nmds_data = high_occ_nmds$nmds_df,
  pcoa_data = high_asv_pcoa$pcoa_df,
  nmds_colors = c("#E7B800", "#FC4E07"),
  pcoa_colors = c("#E7B800", "#FC4E07"),
  pcoa_labels = c("CABBI", "GLBRC"),
  nmds_crop_labels = crop_labels,
  pcoa_crop_labels = crop_labels,
  point_size = 1.25,
  crop_theme = threshold_60_crops,
  brc_theme = threshold_60_brc
)


# high_nmds_pcoa_plots$nmds$crops
# high_nmds_pcoa_plots$nmds$brc
# high_nmds_pcoa_plots$pcoa$crops
# high_nmds_pcoa_plots$pcoa$brc

#--------------------------------------------------------
# NMDS & PCoA: FULL ASV COMMUNITY
#--------------------------------------------------------

thres_core_100_curve <- filtered_phyloseq %>%
  brc_occore(., threshold = 0) %>%
  brc_occ_curve(
    core_summary_list = .,
    color_values = NULL
  ) +
  scale_y_limits +
  scale_fill_manual(
    labels = c("ASVs"),
    values = c("#377EB8")
  ) +
  ggtitle("All ASVs Across Crops and BRCs") +
  title_size +
  guides(fill = guide_legend(title = "Individual ASVs"))


all_asv_nmds_pcoa_plots <- brc_paper_ordinations(
  nmds_data = nmds_result_filtered$nmds_df,
  pcoa_data = pcoa_results_filtered$sample_reads_500_asv_reads_20$pcoa_df,
  nmds_colors = c("#E7B800", "#383961", "#FC4E07", "#1E692D"),
  pcoa_colors = c("#E7B800", "#383961", "#FC4E07", "#1E692D"),
  nmds_labels = c("CABBI", "CBI", "GLBRC", "JBEI"),
  pcoa_labels = c("CABBI", "CBI", "GLBRC", "JBEI"),
  point_size = 1.25,
  crop_theme = list(
    ggtitle("All ASVs across BRC crops"),
    threshold_100_palette
  ),
  brc_theme = list(ggtitle("All ASVs across BRCs"), title_size)
)

# all_asv_nmds_pcoa_plots$pcoa$crops

#--------------------------------------------------------
# ARRANGEMENT
#--------------------------------------------------------

combined_multicores <- ggpubr::ggarrange(
  # First row - 2 plots
  ggpubr::ggarrange(
    thres_core_100_curve,
    thres_core_60_curve,
    core_abun_plot,
    ncol = 3,
    labels = c("A", "D", "G"),
    common.legend = TRUE,
    legend = "bottom"
  ),
  " ",
  # Second row - 4 plots
  ggpubr::ggarrange(
    all_asv_nmds_pcoa_plots$pcoa$crops,
    all_asv_nmds_pcoa_plots$pcoa$brc,
    high_nmds_pcoa_plots$pcoa$crops,
    high_nmds_pcoa_plots$pcoa$brc,
    core_nmds_pcoa_plots$pcoa$crops,
    core_nmds_pcoa_plots$pcoa$brc,
    ncol = 6,
    labels = c("B", "C", "E", "F", "H", "I"),
    common.legend = TRUE,
    legend = "bottom"
  ),

  nrow = 3,
  heights = c(1.5, 0.1, 1),
  common.legend = FALSE,
  legend = "bottom"
)

combined_multicores

## Need final legend arrangement in inkscape

ggsave(
  "data/output/plots/combined_multicores_pcoa.svg",
  plot = combined_multicores,
  dpi = 300,
  width = 375,
  height = 250,
  units = "mm"
)

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
  ) +
    title_size

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


#-------------------------------------------------------------
# Combined Relative Abundance  and ABundance/Occupancy Curves
#-------------------------------------------------------------


box_2 <- ggpubr::ggarrange(
    ggpubr::ggarrange(
        rel_abund_vertical,
        ncol = 1,
        labels = c("A"),
        common.legend = TRUE,
        legend = "bottom"
    ),
    " ",
    # Second row - 4 plots
    ggpubr::ggarrange(
        thres_core_100_curve,
        core_abun_plot,
        ncol = 2,
        labels = c("B", "C"),
        common.legend = TRUE,
        legend = "bottom"
    ),
    
    nrow = 3,
    heights = c(1.5, 0.1, 1),
    common.legend = FALSE,
    legend = "bottom"
)


 box_2

 ggsave(
     "data/output/plots/box_2.png",
     plot = box_2,
     dpi = 300,
     width = 325,
     height = 275,
     units = "mm"
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
