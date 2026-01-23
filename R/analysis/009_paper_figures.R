#########################################################
# CORE MICROBIOME SELECTION - REFACTORED
# Figures publish ready
# For the most part it is ordinations already explored in
# 007_ordinations.R but with different aesthetics and arrangement.
#
# We focus on PCoA for the publication
#
# Project:  Inter-BRC-Core-Microbiome
# Author: Bolívar Aponte Rolón
# Date: 2025-05-13
# Last modified: 2025-01-20 (Refactored)
#
# REFACTORING NOTES:
# - Fixed syntax errors on lines 408-411 and 486-489
# - Consolidated theme definitions into brc_themes.R
# - Uses new unified ordination interface (optional)
# - Cleaner, more maintainable code structure
#########################################################

# SETUP AND DEPENDENCIES ----

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
load(here::here("data/output/bc_core_nmds_007.rda"))
load(here::here("data/output/bc_noncore_nmds_007.rda"))
load(here::here("data/output/high_occ_nmds_007.rda"))
load(here::here("data/output/low_occ_nmds_007.rda"))
load(here::here("data/output/pcoa_results_007.rda"))
load(here::here("data/output/core_summary_lists_old.rda"))
load(here::here("data/output/distance_matrices_007.rda"))
load(here::here("data/output/phyloseq_objects/braycurt_core_006.rda"))
load(here::here("data/output/phyloseq_objects/braycurt_noncore_006.rda"))
load(here::here("data/output/phyloseq_objects/prevalence_core_004.rda"))
load(here::here("data/output/phyloseq_objects/filtered_phyloseq_002.rda"))


# THEMES AND LABELS (now using brc_themes.R) ----
# Note: Most theme definitions are now in brc_themes.R
# We keep only plot-specific customizations here

# Get standard labels
crop_labels <- brc_crop_labels()

# Define plot-specific theme lists using new helpers
core_50_crops <- brc_theme_core50_crops()
core_50_brc <- brc_theme_core50_brc()
threshold_60_core <- brc_theme_threshold60_core()
threshold_60_crops <- brc_theme_threshold60_crops()
threshold_60_brc <- brc_theme_threshold60_brc()
threshold_100_palette <- brc_theme_threshold100_palette()
override_scale <- brc_scale_override_single()


# NMDS & PCoA: BRC BRAY-CURTIS CORE AND ASSOCIATED ----

# Abundance/Occupancy curves
core_abun_plot <- brc_occ_curve(
  core_summary_list = core_summary_lists,
  color_values = NULL
) +
  override_scale +
  geom_hline(yintercept = 0.1688, linetype = 2, linewidth = 1, color = "red") +
  ggtitle("ASVs contributing >2% to Bray-Curtis Dissimilarity") +
  brc_theme_title_size() +
  labs(y = "Occupancy (%)") +
  guides(fill = "none") +
  brc_scale_fill_core()

# NMDS & PCoA plots ----

core_nmds_pcoa_plots <- brc_paper_ordinations(
  nmds_data = bc_core_nmds$nmds_df,
  pcoa_data = pcoa_results$bc_core_asv_pcoa$pcoa_df,
  nmds_colors = c("#E7B800", "#FC4E07"),
  pcoa_colors = c("#E7B800", "#383961", "#FC4E07"),
  nmds_labels = c("CABBI", "GLBRC"),
  pcoa_labels = c("CABBI", "CBI", "GLBRC"),
  point_size = 1.25,
  crop_theme = core_50_crops,
  brc_theme = core_50_brc
)

# Individual plots accessible via:
# core_nmds_pcoa_plots$nmds$crops
# core_nmds_pcoa_plots$nmds$brc
# core_nmds_pcoa_plots$pcoa$crops
# core_nmds_pcoa_plots$pcoa$brc

# NMDS & PCoA: BRC CORES SELECTED BY THRESHOLD ----

thres_core_60_curve <- filtered_phyloseq %>%
  brc_occore(.) %>%
  brc_occ_curve(
    core_summary_list = .,
    color_values = NULL
  ) +
  override_scale +
  labs(y = "Occupancy (%)") +
  geom_hline(yintercept = 0.6, linetype = 2, linewidth = 1, color = "red") +
  threshold_60_core +
  brc_scale_fill_core()


high_nmds_pcoa_plots <- brc_paper_ordinations(
  nmds_data = high_occ_nmds$nmds_df,
  pcoa_data = pcoa_results$high_asv_pcoa$pcoa_df,
  nmds_colors = c("#E7B800", "#FC4E07"),
  pcoa_colors = c("#E7B800", "#FC4E07"),
  pcoa_labels = c("CABBI", "GLBRC"),
  point_size = 1.25,
  crop_theme = threshold_60_crops,
  brc_theme = threshold_60_brc
)

# Individual plots accessible via:
# high_nmds_pcoa_plots$nmds$crops
# high_nmds_pcoa_plots$nmds$brc
# high_nmds_pcoa_plots$pcoa$crops
# high_nmds_pcoa_plots$pcoa$brc

# NMDS & PCoA: FULL ASV COMMUNITY ----

thres_core_100_curve <- filtered_phyloseq %>%
  brc_occore(., threshold = 0) %>%
  brc_occ_curve(
    core_summary_list = .,
    color_values = NULL
  ) +
  brc_scale_y_percent() +
  scale_fill_manual(
    labels = c("ASVs"),
    values = c("#377EB8")
  ) +
  ggtitle("All ASVs Across Crops and BRCs") +
  brc_theme_title_size() +
  labs(y = "Occupancy (%)") +
  guides(fill = "none")


all_asv_nmds_pcoa_plots <- brc_paper_ordinations(
  nmds_data = high_occ_nmds$nmds_df, #PLACEHOLDER - REPLACE WITH ALL ASV NMDS
  pcoa_data = pcoa_results$all_asv_pcoa$pcoa_df,
  nmds_colors = c("#E7B800", "#383961", "#FC4E07", "#1E692D"),
  pcoa_colors = c("#E7B800", "#383961", "#FC4E07", "#1E692D"),
  nmds_labels = c("CABBI", "CBI", "GLBRC", "JBEI"),
  pcoa_labels = c("CABBI", "CBI", "GLBRC", "JBEI"),
  point_size = 1.25,
  crop_theme = list(
    ggtitle("All ASVs across BRC crops"),
    threshold_100_palette
  ),
  brc_theme = list(
    ggtitle("All ASVs across BRCs"),
    brc_theme_title_size()
  )
)

# Individual plot example:
all_asv_nmds_pcoa_plots$pcoa$crops

# ARRANGEMENT ----

combined_multicores <- ggpubr::ggarrange(
  # First row - 3 plots
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
  # Second row - 6 plots
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

# Relative abundance visualization per BRC and Crop ----
# Due to the size of the dataset, the following code needs to be executed in
# machine with >32GB RAM.

# First convert phyloseq to microtable object
micro_obj <- file2meco::phyloseq2meco(filtered_phyloseq)


# Relative Abundance
abund_obj <- trans_abund$new(
  micro_obj,
  taxrank = "phylum",
  delete_taxonomy_prefix = TRUE
)

rel_abund <- abund_obj$plot_bar(
  others_color = "grey70",
  facet = c("brc", "crop"),
  barwidth = 1,
  xtext_keep = FALSE,
  legend_text_italic = FALSE,
)

# Fixing labels using helper function
new_labels <- ggplot2::as_labeller(brc_full_labels())

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
  ggtitle("Relative Abundance of Bacterial Taxa across Samples") +
  brc_theme_title_size() +
  theme(axis.text.x = element_blank())

# Horizontal faceted nested wrapped plot
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


# Clustering threshold analysis by Jae ----
raw <- read_csv(here::here("data/output/sabr/seq-percentage_count_data.csv"))

if (is.na(names(raw)[1]) || names(raw)[1] == "") {
  names(raw)[1] <- "Clustering"
}

# Convert to long format
long_df <- raw %>%
  pivot_longer(
    cols = -1,
    names_to = "Occurrence",
    values_to = "Proportion"
  ) %>%
  mutate(
    Occurrence = gsub("%", "", Occurrence),
    Occurrence = as.numeric(Occurrence),
    Clustering = as.character(.[[1]])
  )

# Define axis and legend order
occ_levels <- c(50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)
clust_levels <- c(85, 87, 90, 93, 95, 97, 100)

plot_df <- long_df %>%
  mutate(
    Occurrence = factor(
      Occurrence,
      levels = occ_levels,
      ordered = TRUE,
      labels = paste0(occ_levels, "%")
    ),
    Clustering = factor(
      Clustering,
      levels = paste0(clust_levels, "%"),
      ordered = TRUE
    )
  )

# Color designation
custom_colors <- c(
  "85%" = "#FF68A1",
  "87%" = "#E68613",
  "90%" = "#0CB702",
  "93%" = "#00A9FF",
  "95%" = "#C77CFF",
  "97%" = "#999999",
  "100%" = "#ABA300"
)

cluster_plot <- ggplot(
  plot_df,
  aes(
    x = Occurrence,
    y = Proportion,
    color = Clustering,
    group = Clustering
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Core OTU Sequence Proportions across Clustering and Occurrence Thresholds",
    x = "Occurrence Threshold (%)",
    y = "Proportion of \nSequences (%)",
    color = "Clustering\nThreshold"
  ) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  scale_x_discrete(labels = function(x) gsub("%", "", x)) +
  theme_bw(base_size = 13) +
  brc_theme_title_size()

cluster_plot

# Combined Relative Abundance, Abundance/Occupancy Curves ----
# & Clustering threshold analysis
box_2 <- ggpubr::ggarrange(
  # First row - Relative abundance
  ggpubr::ggarrange(
    rel_abund_vertical,
    ncol = 1,
    labels = c("A"),
    common.legend = TRUE,
    legend = "bottom"
  ),
  " ",
  # Second row - Occupancy curves
  ggpubr::ggarrange(
    thres_core_100_curve,
    core_abun_plot,
    ncol = 2,
    labels = c("B", "C"),
    common.legend = TRUE,
    legend = "bottom"
  ),
  " ",
  # Third row - Clustering analysis
  ggpubr::ggarrange(
    cluster_plot,
    ncol = 1,
    labels = c("D"),
    common.legend = TRUE,
    legend = "bottom"
  ),

  nrow = 5,
  heights = c(1.6, 0.1, 1.2, 0.1, 1.6),
  common.legend = FALSE,
  legend = "bottom"
)

box_2

ggsave(
  "data/output/plots/box_2.png",
  plot = box_2,
  dpi = 300,
  width = 325,
  height = 300,
  units = "mm",
  bg = "white"
)

#-------------------------------------------------------
# Additional visualizations  ----
#-------------------------------------------------------

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
