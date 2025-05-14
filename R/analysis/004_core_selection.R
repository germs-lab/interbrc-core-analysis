#########################################################
# CORE MICROBIOME SELECTION
# Extract and analyze core microbiome across samples
#
# Project:  Inter-BRC-Core-Microbiome
# Original Author: Brandon Kristy
# Modified by: Bolívar Aponte Rolón
# Date: 2025-02-20
#########################################################

# DESCRIPTION:
# This script identifies the core microbiome across samples using multiple approaches:
# 1. Extract_core() method based on Bray-Curtis dissimilarity
# 2. Threshold-based approach
# 3. Neutral model fitting for abundance-occupancy patterns

#--------------------------------------------------------
# SETUP AND DEPENDENCIES
#--------------------------------------------------------
source("R/utils/000_setup.R")
if (exists("phyloseq")) remove(phyloseq)

#--------------------------------------------------------
# CORE EXTRACTION USING EXTRACT_CORE()
#--------------------------------------------------------
# Ensure minimum sample quality
filtered_phyloseq <- prune_samples(
  sample_sums(filtered_phyloseq) >= 100,
  filtered_phyloseq
)

# Extract core microbiome across all sites (with minimum 2% increase in Bray-Curtis)
braycore_summary <- extract_core(
  physeq,
  Var = "site",
  method = "increase",
  increase_value = 2
)

# Minimum seq depth was ~10,000 reads.

# Save results to avoid recomputation
# save(braycore_summary, file = "data/output/braycore_summary.rda")

#--------------------------------------------------------
# VISUALIZATION OF BRAY-CURTIS AND OCCUPANCY PATTERNS
#--------------------------------------------------------
# Generate Bray-Curtis dissimilarity curve
bray_curtis_curve <- brc_bc_curve(core_summary_list = braycore_summary)


# Generate abundance-occupancy plot
occ_abun_plot <- brc_occ_curve(core_summary_list = braycore_summary)
occ_abun_plot


# Save abundance-occupancy plot
ggsave(
  filename = "bray_curtis_abundance_occupancy.png",
  occ_abun_plot,
  path = "data/output/plots/",
  dpi = 300,
  width = 6,
  height = 4
)

# #--------------------------------------------------------
# # NEUTRAL MODEL FITTING FOR ABUNDANCE-OCCUPANCY
# #--------------------------------------------------------
# # Extract data for model fitting
# taxon <- braycore_summary[[7]]
# spp <- t(braycore_summary[[5]])
# occ_abun <- braycore_summary[[4]]
# names(occ_abun)[names(occ_abun) == "otu"] <- "OTU_ID"

# # Fit neutral community model
# obs.np <- sncm.fit2(spp, taxon, stats = FALSE, pool = NULL)
# sta.np <- sncm.fit2(spp, taxon, stats = TRUE, pool = NULL)

# # Classify taxa based on model predictions
# obs.np$fit_class <- "As predicted"
# obs.np[which(obs.np$freq < obs.np$pred.lwr), "fit_class"] <- "Below prediction"
# obs.np[which(obs.np$freq > obs.np$pred.upr), "fit_class"] <- "Above prediction"
# obs.np[which(is.na(obs.np$freq)), "fit_class"] <- "NA"

# # Combine data for visualization and analysis
# obs.np <- tibble::rownames_to_column(obs.np, "OTU_ID")
# as.data.frame(dplyr::left_join(occ_abun, obs.np, by = "OTU_ID")) -> fit_table

# # Calculate model statistics
# sta.np$above.pred <-
#   sum(obs.np$freq > (obs.np$pred.upr), na.rm = TRUE) / sta.np$Richness
# sta.np$below.pred <-
#   sum(obs.np$freq < (obs.np$pred.lwr), na.rm = TRUE) / sta.np$Richness
# fit_res <- as.data.frame(sta.np)
# rownames(fit_res) <- "Value"

# fit_res
# list_tab <- list(fit_res, fit_table)

# #--------------------------------------------------------
# # NEUTRAL MODEL VISUALIZATION
# #--------------------------------------------------------
# # Prepare data for plotting
# obs1 <- as.data.frame(list_tab[[2]])
# obs1 <- obs1[!is.na(obs1$p), ]
# obs2 <- as.data.frame(list_tab[[1]])

# # Add categories for visualization
# obs1 <- obs1 %>%
#   mutate(fill_fit_class = paste0(fill, ":", fit_class)) %>%
#   mutate(
#     fill_fit_class = case_when(
#       str_detect(fill, "no") ~ "Non Core Taxa",
#       TRUE ~ fill_fit_class
#     )
#   )

# # Plot neutral model fit
# obs1 %>%
#   ggplot(aes(x = log10(otu_rel), y = otu_occ)) +
#   scale_fill_npg(
#     name = "Core Membership: Model Predictions",
#     labels = c(
#       "Core: Above Prediction",
#       "Core: As Predicted",
#       "Core: Below Prediction",
#       "Non-Core Taxa"
#     )
#   ) +
#   geom_point(
#     aes(fill = fill_fit_class),
#     pch = 21,
#     alpha = 0.75,
#     size = 2.2
#   ) +
#   geom_line(
#     color = "red",
#     data = obs1,
#     size = 0.8,
#     aes(y = obs1$freq.pred, x = log10(obs1$p)),
#     alpha = 0.55
#   ) +
#   geom_line(
#     color = "black",
#     lty = "twodash",
#     size = 0.9,
#     data = obs1,
#     aes(y = obs1$pred.upr, x = log10(obs1$p)),
#     alpha = 0.55
#   ) +
#   geom_line(
#     color = "black",
#     lty = "twodash",
#     size = 0.9,
#     data = obs1,
#     aes(y = obs1$pred.lwr, x = log10(obs1$p)),
#     alpha = 0.55
#   ) +
#   labs(x = "Log10(mean abundance)", y = "Occupancy") +
#   annotate(
#     "text",
#     -Inf,
#     Inf,
#     label = paste("italic(R)^2 ==", round(obs2$Rsqr, 3)),
#     parse = TRUE,
#     size = 4.8,
#     hjust = -0.2,
#     vjust = 1.2
#   ) +
#   annotate(
#     "text",
#     -Inf,
#     Inf,
#     label = paste("italic(m) ==", round(obs2$m, 3)),
#     parse = TRUE,
#     size = 4.8,
#     hjust = -0.2,
#     vjust = 3.2
#   ) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(size = 14),
#     title = element_text(size = 14),
#     axis.title.y = element_text(size = 14),
#     strip.text.x = element_text(size = 10),
#     strip.text.y = element_text(size = 14),
#     axis.text.x = element_text(size = 12, color = "black"),
#     axis.text.y = element_text(size = 12, color = "black"),
#     legend.text = element_text(size = 12),
#     legend.title = element_text(size = 14),
#     plot.margin = unit(c(.5, 1, .5, .5), "cm")
#   )

# # Generate summary table of core ASVs
# core_table <- obs1 %>%
#   filter(fill == "core") %>%
#   select(OTU_ID, family, genus, fit_class)

# core_table

#--------------------------------------------------------
# THRESHOLD-BASED CORE SELECTION
#--------------------------------------------------------
# Extract core ASVs based on presence threshold (Jae's method)
prevalence_core <- filter_core(
  filtered_phyloseq,
  threshold = 0.6,
  as = "rows"
)

# Save threshold-based core ASVs
save(
  prevalence_core,
  file = "data/output/phyloseq_objects/prevalence_core.rda"
)

#--------------------------------------------------------
# JBEI-SPECIFIC CORE ANALYSIS
#--------------------------------------------------------
# Clean up JBEI metadata
new_metadata <- drought_jbei %>%
  sample_data() %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(
    across(brc, ~ str_to_lower(.)),
    across(everything(.), ~ as.character(.)),
    new_row = x_sample_id
  ) %>%
  column_to_rownames(., var = "new_row") %>%
  sample_data()

# Update JBEI phyloseq object with cleaned metadata
sample_data(drought_jbei) <- new_metadata

# Filter JBEI samples for quality
drought_jbei <- prune_samples(
  sample_sums(drought_jbei) >= 100,
  drought_jbei
)

drought_jbei <- filter_taxa(
  drought_jbei,
  function(x) {
    sum(x > 100) > (0.00 * length(x)) # Results depend on this cut-off.
  },
  TRUE
)

# Extract JBEI-specific core
jbei_braycore_summary <- extract_core(
  drought_jbei,
  Var = "treatment",
  method = "increase",
  increase_value = 2
)

# Generate JBEI-specific visualizations
bray_curtis_curve <- brc_bc_curve(
  core_summary_list = jbei_braycore_summary,
  max_otus = 100,
  threshold = 1.02
)
bray_curtis_curve

occ_abun_plot <- brc_bc_occ_curve(core_summary_list = jbei_braycore_summary)
occ_abun_plot
