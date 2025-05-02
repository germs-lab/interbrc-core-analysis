## Inter-BRC Core Analysis
## This was part of the `core_microbiome_functions.R` by Brandon Kristy.
## Split by Bolívar Aponte Rolón 2025-02-20

# Setup
source("R/utils/000_setup.R")
if (exists("phyloseq")) remove(phyloseq)


###################################################
### Microbiome Core Selection via extract_core() ###
###################################################

# Check OTU table
filtered_phyloseq <- prune_samples(
  sample_sums(filtered_phyloseq) >= 100,
  filtered_phyloseq
)

# Extract the 'spatial' core microbiome across all sites. The 'Var' in the ExtractCore is 'site'.
core_summary_lists <- extract_core(
  filtered_phyloseq,
  Var = "site",
  method = "increase",
  increase_value = 2
) # Minimum seq depth was ~10,000 reads.

# Save, since it takes a long time.
# save(core_summary_lists, file = "data/output/core_summary_lists.rda")

# Plot Bray-Curtis Dissimilarity Curve:
bray_curtis_curve <- brc_bc_curve(core_summary_list = core_summary_lists)
occ_abun_plot <- brc_bc_occ_curve(core_summary_list = core_summary_lists)

occ_abun_plot

ggsave(
  filename = "abundance_occupancy.png",
  occ_abun_plot,
  path = "data/output/plots/",
  dpi = 300,
  width = 6,
  height = 4
)


# Fit Abundance-Occupancy Distribution to a Neutral Model
# Fit neutral model
taxon <- core_summary_lists[[7]]
spp <- t(core_summary_lists[[5]])
occ_abun <- core_summary_lists[[4]]
names(occ_abun)[names(occ_abun) == "otu"] <- "OTU_ID"

# source community pool
# meta <- core_summary_lists[[6]]

# Fitting model

obs.np <- sncm.fit2(spp, taxon, stats = FALSE, pool = NULL)
sta.np <- sncm.fit2(spp, taxon, stats = TRUE, pool = NULL)

obs.np$fit_class <- "As predicted"
obs.np[which(obs.np$freq < obs.np$pred.lwr), "fit_class"] <-
  "Below prediction"
obs.np[which(obs.np$freq > obs.np$pred.upr), "fit_class"] <-
  "Above prediction"
obs.np[which(is.na(obs.np$freq)), "fit_class"] <- "NA"

obs.np <- tibble::rownames_to_column(obs.np, "OTU_ID")
as.data.frame(dplyr::left_join(occ_abun, obs.np, by = "OTU_ID")) -> fit_table
#
sta.np$above.pred <-
  sum(obs.np$freq > (obs.np$pred.upr), na.rm = TRUE) / sta.np$Richness
sta.np$below.pred <-
  sum(obs.np$freq < (obs.np$pred.lwr), na.rm = TRUE) / sta.np$Richness
fit_res <- as.data.frame(sta.np)
rownames(fit_res) <- "Value"

fit_res
list_tab <- list(fit_res, fit_table)

# Plot Neutral Function
obs1 <- as.data.frame(list_tab[[2]])
obs1 <- obs1[!is.na(obs1$p), ]
obs2 <- as.data.frame(list_tab[[1]])

obs1 <- obs1 %>%
  mutate(fill_fit_class = paste0(fill, ":", fit_class)) %>%
  mutate(
    fill_fit_class = case_when(
      str_detect(fill, "no") ~ "Non Core Taxa",
      TRUE ~ fill_fit_class
    )
  )


obs1 %>%
  ggplot(aes(x = log10(otu_rel), y = otu_occ)) +
  scale_fill_npg(
    name = "Core Membership: Model Predictions",
    labels = c(
      "Core: Above Prediction",
      "Core: As Predicted",
      "Core: Below Prediction",
      "Non-Core Taxa"
    )
  ) +
  geom_point(
    aes(fill = fill_fit_class),
    pch = 21,
    alpha = 0.75,
    size = 2.2
  ) +
  geom_line(
    color = "red",
    data = obs1,
    size = 0.8,
    aes(y = obs1$freq.pred, x = log10(obs1$p)),
    alpha = 0.55
  ) +
  geom_line(
    color = "black",
    lty = "twodash",
    size = 0.9,
    data = obs1,
    aes(y = obs1$pred.upr, x = log10(obs1$p)),
    alpha = 0.55
  ) +
  geom_line(
    color = "black",
    lty = "twodash",
    size = 0.9,
    data = obs1,
    aes(y = obs1$pred.lwr, x = log10(obs1$p)),
    alpha = 0.55
  ) +
  labs(x = "Log10(mean abundance)", y = "Occupancy") +
  annotate(
    "text",
    -Inf,
    Inf,
    label = paste("italic(R)^2 ==", round(obs2$Rsqr, 3)),
    parse = TRUE,
    size = 4.8,
    hjust = -0.2,
    vjust = 1.2
  ) +
  annotate(
    "text",
    -Inf,
    Inf,
    label = paste("italic(m) ==", round(obs2$m, 3)),
    parse = TRUE,
    size = 4.8,
    hjust = -0.2,
    vjust = 3.2
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    title = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.margin = unit(c(.5, 1, .5, .5), "cm")
  )

## Create Table
# Generate a table of core OTUs combined with the fitted model predictions
core_table <- obs1 %>%
  filter(fill == "core") %>%
  select(OTU_ID, family, genus, fit_class)

core_table


##############################################
### Microbiome Core Selection by Threshold ###
##############################################
# Analysis based on Jae's code
# Depending on how your phyloseq object's otu table is structured (e.g., if taxa_are_rows = FALSE ),
# you might have to play with nrow()/ncol() and rowSums/colSums()

# Load phyloseq object
core_asvs_threshold <- filter_core(
  filtered_phyloseq,
  threshold = 0.6,
  as = "rows"
)

save(
  core_asvs_threshold,
  file = "data/output/phyloseq_objects/core_asvs_threshold.rda"
)


########################################################
### Core Selection of JBEI dataset via extract_core() ###
########################################################

# Data set clean up
new_metadata <- drought_jbei %>%
  sample_data() %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(
    across(brc, ~ str_to_lower(.)),
    across(everything(.), ~ as.character(.)),
    new_row = x_sample_id
  ) %>%
  column_to_rownames(., var = "new_row") %>% # Workaround to inserting "sa1" type rownames
  sample_data()

# Update phyloseq object
sample_data(drought_jbei) <- new_metadata

# save(drought_jbei, file = "data/output/phyloseq_objects/jbei/drought_jbei.rda")

# Check OTU table
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

# Extract core
jbei_core_summary_lists <- extract_core(
  drought_jbei,
  Var = "treatment",
  method = "increase",
  increase_value = 2
)

# Plot Bray-Curtis Dissimilarity Curve:
bray_curtis_curve <- brc_bc_curve(
  core_summary_list = jbei_core_summary_lists,
  max_otus = 100,
  threshold = 1.02
)
bray_curtis_curve

occ_abun_plot <- brc_bc_occ_curve(core_summary_list = jbei_core_summary_lists)
occ_abun_plot
