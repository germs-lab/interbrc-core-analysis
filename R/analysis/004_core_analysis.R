## Inter-BRC Core Analysis
## This was part of the `core_microbiome_functions.R` by Brandon Kristy.
## Split by Bolívar Aponte Rolón 2025-02-20


# Setup
source("R/utils/000_setup.R")

###################################################
### Microbiome Core Selection via ExtractCore() ###
###################################################

# Check OTU table
filtered_phyloseq <- prune_samples(sample_sums(filtered_phyloseq) >= 100, filtered_phyloseq)

# Remove major outliers determined by Bray-Curtis beta diversity analysis
# filtered_phyloseq <- subset_samples(
#     filtered_phyloseq,
#     sample_names(filtered_phyloseq) != "X2016.MMPRNT.RHN.G5.r2.Unfert.pr3"
# )

# Extract the 'spatial' core microbiome across all sites. The 'Var' in the ExtractCore is 'site'.
core_summary_lists <- ExtractCore(filtered_phyloseq,
  Var = "site",
  method = "increase",
  increase_value = 2
) # Minimum seq depth was ~10,000 reads.

# Save, since it takes a long time.
# save(core_summary_lists, file = "data/output/core_summary_lists.rda")


# Plot Bray-Curtis Dissimilarity Curve:
max <- 100 # Number of ranked-OTUs to plot

BC_ranked <- core_summary_lists[[2]] %>%
  dplyr::mutate(., rank = factor(.$rank, levels = .$rank)) %>%
  drop_na()

BC_ranked_max <- BC_ranked[1:max, ]


BC_ranked_max %>%
  ggplot(aes(x = rank[1:max], y = proportionBC)) +
  geom_point(
    pch = 21,
    col = "black",
    alpha = 0.85,
    size = 3.5
  ) +
  geom_vline(
    xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),
    lty = 4,
    col = "darkred",
    cex = 0.5
  ) +
  scale_x_discrete(
    limits = BC_ranked$rank[1:max], # making the x axis more readable
    breaks = BC_ranked$rank[1:max][seq(1, length(BC_ranked$rank[1:max]),
      by =
        10
    )]
  ) +
  xlab("Ranked OTUs") +
  ylab("Contribution to Bray-Curtis Similarity") +
  annotate(
    geom = "text",
    x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),
    y = 0.18,
    label = paste("Last 2% increase\n(", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >=
      1.02)])), " OTUs)", sep = ""),
    color = "darkred",
    size = 4,
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    title = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.margin = unit(c(.5, 1, .5, .5), "cm")
  )

# Plot Abundance Occupancy curve
occ_abun_plot <- core_summary_lists[[4]] %>%
  ggplot(aes(
    x = log10(otu_rel),
    y = otu_occ,
    fill = fill
  )) +
  # scale_fill_npg(
  #   name = "Core Membership",
  #   labels = c("Core Taxa", "Non Core Taxa")
  # ) +
  geom_point(pch = 21, alpha = 1, size = 2.5) +
  labs(x = "log10(Mean Relative Abundance)", y = "Occupancy") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12),
    title = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.margin = unit(c(.5, 1, .5, .5), "cm")
  )
occ_abun_plot

ggsave(filename = "abundance_occupancy.png", occ_abun_plot, path = "data/output/plots/", dpi = 300, width = 6, height = 4)


# Fit Abundance-Occupancy Distribution to a Neutral Model
# Fit neutral model
taxon <- core_summary_lists[[7]]
spp <- t(core_summary_lists[[5]])
occ_abun <- core_summary_lists[[4]]
names(occ_abun)[names(occ_abun) == "otu"] <- "OTU_ID"

# source community pool
# meta <- core_summary_lists[[6]]

# fitting model
debugonce(sncm.fit)
obs.np <- sncm.fit(spp, taxon, stats = FALSE, pool = NULL)
sta.np <- sncm.fit(spp, taxon, stats = TRUE, pool = NULL)

obs.np$fit_class <- "As predicted"
obs.np[which(obs.np$freq < obs.np$pred.lwr), "fit_class"] <-
  "Below prediction"
obs.np[which(obs.np$freq > obs.np$pred.upr), "fit_class"] <-
  "Above prediction"
obs.np[which(is.na(obs.np$freq)), "fit_class"] <- "NA"

obs.np <- tibble::rownames_to_column(obs.np, "OTU_ID")
as.data.frame(left_join(occ_abun, obs.np, by = "OTU_ID")) -> fit_table
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
  mutate(fill_fit_class = paste0(fill, ":", fit_class))

obs1 <- obs1 %>%
  mutate(across(
    "fill_fit_class",
    str_replace,
    "no:Below prediction",
    "Non Core Taxa"
  )) %>%
  mutate(across(
    "fill_fit_class",
    str_replace,
    "no:Above prediction",
    "Non Core Taxa"
  )) %>%
  mutate(across(
    "fill_fit_class",
    str_replace,
    "no:As predicted",
    "Non Core Taxa"
  ))


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
core_asvs_threshold <- filter_core(filtered_phyloseq, threshold = 0.6, as = "rows")

physeq_high_occ <- core_asvs_threshold$physeq_high_occ
physeq_low_occ <- core_asvs_threshold$physeq_low_occ

# save(core_asvs_threshold, file = "data/output/phyloseq_objects/core_asvs_threshold.rda")

physeq_high_occ_matrix <- physeq_high_occ@otu_table %>%
    t() %>% # Samples as rows
    as.data.frame() %>%
    .[rowSums(.) > 0, ] %>% # Keep only samples with a non-zero sum. Not all samples have the "core".
    as.matrix()


physeq_low_occ_matrix <- physeq_low_occ@otu_table %>%
    t() %>% # Samples as rows
    as.data.frame() %>%
    .[rowSums(.) > 0, ] %>% # Keep only samples with a non-zero sum. Not all samples have the "core".
    as.matrix()


save(physeq_high_occ_matrix, file = "data/output/high_occ_matrix.rda")
save(physeq_low_occ_matrix, file = "data/output/low_occ_matrix.rda")
