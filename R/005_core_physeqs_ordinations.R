## Inter-BRC Core Analysis
## Core & non_core Phyloseqs and Ordination Plots
## By Bolívar Aponte Rolón

# Setup
source("R/000_setup.R")
load("data/output/core_summary_lists.rda")
load("data/output/phyloseq_objects/filtered_phyloseq.rda")

#################################################
### Core and Non_Core objects and fasta files ###
#################################################

# Subset by "core" and "non-core" names using "core_summary_lists" object

## Name strings to subset
core_asv_strings <- core_summary_lists[[4]] %>%
  dplyr::filter(., .$fill == "core") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

non_core_asv_strings <- core_summary_lists[[4]] %>%
  dplyr::filter(., .$fill == "no") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

# "core and "non_core" communities as matrices and phyloseq objects
## ASV (OTU) Matrices
## Remove samples where row sum 0 through ExtractMatrix()

core_asv_matrix <- ExtractMatrix(filtered_phyloseq, .vec = core_asv_strings) # Useful for manual ordinations

non_core_asv_matrix <- ExtractMatrix(filtered_phyloseq, .vec = non_core_asv_strings)

## Sample names and data
core_sample_strings <- rownames(core_asv_matrix)
non_core_sample_strings <- rownames(non_core_asv_matrix)

#################################################
## Core dimensions is 50 ASVs out of a total of 23473 non-core ASVs in 1813 samples
## When non-zero sums samples are gathered we en up with 50 ASVs present in 1733 samples
## Core ASVs are present in 1733 / 1813 samples (96%)
#################################################

## Phyloseqs
core_brc_phyloseq <- prune_samples(
  sort(sample_names(filtered_phyloseq)) %in% sort(core_sample_strings),
  filtered_phyloseq
) %>%
  prune_taxa(rownames(.@otu_table) %in% core_asv_strings, .)

non_core_brc_phyloseq <- prune_samples(
  sort(sample_names(filtered_phyloseq)) %in% sort(non_core_sample_strings),
  filtered_phyloseq
) %>%
  prune_taxa(rownames(.@otu_table) %in% non_core_asv_strings, .)

## Checking that core is filtered out
rownames(core_brc_phyloseq@otu_table) %in% non_core_asv_strings # Should return FALSE for all 50 ASVs

# Save phyloseqs
save(core_brc_phyloseq, file = "data/output/phyloseq_objects/core_brc_phyloseq.rda")
save(non_core_brc_phyloseq, file = "data/output/phyloseq_objects/non_core_brc_phyloseq.rda")

# Creating Core and non_core Fasta files
## Subsetting fasta based on `subset.fasta` function from
## https://github.com/GuillemSalazar/FastaUtils/blob/master/R/FastaUtils_functions.R

## Core
subset.fasta(
  file = "data/output/rep_asv_seqs.fasta",
  subset = core_asv_strings,
  out = "data/output/core_asv_seqs.fasta"
)

## Non-Core
subset.fasta(
  file = "data/output/rep_asv_seqs.fasta",
  subset = non_core_asv_strings,
  out = "data/output/non_core_asv_seqs.fasta"
)


#################################################
###### Ordinations and Community Analysis #######
#################################################

# Hellinger transformation of matrices
# To be used for NMDS, dbRDA and adonis2
core_hell_matrix <- decostand(t(core_asv_matrix),
  MARGIN = 1,
  method = "hellinger"
) # Now we need samples as columns and ASV are rows

###########

# Ordinations
core_asv_dist <- vegdist(t(core_hell_matrix), method = "bray", upper = FALSE, binary = FALSE, na.rm = TRUE)

## Choosing the number of dimensions
NMDS.scree <- function(x) { # where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress),
    xlim = c(1, 10),
    ylim = c(0, 0.30),
    xlab = "# of Dimensions",
    ylab = "Stress",
    main = "NMDS stress plot"
  )

  for (i in 1:10) {
    points(rep(i + 1, 10), replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

NMDS.scree(core_asv_dist)

NMDS <- metaMDS(as.matrix(core_asv_dist),
  distance = "bray",
  display = c("sites", "species"),
  noshare = TRUE,
  autotransform = FALSE,
  wascores = TRUE,
  tidy = TRUE,
  k = 4, trymax = 500
)
stressplot(NMDS)

# Adding site scores to `NMDS`
sppscores(NMDS) <- core_hell_matrix

## Scores and sample data for NMDS
nmds_scores <- as_tibble(vegan::scores(NMDS)$sites)
core_brc_sample_df <- core_brc_phyloseq@sam_data |>
  data.frame()


# NMDS Aesthetics ####

ggplot(nmds_scores, x = NMDS1, y = NMDS2) +
  geom_point(data = nmds_scores, aes(size = 3, , color = core_brc_sample_df$brc, alpha = 0.5, stroke = 1))



scale_shape_manual(values = c(15:22)) +
  # scale_color_manual(values = met.brewer(name = "Nizami", n = 8, type="discrete")) +
  geom_hline(yintercept = 0, colour = "grey50", size = 0.65) +
  geom_vline(xintercept = 0, colour = "grey50", size = 0.65) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right", legend.title = element_text(),
    legend.text = element_text(face = "italic")
  ) +
  stat_ellipse(
    mapping = NULL, data = nmds_scores, geom = "path", size = 1.3,
    position = "identity", type = "t", linetype = 1, level = 0.95, segments = 51,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE
  )

NM.DS



# Print the plot
print(nmds_plot)
