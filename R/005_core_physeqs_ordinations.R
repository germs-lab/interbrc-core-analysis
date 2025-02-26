## Inter-BRC Core Analysis
## Core & non_core Phyloseqs and Ordination Plots
## By Bolívar Aponte Rolón

# Setup
source("R/000_setup.R")
load("data/output/spatial_core.rda")
load("data/output/phyloseq_objects/filtered_phyloseq.rda")

# Subset by "core" and "non-core" names using "spatiol_core" object

## Name strings to subset
core_names <- spatial_core[[4]] %>%
  dplyr::filter(., .$fill == "core") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

non_core_names <- spatial_core[[4]] %>%
  dplyr::filter(., .$fill == "no") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

# "core and "non_core" communities as matrices and phyloseq objects
## Making new elements for phyloseqs

## ASV (OTU) Matrices
core_asv_matrix <- filtered_phyloseq@otu_table |>
  t() |>
  as.data.frame() |>
  select(contains(core_names)) |>
  as.matrix()

test <- subset(otu_table(filtered_phyloseq), rownames(otu_table(filtered_phyloseq)) %in% core_names)

non_core_asv_matrix <- filtered_phyloseq@otu_table |>
  t() |>
  as.data.frame() |>
  select(contains(non_core_names)) |>
  as.matrix()

## sample data
all_samples <- sample_data(filtered_phyloseq@sam_data) # used for core and non_core

## taxonomy
core_taxa <- filtered_phyloseq@tax_table %>%
  as.data.frame() %>%
  dplyr::filter(rownames(.) %in% core_names) %>%
  as.matrix() %>%
  phyloseq::tax_table()

non_core_taxa <- filtered_phyloseq@tax_table %>%
  as.data.frame() %>%
  dplyr::filter(rownames(.) %in% non_core_names) %>%
  as.matrix() %>%
  phyloseq::tax_table()

## Phyloseqs
core_brc_phyloseq <- phyloseq(
  otu_table(core_asv_matrix, taxa_are_rows = FALSE),
  core_taxa,
  all_samples
)

non_core_brc_phyloseq <- phyloseq(
  otu_table(non_core_asv_matrix, taxa_are_rows = FALSE),
  non_core_taxa,
  all_samples
)

# Save
save(core_brc_phyloseq, file = "data/output/phyloseq_objects/core_brc_phyloseq.rda")
save(non_core_brc_phyloseq, file = "data/output/phyloseq_objects/non_core_brc_phyloseq.rda")


# Hellinger transformation of matrices rarefied to 750 reads
# To be used for dbRDA and adonis2
core_hell_matrix <- decostand(core_asv_matrix, MARGIN = 1, method = "hellinger")
###########

# Ordinations
core_asv_dist <- vegdist(t(core_hell_matrix), method = "bray", binary = FALSE)


# Choosing the number of dimensions
NMDS.scree <- function(x) { # where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10), ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1, 10), replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

# NMDS.scree(x) #4 dimension seem to be appropriate to keep the stress around 0.15. More dimension will complicate the interpretation of results.

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
sppscores(NMDS) <- core_asv_matrix

data.scores <- as_tibble(vegan::scores(NMDS)$sites)
# NMDS Aesthetics ####

ggplot(data.scores, aes(NMDS1, NMDS2)) +
  geom_point(data = data.scores, aes(size = 3, alpha = 0.5, stroke = 1))



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
    mapping = NULL, data = data.scores, geom = "path", size = 1.3,
    position = "identity", type = "t", linetype = 1, level = 0.95, segments = 51,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE
  )

NM.DS



# Print the plot
print(nmds_plot)
