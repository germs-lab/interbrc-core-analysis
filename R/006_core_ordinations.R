## Inter-BRC Core Analysis
## Core & non_core Ordinations
## By Bolívar Aponte Rolón

# Setup
source("R/000_setup.R")
load("data/output/core_summary_lists.rda")
load("data/output/phyloseq_objects/filtered_phyloseq.rda")
load("data/output/core_asv_matrix.rda")
load("data/output/non_core_asv_matrix.rda")


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
set.seed(484035)
nmds_screen_parallel(core_asv_dist, ncores = 32) #  Results: Two dimensions keeps stress below 0.20

NMDS <- metaMDS(as.matrix(core_asv_dist),
  distance = "bray",
  display = c("sites", "species"),
  noshare = TRUE,
  autotransform = FALSE,
  wascores = TRUE,
  tidy = TRUE,
  k = 2,
  trymax = 500,
  parallel = parallel::detectCores()
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
