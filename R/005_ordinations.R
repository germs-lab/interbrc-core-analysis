## Inter-BRC Core Analysis
## Ordination Plots
## By Bolívar Aponte Rolón

# Setup
source("R/000_setup.R")
load("data/output/spatial_core.rda")
load("data/output/phyloseq_objects/filtered_phyloseq.rda")

# Subset by "core" and "non-core" names using "spatiol_core" object

## Names strings to subset
core_names <- spatial_core[[4]] %>%
  dplyr::filter(., .$fill == "core") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

non_core_names <- spatial_core[[4]] %>%
  dplyr::filter(., .$fill == "no") %>%
  column_to_rownames(., var = "otu") %>%
  rownames()

## ASV matrix for "core and "non_core" communities
core_asv_matrix <- filtered_phyloseq@otu_table |> 
  as.data.frame() |>
  select(contains(core_names)) |>
  as.matrix()

non_core_asv_matrix <- filtered_phyloseq@otu_table |> 
  t() |>
  as.data.frame() |>
  select(contains(non_core_names)) |> 
  as.matrix()


# # Hellinger transformation of matrices rarefied to 750 reads
# # To be used for dbRDA and adonis2
# rrfy_hell_matrix <- decostand(as.matrix(t(rep_otu_df(pseq_rrfb))), MARGIN = 1, method = "hellinger")
# ###########

# # Ordinations
# otu_e.dist <- vegdist(t(rabun), method = "bray", binary = FALSE)




# ###### Jae's methods
# ### NMDS
# # calculating Bray-Curtis distance
# bc_dist <- phyloseq::distance(ps, method = "bray")

# # NMDS analysis (k=2 reduces to 2 dimensions)
# nmds <- ordinate(ps, method = "NMDS", distance = "bray", trymax = 100)

# # Visualizing NMDS results
# # Use the "Date" column in the metadata to create four panels by sampling date
# # Convert "Date" to a factor for proper panel separation
# sample_data(ps)$Date <- as.factor(sample_data(ps)$Date)
# colnames(sample_data(ps))

# # Generate NMDS plot
# nmds_plot <- plot_ordination(ps, nmds, color = "Plant", shape = "Location") +
#   geom_point(size = 3, alpha = 0.8) +
#   facet_wrap(~Date, ncol = 2) + # Separate panels by Date (4 panels)
#   labs(
#     title = "NMDS of Microbial Communities by Date",
#     x = "NMDS1",
#     y = "NMDS2"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "right",
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )

# # Print the plot
# print(nmds_plot)
