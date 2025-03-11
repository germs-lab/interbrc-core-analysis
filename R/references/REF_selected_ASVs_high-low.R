# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(magrittr)

# Load phyloseq object
ps <- readRDS(file = "/Users/jaejinlee/Files/Data/2023SABR_amplicon/analysis/ps_object.rds")

# Convert to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Convert OTU table to data frame
asv_table_data <- as.data.frame(otu_table(ps_rel))

# Calculate occurrence of each ASV across samples
asv_sample_counts <- colSums(asv_table_data > 0)

# Get total number of samples
total_samples <- nrow(asv_table_data)

# Define occurrence threshold (60% of samples)
threshold_60 <- total_samples * 0.6

# Select ASVs based on occurrence criteria
high_occurrence_asvs <- names(asv_sample_counts[asv_sample_counts >= threshold_60])
low_occurrence_asvs <- names(asv_sample_counts[asv_sample_counts < threshold_60])

# Print ASV counts
cat("ASVs found in ≥60% samples:", length(high_occurrence_asvs), "\n")
cat("ASVs found in <60% samples:", length(low_occurrence_asvs), "\n")

# Filter phyloseq objects for each category
ps_high_occ <- prune_taxa(high_occurrence_asvs, ps_rel)
ps_low_occ <- prune_taxa(low_occurrence_asvs, ps_rel)
ps_high_occ
ps_low_occ

set.seed(1234)

# User-defined number of ASVs to sample
num_asvs <- 55  # Change this value if needed (e.g., 10, 20, 30, 56, etc.)

# Ensure num_asvs is not greater than available ASVs
num_asvs <- min(num_asvs, length(high_occurrence_asvs), length(low_occurrence_asvs))

# Function to select ASVs while avoiding empty samples
select_valid_asvs_fast <- function(asv_list, ps_obj, num_asvs) {
  asv_table <- as.data.frame(otu_table(ps_obj))
  
  # Prioritize ASVs that are detected in the highest number of samples
  sorted_asvs <- asv_list[order(-colSums(asv_table[, asv_list] > 0))]  # Sort by occurrence frequency
  selected_asvs <- sorted_asvs[1:num_asvs]  # Select the top 'num_asvs' ASVs
  
  return(selected_asvs)
}

# Select ASVs ensuring no empty samples remain
selected_high_occ_asvs <- select_valid_asvs_fast(high_occurrence_asvs, ps_rel, num_asvs)
selected_low_occ_asvs <- select_valid_asvs_fast(low_occurrence_asvs, ps_rel, num_asvs)

# Prune the phyloseq object to include only selected ASVs
ps_high_occ_selected <- prune_taxa(selected_high_occ_asvs, ps_rel)
ps_low_occ_selected <- prune_taxa(selected_low_occ_asvs, ps_rel)
ps_high_occ_selected
ps_low_occ_selected

# Check if any samples are empty after pruning
empty_samples_high <- sum(rowSums(otu_table(ps_high_occ_selected)) == 0)
empty_samples_low <- sum(rowSums(otu_table(ps_low_occ_selected)) == 0)

cat("Empty samples in empty_samples_high:", empty_samples_high, "\n")
cat("Empty samples in ps_low_occ_selected:", empty_samples_low, "\n")


# Perform NMDS on selected ASVs
nmds_high_occ_selected <- ordinate(ps_high_occ_selected, method = "NMDS", distance = "bray", trymax = 100)
nmds_low_occ_selected <- ordinate(ps_low_occ_selected, method = "NMDS", distance = "bray", trymax = 100)

# Plot NMDS results
nmds_plot_high_occ_selected <- plot_ordination(ps_high_occ_selected, nmds_high_occ_selected, color = "Plant", shape = "Location") +
  geom_point(size = 3, alpha = 0.8) +
  #facet_wrap(~Date, ncol = 2) +
  labs(title = "NMDS of Selected Highly Prevalent ASVs", x = "NMDS1", y = "NMDS2") +
  theme_minimal()

nmds_plot_low_occ_selected <- plot_ordination(ps_low_occ_selected, nmds_low_occ_selected, color = "Plant", shape = "Location") +
  geom_point(size = 3, alpha = 0.8) +
  #facet_wrap(~Date, ncol = 2) +
  labs(title = "NMDS of Selected Low Prevalence ASVs", x = "NMDS1", y = "NMDS2") +
  theme_minimal()

# Print NMDS plots
print(nmds_plot_high_occ_selected)
print(nmds_plot_low_occ_selected)

### Procrustes Analysis
# Extract NMDS scores (coordinates)
nmds_high_coords_selected <- scores(nmds_high_occ_selected)
nmds_low_coords_selected <- scores(nmds_low_occ_selected)

# Perform Procrustes analysis
procrustes_result_selected <- procrustes(nmds_high_coords_selected, nmds_low_coords_selected, symmetric = TRUE)

# Statistical significance test (Procrustes-based protest)
protest_result_selected <- protest(nmds_high_coords_selected, nmds_low_coords_selected, permutations = 999)

# Print Procrustes results
print(procrustes_result_selected)
print(protest_result_selected)

# Visualization of Procrustes results
plot(procrustes_result_selected, main = "Procrustes Analysis (Selected ASVs: High vs Low Occurrence)")


### Beta Dispersion Analysis
dispersion_high_selected <- betadisper(bray_dist_high_selected, sample_data(ps_high_occ_selected)$Plant)
dispersion_low_selected <- betadisper(bray_dist_low_selected, sample_data(ps_low_occ_selected)$Plant)

# ANOVA test for dispersion differences
anova(dispersion_high_selected)
anova(dispersion_low_selected)

# Visualization of beta dispersion
boxplot(dispersion_high_selected, main="Beta Dispersion - Selected High Occurrence ASVs")
boxplot(dispersion_low_selected, main="Beta Dispersion - Selected Low Occurrence ASVs")





###### Additional: whole vs selected ASVs

# Perform NMDS on the entire data (without splitting into high/low occurrence)
nmds_all <- ordinate(ps_rel, method = "NMDS", distance = "bray", trymax = 100)

# Plot NMDS results for the entire dataset
nmds_plot_all <- plot_ordination(ps_rel, nmds_all, color = "Plant", shape = "Location") +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "NMDS of All ASVs", x = "NMDS1", y = "NMDS2") +
  theme_minimal()

# Print NMDS plot
print(nmds_plot_all)

### Procrustes Analysis (Whole Data vs. High Occurrence and Low Occurrence)
# Extract NMDS scores for all, high occurrence, and low occurrence
nmds_all_coords <- scores(nmds_all)
nmds_high_coords <- scores(nmds_high_occ_selected)
nmds_low_coords <- scores(nmds_low_occ_selected)

# Perform Procrustes analysis between the entire dataset and high/low occurrence
procrustes_all_vs_high <- procrustes(nmds_all_coords, nmds_high_coords, symmetric = TRUE)
procrustes_all_vs_low <- procrustes(nmds_all_coords, nmds_low_coords, symmetric = TRUE)

# Statistical significance test (Procrustes-based protest)
protest_all_vs_high <- protest(nmds_all_coords, nmds_high_coords, permutations = 999)
protest_all_vs_low <- protest(nmds_all_coords, nmds_low_coords, permutations = 999)

# Print Procrustes results for all vs high and all vs low
print(procrustes_all_vs_high)
print(procrustes_all_vs_low)
print(protest_all_vs_high)
print(protest_all_vs_low)

# Visualization of Procrustes results for all vs high and all vs low
plot(procrustes_all_vs_high, main = "Procrustes Analysis (All vs High Occurrence)")
plot(procrustes_all_vs_low, main = "Procrustes Analysis (All vs Low Occurrence)")


### Beta Dispersion Analysis (Whole Data vs. High Occurrence and Low Occurrence)
dispersion_all <- betadisper(bray_dist_all, sample_data(ps_rel)$Plant)
dispersion_high <- betadisper(bray_dist_high, sample_data(ps_high_occ_selected)$Plant)
dispersion_low <- betadisper(bray_dist_low, sample_data(ps_low_occ_selected)$Plant)

# ANOVA test for dispersion differences
anova(dispersion_all)
anova(dispersion_high)
anova(dispersion_low)

# Visualization of beta dispersion for all, high occurrence, and low occurrence
boxplot(dispersion_all, main="Beta Dispersion - All ASVs")
boxplot(dispersion_high, main="Beta Dispersion - High Occurrence ASVs")
boxplot(dispersion_low, main="Beta Dispersion - Low Occurrence ASVs")




###### Negative control test
nmds_high_occ_selected
nmds_low_occ_selected

# Extract the 'sites' matrix (coordinates)
nmds_low_sites <- nmds_low_coords_selected$sites

# Convert to data frame
nmds_low_sites_df <- as.data.frame(nmds_low_sites)

# Randomly shuffle the NMDS1 and NMDS2 coordinates
set.seed(123)  # Set seed for reproducibility
shuffled_nmds1 <- sample(nmds_low_sites_df$NMDS1)
shuffled_nmds2 <- sample(nmds_low_sites_df$NMDS2)

# Create a new object for the negative control with shuffled coordinates
nmds_negative_ctrl <- list(
  sites = data.frame(NMDS1 = shuffled_nmds1, NMDS2 = shuffled_nmds2),
  species = nmds_low_coords_selected$species  # Keep the species part the same
)

# Ensure sample names match between the negative control and ps_rel
rownames(nmds_negative_ctrl$sites) <- sample_names(ps_rel)

# Plot NMDS results for the negative control with Plant and Location
nmds_plot_negative_ctrl <- plot_ordination(ps_rel, nmds_negative_ctrl, color = "Plant", shape = "Location") +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "NMDS of Negative Control (Shuffled Coordinates)", x = "NMDS1", y = "NMDS2") +
  theme_minimal()

# Print the plot
print(nmds_plot_negative_ctrl)


# Perform Procrustes analysis between high occurrence and negative control
procrustes_result_negative_ctrl <- procrustes(scores(nmds_high_occ_selected), nmds_negative_ctrl$sites, symmetric = TRUE)

# Perform statistical significance test (Procrustes-based protest)
protest_result_negative_ctrl <- protest(scores(nmds_high_occ_selected), nmds_negative_ctrl$sites, permutations = 999)

# Print results
print(procrustes_result_negative_ctrl)
print(protest_result_negative_ctrl)

# visualize the Procrustes analysis
plot(procrustes_result_negative_ctrl, main = "Procrustes Analysis (High Occurrence vs Negative Control)")

