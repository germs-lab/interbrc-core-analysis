setwd("~/Dropbox/Projects/GLBRC/Inter_BRC_workshop/qiime_output/")

library(phyloseq)
library(tidyverse)

taxonomy <- read.delim("~/Library/CloudStorage/Dropbox/Projects/GLBRC/Inter_BRC_workshop/qiime_output/taxonomy.tsv", comment.char="#")
metadata <- read.delim("~/Library/CloudStorage/Dropbox/Projects/GLBRC/Inter_BRC_workshop/qiime_output/metadata.tsv", row.names =1)
feature.table.modified <- read.delim("~/Library/CloudStorage/Dropbox/Projects/GLBRC/Inter_BRC_workshop/qiime_output/feature-table-modified.tsv")

feature.table.modified <- feature.table.modified %>% remove_rownames %>% column_to_rownames(var="OTU.ID")


ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

taxonomy_modified <- taxonomy %>%
  mutate_at('Taxon',str_replace_all, "[a-z]__","") %>%
  separate(Taxon, sep = ';', into=ranks,remove = TRUE) %>%
  column_to_rownames(var = "Feature.ID") %>%
  as.matrix()


# Import 
TAX <- tax_table(taxonomy_modified)
feature.table.modified <- as.matrix(feature.table.modified)
OTU <- otu_table(feature.table.modified, taxa_are_rows = T)
metadata <- sample_data(metadata)

# Merge
phyloseq <- phyloseq(OTU, TAX, metadata)

# Phyloseq Filtering 
# Remove singletons
phyloseq_removed_singletons <- prune_taxa(taxa_sums(phyloseq) > 1, phyloseq)
tax.remove <- ntaxa(phyloseq) - ntaxa(phyloseq_removed_singletons) #1,675 ASVs
# Remove plant contaminants
phyloseq_removed_singletons_plant <- subset_taxa(phyloseq_removed_singletons, kingdom != "Eukaryota" & kingdom != "Unassigned")



n.filtered <- ntaxa(phyloseq_removed_singletons)- ntaxa(phyloseq_removed_singletons_plant) # 1,989 ASVs removed

# remove low sequence coverage samples
## all samples 
sorted <- sort(sample_sums(phyloseq_removed_singletons_plant))  
min <- min(sample_sums(phyloseq_removed_singletons_plant)) #0
max <- max(sample_sums(phyloseq_removed_singletons_plant)) # 237,153
mean <- mean(sample_sums(phyloseq_removed_singletons_plant)) # 19,051.85
median <- median(sample_sums(phyloseq_removed_singletons_plant)) #16934

# Remove samples with low number of reads (<100) - there's a couple of empty samples 
phyloseq_removed_singletons_plant <- prune_samples(sample_sums(phyloseq_removed_singletons_plant)>=100, phyloseq_removed_singletons_plant)



## Filter by read abundance: each ASV must have a minimum of 20 reads:

# Read abundance 
# ASV msut have at least 20 reads 
filtered_phyloseq <- filter_taxa(phyloseq_removed_singletons_plant, function(x) sum(x > 20) > (0.00*length(x)), TRUE)


# Export FILTERED and UNFILTERED phyloseq objects for analytics team: 

# Filtered object 
save(filtered_phyloseq, file = "filtered_phyloseq.rda")

# Unfiltered object 
save(unfiltered_phyloseq, file = "unfiltered_phyloseq.rda")





