# Inter-BRC-Core-Analysis
# Summary Stats and Column rename
# Bolívar Aponte Rolón
# 2025-02-18


# Setup
source("R/000_setup.R")


## Source files
load(
    "data/output/phyloseq_objects/filtered_phyloseq.rda"
)
load(
    "data/output/phyloseq_objects/unfiltered_phyloseq.rda"
)


## Summary stats

metagMisc::phyloseq_summary(filtered_phyloseq, more_stats = F, long = F)

percent_phyla_samples <- metagMisc::phyloseq_ntaxa_by_tax(
    filtered_phyloseq,
    TaxRank = "phylum",
    relative = F,
    add_meta_data = F
) |>
    as.data.frame() |>
    mutate(sum = sum(N.OTU)) |>
    group_by(phylum) |>
    summarise(occurance_in_samples = n())


# Saving phyloseq as DF

metagMisc::phyloseq_to_df(
    filtered_phyloseq,
    addtax = T,
    addtot = F,
    addmaxrank = F,
    sorting = "abundance"
) %>%
    # rename(ASV = OTU) %>%
    write.csv(.,
              file.path("data/output/phyloseq_objects/filtered_phyloseq_df.csv"))

