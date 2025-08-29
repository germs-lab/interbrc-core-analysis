#########################################################
# Clustering Threshold Analyiss for SABR data
#
# Project:  Inter-BRC-Core-Microbiome
# Authors: Jaejin Lee & Bolívar Aponte Rolón
# Last modified: 2025-08-29
#########################################################

#--------------------------------------------------------
# Setup
#--------------------------------------------------------
library(phyloseq)
library(dplyr)
library(tidyr)
library(readr)
library(here)

# basic settings
similarities <- c(85, 87, 90, 93, 95, 97)
occ_cutoffs <- c(50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)

rds_in <- here::here("data/input/sabr/all_ps_SABR.rds")
out_dir <- here::here("data/output/sabr")

ps_list <- readRDS(rds_in)
stopifnot(all(paste0("ps_", similarities) %in% names(ps_list)))

for (nm in names(ps_list)) {
  assign(nm, ps_list[[nm]], envir = .GlobalEnv)
}

# calculating: (1) Core OTU#, (2) Core OTU seq percentage(%)
results <- list()

for (sim in similarities) {
  ps <- get(paste0("ps_", sim))
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(otu_table(ps))) {
    otu_mat <- t(otu_mat)
  } # 행=OTU

  total_counts_all <- sum(otu_mat)
  n_all <- nsamples(ps)
  if (n_all == 0 || total_counts_all == 0) {
    next
  }

  occ_counts <- rowSums(otu_mat > 0)

  for (occ in occ_cutoffs) {
    min_occ <- ceiling((occ / 100) * n_all)
    core_otus <- names(occ_counts[occ_counts >= min_occ])

    core_counts <- if (length(core_otus) > 0) {
      sum(otu_mat[core_otus, , drop = FALSE])
    } else {
      0
    }
    percent_core <- round(100 * core_counts / total_counts_all, 2)

    results[[length(results) + 1]] <- data.frame(
      Similarity = sim,
      Occurrence = occ,
      Core_OTUs = length(core_otus),
      Core_Counts = core_counts,
      Total_Counts = total_counts_all,
      Percent_Core = percent_core,
      stringsAsFactors = FALSE
    )
  }
}

summary_long <- bind_rows(results) %>%
  arrange(Similarity, desc(Occurrence))

# Wide tables
otu_counts_wide <- summary_long %>%
  select(Similarity, Occurrence, Core_OTUs) %>%
  pivot_wider(names_from = Occurrence, values_from = Core_OTUs) %>%
  arrange(Similarity)

seq_pct_wide <- summary_long %>%
  select(Similarity, Occurrence, Percent_Core) %>%
  pivot_wider(names_from = Occurrence, values_from = Percent_Core) %>%
  arrange(Similarity)

# Save
write_csv(summary_long, file.path(out_dir, "core_summary_long.csv"))
write_csv(otu_counts_wide, file.path(out_dir, "core_otu_counts.csv"))
write_csv(seq_pct_wide, file.path(out_dir, "core_seq_percent.csv"))
