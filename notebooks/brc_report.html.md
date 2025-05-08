---
title: "Inter-BRC Core Microbiome"
subtitle: "Determining a microbial 'core' in jbei samples"
author: "Bolívar Aponte Rolón"
date: "2025-04-21"
date-modified: "2025-04-23"
date-format: "long"
format:
  html:
    toc: true
    toc-location: left
    toc-depth: 2
    number-sections: true
    number-depth: 1
    theme: lumen
    highlight-style: github
    code-overflow: wrap
    code-fold: false
    code-copy: true
    code-link: false
    code-tools: false
    code-block-border-left: "#0C3823"
    code-block-bg: "#eeeeee"
    fig-cap-location: margin
    linestretch: 1.25
    fontsize: 14pt
    embed-resources: true
    css: styles.css
execute:
  echo: true
  warning: false
  message: false
  keep-md: true
editor:
  markdown:
    wrap: 72
    canonical: true
params:
  BRC: "BRC"
---




# Objective
Analyses for determining a microbial "core" using ASV contribution to Bray-Curtis dissimilarity.

# Setup


::: {.cell}

```{.r .cell-code}
# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

source(here::here("R/utils/000_setup.R"))
if (exists("phyloseq")) remove(phyloseq)

# # Check if filtered_phyloseq exists, if not, try to load it directly
# # Needed when rendering directly from terminal
# if (!exists("filtered_phyloseq")) {
#   message(
#     "filtered_phyloseq not loaded by setup script, trying to load it directly..."
#   )

#   # Try to load the filtered_phyloseq object directly
#   filtered_phyloseq_path <- here::here(
#     "data/output/phyloseq_objects/filtered_phyloseq.rda"
#   )
#   if (file.exists(filtered_phyloseq_path)) {
#     load(filtered_phyloseq_path)
#     message("Successfully loaded filtered_phyloseq directly from file")
#   } else {
#     stop("Could not find filtered_phyloseq at path: ", filtered_phyloseq_path)
#   }
# }

## Column name clean-up
colnames(sample_data(filtered_phyloseq)) <- filtered_phyloseq@sam_data |>
  janitor::clean_names() |>
  colnames()
```
:::



# Microbiome Core Selection via `extract_core()`


::: {.cell}

```{.r .cell-code}
# Check OTU table
# filtered_phyloseq <- drought_jbei
filtered_phyloseq <- phyloseq::prune_samples(
  phyloseq::sample_sums(filtered_phyloseq) >= 100,
  filtered_phyloseq
)

filtered_phyloseq <- phyloseq::subset_samples(
  filtered_phyloseq,
  brc == "jbei"
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

ggsave(
  filename = str_glue("{params$BRC}_bray_curtis_curve.png"),
  bray_curtis_curve,
  path = here::here(str_glue("data/output/plots/{params$BRC}/")),
  dpi = 300,
  width = 6,
  height = 4
)

ggsave(
  filename = str_glue("{params$BRC}_abundance_occupancy.png"),
  occ_abun_plot,
  path = here::here(str_glue("data/output/plots/{params$BRC}/")),
  dpi = 300,
  width = 6,
  height = 4
)

# Print
bray_curtis_curve +
  ggtitle(str_glue("{params$BRC} Core Microbial Community ranked"))
```

::: {.cell-output-display}
![](brc_report_files/figure-html/process-data-1.png){width=672}
:::

```{.r .cell-code}
occ_abun_plot + ggtitle(str_glue("{params$BRC} Core Microbial Occurance Curve"))
```

::: {.cell-output-display}
![](brc_report_files/figure-html/process-data-2.png){width=672}
:::
:::




# Microbiome Core Selection by Threshold

Analysis based on Jae's code.
Depending on how your phyloseq object's otu table is structured (e.g., if taxa_are_rows = FALSE ) you might have to play with nrow()/ncol() and rowSums/colSums()



::: {.cell}

```{.r .cell-code}
# core_asvs_threshold <- filter_core(
#   filtered_phyloseq,
#   threshold = 0.6,
#   as = "rows"
# )

# save(
#   core_asvs_threshold,
#   file = here::here(
#     "data/output/phyloseq_objects/{params$BRC}/{params$BRC}_high_threshold.rda"
#   )
# )
```
:::