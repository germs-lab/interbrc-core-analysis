---
title: "Inter-BRC Core Microbiome"
subtitle: "Determining a microbial 'core' in JBEI samples"
author: "Bolívar Aponte Rolón"
date: "2025-04-21"
date-modified: "2025-04-28"
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
    code-fold: true
    code-summary: "See code"
    code-copy: true
    code-link: false
    code-tools: false
    code-block-border-left: "#0C3823"
    code-block-bg: "#eeeeee"
    fig-cap-location: margin
    fig-height: 7
    fig-width: 7
    linestretch: 1.25
    fontsize: 14pt
    embed-resources: true
    #css: styles.css
    mermaid:
      theme: neutral
execute:
  echo: false
  warning: false
  message: false
  keep-md: true
editor:
  markdown:
    wrap: 72
    canonical: true
---




# Objective
Analyses for determining a microbial "core" using ASV contribution to Bray-Curtis dissimilarity.

# Setup


::: {.cell}

:::



# Microbiome Core Selection via `extract_core()` in all BRCs

To obtain the microbial core determined by ASV's contribution to Bray-Curtis dissimilarity, we:

- Pruned samples
  - Used `prune_samples` and `sample_sums` to remove any samples that had fewer than 100 total reads

- Filtered ASVs based on read counts
  - Used `filter_taxa` function to keep only ASVs that had:
    - More than 20 reads in at least one sample

All these steps were performed in `002_phyloseq.R` not shown here.

After which we:
- Extracted the core microbiome across different sites
  - Parameters:
    - Variable of interest was "site"
    - Used an "increase" method
    - With an increase value of 2




::::::{.cell layout-align="default"}

:::::{.cell-output-display}

::::{}
`<figure class=''>`{=html}

:::{}

<pre class="mermaid mermaid-js">flowchart LR
    A[All BRCs Raw Data] --&gt; B[DADA2 Pipeline]
    B --&gt; C[filtered_phyloseq&lt;br/&gt;All BRCs Dataset]
    C --&gt;|extract_core| D[BC-core]
    subgraph &quot;Initial Processing&quot;
        A
        B
    end
    subgraph &quot;Filtering Steps&quot;
        C
        style C fill:#fff,stroke:#333
    end
    subgraph &quot;All BRC core&quot;
        D
        style D fill:#fff,stroke:#333
    end
</pre>
:::
`</figure>`{=html}
::::
:::::
::::::


::: {.cell}
::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/process-BRC-data-1.png){width=672}
:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/process-BRC-data-2.png){width=672}
:::
:::




# Microbiome Core Selection via `extract_core()` in JBEI dataset

Now let's looks at the JBEI dataset subsetted from the filtered_phyloseq which contains all BRCs and was submitted to DADA2 pipeline and filtered as one data set. *The result is a phyloseq object with 59,961 ASvs and 16 samples.*




::::::{.cell layout-align="default"}

:::::{.cell-output-display}

::::{}
`<figure class=''>`{=html}

:::{}

<pre class="mermaid mermaid-js">flowchart LR
    A[All BRCs Raw Data] --&gt; B[DADA2 Pipeline]
    B --&gt; C[filtered_phyloseq&lt;br/&gt;All BRCs Dataset]
    C --&gt;|subset_samples&lt;br/&gt;brc == &#39;jbei&#39;| D[JBEI Dataset]
    subgraph &quot;Initial Processing&quot;
        A
        B
    end
    subgraph &quot;Filtering Steps&quot;
        C
        style C fill:#fff,stroke:#333
    end
    subgraph &quot;Final Dataset&quot;
        D
        style D fill:#fff,stroke:#333
    end
</pre>
:::
`</figure>`{=html}
::::
:::::
::::::


::: {.cell}
::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/process-JBEI-data-1.png){width=672}
:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/process-JBEI-data-2.png){width=672}
:::
:::



Notice the odd pattern! Not a smotth curve. It should look like the previous plot. The data is structure the same way. 

Let's look at the JBEI dataset filtered and prepared by Yen.

# JBEI core microbiome (Yen's version)

Here we cleaned up the data set a little bit to make sure that the data type were the same as previous data set.Ruling out data type as an issue.

We then filtered for minimum sample reads and ASV read counts in the same manner. *The result is a phyloseq object with 2,213 ASvs and 16 samples.*
We have about ~10 core ASVs.



::: {.cell}
::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/process-JBE-YEN-data-1.png){width=672}
:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/process-JBE-YEN-data-2.png){width=672}
:::
:::



Looking at the data closely, Phillip and I noticed that there ar many repeated read counts at low levels <100 reads. For example, we saw many ASV that had only 21 reads. They also formed a diagonal pattern. Much like a comparison matrix.

Does increasing the filtering parameters help? Let's change the minimum ASV reads counts.



::: {.cell}
::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-7-1.png){width=672}
:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-7-2.png){width=672}
:::
:::



*The result is a phyloseq object with 473 ASVs and 16 samples.* We have about ~10 core ASVs.

It still doesn't look great.

# Analyzing read patterns for JBEI
## Analysis Across Filtering Levels



::: {.cell}

:::

::: {.cell}
::: {.cell-output .cell-output-stdout}

```

###  Yen's original data 
```


:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-9-1.png){width=672}
:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-9-2.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```

###  More than 20 reads per ASV 
```


:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-9-3.png){width=672}
:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-9-4.png){width=672}
:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-9-5.png){width=672}
:::

::: {.cell-output .cell-output-stdout}

```

###  More than 100 reads per ASV 
```


:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-9-6.png){width=672}
:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-9-7.png){width=672}
:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-9-8.png){width=672}
:::

::: {.cell-output-display}
![](jbei_brc_report_files/figure-html/unnamed-chunk-9-9.png){width=672}
:::
:::




Concerning Patterns Identified:

1. High frequency of exactly 21 reads in multiple ASVs
2. Diagonal pattern formation in read count distribution
3. "matrix-like" pattern at low read counts (<100)

These patterns strongly suggest:

* Potential technical artifacts in the sequencing or processing pipeline
* Possible cross-contamination between samples
* Risk of false positive ASVs
* ***Selection of "core" with this method is not adequate due to sample size***

An ASV present with 21 reads in one sample contribution is proportionally greater. Bray-Curstis dissimilarity is sensitive to this. 


## PERMANOVA (dbRDA)

Code from Yen. Use `drought_JBEI.rds` with transformed sample counts (relative abundance)



::: {.cell}

```{.r .cell-code}
#|
# import and compute bray-curtis distance
ps_jbei <- readRDS(here::here(
  "data/output/phyloseq_objects/jbei/drought_JBEI.rds"
))

ps_rel <- transform_sample_counts(ps_jbei, function(x) x / sum(x))
BC.dist <- vegan::vegdist(t(otu_table(ps_rel)), method = "bray")

## permanova looking at harvest and treatment effect
vegan::adonis2(
  BC.dist ~ Harvest * Treatment,
  data = data.frame(sample_data(ps_rel)),
  nperm = 1000,
  method = "bray",
  by = "terms"
)
```

::: {.cell-output .cell-output-stdout}

```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = BC.dist ~ Harvest * Treatment, data = data.frame(sample_data(ps_rel)), method = "bray", by = "terms", nperm = 1000)
                  Df SumOfSqs      R2      F Pr(>F)    
Harvest            1  0.20456 0.28493 5.7995  0.001 ***
Treatment          1  0.04793 0.06676 1.3587  0.153    
Harvest:Treatment  1  0.04219 0.05876 1.1960  0.206    
Residual          12  0.42327 0.58956                  
Total             15  0.71794 1.00000                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


:::

```{.r .cell-code}
## remove harvest effect
vegan::adonis2(
  BC.dist ~ Treatment,
  data = data.frame(sample_data(ps_rel)),
  nperm = 1000,
  method = "bray",
  by = "terms",
  strata = sample_data(ps_rel)$Harvest
)
```

::: {.cell-output .cell-output-stdout}

```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Blocks:  strata 
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = BC.dist ~ Treatment, data = data.frame(sample_data(ps_rel)), method = "bray", by = "terms", strata = sample_data(ps_rel)$Harvest, nperm = 1000)
          Df SumOfSqs      R2      F Pr(>F)  
Treatment  1  0.04793 0.06676 1.0014   0.05 *
Residual  14  0.67001 0.93324                
Total     15  0.71794 1.00000                
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


:::
:::




## Applying `adonis2 analysis to filtered datasets


::: {.cell}

```{.r .cell-code}
# Phyloseq objects to analyze
ps_list <- list(
  filtered_jbei = jbei_phyloseq,
  filtered_20 = drought_jbei_filtered_20,
  filtered_100 = drought_jbei_filtered_100
)

formulas <- list(
  full_model = "harvest * treatment",
  treatment_only = "treatment"
)
# Strata
strata_vars <- list(
  full_model = NULL,
  treatment_only = "harvest"
)

# Replicating Yen's results
run_adonis2_analysis(ps_jbei, formula = "Harvest * Treatment")
```

::: {.cell-output .cell-output-stdout}

```
$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
                  Df SumOfSqs      R2      F Pr(>F)    
Harvest            1  0.20456 0.28493 5.7995  0.001 ***
Treatment          1  0.04793 0.06676 1.3587  0.133    
Harvest:Treatment  1  0.04219 0.05876 1.1960  0.196    
Residual          12  0.42327 0.58956                  
Total             15  0.71794 1.00000                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


:::

```{.r .cell-code}
run_adonis2_analysis(ps_jbei, formula = "Treatment", strata = "Harvest")
```

::: {.cell-output .cell-output-stdout}

```
$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)
Treatment  1  0.04793 0.06676 1.0014  0.354
Residual  14  0.67001 0.93324              
Total     15  0.71794 1.00000              

$treatment_model_strata
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Blocks:  strata 
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", strata = strata_var, nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)  
Treatment  1  0.04793 0.06676 1.0014  0.058 .
Residual  14  0.67001 0.93324                
Total     15  0.71794 1.00000                
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


:::
:::




## PERMANOVA on filtered datasets



::: {.cell}

```{.r .cell-code}
#|
results_all <- map(ps_list, function(ps) {
  map(formulas, ~ run_adonis2_analysis(ps, formula = .x))
})


results_with_strata <- map(ps_list, function(ps) {
  map2(formulas, strata_vars, function(form, strat) {
    run_adonis2_analysis(ps, formula = form, strata = strat)
  })
})

# # Access results
results_all
```

::: {.cell-output .cell-output-stdout}

```
$filtered_jbei
$filtered_jbei$full_model
$filtered_jbei$full_model$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
                  Df SumOfSqs      R2      F Pr(>F)
harvest            1   0.3885 0.05375 0.8067  0.911
treatment          1   0.4920 0.06807 1.0215  0.342
harvest:treatment  1   0.5676 0.07853 1.1785  0.145
Residual          12   5.7793 0.79964              
Total             15   7.2274 1.00000              


$filtered_jbei$treatment_only
$filtered_jbei$treatment_only$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)
treatment  1   0.4920 0.06807 1.0226  0.349
Residual  14   6.7354 0.93193              
Total     15   7.2274 1.00000              



$filtered_20
$filtered_20$full_model
$filtered_20$full_model$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
                  Df SumOfSqs      R2      F Pr(>F)    
harvest            1  0.19525 0.30142 6.3055  0.001 ***
treatment          1  0.04307 0.06649 1.3909  0.136    
harvest:treatment  1  0.03788 0.05847 1.2232  0.221    
Residual          12  0.37158 0.57362                  
Total             15  0.64778 1.00000                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


$filtered_20$treatment_only
$filtered_20$treatment_only$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)
treatment  1  0.04307 0.06649 0.9972  0.377
Residual  14  0.60471 0.93351              
Total     15  0.64778 1.00000              



$filtered_100
$filtered_100$full_model
$filtered_100$full_model$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
                  Df SumOfSqs      R2      F Pr(>F)    
harvest            1  0.19770 0.37520 9.0381  0.001 ***
treatment          1  0.03611 0.06852 1.6506  0.116    
harvest:treatment  1  0.03063 0.05813 1.4002  0.175    
Residual          12  0.26248 0.49815                  
Total             15  0.52691 1.00000                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


$filtered_100$treatment_only
$filtered_100$treatment_only$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)
treatment  1  0.03611 0.06852 1.0299  0.355
Residual  14  0.49081 0.93148              
Total     15  0.52691 1.00000              
```


:::

```{.r .cell-code}
results_with_strata
```

::: {.cell-output .cell-output-stdout}

```
$filtered_jbei
$filtered_jbei$full_model
$filtered_jbei$full_model$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
                  Df SumOfSqs      R2      F Pr(>F)
harvest            1   0.3885 0.05375 0.8067  0.914
treatment          1   0.4920 0.06807 1.0215  0.347
harvest:treatment  1   0.5676 0.07853 1.1785  0.158
Residual          12   5.7793 0.79964              
Total             15   7.2274 1.00000              


$filtered_jbei$treatment_only
$filtered_jbei$treatment_only$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)
treatment  1   0.4920 0.06807 1.0226  0.374
Residual  14   6.7354 0.93193              
Total     15   7.2274 1.00000              

$filtered_jbei$treatment_only$treatment_model_strata
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Blocks:  strata 
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", strata = strata_var, nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)
treatment  1   0.4920 0.06807 1.0226  0.362
Residual  14   6.7354 0.93193              
Total     15   7.2274 1.00000              



$filtered_20
$filtered_20$full_model
$filtered_20$full_model$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
                  Df SumOfSqs      R2      F Pr(>F)    
harvest            1  0.19525 0.30142 6.3055  0.001 ***
treatment          1  0.04307 0.06649 1.3909  0.141    
harvest:treatment  1  0.03788 0.05847 1.2232  0.204    
Residual          12  0.37158 0.57362                  
Total             15  0.64778 1.00000                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


$filtered_20$treatment_only
$filtered_20$treatment_only$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)
treatment  1  0.04307 0.06649 0.9972  0.396
Residual  14  0.60471 0.93351              
Total     15  0.64778 1.00000              

$filtered_20$treatment_only$treatment_model_strata
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Blocks:  strata 
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", strata = strata_var, nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)  
treatment  1  0.04307 0.06649 0.9972  0.067 .
Residual  14  0.60471 0.93351                
Total     15  0.64778 1.00000                
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



$filtered_100
$filtered_100$full_model
$filtered_100$full_model$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
                  Df SumOfSqs      R2      F Pr(>F)    
harvest            1  0.19770 0.37520 9.0381  0.001 ***
treatment          1  0.03611 0.06852 1.6506  0.117    
harvest:treatment  1  0.03063 0.05813 1.4002  0.168    
Residual          12  0.26248 0.49815                  
Total             15  0.52691 1.00000                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


$filtered_100$treatment_only
$filtered_100$treatment_only$full_model
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)
treatment  1  0.03611 0.06852 1.0299  0.349
Residual  14  0.49081 0.93148              
Total     15  0.52691 1.00000              

$filtered_100$treatment_only$treatment_model_strata
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Blocks:  strata 
Permutation: free
Number of permutations: 999

vegan::adonis2(formula = formula, data = sample_data_df, method = method, by = "terms", strata = strata_var, nperm = nperm)
          Df SumOfSqs      R2      F Pr(>F)  
treatment  1  0.03611 0.06852 1.0299  0.031 *
Residual  14  0.49081 0.93148                
Total     15  0.52691 1.00000                
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


:::
:::