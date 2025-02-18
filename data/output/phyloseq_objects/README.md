# inter_brc_core_microbiome

**Core Microbiome Selection Summary**

The team evaluated multiple filtering strategies for our merged amplicon dataset, which originally contained over 34 million features and 218,608 ASVs. Initial filtering had already removed plant/mitochondrial DNA and singletons, resulting in 214,944 ASVs.

**Filtering Options Considered:**

1. **5% Sample Prevalence:**
    - Result: 1,849 ASVs
    - Note: Too restrictive given the diversity of crops and compartments.
2. **20 Reads Minimum:**
    - Result: ~58,961 ASVs
    - Proposal: Adequate representation (observing an ASV at least 20 times is considered sufficient).
3. **150 Reads Minimum:**
    - Result: 2,117 ASVs
    - Note: Also too strict.
4. **5% Prevalence + 20 Reads Minimum:**
    - Result: 1,086 ASVs
    - Note: Combined filtering further reduced the dataset, deemed overly punitive.

**Final Decision:**

- **20 Reads Minimum Filtering** was selected as the optimal approach. This threshold preserves a robust dataset of approximately 59K ASVs (with a phyloseq object size of ~0.9GB), ensuring sufficient data for downstream analyses while accommodating the heterogeneous nature of the samples.
- The approach avoids the pitfalls of prevalence-based filtering, which could disproportionately exclude data from less-sampled projects.

**Additional Resources:**

- **Core Microbiome Functions:**  
    Functions adapted from [Shade and Stopnisek (2019)](https://doi.org/10.1016/j.mib.2019.09.008) have been provided in the `core_microbiome_functions.R` script, available in the phyloseq folder.

To load the filtered phyloseq object in R, use:

```r
load("path/to/directory/filtered_phyloseq.rda")
```

### Removal of absolute paths
*Changes by Bolívar Aponte Rolón 2025-02-14*

Absolute paths were removed from:
```
R/
└── phyloseq_import.R
```
