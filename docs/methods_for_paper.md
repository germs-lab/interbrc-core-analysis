

### Statistical method
Post-DADA2 pipeline and harmonization (this is something to ask Brandon)
 
We implemented Shade and Stopnisek’s [@shade2019] method for core microbiome identification with the "BRCore" package [@aponterolon2026] which uses abundance-occupancy distributions and defines a species core based on its contribution to beta diversity (e.g. Bray-Curtis dissimilarity) when compared to the entire dataset. We used a feedstock type as the prioritizing variable for occurrence detection. Briefly, we used the "last percent Bray-Curtis increase" parameter, which selects the final rank for which adding one more taxon still increases the explained Bray-Curtis similarity by at least `p`, where `p` is your chosen percent threshold (See further details in [@shade2019]). 

Additionally, we selected were 100 sample reads and 20 ASV reads to balance data retention while suppressing low-abundance noise. From this, we selected a microbial core based on a 60% occupancy across samples and feedstock types. We then used a custom wrapper function that employs `decostand()`, `vegdist()`, and `metaMDS()` from the `vegan` R package [@oksanen2022] to perform a Hellinger transformation on the selected cores' abundance matrices, calculate the Bray-Curtis dissimilarity matrices, and finally calculate a Non-Metric Multidimensional Scaling (NMDS) analysis. We optimized the NMDS setting parameters to `k`= 2 dimension, `trymax`= 100, and `maxit` = 999  [@oksanen2022].  We also utilized a custom function wrapping `wcmdscale` for Principal Coordinate Analysis (PCoA) and optimized setting `k` = 2 [@oksanen2022].  

Lastly, on a subset of data focuse on samples from SABR, we performed a clustering-threshold analysis to quantify how the proportion of sequences attributed to the “core” changes with OTU clustering stringency and the core-occurrence criterion. We determined occupancy thresholds in 5% increment starting t at 50% occupancy and  85%  to 100% clustering thresholds.  All analyses were performed using  *R* [v.4.4.1; @rcoreteam2024]. Visualization were generated with *ggplot2* package [@wickham2016]



Results


## Methods

### Data preprocessing
To reduce low-abundance noise while retaining sufficient signal we selected a threshold of 100 reads per sample and 20 reads per ASV. All downstream analyses were performed on data filtered to these thresholds.

### Core microbiome identification 

#### Beta Diversity & Occupancy-based core
We identified core taxa following the abundance–occupancy framework of Shade and Stopnisek [@shade2019], implemented with the BRCore package [@aponterolon2026]. Briefly, BRCore ranks taxa by their contribution to Bray–Curtis similarity and applies the “last percent Bray–Curtis increase” criterion: taxa are included up to the final rank at which adding one more taxon still increases the explained Bray–Curtis similarity by at least `p`. We set p = 2%, consistent with prior applications of the method. Occurrence detection was prioritized by feedstock (crop) category to ensure core membership reflects consistency across feedstock types.

Complementary to the Bray–Curtis contribution approach, we defined a threshold-based core using occupancy across samples and feedstock types. Taxa present in at least 60% of samples (≥60% occupancy) were considered members of the occupancy-based core. We use both definitions in the manuscript: the Bray–Curtis contribution “core” to capture taxa most influential to community dissimilarity, and the occupancy “core” to capture widespread taxa across feedstocks.

### Dissimilarity and ordination analyses
For ordination analyses, we applied the following workflow using functions from the vegan R package [@oksanen2022]: (1) Hellinger transformation of abundance matrices with `decostand(method = "hellinger")`; (2)  Bray–Curtis dissimilarity computation with `vegdist(method = "bray")`; (3) Non-metric multidimensional scaling (NMDS) with `metaMDS(k = 2, trymax = 100)`. We set the dimensionality to k = 2 and increased the optimization effort with `trymax = 100`; and allowed up to ~1,000 iterations (`maxit` ≈ 999); (4) Principal Coordinates Analysis (PCoA) using a custom wrapper around `wcmdscale` with k = 2.

Ordination scores were visualized by feedstock and BRC, using consistent aesthetics and legends across figures to facilitate comparison. Note that ordination axes display ordination scores (not percentages).

### Clustering-threshold analysis (SABR subset)
On the SABR subset, we assessed how the fraction of sequences attributed to the “core” changes under increasing OTU clustering stringency and occupancy criteria. We examined occupancy thresholds in 5% increments, starting at 50%. Clustering thresholds ranged from 85% to 100%. These analyses quantify the sensitivity of core membership to clustering granularity and to the occupancy cutoff.

### Software and reproducibility
All analyses were conducted in R v4.4.1 [@rcoreteam2024]. Visualizations were produced with ggplot2 [@wickham2016], with standardized themes for comparability across figures. Core identification and supporting utilities were implemented with BRCore [@aponterolon2026] and vegan [@oksanen2022].


## Results

Th selection 100 reads per sample and 20 reads per ASV resulted in retention of 1813 samples and 59961 AVS across BRCs. Through our beta diversity core analysis we found 50 "core" ASVs that contributed at least 2% to beta diversity across feedstocks and BRCs. Meanwhile, with out abundance-occupancy threshold approach we found only 12 ASVs present across >= 60% of samples. 


