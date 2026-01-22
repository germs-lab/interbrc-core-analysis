# Refactoring Work Plan: BRC Core Analysis Functions

**Date Created:** 2025-01-20  
**Purpose:** Document the current state of `brc_*` functions and outline consolidation opportunities

---

## Current Function Inventory

### 1. Ordination Functions

#### `brc_nmds.R` - NMDS Ordination Calculation
- **Purpose:** Performs Hellinger transformation, Bray-Curtis distance calculation, and NMDS ordination
- **Key Parameters:**
  - `asv_matrix`: ASV count matrix (samples x ASVs)
  - `physeq`: phyloseq object with metadata
  - `k`: Number of dimensions (default: 2)
  - `trymax`: Maximum random starts (default: 100)
  - `ncores`: Parallel processing cores
- **Returns:** List with `nmds_scores`, `nmds_df`, `ordi`
- **Features:**
  - Input validation
  - NA/NaN checking
  - Convergence warnings
  - Species scores integration

#### `brc_pcoa.R` - PCoA Ordination Calculation
- **Purpose:** Principal Coordinates Analysis on distance matrices
- **Key Parameters:**
  - `asv_matrix`: Matrix or distance object
  - `physeq`: phyloseq object
  - `k`: Number of dimensions (default: 2)
  - `ncores`: Parallel processing cores
- **Returns:** List with `pcoa_df`, `ordi`, `distance_matrix`
- **Features:**
  - Accepts distance matrices directly
  - Uses `vegan::wcmdscale`
  - Simpler than NMDS (deterministic)

#### `brc_gg_ordi.R` - Standard Ordination Plot Visualization
- **Purpose:** Standardized ordination plotting with ggplot2
- **Key Parameters:**
  - `.data`: Data frame with ordination coordinates
  - `ordi`: Method type ("NMDS" or "PCoA")
  - `.color`, `.shape`: Aesthetic mappings
  - `color_values`, `.labels`: Custom styling
- **Returns:** ggplot object
- **Features:**
  - Automatic ellipse drawing
  - Reference lines at x=0, y=0
  - Flexible legend handling
  - NA filtering

#### `brc_flex_ordi.R` - Flexible Ordination Visualization
- **Purpose:** Plot dbRDA, capscale, PCoA, NMDS with biplot arrows
- **Key Parameters:**
  - `ordination_result`: vegan ordination object
  - `metadata_df`: Sample metadata
  - `color_var`, `shape_var`: Visual grouping
  - `biplot`: Show constraint arrows
  - `ellipse`: Draw confidence ellipses
- **Returns:** ggplot object
- **Features:**
  - Works with constrained ordinations (dbRDA, RDA)
  - Biplot arrows for constraints
  - Automatic variance explained labels

#### `brc_paper_ordinations.R` - Paper-Ready Ordination Plots
- **Purpose:** Generate matching NMDS and PCoA plots for publications
- **Key Parameters:**
  - `nmds_data`, `pcoa_data`: Pre-calculated ordination results
  - Separate colors/labels for each method
  - `crop_theme`, `brc_theme`: Custom theme lists
- **Returns:** Nested list with `nmds` and `pcoa` plots (crops and BRC)
- **Features:**
  - Creates 4 plots: NMDS crops/BRC, PCoA crops/BRC
  - Uses `brc_gg_ordi()` internally
  - Consistent styling across panels

---

### 2. Core Microbiome Analysis Functions

#### `brc_occore.R` - Core Microbiome Occupancy Calculation
- **Purpose:** Calculate core microbiome based on occupancy thresholds
- **Key Parameters:**
  - `physeq`: phyloseq object
  - `threshold`: Occupancy threshold (default varies)
- **Returns:** Core microbiome summary data
- **Features:**
  - Threshold-based filtering
  - Occupancy-abundance calculations

#### `brc_occ_curve.R` - Occupancy-Abundance Curve Plotting
- **Purpose:** Visualize occupancy-abundance relationships
- **Key Parameters:**
  - `core_summary_list`: Core summary data (list or data frame)
  - `point_size`, `point_alpha`, `point_shape`: Aesthetics
  - `color_values`: Custom colors
- **Returns:** ggplot object
- **Features:**
  - Flexible input types (list or data frame)
  - log10 transformation of abundance
  - Custom theme support

#### `brc_bc_curve.R` - Bray-Curtis Dissimilarity Curve Plotting
- **Purpose:** Plot Bray-Curtis dissimilarity curves
- **Key Parameters:** (Need to examine file for details)
- **Returns:** ggplot object
- **Features:** Similar to occupancy curves but for dissimilarity

---

### 3. Analysis Pipeline Functions

#### `brc_run_adonis2_analysis.R` - PERMANOVA Analysis
- **Purpose:** Run PERMANOVA (adonis2) on distance matrices
- **Key Parameters:**
  - Distance matrix
  - Formula for analysis
  - Permutation settings
- **Returns:** PERMANOVA results
- **Features:**
  - Wrapper for `vegan::adonis2`
  - Standardized output format

#### `brc_analyze_dataset.R` - Dataset Analysis Wrapper
- **Purpose:** High-level wrapper for dataset analysis
- **Key Parameters:** (Need to examine)
- **Returns:** Analysis results
- **Features:** Combines multiple analysis steps

#### `brc_analyze_read_patterns.R` - Read Pattern Analysis
- **Purpose:** Analyze read count patterns across samples
- **Key Parameters:** (Need to examine)
- **Returns:** Pattern analysis results

#### `brc_core_select_pipeline.R` - Core Selection Pipeline Wrapper
- **Purpose:** Complete pipeline for core microbiome selection
- **Key Parameters:** (Need to examine)
- **Returns:** Selected core microbiome
- **Features:** End-to-end workflow

---

### 4. Phyloseq Utility Functions

#### `brc_filter_phyloseq.R` - Phyloseq Filtering Utilities
- **Purpose:** Filter phyloseq objects by various criteria
- **Key Parameters:**
  - `physeq`: phyloseq object
  - Various filtering thresholds
- **Returns:** Filtered phyloseq object
- **Features:**
  - Multiple filtering strategies
  - Preserves phyloseq structure

---

### 5. Plotting and Output Functions

#### `brc_ggsave.R` - Batch Plot Saving Utility
- **Purpose:** Save multiple plots with consistent settings
- **Key Parameters:**
  - `objects`: Named list of ggplot objects
  - `save_path`: Directory path for saving
- **Returns:** NULL (saves files)
- **Features:**
  - Batch processing with `purrr::walk2`
  - Consistent dimensions (200x200mm)
  - 300 DPI PNG output

#### `brc_plot_asv_presence_threshold.R` - ASV Presence Threshold Plotting
- **Purpose:** Visualize ASV presence across thresholds
- **Key Parameters:** (Need to examine)
- **Returns:** ggplot object

#### `brc_plot_nonzero_stats.R` - Non-Zero Statistics Plotting
- **Purpose:** Plot statistics for non-zero ASV counts
- **Key Parameters:** (Need to examine)
- **Returns:** ggplot object

---

## Consolidation Opportunities

### 1. **Ordination Functions - HIGH PRIORITY**

**Problem:** Four separate functions with overlapping functionality
- `brc_nmds()` and `brc_pcoa()` both calculate ordinations
- `brc_gg_ordi()` and `brc_flex_ordi()` both visualize ordinations
- `brc_paper_ordinations()` wraps `brc_gg_ordi()`

**Proposed Solution:** Create unified interfaces

#### Option A: Consolidate Calculation Functions
```r
brc_ordination <- function(
  data,
  physeq,
  method = c("NMDS", "PCoA", "dbRDA", "capscale"),
  distance = "bray",
  transform = "hellinger",
  ...
)
```
- Single function with method parameter
- Maintains backward compatibility with aliases
- Reduces code duplication
- Standardizes parameter names

#### Option B: Keep Separate but Harmonize
- Standardize parameter names across functions
- Share common helper functions
- Ensure consistent return structure
- Better documentation

**Recommendation:** Option A - full consolidation
- Reduces maintenance burden
- Easier to add new methods
- More intuitive API

### 2. **Theme and Styling Functions - MEDIUM PRIORITY**

**Problem:** Theme definitions scattered across analysis scripts
- Title sizes defined inline (009_paper_figures.R lines 64-72)
- Scale helpers duplicated
- Label transformations repeated

**Proposed Solution:** Create `brc_themes.R`
```r
# Standard theme components
brc_theme_title_size()
brc_theme_paper()
brc_scale_y_percent()
brc_crop_labels()
```

**Benefits:**
- Consistent styling across plots
- Single source of truth
- Easy to update globally
- Reduces code in analysis scripts

### 3. **Plot Saving Functions - LOW PRIORITY**

**Problem:** Multiple ggsave calls with similar settings
- `brc_ggsave()` exists but limited
- Many manual ggsave calls in scripts

**Proposed Solution:** Enhance `brc_ggsave.R`
- Add preset configurations (paper, presentation, web)
- Support different file formats
- Better error handling
- Progress indicators for batch saves

### 4. **Analysis Pipeline Functions - LOW PRIORITY**

**Problem:** Pipeline functions may have overlapping steps
- Need to audit: `brc_analyze_dataset.R`, `brc_core_select_pipeline.R`
- Potential code duplication

**Proposed Solution:** 
- Document interdependencies
- Create modular helper functions
- Standardize input/output formats

---

## Optimization Recommendations

### Code Quality Improvements

1. **Roxygen Documentation**
   - ✅ Already present in most functions
   - 📝 Need to complete for: `brc_bc_curve.R`, pipeline functions
   - Add `@examples` sections
   - Document return structure more explicitly

2. **Input Validation**
   - ✅ Good validation in `brc_nmds`, `brc_pcoa`, `brc_gg_ordi`
   - 📝 Add validation to plotting functions
   - Standardize error messages (use `cli` package)

3. **Error Handling**
   - ✅ Using `cli::cli_abort()` in newer functions
   - 📝 Migrate older functions to `cli`
   - Add informative warnings for edge cases

4. **Parameter Naming Consistency**
   - ❌ Inconsistent: `.data` vs `data`, `physeq` vs `phyloseq_object`
   - 📝 Standardize to: `data`, `physeq`, `color`, `shape`
   - Use `.` prefix for NSE parameters only

5. **Return Structure Standardization**
   - ✅ Most functions return lists with named elements
   - 📝 Document return structure in roxygen
   - Consider S3 classes for complex returns

### Performance Improvements

1. **Parallel Processing**
   - ✅ Already implemented in `brc_nmds` and `brc_pcoa`
   - 📝 Use helper: `get_available_cores()` consistently
   - Add progress bars for long-running operations

2. **Memory Optimization**
   - 📝 Consider chunked processing for large datasets
   - Test with high-memory datasets (>32GB requirement noted in script)

3. **Code Efficiency**
   - 📝 Profile functions to identify bottlenecks
   - Vectorize where possible
   - Reduce intermediate copies

---

## Implementation Roadmap

### Phase 1: Foundation (Week 1)
- [x] Create work plan document
- [ ] Create `brc_themes.R` with reusable theme components
- [ ] Create `brc_ordination.R` with unified interface
- [ ] Add comprehensive tests

### Phase 2: Script Refactoring (Week 2)
- [ ] Refactor `009_paper_figures.R` to use new functions
- [ ] Fix syntax errors (lines 408-411, 486-489)
- [ ] Validate output equivalence

### Phase 3: Documentation (Week 3)
- [ ] Complete roxygen documentation
- [ ] Add usage vignettes
- [ ] Create function cross-reference

### Phase 4: Extended Consolidation (Future)
- [ ] Audit pipeline functions for consolidation
- [ ] Enhance `brc_ggsave` with presets
- [ ] Performance profiling and optimization

---

## Testing Strategy

### Unit Tests
- Test each function in isolation
- Validate input checking
- Test error conditions

### Integration Tests
- Test function combinations
- Validate workflow outputs
- Compare refactored vs original

### Regression Tests
- Ensure refactored code produces identical results
- Use snapshot testing for plots
- Numerical tolerance for floating-point comparisons

---

## Backward Compatibility

### Deprecation Strategy
1. Keep old function names as aliases
2. Add deprecation warnings with `lifecycle` package
3. Document migration path
4. Provide codemod scripts for bulk updates

### Example:
```r
#' @export
#' @rdname brc_ordination
brc_nmds <- function(...) {
  lifecycle::deprecate_warn("2.0.0", "brc_nmds()", "brc_ordination()")
  brc_ordination(..., method = "NMDS")
}
```

---

## Success Metrics

- [ ] Reduced lines of code by 30%
- [ ] All functions have >80% code coverage
- [ ] Documentation coverage at 100%
- [ ] No breaking changes to existing analysis scripts
- [ ] Faster execution time for ordinations (measure baseline first)

---

## Notes and Considerations

### Dependencies
Current key dependencies:
- `phyloseq`: Core data structure
- `vegan`: Ordination and diversity
- `ggplot2`: Visualization
- `dplyr`/`tidyr`: Data manipulation
- `cli`: Error messaging
- `rlang`: NSE support

### Memory Requirements
- `009_paper_figures.R` notes >32GB RAM requirement for relative abundance section
- Consider streaming/chunked processing for large datasets

### Code Style
- Following tidyverse style guide
- Using `purrr` for functional programming
- Embracing NSE where appropriate (ggplot2 aesthetics)

---

**Last Updated:** 2025-01-20  
**Contributors:** Bolívar Aponte Rolón (original code), AI Assistant (refactoring plan)
