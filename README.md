# Inter-BRC Core Analysis Singularity Container

This Docker/Singularity container enables reproducible execution of microbiome core analysis workflows (scripts 004-007) in HPC environments. These scripts can also be deployed individually on an HPC or powerful local machine. 

## Purpose & Scope

### What This Container Does

- Reproducible execution of microbiome analysis workflows
- Isolation of R package dependencies using `renv`
- Compatibility with HPC environments that prefer Singularity
- Focused analysis on core microbiome selection and ordination analyses

### Container Scope

This container is **specifically designed** for scripts 004-007 series and includes:

- Core microbiome extraction algorithms
- Ordination analysis functions
- Required phyloseq objects and ASV matrices
- Custom BRC analysis functions

### Supported Scripts

- `004_core_selection.R` & `004_core_selection_HPC.R` - Core microbiome identification using Bray-Curtis dissimilarity
- `005_core_selection_per_BRC.R` - Identifies the core microbiome across samples for individual BRCs
- `006_save_main_physeqs.R` - Export phyloseq objects and corresponding FASTA files
- `007_ordinations.R` - Community ordination analysis (NMDS, PCoA, dbRDA)
- `007_ordinations_full.R` - Full dataset NMDS analysis for HPC environments

## What's Included

- R 4.4+ with bioinformatics packages
- Pre-computed phyloseq objects and ASV matrices
- Custom BRC analysis functions
- All required system dependencies (cmake, git, libcurl4-openssl-dev, libssl-dev, libxml2-dev, etc.)
- R Package Management via `renv` for reproducible package management

## Quick Start

1. **Download container:**

```bash
singularity pull docker://ghcr.io/germs-lab/interbrc-core-analysis/interbrc-lite-container:v4
```
    
2. **Run analysis:**

```bash
singularity exec --no-home --pwd /opt/interbrc-core-analysis interbrc-lite_v4.sif Rscript "R/analysis/004_core_selection.R"
```
    

## Usage Instructions

Docker and Singularity differ in their file structure and library path management. Singularity automatically binds host directories, causing R to look for packages on host system. To run analyses in an isolated environment, we avoid mounting local machine or HPC home directories.

Here we bind host file system to the containers file system to ensure that the analyses results save to the host system. See documentation on using `--bind` and ``no-home` flags [here](https://docs.sylabs.io/guides/4.3/user-guide/bind_paths_and_mounts.html#using-bind-or-mount-with-the-writable-flag).

### Example 1: Core Selection Analysis

```bash
singularity exec -bind /path/to/your/output:/opt/interbrc-core-analysis/data/output \
 --no-home --pwd /opt/interbrc-core-analysis interbrc-lite_v4.sif Rscript "R/analysis/004_core_selection.R"
```

### Example 2: Full Ordination Analysis (HPC)

```bash
singularity exec -bind /path/to/your/output:/opt/interbrc-core-analysis/data/output \
--no-home --pwd /opt/interbrc-core-analysis interbrc-lite_v4.sif Rscript "R/analysis/007_ordinations_full.R"
```

### Example 3: Multiple Bind Mounts
```bash
singularity exec \
  --bind /path/to/your/output:/opt/interbrc-core-analysis/data/output \
  --bind /path/to/your/plots:/opt/interbrc-core-analysis/data/output/plots \
  --no-home --pwd /opt/interbrc-core-analysis \
  interbrc-lite_v4.sif Rscript "R/analysis/007_ordinations_full.R"
```
### Example 4 HPC Example
```bash
# On HPC, typically you'd bind your work directory
singularity exec --bind $PWD/output:/opt/interbrc-core-analysis/data/output \
  --no-home --pwd /opt/interbrc-core-analysis \
  interbrc-lite_v4.sif Rscript "R/analysis/007_ordinations_full.R"
```  
Key Point: Without bind mounts, all output stays inside the read-only container and is lost when the container stops running.

## Build Process

_**You can build the Docker container then convert to Singularity (preferable if you want to modify components) or build your Singularity container directly.**_

### Method 1: Pre-built Container (Recommended for Users)

**Direct download:**

```bash
singularity pull docker://ghcr.io/germs-lab/interbrc-core-analysis/interbrc-lite-container:v4
```

**With authentication (if required):**

```bash
export SINGULARITY_DOCKER_USERNAME=username
export SINGULARITY_DOCKER_PASSWORD=password
singularity build interbrc-lite_v4.sif docker://ghcr.io/germs-lab/interbrc-core-analysis/interbrc-lite-container:v4
```

### Method 2: Docker → Singularity (For Customization)

**Local build:**

```bash
# Clone repository first
docker build -t interbrc-lite-container:v4 .
# Convert to Singularity
singularity build interbrc-lite_v4.sif docker-daemon://interbrc-lite-container:v4
```

**Remote build:**

```bash
docker pull ghcr.io/germs-lab/interbrc-core-analysis/interbrc-lite-container:v4
singularity build interbrc-lite_v4.sif docker-daemon://interbrc-lite-container:v4
```

### Method 3: Direct Singularity Build

_**Remember to login to Singularity (see [Singularity remote login](https://docs.sylabs.io/guides/4.3/user-guide/cli/singularity_remote_login.html) & [build environment](https://docs.sylabs.io/guides/4.3/user-guide/build_env.html#docker))**_

```bash
export SINGULARITY_DOCKER_USERNAME=username
export SINGULARITY_DOCKER_PASSWORD=password
singularity build interbrc-lite_v4.sif docker://ghcr.io/germs-lab/interbrc-core-analysis/interbrc-lite-container:v4
```

## Container Structure


```
/opt/interbrc-core-analysis/
├── Dockerfile
├── renv.lock
├── data/output/
│   ├── asv_matrices.rda
│   └── phyloseq_objects/filtered_phyloseq.rda
├── R/
│   ├── analysis/
│   ├── functions/
│   ├── references/
│   └── utils/
└── renv/
/opt/renvcache/          # renv package cache
/opt/Rlibsymlinks/       # R library symlinks
```

## Troubleshooting

### Library Path Issues

If R cannot find packages, ensure you're using isolation flags:

```bash
singularity exec --no-home --pwd /opt/interbrc-core-analysis container.sif Rscript script.R
```

### Environment Variables (if needed)

```bash
SINGULARITYENV_RENV_PATHS_CACHE=/opt/renvcache \
SINGULARITYENV_R_LIBS=/opt/Rlibsymlinks \
singularity exec --no-home --pwd /opt/interbrc-core-analysis container.sif Rscript script.R
```

### Authentication Issues

For private repositories, ensure proper Docker credentials are set before building.

## Important Notes

- Container is read-only once built as Singularity .sif file
- Designed for specific analysis workflows, not general-purpose R environment
- Requires adequate computational resources for ordination analyses
- Best suited for HPC environments with Singularity support
- Docker container runs isolated with its own file system
- The file system from Docker to Singularity requires proper `renv` cache and library path management

---

_Last updated: 2025-07-30_ _Maintained by: @jibarozzo_