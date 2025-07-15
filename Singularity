Bootstrap: docker
From: rocker/r-ubuntu:22.04

%files
    # Copy any external files you need into the container here, if required
    data/output/asv_matrices.rda /opt/interbrc-core-analysis/data/output/asv_matrices.rda
    data/output/phyloseq_objects/filtered_phyloseq.rda /opt/interbrc-core-analysis/data/output/phyloseq_objects/filtered_phyloseq.rda
    renv/ /opt/interbrc-core-analysis/renv/
    R/ /opt/interbrc-core-analysis/R/
    renv.lock /opt/interbrc-core-analysis/

%post
    # Install system dependencies
    apt-get update && apt-get install -y --no-install-recommends \
        cmake \
        git \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libglpk-dev \
        libx11-dev \
        pandoc \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

    # Create project directory
    mkdir -p /opt/interbrc-core-analysis/
    
    # Set working directory
    cd /opt/interbrc-core-analysis

    # Install R package manager and use renv for reproducibility
    Rscript -e "install.packages('renv', repos='https://cloud.r-project.org/')"
    Rscript -e "renv::install('pak')"
    Rscript -e "options(renv.config.pak.enabled = TRUE)"
    Rscript -e "renv::restore()"
    ldconfig

    # Set permissions (optional, adjust as needed)
    useradd -m r-user && \
    chown -R r-user:r-user /opt/interbrc-core-analysis


%labels
    Author jibarozzo
    Version v1.0
    Description Singularity image built from Dockerfile instructions for interbrc-core-analysis

%runscript
    exec R --no-save