# syntax=docker/dockerfile:1

# Template image
FROM rocker/r-ubuntu:22.04 

# Install system dependencies 
RUN apt-get update && apt-get install -y --no-install-recommends \
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
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create a directory for app installation
# Avoid installing to /root, /home, or /tmp
RUN mkdir -p /opt/interbrc-core-analysis

# Set working directory
WORKDIR /opt/interbrc-core-analysis

# Copy only the files needed for package restoration
#COPY renv.lock .

# Copy only necessary project files (lite version)
COPY Dockerfile .

# For full reproducibilty copy all project files
# COPY . /opt/interbrc-core-analysis/
COPY data/output/asv_matrices.rda /opt/interbrc-core-analysis/data/output/asv_matrices.rda
COPY data/output/phyloseq_objects/filtered_phyloseq.rda /opt/interbrc-core-analysis/data/phyloseq_objects/filtered_phyloseq.rda
COPY renv/ /opt/interbrc-core-analysis/renv/
COPY R/ /opt/interbrc-core-analysis/R/


# Install R packages globally (not in /root or /home)
# Using pak
#RUN Rscript -e "install.packages('pak'); pak::pkg_install(c('conflicted', 'styler', 'phyloseq', 'vegan', 'tidyverse', 'minpack.lm', 'Hmisc', 'stats4', 'germs-lab/BRCore@b391575', 'furrr', 'parallelly', 'doParallel', 'future', 'here'))"

# Install R packages using renv for better reproducibility
RUN Rscript -e "install.packages('renv', repos='https://cloud.r-project.org/')" \
    && Rscript -e "renv::init(bare=TRUE)" \
    && Rscript -e "options(renv.config.pak.enabled = TRUE)" \
    && Rscript -e "renv::install('pak')" \
    && Rscript -e "renv::install(c('conflicted', 'styler', 'phyloseq', 'vegan', 'tidyverse', 'minpack.lm', 'Hmisc', 'stats4', 'germs-lab/BRCore@b391575', 'furrr', 'parallelly', 'doParallel', 'future', 'here'))" 

# Run ldconfig to update the shared library cache
RUN ldconfig

# Set R library path explicitly to avoid relying on user's home
ENV R_LIBS_USER=/opt/interbrc-core-analysis/R/library

# Optional: Set a non-root user for better security
RUN useradd -m r-user && \
    chown -R r-user:r-user /opt/interbrc-core-analysis
USER r-user

# Default command when the container runs
CMD ["R", "--no-save"]