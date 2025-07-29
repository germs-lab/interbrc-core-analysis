# syntax=docker/dockerfile:1
# see https://jnicolaus.com/tutorials/2022-04-23-renv-docker-singularity/
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
    libx11-dev \
    pandoc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create a directory for app installation and renv cache
RUN mkdir -p /opt/interbrc-core-analysis \
    && mkdir -p /opt/renvcache \
    && mkdir -p /opt/Rlibsymlinks 

# Set renv environment variables in Renviron.site
RUN echo "RENV_PATHS_CACHE='/opt/renvcache'" >> $(R RHOME)/etc/Renviron.site


# Set working directory
WORKDIR /opt/interbrc-core-analysis

# Copy only the files needed for package restoration
COPY renv.lock .
COPY Dockerfile .
COPY data/output/asv_matrices.rda /opt/interbrc-core-analysis/data/output/asv_matrices.rda
COPY data/output/phyloseq_objects/filtered_phyloseq.rda /opt/interbrc-core-analysis/data/output/phyloseq_objects/filtered_phyloseq.rda
COPY renv/ /opt/interbrc-core-analysis/renv/
COPY R/ /opt/interbrc-core-analysis/R/


# Install R packages using renv for better reproducibility
RUN Rscript -e "install.packages('renv', repos='https://cloud.r-project.org/')" 
RUN cd /opt/interbrc-core-analysis \
    && Rscript -e "options(renv.config.pak.enabled = TRUE)" \
    && Rscript -e "renv::install('pak@0.9.0')" \
    && Rscript -e "renv::restore()" \
    && mv $(Rscript -e "cat(paste0(renv::paths\$library()))") /opt/Rlibsymlinks \
    && echo "R_LIBS=/opt/Rlibsymlinks" >> $(R RHOME)/etc/Renviron.site \
    ldconfig 
    # Run ldconfig to update the shared library cache

# Set permissions for the non-root user
RUN useradd -m r-user \
    && chown -R r-user:r-user /opt/interbrc-core-analysis \
    && chown -R r-user:r-user /opt/renvcache \
    && chown -R r-user:r-user /opt/Rlibsymlinks

USER r-user

# Default command when the container runs
CMD ["R", "--no-save"]