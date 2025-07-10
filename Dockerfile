# syntax=docker/dockerfile:1

# Template image
FROM rocker/r-ubuntu:22.04 

# Create and set working directory
WORKDIR /workspace

# Copy only the files needed for package restoration
COPY renv.lock ./

# Install system dependencies
RUN apt-get update && \
    apt-get install -y \
        cmake \
        git \
        libcurl4-openssl-dev \
        libssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libgit2-dev \
        libglpk-dev \
        libx11-dev \
        libxml2-dev \
        libjpeg-turbo8-dev \
        # Add any other system dependencies your R packages need
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

# Install R packages using pak
RUN Rscript -e "install.packages('pak'); pak::pkg_install(c('conflicted', 'styler', 'phyloseq', 'vegan', 'tidyverse', 'minpack.lm', 'Hmisc', 'stats4', 'vmikk/metagMisc', 'germs-lab/BRCore@b391575', 'furrr', 'parallelly', 'doParallel', 'future'))"

# Alternative: Use renv for better reproducibility
# RUN Rscript -e "install.packages('renv', repos='https://cloud.r-project.org/')" && \
#     Rscript -e "renv::restore()"

# Copy the rest of the project files
COPY . .

# Optional: Set a non-root user for better security
RUN useradd -m r-user && \
    chown -R r-user:r-user /workspace
USER r-user

# Set default command
CMD ["R", "--no-save"]