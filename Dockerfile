# syntax=docker/dockerfile:1

# Template image
FROM rocker/r-ubuntu:22.04 

# Create and set working directory
WORKDIR /workspace

# Copy only the files needed for package restoration
COPY renv.lock ./

# Install renv and restore packages
# System dependencies are handled by `pak` package
RUN Rscript -e "renv_1.1.1 <- 'https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.1.tar.gz'; install.packages(renv_1.1.1, repos=NULL, type = 'source')" && \
    Rscript -e "options(renv.consent = TRUE, renv.config.pak.enabled = TRUE); renv::restore(lockfile = 'renv.lock', clean = TRUE)"


# Copy the rest of the project files
COPY . .

# Optional: Set a non-root user for better security
RUN useradd -m r-user \
    && chown -R r-user:r-user /workspace
USER r-user

# Set default command
CMD ["R"]