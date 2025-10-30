FROM mambaorg/micromamba:1.5.8-focal

WORKDIR /app

# System libraries needed for R & Python builds
USER root
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Copy and build environment as root to avoid permission errors
COPY env/environment.yml /app/env/environment.yml
USER root
RUN micromamba create -y -n gene_nlp -f /app/env/environment.yml && \
    micromamba clean --all --yes
USER $MAMBA_USER

# Pre-install R packages
COPY R/install.R /app/R/install.R

USER root
RUN micromamba run -n gene_nlp bash -c "Rscript /app/R/install.R"
USER $MAMBA_USER

USER root
RUN echo '.libPaths("/opt/conda/envs/gene_nlp/lib/R/library")' | tee -a /opt/conda/envs/gene_nlp/lib/R/etc/Rprofile.site > /dev/null
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
RUN chmod +x /app/scripts/run.sh

CMD micromamba run -n gene_nlp ./scripts/run.sh
