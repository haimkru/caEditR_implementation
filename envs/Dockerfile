# ------------------------------------------------------------------------------
# envs/Dockerfile
#
# Description: Docker image used by every caEditR_implementation WDL task except the
#   tensorQTL tasks (those use gcr.io/broad-cga-francois-gtex/tensorqtl).
#   Published by .github/workflows/docker-image.yml as
#   ghcr.io/haimkru/caeditr_implementation:<tag>.
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   1. micromamba base image; install envs/environment.yml.
#   2. git clone caEditR at a pinned commit (CAEDITR_COMMIT), check the commit,
#      R CMD INSTALL it.
#   3. Copy the pipeline's R/ helpers and scripts/ to /opt/caeditr_implementation.
#   4. Smoke checks: required R packages load, tools are on PATH.
#
# Build locally (Terra machines are linux/amd64):
#   docker buildx build --platform linux/amd64 -f envs/Dockerfile -t caeditr_implementation:local --load .
# ------------------------------------------------------------------------------
FROM mambaorg/micromamba:2.9.0-ubuntu22.04

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH
ENV R_LIBS_USER=/opt/conda/lib/R/library

# caEditR 0.99.3 (github.com/haimkru/caEditR, branch bioconductor); bump both
# values together when a new caEditR version should be used.
ARG CAEDITR_REPOSITORY=https://github.com/haimkru/caEditR.git
ARG CAEDITR_COMMIT=812bf93c9debf5fceddc5105cea1c1cf0a4a26c9
ARG CAEDITR_VERSION=0.99.3
RUN git clone --quiet "$CAEDITR_REPOSITORY" /tmp/caEditR && \
    git -C /tmp/caEditR checkout --quiet "$CAEDITR_COMMIT" && \
    test "$(git -C /tmp/caEditR rev-parse HEAD)" = "$CAEDITR_COMMIT" && \
    R CMD INSTALL --no-test-load /tmp/caEditR && \
    rm -rf /tmp/caEditR && \
    Rscript -e "stopifnot(as.character(packageVersion('caEditR')) == '$CAEDITR_VERSION')"

COPY --chown=$MAMBA_USER:$MAMBA_USER R /opt/caeditr_implementation/R
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts /opt/caeditr_implementation/scripts

RUN Rscript -e 'for (p in c("caEditR", "TCA", "data.table", "optparse", "nnls", "irlba", "qvalue")) \
      suppressPackageStartupMessages(library(p, character.only = TRUE))' && \
    Rscript -e 'stopifnot(as.character(packageVersion("TCA")) == "1.2.1")' && \
    bgzip --version > /dev/null && tabix --version > /dev/null && plink2 --version && \
    micromamba env export --name base > /opt/caeditr_implementation/conda_environment_lock.yml

WORKDIR /workspace
