# ------------------------------------------------------------------------------
# envs/tensorqtl_cpu/Dockerfile
#
# Description: CPU-only tensorQTL image for LOCAL and CI test runs of the
#   tensorQTL tasks (set tensorqtl_docker_image to it). Terra runs keep the AoU
#   default gcr.io/broad-cga-francois-gtex/tensorqtl:latest (GPU, 22 GB), whose
#   versions are pinned here: tensorQTL 1.0.9, torch 2.1.2, Python 3.10,
#   numpy 1.26.3, pandas 2.1.4, pyarrow 14.0.2, scipy 1.11.4; Pgenlib 0.95.1
#   (the Broad 0.90.1 does not compile on arm64; same reader API).
#   Why: the Broad image cannot run on Apple Silicon (PyTorch segfaults under
#   qemu x86 emulation) and does not fit on a GitHub runner. Without rpy2,
#   tensorQTL skips its own q-values (they are not used by the pipeline).
#   Builds natively for linux/amd64 and linux/arm64.
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code: python:3.10-slim + compilers (Pgenlib builds from source on
#   arm64) -> pip install pinned packages -> remove compilers -> import check.
#
# Build: docker buildx build -f envs/tensorqtl_cpu/Dockerfile -t tensorqtl_cpu:local --load .
# ------------------------------------------------------------------------------
FROM python:3.10-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential zlib1g-dev && \
    pip install --no-cache-dir \
      "numpy==1.26.3" "pandas==2.1.4" "scipy==1.11.4" "pyarrow==14.0.2" \
      "Pgenlib==0.95.1" "pandas-plink" && \
    pip install --no-cache-dir --extra-index-url https://download.pytorch.org/whl/cpu "torch==2.1.2" && \
    pip install --no-cache-dir --no-deps "tensorqtl==1.0.9" && \
    pip install --no-cache-dir "qtl" && \
    apt-get purge -y build-essential zlib1g-dev && apt-get autoremove -y && rm -rf /var/lib/apt/lists/* && \
    python3 -c "import tensorqtl, torch, pgenlib; from tensorqtl import cis, genotypeio, pgen; print(torch.__version__)"

CMD ["bash"]
