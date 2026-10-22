# caEditR_implementation

[![Docker Image CI](https://github.com/haimkru/caEditR_implementation/actions/workflows/docker-image.yml/badge.svg)](https://github.com/haimkru/caEditR_implementation/actions/workflows/docker-image.yml)
[![Pipeline CI](https://github.com/haimkru/caEditR_implementation/actions/workflows/pipeline-ci.yml/badge.svg)](https://github.com/haimkru/caEditR_implementation/actions/workflows/pipeline-ci.yml)

Terra (WDL 1.0) pipelines for **cell-type-specific RNA-editing QTLs (edQTLs)**
in All of Us whole-blood RNA-seq, built on
[caEditR](https://github.com/haimkru/caEditR) (version 0.99.3). Three approaches:

1. **Bulk edQTLs** - tensorQTL on the quality-controlled bulk editing ratios
   (standard edQTL pipeline, Li et al. 2022).
2. **tensorQTL on cell-type matrices** - per-sample cell-type editing from
   `caNRD_edit()` (caEditR) and `TCA_Like()` (CRAN TCA), each mapped with tensorQTL.
3. **caNRD-editQTL** - end-to-end joint test on the bulk ratios
   (`caNRD_editQTL()`, scan engine + exact lead-variant re-fit), which avoids
   the cross-cell-type leakage of approach 2.

## Workflows

| Dockstore name | WDL | Purpose |
| --- | --- | --- |
| `editing_edqtl_pipeline` | [workflows/editing_edqtl_pipeline.wdl](workflows/editing_edqtl_pipeline.wdl) | All steps in one submission |
| `prepare_editing_phenotypes` | [workflows/prepare_editing_phenotypes.wdl](workflows/prepare_editing_phenotypes.wdl) | Index inputs, cohort, caNRD-edit + TCA per chunk, merge, tensorQTL phenotype sets |
| `tensorqtl_edqtl` | [workflows/tensorqtl_edqtl.wdl](workflows/tensorqtl_edqtl.wdl) | tensorQTL cis nominal (optionally permutations) on every phenotype set |
| `canrd_editqtl` | [workflows/canrd_editqtl.wdl](workflows/canrd_editqtl.wdl) | caNRD-editQTL cis scan per chunk + aggregation |
| `summarize_edqtl` | [workflows/summarize_edqtl.wdl](workflows/summarize_edqtl.wdl) | Side-by-side comparison of the three approaches |

Tasks are in `workflows/tasks/`. The tensorQTL tasks are copies of
[AoU-Multiomics-Analysis/tensorQTL](https://github.com/AoU-Multiomics-Analysis/tensorQTL)
(see the header of `workflows/tasks/tensorqtl.wdl`).

## Inputs (formats: [docs/data-dictionary.md](docs/data-dictionary.md))

Editing ratio + coverage BED matrices (`site_id` = `chrom:pos:strand`),
linear-CPM expression BED, GTF, CellTypeDeconvolution
`proportions_tca_weights.tsv`, PLINK2 genotypes, covariates TSV
(`sample_id`, genotype PCs, sex, age), optional sample list.

## Scaling: chunks

The work is split into **chunks**, and each chunk runs on its own Terra VM at
the same time as the others (a WDL scatter).

- **A chunk is a number of editing sites, not base pairs.** The sites are
  sorted by position, and each chunk is the next `sites_per_chunk` consecutive
  sites on one chromosome (default 500). Its genomic span is whatever those
  sites cover: a few Mb to ~15 Mb at AoU density. A chunk never crosses a
  chromosome and never splits two sites at the same position, so the last
  chunk of a chromosome is usually smaller.
- **What a chunk VM reads.** It reads only its sites' rows from the bgzipped,
  tabix-indexed editing-ratio and coverage matrices. For caNRD-editQTL it also
  extracts, with `plink2`, only the variants within +/- 1 Mb of those sites
  (first site - 1 Mb to last site + 1 Mb).
- **Chunking does not change the results.** Each site is still tested against
  exactly the variants within +/- 1 Mb of that site (`cis_window`), the same
  window tensorQTL uses. The one exception: TCA fits its shared noise term
  within each chunk.
- **Which stages are chunked.** DeconvolveChunk (caNRD-edit + TCA) and
  CanrdEditqtlChunk. tensorQTL is not chunked: it runs one job per phenotype
  set (bulk, or one method x cell type) over the whole genome.
- **Why it scales.** Memory and time per chunk depend on the chunk size and
  the number of samples, not on the total number of sites or variants. 10k
  samples x 100k sites x 1M variants is ~200 chunks per chunked stage, each
  needing ~4 GB (deconvolution) or ~12 GB (caNRD-editQTL), measured at 10k
  samples. Fewer, larger chunks mean fewer VM start-ups; smaller chunks mean
  less memory per VM. The pilot uses 5000 sites per chunk.

Resources, time and cost:
[docs/TERRA_RUN_GUIDE.md](docs/TERRA_RUN_GUIDE.md#5-scale-and-resources-10k-samples-100k-sites-1m-variants).

## Repository layout

```
R/            shared R helpers (I/O, cohort, site QC, chunks, INT/PCs, genotypes, statistics)
scripts/      one command-line script per WDL task (R or bash)
workflows/    WDL entry points; workflows/tasks/ task definitions
envs/         Dockerfile + conda environment of ghcr.io/haimkru/caeditr_implementation
tests/        R unit tests; fixtures/sim_small = small simulated cohort for CI
tools/        lint and WDL checks (logging fields, Terra workflow-scope file writes)
docs/         Terra run guide, data dictionary
examples/     Terra data-table template (TSV)
```

## Docker image

`envs/Dockerfile` (micromamba, R 4.5.3, caEditR from GitHub at a pinned
commit, TCA 1.2.1, htslib, plink2). Published by `.github/workflows/docker-image.yml`
as `ghcr.io/haimkru/caeditr_implementation:sha-<commit>` / `latest`.
Build locally: `docker buildx build --platform linux/amd64 -f envs/Dockerfile -t caeditr_implementation:local --load .`

`envs/tensorqtl_cpu/Dockerfile`: CPU tensorQTL 1.0.9 image for local and CI
runs of the tensorQTL tasks (`tensorqtl_docker_image`); Terra uses the Broad
GPU image `gcr.io/broad-cga-francois-gtex/tensorqtl:latest`.

Author: Haim Krupkin+Claude. Contact: hkrupkin@stanford.edu
