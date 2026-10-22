# caEditR_edQTL

<!-- workflow-badges:start -->
<!-- workflow-badges:end -->

Terra (WDL 1.0) pipelines for **cell-type-specific RNA-editing QTLs (edQTLs)**
in All of Us whole-blood RNA-seq, built on
[caEditR](https://github.com/haimkru/caEditR) (pinned: 0.99.3, commit
`812bf93`). Three approaches, as in the caEditR manuscript (Methods 4.6):

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

## Scaling

Sites are split into genomic chunks (`sites_per_chunk`, default 500) that are
scattered; every chunk reads only its region of the bgzipped, tabix-indexed
editing matrices and extracts only its cis variants from the PLINK2 files with
`plink2`. Memory per shard is independent of the number of sites and variants;
10k samples x 100k sites x 1M variants = ~200 shards per stage. See
[docs/TERRA_RUN_GUIDE.md](docs/TERRA_RUN_GUIDE.md#5-scale-and-resources-10k-samples-100k-sites-1m-variants).

## Repository layout

```
R/            shared R helpers (I/O, cohort, site QC, chunks, INT/PCs, genotypes, statistics)
scripts/      one command-line script per WDL task (R or bash)
workflows/    WDL entry points; workflows/tasks/ task definitions
envs/         Dockerfile + conda environment of ghcr.io/aou-multiomics-analysis/caeditr_edqtl
tests/        R unit tests; fixtures/sim_small = small simulated cohort for CI
tools/        lint and WDL checks (logging fields, Terra workflow-scope file writes)
docs/         Terra run guide, data dictionary
examples/     Terra data-table template (TSV)
```

## Running

- Terra: [docs/TERRA_RUN_GUIDE.md](docs/TERRA_RUN_GUIDE.md).
- Locally (miniwdl + Docker, same images): `miniwdl run workflows/<wdl> name=value ...`;
  the development project drives this with one script per step.
- CI (`.github/workflows/pipeline-ci.yml`): miniwdl + womtool checks, R lint and
  tests in the image, end-to-end run on `tests/fixtures/sim_small`
  (tensorQTL with the CPU image; the 22 GB GPU image does not fit a GitHub runner).

## Docker image

`envs/Dockerfile` (micromamba, R 4.5.3, caEditR from GitHub at a pinned
commit, TCA 1.2.1, htslib, plink2). Published by `.github/workflows/docker-image.yml`
as `ghcr.io/aou-multiomics-analysis/caeditr_edqtl:sha-<commit>` / `latest`.
Build locally: `docker buildx build --platform linux/amd64 -f envs/Dockerfile -t caeditr_edqtl:local --load .`

`envs/tensorqtl_cpu/Dockerfile`: CPU tensorQTL 1.0.9 image for local and CI
runs of the tensorQTL tasks (`tensorqtl_docker_image`); Terra uses the Broad
GPU image `gcr.io/broad-cga-francois-gtex/tensorqtl:latest`.

Author: Haim Krupkin+Claude. Contact: hkrupkin@stanford.edu
