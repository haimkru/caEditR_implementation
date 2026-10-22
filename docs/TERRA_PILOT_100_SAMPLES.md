# Terra pilot: 100 AoU samples, step by step

This walks through one real run of the whole pipeline on **100 AoU samples**,
using `editing_edqtl_pipeline` as a single Terra submission. It checks the
setup, the files and the cost before the full run. Expect about 1-2 hours and
well under $1 of compute. Power is low at 100 samples, so read the results as
a technical check, not as discoveries.

## 0. One-time setup (only the first time)

Follow `TERRA_RUN_GUIDE.md` section 0:

1. Push the repo.
2. Let `Docker Image CI` build the image. Make the ghcr.io package public and
   copy the image digest `ghcr.io/aou-multiomics-analysis/caeditr_edqtl@sha256:...`.
3. Register the repo on Dockstore.
4. Import `editing_edqtl_pipeline` into the AoU workspace.

Strongly recommended: first run the same workflow on the small simulated
fixture (`tests/fixtures/sim_small`, copied to the bucket). It takes minutes
and costs cents. Use the fixture row of `examples/terra_data_table.tsv` and
the settings in `TERRA_RUN_GUIDE.md` section 4b.

## 1. Inputs in the workspace bucket

You need these files in `gs://<WORKSPACE_BUCKET>/...`. Formats are in
`data-dictionary.md`. All of them can contain every sample; the sample list
picks the 100.

| Input | File |
| --- | --- |
| `editing_ratios` | editing ratio BED matrix (`#chr start end site_id <research_ids>`) |
| `editing_coverage` | matching read-depth BED matrix |
| `expression` | linear-CPM gene BED |
| `gtf` | the GTF used for expression |
| `proportions` | CellTypeDeconvolution `proportions_tca_weights.tsv` |
| `covariates` | `sample_id` + genotype PCs + sex + age |
| `plink_pgen/pvar/psam` | prepare_QTL PLINK2 genotypes |

## 2. Make the 100-sample list (Terra Jupyter / RStudio, R)

Pick 100 samples that are present in **every** input, so none are dropped:

```r
bucket <- Sys.getenv("WORKSPACE_BUCKET")
header <- function(path) strsplit(system(paste("gsutil cat", path, "| gzip -dcf | head -n 1"),
                                         intern = TRUE), "\t")[[1]]
editing <- header(paste0(bucket, "/editing/editing_ratios.bed.gz"))[-(1:4)]
coverage <- header(paste0(bucket, "/editing/editing_coverage.bed.gz"))[-(1:4)]
expr <- header(paste0(bucket, "/expression/expression.cpm.bed.gz"))[-(1:4)]
props <- read.delim(pipe(paste("gsutil cat", paste0(bucket, "/deconvolution/proportions_tca_weights.tsv"))),
                    colClasses = "character")$sample_id
covs <- read.delim(pipe(paste("gsutil cat", paste0(bucket, "/covariates/additional_covariates.tsv"))),
                   colClasses = "character")
covs <- covs$sample_id[stats::complete.cases(covs)]
psam <- read.delim(pipe(paste("gsutil cat", paste0(bucket, "/genotypes/aou.psam"))),
                   colClasses = "character", check.names = FALSE)[[1]]
shared <- Reduce(intersect, list(editing, coverage, expr, props, covs, psam))
cat(length(shared), "samples are in every input\n")
set.seed(1)
writeLines(sample(shared, 100), "sample_list_pilot100.txt")
system(paste("gsutil cp sample_list_pilot100.txt", paste0(bucket, "/sample_lists/")))
```

Adjust the paths to your files. If `length(shared)` is far below the size you
expect, the sample IDs do not match between inputs; fix that before going on.

## 3. Data table (TSV upload)

Edit `examples/terra_data_table_pilot100.tsv` (replace `YOUR_BUCKET` and the
file names). Then upload it: **Data -> Import data -> Upload TSV**. This
creates table `edqtl_run` with row `aou_pilot100`.

## 4. Configure the workflow

Go to **Workflows -> editing_edqtl_pipeline** and select "Run workflow(s)
with inputs defined by data table". Set the root entity type to `edqtl_run`
and select row `aou_pilot100`.

Inputs. Leave everything not listed at its default.

| Input | Value |
| --- | --- |
| `editing_ratios` | `this.editing_ratios` |
| `editing_coverage` | `this.editing_coverage` |
| `expression` | `this.expression` |
| `gtf` | `this.gtf` |
| `proportions` | `this.proportions` |
| `covariates` | `this.covariates` |
| `plink_pgen` / `plink_pvar` / `plink_psam` | `this.plink_pgen` / `this.plink_pvar` / `this.plink_psam` |
| `sample_list` | `this.sample_list` |
| `output_prefix` | `this.output_prefix` |
| `docker_image` | `"ghcr.io/aou-multiomics-analysis/caeditr_edqtl@sha256:<digest>"` |
| `sites_per_chunk` | `5000` (at 100 samples each shard takes seconds; fewer VMs = less start-up cost) |
| `n_editing_pcs` | `5` (15 is for the full cohort) |
| `run_tensorqtl` | `true` |
| `run_tensorqtl_nominal` | `true` (default) |
| `run_tensorqtl_permutations` | `false` (default) |
| `run_canrd_editqtl` | `true` |

The rest keeps the AoU defaults:

- `cis_window = 1000000`, shared by tensorQTL and caNRD-editQTL;
- `min_coverage = 20` and `min_covered_samples = 60`;
- `preemptible_attempts = 2`;
- the tensorQTL GPU machine `g2-standard-16` + `nvidia-l4` in `us-central1-c`.

Outputs: click **Use defaults**, so the results are written back to the row.
Keep **call caching on**, then press **Run analysis**.

## 5. While it runs

Open the submission in Job Manager. The order is:

1. IndexEditingRatios / IndexEditingCoverage;
2. PrepareCohort;
3. DeconvolveChunk (about 20 shards);
4. MergeChunkOutputs;
5. PreparePhenotypes (bulk + 2 x cell types);
6. tensorQTL (nominal) and caNRD-editQTL;
7. summary.

Every task log has `stage=... start_time=... completion_time=...` lines.
Errors from the pipeline itself start with `Input error:`.

## 6. Check the outputs (in this order)

1. **`sample_alignment_report.tsv`**: 100 rows with `kept = TRUE`. Any
   `FALSE` column shows which input lacks a sample.
2. **`cohort_summary.tsv`**: shows samples, cell types, sites, sites with a
   host gene, and chunks. Most sites should have a host gene; if not, check
   that the chromosome names and GTF match.
3. **`site_qc`**: shows how many sites pass coverage QC (`pass = TRUE`).
4. **`phenotype_inventory.tsv`** and **`deconvolution_summary.tsv`**: one
   bulk set plus caNRD-edit and TCA sets per cell type. The caNRD-edit
   missing fraction is high for rare cell types (expected: RNA share < 10%).
5. **`tensorqtl_summary`** and **`canrd_editqtl_summary`**: sites tested and
   counts per cell type.
6. **`edsite_overview`**: the three approaches side by side.
7. **Cost**: the submission's cost in Terra. Use it with the per-stage table
   in `TERRA_RUN_GUIDE.md` section 5 to extrapolate to the full run.

## 7. Common problems

| Message | Fix |
| --- | --- |
| `only N samples are shared by all inputs (need >= 30)` | IDs differ between inputs (e.g. `research_id` vs another ID); see `sample_alignment_report.tsv` |
| `editing ratios and coverage must list the same sites` | the two BEDs have different site rows; make them from the same site list |
| `... must start with columns #chr, start, end, site_id` | wrong header; see `data-dictionary.md` |
| `proportion rows must sum to 1` | use `proportions_tca_weights.tsv`, not the 22-column LM22 table |
| tensorQTL task waits for a long time | GPU quota / availability; see `TERRA_RUN_GUIDE.md` section 0 for GPU options |
| a task is killed for memory | raise its memory input (e.g. `canrd_chunk_memory`) or lower `sites_per_chunk` |

## 8. Then the full run

Use the same workflow and data-table layout without the 100-sample list (or
with the full cohort list). Change `sites_per_chunk` back to `500` and
`n_editing_pcs` to `15`. Settings are in `TERRA_RUN_GUIDE.md` section 4b;
time and cost are in section 5.
