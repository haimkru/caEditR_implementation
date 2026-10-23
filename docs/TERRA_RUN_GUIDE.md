# Running caEditR_implementation on All of Us (Terra)

First real run: follow [TERRA_PILOT_100_SAMPLES.md](TERRA_PILOT_100_SAMPLES.md)
(100 AoU samples, step by step). This guide is the full reference.

Who this is for: whoever runs the pipeline in the AoU Researcher Workbench /
Terra. The same WDLs and Docker images were run locally on simulated data
with miniwdl + Docker (`code/pipeline_steps/step1-5` in the development
project); on Terra only the input paths and the resource settings change.

## 0. One-time setup

1. **Push this repository** to GitHub as `haimkru/caEditR_implementation`.
2. **Docker image.** The `Docker Image CI` GitHub Action builds `envs/Dockerfile`
   (linux/amd64) and pushes `ghcr.io/haimkru/caeditr_implementation`
   with tags `sha-<commit>` and `latest`. Make the package public (GitHub ->
   Packages -> caeditr_implementation -> settings), like the other AoU-Multiomics images.
   Copy the digest printed at the end of the job
   (`ghcr.io/haimkru/caeditr_implementation@sha256:...`) and use it as
   `docker_image` in every run, so a run is reproducible.
   tensorQTL tasks use the AoU default `gcr.io/broad-cga-francois-gtex/tensorqtl:latest`
   (GPU). `envs/tensorqtl_cpu/Dockerfile` is a CPU twin for local tests and CI
   (`tensorqtl_docker_image` input); it is not needed on Terra.
3. **Dockstore.** Log in to dockstore.org with GitHub, install the Dockstore
   GitHub App on the repo (https://github.com/apps/dockstore), then push any
   commit: Dockstore registers the workflows from `.dockstore.yml` on push.
   Publish them under My Workflows. `.dockstore.yml` defines:

   | Dockstore name | WDL | What it does |
   | --- | --- | --- |
   | `editing_edqtl_pipeline` | `workflows/editing_edqtl_pipeline.wdl` | Everything in one submission |
   | `prepare_editing_phenotypes` | `workflows/prepare_editing_phenotypes.wdl` | Step 2: caNRD-edit + TCA cell-type matrices, tensorQTL phenotype sets |
   | `tensorqtl_edqtl` | `workflows/tensorqtl_edqtl.wdl` | Step 3: tensorQTL on bulk, caNRD-edit and TCA sets |
   | `canrd_editqtl` | `workflows/canrd_editqtl.wdl` | Step 4: end-to-end caNRD-editQTL |
   | `summarize_edqtl` | `workflows/summarize_edqtl.wdl` | Step 5a: side-by-side comparison |

   Import each into the AoU workspace ("Launch with Terra" from Dockstore).

## 1. Inputs to place in the workspace bucket

See `docs/data-dictionary.md` for exact formats.

| WDL input | Where it comes from on AoU |
| --- | --- |
| `editing_ratios` | Editing-level matrix (sites x samples, BED layout, `site_id` = `chrom:pos:strand`) |
| `editing_coverage` | Matching read-depth matrix (same sites and samples) |
| `expression` | Linear-CPM gene BED (the `expression` input of CellTypeDeconvolution / prepare_QTL) |
| `gtf` | The GTF used for expression (GENCODE) |
| `proportions` | `tca_weights` output of CellTypeDeconvolution (`proportions_tca_weights.tsv`) |
| `covariates` | `sample_id` + genotype PCs (calculateGenotypePCs), sex, age — the prepare_QTL `AdditionalCovariates` file |
| `plink_pgen/pvar/psam` | prepare_QTL `prepare_VCF` / `convertVCF2Plink` outputs |
| `sample_list` | Optional (e.g. the eQTL sample list) |

If the editing data exist as caEditR CSV matrices instead of BEDs, convert
them first, in a terminal of the image (`docker run ghcr.io/haimkru/caeditr_implementation@sha256:...`):
`Rscript /opt/caeditr_implementation/scripts/convert_caeditr_csv_to_bed.R --input-csv ratios.csv --output-bed editing_ratios.bed.gz`.

## 2. Data table

Terra data tables are uploaded as TSV. Edit `examples/terra_data_table.tsv`
(one row per run, `gs://` paths) and upload it: Data -> Import data ->
Upload TSV. The first column name `entity:edqtl_run_id` creates a table called
`edqtl_run`.

## 3. Option A - one submission

Workflow `editing_edqtl_pipeline`, "Run workflow(s) with inputs defined by
data table", root entity `edqtl_run`. Set the inputs to `this.<column>`
(`this.editing_ratios`, ... , `this.output_prefix`), `docker_image` to the
digest, and "Use defaults" for the rest. Outputs: "Use defaults" writes them
back to the row.

## 4. Option B - step by step (same order as the local step files)

Write every workflow's outputs to the data table ("Outputs" tab -> use
defaults), so the next workflow reads them as `this.<output>`:

1. `prepare_editing_phenotypes`: inputs `this.editing_ratios`, `this.editing_coverage`,
   `this.expression`, `this.gtf`, `this.proportions`, `this.covariates`,
   `this.plink_psam`, `this.sample_list`, `this.output_prefix`.
2. `tensorqtl_edqtl`: `phenotype_labels = this.phenotype_labels`,
   `phenotype_n_sites = this.phenotype_n_sites`, `phenotype_int_beds = this.phenotype_int_beds`,
   `phenotype_covariates = this.phenotype_covariates`, PLINK files, `output_prefix`.
   Default: nominal mode only (all pairs, raw p-values); `run_permutations = true`
   adds tensorQTL's empirical per-site p-values at extra GPU cost.
3. `canrd_editqtl`: `chunk_ids = this.chunk_ids`, `chunks = this.chunks`,
   `site_chunks = this.site_chunks`, `editing_ratios = this.indexed_editing_ratios`,
   `editing_ratios_index = this.indexed_editing_ratios_index`,
   `editing_coverage = this.indexed_editing_coverage`,
   `editing_coverage_index = this.indexed_editing_coverage_index`,
   `site_theta = this.site_theta`, `proportions = this.proportions_aligned`,
   `cohort_samples = this.cohort_samples`, `covariates = this.covariates_aligned`,
   `editing_pcs = this.bulk_editing_pcs`, PLINK files, `output_prefix`.
4. `summarize_edqtl`: `tensorqtl_all_analyses = this.tensorqtl_all_analyses`,
   `canrd_lead_site_celltype = this.lead_site_celltype`, `celltypes = this.celltypes`.

Mapping from the local parameter CSV (`code/pipeline_steps/parameters_*.csv`):
rows with `step = 2/3/4` are the inputs of workflows 1/2/3 above, with the
same names; for AoU use the defaults written in the description column
(e.g. `sites_per_chunk = 500`, `n_editing_pcs = 15`, `min_covered_samples = 60`).

## 4b. Run order and settings: fixture -> 100-sample pilot -> full run

| Setting | Fixture (`tests/fixtures/sim_small`) | Pilot (100 AoU samples) | Full AoU run |
| --- | --- | --- | --- |
| `docker_image` | ghcr.io `@sha256:` digest | same | same |
| `sample_list` | fixture `sample_list.txt` | 100 research_ids | full cohort list, or unset = every sample present in all inputs |
| `sites_per_chunk` | 10 | 5000 (seconds per shard; fewer VMs = less start-up cost) | **500** (default, ~200 shards) |
| `n_editing_pcs` | 3 | 5 (15 PCs + covariates is too many for 100 samples) | **15** (default) |
| `min_covered_samples` | 30 | 60 (default) | 60 (default); consider also `min_covered_fraction` (e.g. 0.1) because 60 is only 0.6 % of 10k samples |
| `cis_window` | 1000000 | 1000000 | 1000000 (same for tensorQTL and caNRD-editQTL) |
| tensorQTL | nominal (default) | nominal | nominal (GPU g2-standard-16 + L4) |
| `preemptible_attempts` | 2 | 2 | 2 |

Expected: fixture minutes and cents; pilot well under $1 (the two index tasks
and input downloads stay full size); full run ~4-6 h wall time if the L4 GPU
quota allows the 13 tensorQTL jobs in parallel (~8-10 h with 1 GPU), ~150
vCPU-hours + 4-9 GPU-hours, ~$4-6 preemptible (~$13-19 on demand).
Check `sample_alignment_report.tsv` after PrepareCohort: the cohort is the
intersection of all inputs (and the sample list); it must have >= 30 samples.

## 5. Scale and resources (10k samples, 100k sites, 1M variants)

Measured (not guessed) with `code/analysis/profile_pipeline_resources.sh` on a
simulated cohort with **10,000 samples** and AoU-like cis density (710 cis
variants per site in the +/- 1 Mb window), native amd64-class CPU, then
projected linearly in the number of sites (every stage is linear in sites by
design):

| Task | Shards | Projected peak memory | Time per shard | CPU-hours | Default VM |
| --- | --- | --- | --- | --- | --- |
| IndexBedMatrix | 2 | 0.1 GB (+ sort buffer) | ~9 min | 0.3 | 4 CPU, 16 GB, 200 GB disk |
| PrepareCohort | 1 | <= 24 GB (upper bound) | ~9 min | 0.2 | 4 CPU, 32 GB |
| DeconvolveChunk | 200 | 4.2 GB | ~1.2 min | 4 | 1 CPU, 6 GB, 100 GB disk |
| MergeChunkOutputs | 1 | streams text | minutes | <1 | 4 CPU, 16 GB, 250 GB disk |
| PreparePhenotypes | 1 + 2 x cell types | 2.1 GB (streams site blocks; Gram matrix 10k x 10k) | ~60 min | ~1 per set | 2 CPU, 8 GB |
| CanrdEditqtlChunk | 200 | 12.4 GB | ~3.8 min | 13 | 4 CPU, 24 GB, 150 GB disk |
| tensorQTL cis nominal | 1 + 2 x cell types | AoU tensorQTL defaults | | GPU | g2-standard-16 + L4 |

Total ~20 CPU-hours + the GPU tensorQTL runs; wall time is set by the scatter
width (200 shards in parallel -> roughly 1-2 h end to end, quota permitting).
If a task is killed for memory, raise its `*_memory` input or lower
`sites_per_chunk` (memory and time per shard are proportional to it).
`preemptible_attempts = 2` and `max_retries = 1` are the defaults; Cromwell
call caching lets a failed run resume.

### Compute cost (estimate; check current GCP prices in your region)

Assumptions: GCP list prices of roughly $0.032 per vCPU-hour + $0.0042 per
GB-hour on demand (preemptible/spot ~60-70 % less), ~3 min of VM start-up and
file localization per shard, g2-standard-16 + L4 at ~$1.2/h on demand.

| Part | VM-hours | On demand | Preemptible |
| --- | --- | --- | --- |
| DeconvolveChunk (200 x 1 CPU / 6 GB) | ~14 | ~$1 | ~$0.3 |
| CanrdEditqtlChunk (200 x 4 CPU / 24 GB) | ~23 | ~$5 | ~$1.5 |
| PreparePhenotypes (13 x 2 CPU / 8 GB) | ~14 | ~$1.5 | ~$0.5 |
| Index / cohort / merge / aggregate | ~2 | <$1 | <$0.5 |
| tensorQTL cis nominal (13 sets, GPU; one pass, no permutations) | ~2-4 | ~$2-5 | ~$1-2 |
| **Total compute** | | **~$10-15** | **~$3-5** |

Storage of the outputs (~100-150 GB incl. caNRD-editQTL all pairs) is ~$3-5
per month. How to keep it low:

- Keep `preemptible_attempts = 2` (default) and `tensorqtl_num_preempt = 1`.
- Keep `run_tensorqtl_permutations = false` (default): permutations (10,000
  per site) multiply the GPU time and only add tensorQTL's own empirical
  p-values, which are not used.
- Restrict `deconvolution_methods` (e.g. only `canrd_edit`) to halve the
  tensorQTL runs if TCA is not needed.
- Set `canrd_write_all_pairs = false` if only the lead / site tables are needed
  (smaller outputs, less storage).
- Keep Terra call caching on: a resubmission reuses every finished task.
- Run the small fixture first (cents) to catch setup errors before paying for
  the full run.

## 6. Figures

In a Terra RStudio / Jupyter (R) session, copy the finalized tables into
`results/finalized_tables/<run>/`, e.g. in R:
`system("gsutil -u $GOOGLE_PROJECT cp gs://<bucket>/submissions/<...>/*.tsv* results/finalized_tables/<run>/")`,
then run, in R,
`system('Rscript code/analysis/figures/make_figures.R --run-name <run> --prefix <prefix> --truth-dir ""')`
(fig1 = edSite counts; figures 2-4 need simulation truth). Needs `ggplot2`,
`svglite`, `data.table`, `optparse`.

## 7. Checklist before submitting

1. `editing_ratios` and `editing_coverage` list the same sites and samples.
2. Chromosome names agree between the BEDs, the GTF and the `.pvar` (`chr1`).
3. `gene_id`s in the GTF match the expression BED (including version suffixes).
4. The proportions are the CellTypeDeconvolution `tca_weights` for these samples.
5. Covariates have no constant columns and no names starting with `editing_PC`.
6. `docker_image` is a digest; record the repository commit.

## 8. What has and has not been tested

- Passed: `miniwdl check` and Cromwell `womtool validate` (v87) of every WDL;
  logging / workflow-scope file-write checks; R unit tests.
- Passed: complete local runs of `prepare_editing_phenotypes`, `tensorqtl_edqtl`,
  `canrd_editqtl`, `summarize_edqtl` with miniwdl in Docker, on simulated
  data (100 x 100 x 100 and 1,000 samples x 2,000 sites x 20,000 variants),
  with both images built natively for the Apple-Silicon test machine from the
  same Dockerfiles (x86 emulation crashes R arrow and PyTorch); tensorQTL with
  the CPU image `envs/tensorqtl_cpu` (same tensorQTL 1.0.9 as the Broad image).
  The linux/amd64 images Terra uses are built and run end to end by the
  GitHub Actions `Pipeline CI` (pending the first push).
- Not yet done: a Terra run (Cromwell on Google Batch, GPU tensorQTL, cloud
  file localization). Run the small simulated fixture
  (`tests/fixtures/sim_small`) on Terra once before the full AoU run.
