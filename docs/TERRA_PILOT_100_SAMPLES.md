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
   copy the image digest `ghcr.io/haimkru/caeditr_implementation@sha256:...`.
3. Register the repo on Dockstore.
4. Import `editing_edqtl_pipeline` into the AoU workspace.

Strongly recommended: first run the same workflow on the small simulated
fixture. It takes minutes and costs cents. Copy the fixture to the bucket in R:

```r
bucket <- Sys.getenv("WORKSPACE_BUCKET")
system("git clone https://github.com/haimkru/caEditR_implementation.git")
system(paste0("gsutil -m cp caEditR_implementation/tests/fixtures/sim_small/* ", bucket, "/caeditr_test/"))
```

Then use the `sim_small_fixture` row of `examples/terra_data_table.tsv` (with
`gs://YOUR_BUCKET/caeditr_edqtl_test/` replaced by your
`$WORKSPACE_BUCKET/caeditr_test/`) and the fixture settings in
`TERRA_RUN_GUIDE.md` section 4b.

## 1. The 7 inputs and where they come from on AoU

| # | Input (WDL name) | Where it comes from on AoU | Typical file |
| --- | --- | --- | --- |
| 1 | genotypes (`plink_pgen`, `plink_pvar`, `plink_psam`) | prepare_QTL `prepare_VCF` (the genotypes of the 9k eQTL tensorQTL runs) | `<prefix>.pgen/.pvar/.psam` under `submissions/.../call-plink2/` |
| 2 | RNA editing ratio per site (`editing_ratios`) | your editing quantification (not produced by any AoU-Multiomics repo) | BED matrix, see `data-dictionary.md` |
| 3 | coverage per site (`editing_coverage`) | same quantification as #2 | BED matrix, same sites and samples |
| 4 | expression (`expression`) | the linear-CPM BED given to the cell-type deconvolution (its `expression` input, or its `filtered_expression` output); prepare_eQTL `.expression.raw.bed.gz` from a count GCT is TMM CPM too | `*.bed.gz` |
| 5 | GENCODE annotation (`gtf`) | the notebook disk: `/home/jupyter/multiomics_qtls/edit/gencode.v48.annotation.gtf.gz`. Copy it to the bucket, because workflows only read `gs://` files | `gencode.v48.annotation.gtf.gz` |
| 6 | covariates (`covariates`): genotype PCs, sex, age | genotype PCs: `prepare_VCF` / `calculateGenotypePCs` (= the eQTL `AdditionalCovariates`); add sex and age from the AoU demographics | `<prefix>_genetic_PCs.tsv` (+ sex, age columns) |
| 7 | cell-type proportions (`proportions`) | the `tca_weights` output of the deconvolution / prepare_cell_type_eQTL run | `tca_weights.tsv` or `proportions_tca_weights.tsv` |

The **sample list** of the AoU eQTL cohort is
`gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/eQTL_9k_set/sample_lists/COMB_eQTL_9k_sample_list.tsv`
(a `research_id` column). The ID column may be called `research_id` or
`sample_id`; the pipeline accepts both.

Tip: the exact genotype and covariate files of the 9k eQTL runs are listed
in Terra **Job History -> that tensorQTL submission -> Inputs**
(`plink_pgen`, `covariates`, ...).

Find candidates for each input in your bucket (R; this includes workflow outputs under `submissions/`):

```r
project <- Sys.getenv("GOOGLE_PROJECT")
bucket  <- Sys.getenv("WORKSPACE_BUCKET")
all_files <- system(paste("gsutil -u", project, "ls -r", bucket), intern = TRUE)
find <- function(pattern) all_files[grepl(pattern, all_files, ignore.case = TRUE)]
cat("\n== 1 genotypes ==\n");     print(find("\\.(pgen|pvar|psam)$"))
cat("\n== 2-3 editing ==\n");     print(find("edit|redi|coverage|a2g|a_to_g"))
cat("\n== 4 expression ==\n");    print(find("expression.*bed\\.gz$|cpm.*bed\\.gz$|gene_reads\\.gct"))
cat("\n== 5 gtf ==\n");           print(find("\\.gtf"))
cat("\n== 6 covariates ==\n");    print(find("genetic_PCs|covariate"))
cat("\n== 7 proportions ==\n");   print(find("tca_weights"))
cat("\n== sample lists ==\n");    print(find("sample_list"))
# copy the GENCODE GTF from the notebook disk to the bucket (workflows need gs:// files)
system(paste0("gsutil -u ", project, " cp /home/jupyter/multiomics_qtls/edit/gencode.v48.annotation.gtf.gz ",
              bucket, "/reference/"))
```

Use the paths it prints in step 2 below and in the data table (step 3).

## 2. Make the 100-sample list (all in R: Terra RStudio or a Jupyter R notebook)

Everything in this step is R; paste it into the R console or an R notebook
cell (not a terminal).

**2a. Find the real paths of your files:**

```r
project <- Sys.getenv("GOOGLE_PROJECT")    # billing project (needed for requester-pays buckets)
bucket  <- Sys.getenv("WORKSPACE_BUCKET")  # or another gs:// bucket that holds the inputs
files <- system(paste("gsutil -u", project, "ls -r", bucket), intern = TRUE)
files <- files[grepl("bed\\.gz$|\\.tsv$|\\.psam$|\\.gtf", files) & !grepl("/submissions/", files)]
print(files)
```

**2b. Check the files and pick the 100 samples.** Put the paths printed above into `paths`. The code checks every file, then
picks 100 samples present in **every** input, so none are dropped:

```r
project <- Sys.getenv("GOOGLE_PROJECT")    # billing project (needed for requester-pays buckets)
bucket  <- Sys.getenv("WORKSPACE_BUCKET")
paths <- c(                                 # <- edit to your files
  editing_ratios   = paste0(bucket, "/editing/editing_ratios.bed.gz"),
  editing_coverage = paste0(bucket, "/editing/editing_coverage.bed.gz"),
  expression       = paste0(bucket, "/expression/expression.cpm.bed.gz"),
  proportions      = paste0(bucket, "/deconvolution/proportions_tca_weights.tsv"),
  covariates       = paste0(bucket, "/covariates/additional_covariates.tsv"),
  psam             = paste0(bucket, "/genotypes/aou.psam")
)

# stop with a clear message if a file is missing or unreadable
gs_check <- function(path) {
  status <- system(paste("gsutil -u", project, "-q stat", shQuote(path)))
  if (status != 0) stop("Not found / not readable: ", path,
                        "\n  list the folder with: gsutil -u ", project, " ls ", dirname(path))
}
for (path in paths) gs_check(path)

gs_stream <- function(path) paste("gsutil -u", project, "cat", shQuote(path), "| gzip -dcf")
header <- function(path) {
  line <- system(paste(gs_stream(path), "| head -n 1"), intern = TRUE)
  if (length(line) == 0) stop("Empty file: ", path)
  strsplit(line, "\t")[[1]]
}
read_tsv_gs <- function(path) read.delim(pipe(gs_stream(path)), colClasses = "character", check.names = FALSE)

editing  <- header(paths[["editing_ratios"]])[-(1:4)]
coverage <- header(paths[["editing_coverage"]])[-(1:4)]
expr     <- header(paths[["expression"]])[-(1:4)]
props    <- read_tsv_gs(paths[["proportions"]])$sample_id
covs     <- read_tsv_gs(paths[["covariates"]])
covs     <- covs$sample_id[stats::complete.cases(covs)]
psam     <- read_tsv_gs(paths[["psam"]])[[1]]          # first column = IID

counts <- c(editing = length(editing), coverage = length(coverage), expression = length(expr),
            proportions = length(props), covariates = length(covs), genotypes = length(psam))
print(counts)
shared <- Reduce(intersect, list(editing, coverage, expr, props, covs, psam))
cat(length(shared), "samples are in every input\n")

set.seed(1)
writeLines(sample(shared, 100), "sample_list_pilot100.txt")
system(paste("gsutil -u", project, "cp sample_list_pilot100.txt", paste0(bucket, "/sample_lists/")))
```

If `length(shared)` is far below the counts printed above, the sample IDs do
not match between inputs (e.g. one file uses `research_id`, another a
different ID); fix that before going on.

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
| `docker_image` | `"ghcr.io/haimkru/caeditr_implementation@sha256:d7758d4b4fbee8f8d18ac8b70ead6e5d772b7c43f7056973d1d9562c9d9496af"` (current image, built from commit 642f665; use the newest digest from the `Docker Image CI` log after any change to `R/`, `scripts/` or `envs/`) |
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
