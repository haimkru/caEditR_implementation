#!/usr/bin/env Rscript
################################################################################
# scripts/prepare_cohort.R
#
# Description:
#   First stage of the caEditR_edQTL pipeline (WDL task PrepareCohort).
#   Defines the analysis cohort and every per-site / per-sample input that
#   later stages share, WITHOUT loading the large editing matrices:
#     * cohort samples = samples present in all inputs (+ optional list);
#     * site -> host gene mapping (GTF gene bodies);
#     * theta = per-cell-type host-gene expression, estimated by NNLS of the
#       bulk CPM on the cell-type proportions (caEditR::estimate_theta_nnls);
#     * genomic chunks of sites for the Terra scatter.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   1. Read sample IDs from: editing-ratio BED header, coverage BED header,
#      expression BED header, proportions (sample_id), PLINK2 psam (IID),
#      covariates (sample_id, complete numeric rows only), optional list.
#   2. Cohort = intersection, in editing-BED order; write the report.
#   3. Read the 4 coordinate columns of both editing BEDs; they must list the
#      same sites in the same order (both come from index_bed_matrix.sh).
#   4. Read GTF genes; find genes overlapping any site; read only those genes
#      from the expression BED (cohort samples only).
#   5. Map each site to one host gene (same strand preferred, then highest
#      mean CPM).
#   6. theta[gene, cell type] = NNLS(CPM ~ proportions), floored at
#      --theta-floor; site theta = host-gene theta; sites without a host gene
#      get theta = 1 in every cell type (RNA share = proportion).
#   7. Chunk the sites (<= --sites-per-chunk per chunk, one chromosome each).
#   8. Write everything below.
#
# Inputs:
#   --editing-ratios     sorted, bgzipped, tabix-indexed BED (site_id; ratios)
#   --editing-coverage   same layout, read depth (A+G reads)
#   --expression         BED(.gz) of linear CPM, 4th column gene_id
#   --gtf                GTF(.gz) with gene records (gene_id, gene_name)
#   --proportions        TSV sample_id + one column per cell group
#                        (CellTypeDeconvolution proportions_tca_weights.tsv)
#   --psam               PLINK2 .psam of the genotypes
#   --covariates         TSV sample_id + numeric covariates (genotype PCs, sex, age)
#   --sample-list        optional, one sample ID per line
#   --sites-per-chunk    sites per scatter shard (default 500)
#   --theta-floor        NNLS floor = "not expressed" (default 0.001)
#   --output-dir         output directory
# Outputs (in --output-dir):
#   cohort_samples.txt, sample_alignment_report.tsv, celltypes.tsv,
#   proportions_aligned.tsv, covariates_aligned.tsv, site_to_gene.tsv,
#   gene_theta.tsv, site_theta.tsv, chunks.tsv, site_chunks.tsv,
#   chunk_ids.txt, cohort_summary.tsv
################################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

## ---- load the shared helpers in <repo>/R -----------------------------------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
for (helper in list.files(file.path(script_dir, "..", "R"), pattern = "[.]R$", full.names = TRUE)) {
  source(helper)
}

## ---- 0. options ---------------------------------------------------------------
option_list <- list(
  make_option("--editing-ratios", type = "character"),
  make_option("--editing-coverage", type = "character"),
  make_option("--expression", type = "character"),
  make_option("--gtf", type = "character"),
  make_option("--proportions", type = "character"),
  make_option("--psam", type = "character"),
  make_option("--covariates", type = "character"),
  make_option("--sample-list", type = "character", default = ""),
  make_option("--sites-per-chunk", type = "integer", default = 500L),
  make_option("--theta-floor", type = "double", default = 1e-3),
  make_option("--output-dir", type = "character", default = "outputs")
)
opt <- clean_option_names(parse_args(OptionParser(option_list = option_list)))
require_options(opt, c("editing_ratios", "editing_coverage", "expression", "gtf",
                       "proportions", "psam", "covariates"))
for (name in c("editing_ratios", "editing_coverage", "expression", "gtf", "proportions",
               "psam", "covariates")) {
  require_readable_file(opt[[name]], name)
}
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(name) file.path(opt$output_dir, name)
log_info("prepare_cohort.R (caEditR_edQTL ", pipeline_version(), ") starting")

## ---- 1. sample IDs of every input -------------------------------------------------
editing_columns <- read_bed_header(opt$editing_ratios)
editing_samples <- check_bed_columns(editing_columns, "site_id", opt$editing_ratios)
coverage_columns <- read_bed_header(opt$editing_coverage)
coverage_samples <- check_bed_columns(coverage_columns, "site_id", opt$editing_coverage)
expression_columns <- read_bed_header(opt$expression)
expression_samples <- check_bed_columns(expression_columns, "gene_id", opt$expression)

proportions_dt <- read_sample_table(opt$proportions, "sample_id")
celltypes <- setdiff(names(proportions_dt), "sample_id")
if (length(celltypes) < 2) {
  stop_input("proportions must have at least two cell-type columns")
}
proportions <- sample_table_to_matrix(proportions_dt, "sample_id")
if (any(!is.finite(proportions)) || any(proportions < 0)) {
  stop_input("proportions must be finite and non-negative")
}
row_sum_error <- max(abs(rowSums(proportions) - 1))
if (row_sum_error > 1e-4) {
  stop_input("proportion rows must sum to 1 (max error ", signif(row_sum_error, 3), ")")
}

genotype_samples <- read_psam_sample_ids(opt$psam)

covariates_dt <- read_sample_table(opt$covariates, "sample_id")
covariates <- sample_table_to_matrix(covariates_dt, "sample_id")
if (ncol(covariates) == 0) {
  stop_input("covariates must have at least one numeric column besides sample_id")
}
complete_covariate_samples <- rownames(covariates)[stats::complete.cases(covariates)]

sample_list <- NULL
if (nzchar(opt$sample_list)) {
  require_readable_file(opt$sample_list, "sample_list")
  sample_list <- read_sample_list(opt$sample_list)
}

## ---- 2. cohort --------------------------------------------------------------------
alignment <- align_cohort_samples(
  list(editing_ratios = editing_samples,
       editing_coverage = coverage_samples,
       expression = expression_samples,
       proportions = rownames(proportions),
       genotypes = genotype_samples,
       covariates_complete = complete_covariate_samples),
  sample_list = sample_list
)
cohort <- alignment$samples
log_info("Cohort: ", length(cohort), " samples present in every input")
if (length(cohort) < 30) {
  stop_input("only ", length(cohort), " samples are shared by all inputs (need >= 30); ",
             "see sample_alignment_report.tsv")
}
writeLines(cohort, out("cohort_samples.txt"))
write_tsv(alignment$report, out("sample_alignment_report.tsv"))

celltype_table <- data.table(celltype = celltypes, slug = slugify_celltype(celltypes))
if (anyDuplicated(celltype_table$slug) > 0) {
  stop_input("two cell types have the same file-name slug: ",
             paste(celltype_table$celltype[duplicated(celltype_table$slug)], collapse = ", "))
}
write_tsv(celltype_table, out("celltypes.tsv"))

proportions_cohort <- proportions[cohort, , drop = FALSE]
write_tsv(data.table(sample_id = cohort, as.data.table(proportions_cohort)),
          out("proportions_aligned.tsv"))
covariates_cohort <- covariates[cohort, , drop = FALSE]
write_tsv(data.table(sample_id = cohort, as.data.table(covariates_cohort)),
          out("covariates_aligned.tsv"))

## ---- 3. sites (coordinates only) ----------------------------------------------------
read_site_coordinates <- function(path) {
  fread(cmd = stream_command(path), sep = "\t", header = TRUE, select = 1:4,
        colClasses = list(character = c(1, 4), integer = c(2, 3)), showProgress = FALSE)
}
site_meta <- read_site_coordinates(opt$editing_ratios)
coverage_meta <- read_site_coordinates(opt$editing_coverage)
if (!identical(site_meta$site_id, coverage_meta$site_id)) {
  stop_input("editing ratios and coverage must list the same sites in the same order ",
             "(index both with scripts/index_bed_matrix.sh)")
}
if (anyDuplicated(site_meta$site_id) > 0) {
  stop_input("duplicated site_id in editing ratios, e.g. ", site_meta$site_id[duplicated(site_meta$site_id)][1])
}
log_info("Sites: ", nrow(site_meta))

## ---- 4. candidate host genes and their expression ----------------------------------------
genes <- read_gtf_genes(opt$gtf)
log_info("GTF genes: ", nrow(genes))
site_points <- data.table(chr = site_meta[["#chr"]], start = site_meta$end, end = site_meta$end)
gene_ranges <- copy(genes)
setkey(gene_ranges, chr, start, end)
candidate_hits <- foverlaps(site_points, gene_ranges, type = "within", nomatch = NULL)
candidate_genes <- unique(candidate_hits$gene_id)
log_info("Genes overlapping at least one site: ", length(candidate_genes))

# read only candidate genes from the expression BED (awk filter on column 4)
id_file <- tempfile(fileext = ".txt")
writeLines(candidate_genes, id_file)
expression_command <- paste(
  stream_command(opt$expression),
  "| awk -F'\\t' 'NR == FNR { keep[$1]; next } FNR == 1 || ($4 in keep)'",
  shQuote(id_file), "-"
)
expression_dt <- fread(cmd = expression_command, sep = "\t", header = TRUE,
                       select = c("#chr", "start", "end", "gene_id", cohort),
                       colClasses = list(character = c("#chr", "gene_id")),
                       check.names = FALSE, showProgress = FALSE)
expression <- as.matrix(expression_dt[, cohort, with = FALSE])
storage.mode(expression) <- "double"
rownames(expression) <- expression_dt$gene_id
if (any(!is.finite(expression)) || any(expression < 0)) {
  stop_input("expression must be finite, non-negative linear CPM")
}
gene_mean_expression <- rowMeans(expression)
log_info("Expressed candidate host genes read: ", nrow(expression))

## ---- 5. site -> host gene -----------------------------------------------------------------
site_to_gene <- map_sites_to_host_genes(site_meta, genes, gene_mean_expression)
site_to_gene <- cbind(site_meta[, .(chr = `#chr`, bed_start = start, bed_end = end)], site_to_gene)
write_tsv(site_to_gene, out("site_to_gene.tsv"))
n_mapped <- sum(!is.na(site_to_gene$host_gene_id))
log_info("Sites mapped to a host gene: ", n_mapped, " / ", nrow(site_to_gene))

## ---- 6. theta -----------------------------------------------------------------------------
host_genes <- unique(stats::na.omit(site_to_gene$host_gene_id))
if (length(host_genes) > 0) {
  gene_theta <- caEditR::estimate_theta_nnls(
    expression[host_genes, , drop = FALSE],
    proportions_cohort,
    floor = opt$theta_floor
  )
} else {
  gene_theta <- matrix(numeric(0), 0, length(celltypes), dimnames = list(NULL, celltypes))
}
write_tsv(data.table(gene_id = rownames(gene_theta), as.data.table(gene_theta)), out("gene_theta.tsv"))

site_theta <- matrix(1, nrow(site_meta), length(celltypes),
                     dimnames = list(site_meta$site_id, celltypes))
mapped <- !is.na(site_to_gene$host_gene_id)
site_theta[mapped, ] <- gene_theta[site_to_gene$host_gene_id[mapped], , drop = FALSE]
write_tsv(data.table(site_id = site_meta$site_id, as.data.table(site_theta),
                     theta_source = ifelse(mapped, "host_gene_nnls", "no_host_gene_uniform")),
          out("site_theta.tsv"))

## ---- 7. chunks -----------------------------------------------------------------------------
chunking <- define_site_chunks(site_meta, opt$sites_per_chunk)
write_tsv(chunking$chunks, out("chunks.tsv"))
write_tsv(chunking$site_chunks, out("site_chunks.tsv"))
writeLines(chunking$chunks$chunk_id, out("chunk_ids.txt"))
log_info("Chunks: ", nrow(chunking$chunks), " (<= ", opt$sites_per_chunk, " sites each)")

## ---- 8. summary ------------------------------------------------------------------------------
summary_table <- data.table(
  field = c("pipeline_version", "n_samples_cohort", "n_celltypes", "celltypes", "n_sites",
          "n_sites_with_host_gene", "n_host_genes", "n_chunks", "sites_per_chunk",
          "theta_floor", "n_covariates"),
  value = c(pipeline_version(), length(cohort), length(celltypes), paste(celltypes, collapse = ";"),
            nrow(site_meta), n_mapped, length(host_genes), nrow(chunking$chunks),
            opt$sites_per_chunk, opt$theta_floor, ncol(covariates))
)
write_tsv(summary_table, out("cohort_summary.tsv"))
log_info("prepare_cohort.R finished")
