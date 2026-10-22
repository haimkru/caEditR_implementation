#!/usr/bin/env Rscript
################################################################################
# scripts/deconvolve_chunk.R
#
# Description:
#   Scatter stage of the caEditR_edQTL pipeline (WDL task DeconvolveChunk).
#   For one genomic chunk of editing sites:
#     * per-site coverage QC (same rule in every stage);
#     * writes the QC-passing BULK ratios (input of the standard edQTL scan);
#     * caNRD-edit (caEditR::caNRD_edit, maximum-likelihood estimator): one
#       sites x samples matrix per cell type;
#     * TCA (caEditR::TCA_Like = CRAN TCA): one matrix per cell type.
#   Chunks are merged genome-wide afterwards (merge_chunk_outputs.sh).
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   1. Load cohort samples, proportions and the chunk's ratios / coverage /
#      theta (R/chunk_data.R, tabix region read).
#   2. Site QC: >= --min-covered-samples samples with >= --min-coverage
#      reads (and >= --min-covered-fraction), non-zero sd; write the table.
#   3. Bulk: write QC-passing ratios (NA where a sample has no reads).
#   4. caNRD-edit on QC-passing sites with coverage and theta
#      (theta_floor = --theta-floor); write per-cell-type matrices and the
#      per-site diagnostics (status, condition number, excluded cell types).
#   5. TCA on QC-passing sites; TCA needs complete data, so missing ratios
#      are replaced by the site mean (recorded in the summary); write
#      per-cell-type matrices.
#   6. Write a small per-chunk summary (method, cell type, sites, missing).
#
# Inputs:
#   --chunk-id, --chunks, --site-chunks        chunk definitions (prepare_cohort.R)
#   --editing-ratios, --editing-coverage       bgzipped + tabix BEDs
#   --site-theta, --proportions, --samples     prepare_cohort.R outputs
#   --min-coverage (20) --min-covered-samples (60) --min-covered-fraction (0)
#   --theta-floor (0.001) --canrd-min-coverage (10)
#   --methods "canrd_edit,tca"                 which deconvolutions to run
#   --output-dir
# Outputs (in --output-dir):
#   beds/bulk.all.<chunk>.bed.gz
#   beds/canrd_edit.<celltype_slug>.<chunk>.bed.gz
#   beds/tca.<celltype_slug>.<chunk>.bed.gz
#   reports/<chunk>.site_qc.tsv
#   reports/<chunk>.canrd_diagnostics.tsv
#   reports/<chunk>.deconvolution_summary.tsv
################################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(caEditR)
})

## ---- load the shared helpers in <repo>/R -----------------------------------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
for (helper in list.files(file.path(script_dir, "..", "R"), pattern = "[.]R$", full.names = TRUE)) {
  source(helper)
}

## ---- 0. options ---------------------------------------------------------------
option_list <- list(
  make_option("--chunk-id", type = "character"),
  make_option("--chunks", type = "character"),
  make_option("--site-chunks", type = "character"),
  make_option("--editing-ratios", type = "character"),
  make_option("--editing-coverage", type = "character"),
  make_option("--site-theta", type = "character"),
  make_option("--proportions", type = "character"),
  make_option("--samples", type = "character"),
  make_option("--min-coverage", type = "double", default = 20),
  make_option("--min-covered-samples", type = "integer", default = 60L),
  make_option("--min-covered-fraction", type = "double", default = 0),
  make_option("--theta-floor", type = "double", default = 1e-3),
  make_option("--canrd-min-coverage", type = "double", default = 10),
  make_option("--methods", type = "character", default = "canrd_edit,tca"),
  make_option("--output-dir", type = "character", default = "outputs")
)
opt <- clean_option_names(parse_args(OptionParser(option_list = option_list)))
require_options(opt, c("chunk_id", "chunks", "site_chunks", "editing_ratios", "editing_coverage",
                       "site_theta", "proportions", "samples"))
methods <- split_comma_option(opt$methods)
unknown_methods <- setdiff(methods, c("canrd_edit", "tca"))
if (length(unknown_methods) > 0) {
  stop_input("--methods may only contain canrd_edit and tca; got ", paste(unknown_methods, collapse = ", "))
}
bed_dir <- file.path(opt$output_dir, "beds")
report_dir <- file.path(opt$output_dir, "reports")
dir.create(bed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
chunk_id <- opt$chunk_id
log_info("deconvolve_chunk.R (caEditR_edQTL ", pipeline_version(), ") chunk ", chunk_id)

## ---- 1. inputs ------------------------------------------------------------------
cohort <- load_cohort_tables(opt$samples, opt$proportions)
celltypes <- cohort$celltypes
slugs <- slugify_celltype(celltypes)
chunk <- load_chunk_editing(chunk_id, opt$chunks, opt$site_chunks, opt$editing_ratios,
                            opt$editing_coverage, opt$site_theta, cohort$samples, celltypes)

## ---- 2. site QC ----------------------------------------------------------------
qc <- site_coverage_qc(chunk$ratio, chunk$coverage,
                       min_coverage = opt$min_coverage,
                       min_covered_samples = opt$min_covered_samples,
                       min_covered_fraction = opt$min_covered_fraction)
qc[, chunk_id := chunk_id]
write_tsv(qc, file.path(report_dir, paste0(chunk_id, ".site_qc.tsv")))
pass_sites <- qc[pass == TRUE, site_id]
log_info("Sites passing QC: ", length(pass_sites), " / ", nrow(qc))

keep <- match(pass_sites, rownames(chunk$ratio))
meta <- chunk$meta[keep]
ratio <- chunk$ratio[keep, , drop = FALSE]
coverage <- chunk$coverage[keep, , drop = FALSE]
theta <- chunk$theta[keep, , drop = FALSE]

summary_rows <- list()
add_summary <- function(method, celltype, values, note = "") {
  summary_rows[[length(summary_rows) + 1]] <<- data.table(
    chunk_id = chunk_id, method = method, celltype = celltype,
    n_sites = nrow(values), n_samples = ncol(values),
    fraction_missing = if (length(values)) mean(is.na(values)) else NA_real_,
    n_sites_all_missing = sum(rowSums(!is.na(values)) == 0),
    note = note
  )
}
chunk_bed <- function(method, slug) {
  file.path(bed_dir, paste0(method, ".", slug, ".", chunk_id, ".bed.gz"))
}

## ---- 3. bulk ratios ------------------------------------------------------------------
write_bed_matrix(meta, ratio, chunk_bed("bulk", "all"), index = FALSE)
add_summary("bulk", "all", ratio)

## ---- 4. caNRD-edit -----------------------------------------------------------------------
if ("canrd_edit" %in% methods) {
  if (length(pass_sites) > 0) {
    canrd <- suppressMessages(caNRD_edit(
      bulk_editing = ratio,
      coverage = coverage,
      proportions = cohort$proportions,
      theta = theta,
      theta_floor = opt$theta_floor,
      min_coverage = opt$canrd_min_coverage,
      estimator = "ml"
    ))
    diagnostics <- as.data.table(canrd$diagnostics)
    diagnostics[, chunk_id := chunk_id]
    for (h in seq_along(celltypes)) {
      values <- canrd$deconvolved[[celltypes[h]]][pass_sites, , drop = FALSE]
      write_bed_matrix(meta, values, chunk_bed("canrd_edit", slugs[h]), index = FALSE)
      add_summary("canrd_edit", celltypes[h], values)
    }
  } else {
    diagnostics <- data.table(site_id = character(0), chunk_id = character(0))
    for (h in seq_along(celltypes)) {
      write_bed_matrix(meta, ratio, chunk_bed("canrd_edit", slugs[h]), index = FALSE)
      add_summary("canrd_edit", celltypes[h], ratio, "no QC-passing sites")
    }
  }
  write_tsv(diagnostics, file.path(report_dir, paste0(chunk_id, ".canrd_diagnostics.tsv")))
}

## ---- 5. TCA --------------------------------------------------------------------------
if ("tca" %in% methods) {
  if (length(pass_sites) > 0) {
    # TCA needs complete data: impute missing ratios with the site mean
    n_imputed <- sum(is.na(ratio))
    ratio_complete <- impute_row_means(ratio)
    tca <- TCA_Like(ratio_complete, cohort$proportions)
    for (h in seq_along(celltypes)) {
      values <- tca$deconvolved[[celltypes[h]]][pass_sites, , drop = FALSE]
      write_bed_matrix(meta, values, chunk_bed("tca", slugs[h]), index = FALSE)
      add_summary("tca", celltypes[h], values, paste0("imputed_bulk_values=", n_imputed))
    }
  } else {
    for (h in seq_along(celltypes)) {
      write_bed_matrix(meta, ratio, chunk_bed("tca", slugs[h]), index = FALSE)
      add_summary("tca", celltypes[h], ratio, "no QC-passing sites")
    }
  }
}

## ---- 6. summary ------------------------------------------------------------------------
write_tsv(rbindlist(summary_rows), file.path(report_dir, paste0(chunk_id, ".deconvolution_summary.tsv")))
log_info("deconvolve_chunk.R finished chunk ", chunk_id)
