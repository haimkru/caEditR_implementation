#!/usr/bin/env Rscript
################################################################################
# scripts/canrd_editqtl_chunk.R
#
# Description:
#   Scatter stage of the END-TO-END cell-type edQTL scan (WDL task
#   CanrdEditqtlChunk): caEditR::caNRD_editQTL() tests, directly on the BULK
#   editing ratios, whether each cis variant changes editing in each cell
#   type (genotype x RNA-share model, caNRD variance model, read-depth
#   noise). One chunk of sites per task; all cis variants of those sites.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   1. Load cohort, proportions, the chunk's ratios / coverage / theta
#      (same loader and same site QC as deconvolve_chunk.R).
#   2. Covariates = supplied covariates (+ optionally the first
#      --n-editing-pcs bulk editing PCs); constant columns dropped.
#   3. Genotypes: plink2 extracts ALT dosages of variants with MAF >=
#      --maf in [first site - window, last site + window] for the cohort.
#   4. Pairs: every variant within +/- --cis-window of each QC-passing site.
#   5. caNRD_editQTL(engine = "scan"): per site the variance components are
#      fitted once without genotype, every variant is tested by GLS at that
#      variance, and the lead variant (smallest site-level p) is re-fitted
#      exactly with likelihood-ratio tests (refine = "lead").
#   6. Write: all pairs (one row per site x variant x cell type), lead
#      variants (exact re-fit rows), per-site summary (number of variants
#      tested -> used for the multiple-testing correction), per site x cell
#      type minimum p, scan diagnostics and site QC.
#
# Inputs:
#   --chunk-id, --chunks, --site-chunks, --editing-ratios, --editing-coverage,
#   --site-theta, --proportions, --samples      (as deconvolve_chunk.R)
#   --covariates          covariates_aligned.tsv (sample_id + numeric)
#   --editing-pcs         optional bulk <label>.phenotype_PCs.tsv
#   --n-editing-pcs       how many of those PCs to add (default 0)
#   --pfile-prefix        PLINK2 prefix (.pgen/.pvar/.psam)
#   --cis-window (1e6) --maf (0.01) --min-minor-allele-samples (10)
#   --min-samples (30) --min-mean-phi (0.1) --theta-floor (0.001)
#   --min-coverage (20) --min-covered-samples (60) --min-covered-fraction (0)
#   --refine ("lead" or a p-value threshold)
#   --write-all-pairs (true)
#   --output-dir
# Outputs (in --output-dir, prefix <chunk>.canrd_editqtl):
#   .all_pairs.tsv.gz, .lead_variants.tsv, .site_summary.tsv,
#   .site_celltype_min_p.tsv, .scan_info.tsv, .site_qc.tsv
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
  make_option("--covariates", type = "character"),
  make_option("--editing-pcs", type = "character", default = ""),
  make_option("--n-editing-pcs", type = "integer", default = 0L),
  make_option("--pfile-prefix", type = "character"),
  make_option("--cis-window", type = "double", default = 1e6),
  make_option("--maf", type = "double", default = 0.01),
  make_option("--min-minor-allele-samples", type = "integer", default = 10L),
  make_option("--min-samples", type = "integer", default = 30L),
  make_option("--min-mean-phi", type = "double", default = 0.10),
  make_option("--theta-floor", type = "double", default = 1e-3),
  make_option("--min-coverage", type = "double", default = 20),
  make_option("--min-covered-samples", type = "integer", default = 60L),
  make_option("--min-covered-fraction", type = "double", default = 0),
  make_option("--refine", type = "character", default = "lead"),
  make_option("--write-all-pairs", type = "character", default = "true"),
  make_option("--output-dir", type = "character", default = "outputs")
)
opt <- clean_option_names(parse_args(OptionParser(option_list = option_list)))
require_options(opt, c("chunk_id", "chunks", "site_chunks", "editing_ratios", "editing_coverage",
                       "site_theta", "proportions", "samples", "covariates", "pfile_prefix"))
for (extension in c(".pgen", ".pvar", ".psam")) {
  require_readable_file(paste0(opt$pfile_prefix, extension), paste0("genotypes", extension))
}
refine <- if (grepl("^[0-9.eE-]+$", opt$refine)) as.numeric(opt$refine) else opt$refine
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
chunk_id <- opt$chunk_id
out <- function(suffix) file.path(opt$output_dir, paste0(chunk_id, ".canrd_editqtl", suffix))
log_info("canrd_editqtl_chunk.R (caEditR_edQTL ", pipeline_version(), ") chunk ", chunk_id)

## ---- 1. editing data and site QC --------------------------------------------------
cohort <- load_cohort_tables(opt$samples, opt$proportions)
celltypes <- cohort$celltypes
chunk <- load_chunk_editing(chunk_id, opt$chunks, opt$site_chunks, opt$editing_ratios,
                            opt$editing_coverage, opt$site_theta, cohort$samples, celltypes)
qc <- site_coverage_qc(chunk$ratio, chunk$coverage,
                       min_coverage = opt$min_coverage,
                       min_covered_samples = opt$min_covered_samples,
                       min_covered_fraction = opt$min_covered_fraction)
qc[, chunk_id := chunk_id]
write_tsv(qc, out(".site_qc.tsv"))
pass_sites <- qc[pass == TRUE, site_id]
log_info("Sites passing QC: ", length(pass_sites), " / ", nrow(qc))
ratio <- chunk$ratio[pass_sites, , drop = FALSE]
coverage <- chunk$coverage[pass_sites, , drop = FALSE]
theta <- chunk$theta[pass_sites, , drop = FALSE]
sites <- parse_site_ids(pass_sites)
sites[, chr := chunk$meta[match(pass_sites, chunk$meta$site_id)][["#chr"]]]

## ---- 2. covariates ------------------------------------------------------------------
covariates <- sample_table_to_matrix(read_sample_table(opt$covariates, "sample_id"), "sample_id")
covariates <- covariates[cohort$samples, , drop = FALSE]
if (opt$n_editing_pcs > 0) {
  require_readable_file(opt$editing_pcs, "editing_pcs")
  pcs <- sample_table_to_matrix(read_sample_table(opt$editing_pcs, "sample_id"), "sample_id")
  n_use <- min(opt$n_editing_pcs, ncol(pcs))
  covariates <- cbind(covariates, pcs[cohort$samples, seq_len(n_use), drop = FALSE])
}
covariate_sd <- apply(covariates, 2, stats::sd)
covariates <- covariates[, is.finite(covariate_sd) & covariate_sd > 0, drop = FALSE]
log_info("Covariates used: ", paste(colnames(covariates), collapse = ", "))

## ---- 3-4. genotypes and cis pairs ----------------------------------------------------
pairs <- data.table(site_id = character(0), variant_id = character(0), distance = numeric(0))
dosages <- NULL
if (nrow(sites) > 0) {
  genotype_data <- extract_cis_dosages(
    pfile_prefix = opt$pfile_prefix,
    samples = cohort$samples,
    chr = sites$chr[1],
    from_bp = min(sites$pos) - opt$cis_window,
    to_bp = max(sites$pos) + opt$cis_window,
    maf = opt$maf,
    work_dir = file.path(opt$output_dir, "tmp_genotypes")
  )
  dosages <- genotype_data$dosages
  pairs <- build_cis_pairs(sites[, .(site_id, chr, pos)], genotype_data$variants, opt$cis_window)
}
log_info("Cis variants: ", if (is.null(dosages)) 0 else nrow(dosages), "; pairs: ", nrow(pairs))

## ---- 5. caNRD-editQTL scan --------------------------------------------------------------
pair_columns <- c("site_id", "variant_id", "celltype", "status", "beta", "se", "p", "p_site",
                  "mean_phi", "vif", "weakly_identifiable", "n_samples", "n_minor_allele_samples",
                  "n_celltypes", "refined")
if (nrow(pairs) > 0) {
  result <- caNRD_editQTL(
    bulk_editing = ratio,
    genotypes = dosages,
    proportions = cohort$proportions,
    theta = theta,
    theta_floor = opt$theta_floor,
    pairs = as.data.frame(pairs[, .(site_id, variant_id)]),
    coverage = coverage,
    covariates = if (ncol(covariates) > 0) covariates else NULL,
    min_mean_phi = opt$min_mean_phi,
    min_samples = opt$min_samples,
    min_minor_allele_samples = opt$min_minor_allele_samples,
    engine = "scan",
    refine = refine,
    impute_genotypes = "mean"
  )
  scan_info <- attr(result, "scan_info")
  result <- as.data.table(result)
  result[, distance := pairs$distance[match(paste(site_id, variant_id), paste(pairs$site_id, pairs$variant_id))]]
} else {
  result <- data.table(matrix(nrow = 0, ncol = length(pair_columns), dimnames = list(NULL, pair_columns)))
  result[, distance := numeric(0)]
  scan_info <- NULL
}
if (!"refined" %in% names(result)) {
  result[, refined := FALSE]
}

## ---- 6. outputs ------------------------------------------------------------------------------
keep_columns <- intersect(c(pair_columns, "distance"), names(result))
if (as_flag(opt$write_all_pairs)) {
  write_tsv(result[, keep_columns, with = FALSE], out(".all_pairs.tsv.gz"))
}
lead <- result[refined == TRUE]
lead[, chunk_id := chunk_id]
write_tsv(lead, out(".lead_variants.tsv"))

tested <- result[status == "tested"]
site_summary <- data.table(site_id = pass_sites)
site_summary[, n_variants_in_window := pairs[, .N, by = site_id][match(pass_sites, site_id), N]]
site_summary[is.na(n_variants_in_window), n_variants_in_window := 0L]
site_summary[, n_variants_tested := tested[, uniqueN(variant_id), by = site_id][match(pass_sites, site_id), V1]]
site_summary[is.na(n_variants_tested), n_variants_tested := 0L]
lead_site <- unique(lead[, .(site_id, lead_variant_id = variant_id, p_site_lead = p_site,
                             lead_distance = distance)], by = "site_id")
site_summary <- merge(site_summary, lead_site, by = "site_id", all.x = TRUE, sort = FALSE)
site_summary[, chunk_id := chunk_id]
write_tsv(site_summary, out(".site_summary.tsv"))

if (nrow(tested) > 0) {
  min_p <- tested[order(p), .(min_p_scan = p[1], variant_at_min_p = variant_id[1],
                              beta_at_min_p = beta[1], n_variants_tested = uniqueN(variant_id)),
                  by = .(site_id, celltype)]
} else {
  min_p <- data.table(site_id = character(0), celltype = character(0), min_p_scan = numeric(0),
                      variant_at_min_p = character(0), beta_at_min_p = numeric(0),
                      n_variants_tested = integer(0))
}
min_p[, chunk_id := chunk_id]
write_tsv(min_p, out(".site_celltype_min_p.tsv"))

if (!is.null(scan_info)) {
  write_tsv(as.data.table(scan_info), out(".scan_info.tsv"))
} else {
  write_tsv(data.table(site_id = character(0)), out(".scan_info.tsv"))
}
unlink(file.path(opt$output_dir, "tmp_genotypes"), recursive = TRUE)
log_info("canrd_editqtl_chunk.R finished chunk ", chunk_id, ": ", nrow(result), " result rows, ",
         nrow(lead_site), " lead variants")
