#!/usr/bin/env Rscript
################################################################################
# scripts/aggregate_tensorqtl.R
#
# Description:
#   Gather stage of the "standard" tensorQTL edQTL scans (WDL task
#   AggregateTensorqtl). One phenotype set = bulk, caNRD-edit x cell type or
#   TCA x cell type. Default input = tensorQTL cis NOMINAL results (every site
#   x variant pair, parquet per chromosome); the cheaper default, since the
#   multiple-testing correction is done by the user, not by tensorQTL.
#   Optional input = tensorQTL cis PERMUTATION results (cis_qtl.txt.gz).
#   Output: one row per site and phenotype set with the lead variant, the
#   number of variants tested and the raw p-value / slope / SE.
#
# Version:      0.2.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   Nominal (--cis-nominal-list):
#     1. Group the parquet files by phenotype set (file name
#        <prefix>.<analysis>.<celltype>.cis_qtl_pairs.<chr>.parquet).
#     2. Read one chromosome file at a time (memory-safe at AoU size); per
#        site keep the variant with the smallest pval_nominal and count the
#        variants tested.
#     3. Convenience columns: p_for_q = min(1, lead p x variants tested)
#        (Bonferroni over the site's cis variants), q = Storey q across sites
#        of the set, edsite = q <= --fdr. Raw p-values stay in the table.
#   Permutation (--cis-permutation-list): read the tables; p_for_q =
#     pval_beta (pval_perm if beta is unavailable); q, edsite as above.
#   Write the combined table and counts per set.
#
# Inputs:
#   --cis-nominal-list      text file, one tensorQTL nominal parquet per line
#   --cis-permutation-list  text file, one cis_qtl.txt.gz per line (optional)
#   --fdr (0.05) --qvalue-lambda (optional) --prefix --output-dir
# Outputs (in --output-dir):
#   <prefix>.tensorqtl_cis.all_analyses.tsv.gz   one row per site x set
#   <prefix>.tensorqtl_cis.summary.tsv           counts per set
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
  make_option("--cis-nominal-list", type = "character", default = ""),
  make_option("--cis-permutation-list", type = "character", default = ""),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--qvalue-lambda", type = "character", default = ""),
  make_option("--prefix", type = "character", default = "cohort"),
  make_option("--output-dir", type = "character", default = "outputs")
)
opt <- clean_option_names(parse_args(OptionParser(option_list = option_list)))
qvalue_lambda <- if (nzchar(opt$qvalue_lambda)) as.numeric(opt$qvalue_lambda) else NULL
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(suffix) file.path(opt$output_dir, paste0(opt$prefix, ".tensorqtl_cis", suffix))
log_info("aggregate_tensorqtl.R (caEditR_implementation ", pipeline_version(), ")")

read_path_list <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(character(0))
  paths <- readLines(path, warn = FALSE)
  paths[nzchar(paths)]
}
# "<prefix>.<analysis>.<celltype>" -> analysis and cell-type slug
set_pattern <- "^(.*)\\.(bulk|canrd_edit|tca)\\.([A-Za-z0-9_]+)$"
parse_set <- function(label) {
  if (!grepl(set_pattern, label)) {
    stop_input("cannot parse analysis / cell type from ", label,
               " (expected <prefix>.<bulk|canrd_edit|tca>.<celltype>)")
  }
  list(analysis = sub(set_pattern, "\\2", label), celltype_slug = sub(set_pattern, "\\3", label))
}

results <- list()

## ---- nominal results ----------------------------------------------------------------
nominal_files <- read_path_list(opt$cis_nominal_list)
if (length(nominal_files) > 0) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("the arrow R package is needed to read tensorQTL parquet files", call. = FALSE)
  }
  nominal_label <- sub("\\.cis_qtl_pairs\\..*\\.parquet$", "", basename(nominal_files))
  for (label in unique(nominal_label)) {
    set <- parse_set(label)
    lead_by_chromosome <- lapply(nominal_files[nominal_label == label], function(path) {
      pairs <- as.data.table(arrow::read_parquet(path))
      if (nrow(pairs) == 0) return(NULL)
      pairs[, n_variants_tested := .N, by = phenotype_id]
      lead <- pairs[pairs[, .I[which.min(pval_nominal)], by = phenotype_id]$V1]
      rm(pairs)
      lead
    })
    lead <- rbindlist(lead_by_chromosome)
    if (nrow(lead) == 0) next
    setnames(lead, "phenotype_id", "site_id")
    lead[, `:=`(analysis = set$analysis, celltype_slug = set$celltype_slug, source = "cis_nominal",
                p_for_q = bonferroni_within(pval_nominal, n_variants_tested))]
    results[[length(results) + 1]] <- lead
    log_info("Nominal set ", label, ": ", nrow(lead), " sites")
  }
}

## ---- permutation results (optional) -------------------------------------------------------
permutation_files <- read_path_list(opt$cis_permutation_list)
for (path in permutation_files) {
  label <- sub("(\\.INT)?\\.cis_qtl\\.txt\\.gz$", "", basename(path))
  set <- parse_set(label)
  dt <- fread(cmd = stream_command(path), sep = "\t", header = TRUE,
              colClasses = list(character = c("phenotype_id", "variant_id")))
  if ("qval" %in% names(dt)) setnames(dt, "qval", "tensorqtl_qval")
  if ("pval_nominal_threshold" %in% names(dt)) setnames(dt, "pval_nominal_threshold", "tensorqtl_pval_nominal_threshold")
  setnames(dt, "phenotype_id", "site_id")
  if ("num_var" %in% names(dt)) dt[, n_variants_tested := num_var]
  p_beta_ok <- any(is.finite(dt$pval_beta))
  dt[, `:=`(analysis = set$analysis, celltype_slug = set$celltype_slug, source = "cis_permutations",
            p_for_q = if (p_beta_ok) pval_beta else pval_perm)]
  results[[length(results) + 1]] <- dt
  log_info("Permutation set ", label, ": ", nrow(dt), " sites")
}

if (length(results) == 0) {
  stop_input("no tensorQTL results given (--cis-nominal-list / --cis-permutation-list)")
}

## ---- q-values (convenience) and output -------------------------------------------------------
all_results <- rbindlist(results, fill = TRUE)
all_results[, q := storey_qvalue(p_for_q, qvalue_lambda), by = .(source, analysis, celltype_slug)]
all_results[, edsite := !is.na(q) & q <= opt$fdr]
setcolorder(all_results, intersect(c("source", "analysis", "celltype_slug", "site_id", "variant_id",
                                     "n_variants_tested", "pval_nominal", "slope", "slope_se",
                                     "p_for_q", "q", "edsite"), names(all_results)))
write_tsv(all_results, out(".all_analyses.tsv.gz"))

summary_table <- all_results[, .(n_sites_tested = .N,
                                 n_edsites = sum(edsite),
                                 median_lead_pval_nominal = stats::median(pval_nominal, na.rm = TRUE)),
                             by = .(source, analysis, celltype_slug)]
summary_table[, fdr := opt$fdr]
setorder(summary_table, source, analysis, celltype_slug)
write_tsv(summary_table, out(".summary.tsv"))
log_info("aggregate_tensorqtl.R finished")
