#!/usr/bin/env Rscript
################################################################################
# scripts/summarize_edqtl_results.R
#
# Description:
#   Final stage (WDL task SummarizeEdqtlResults): puts the three edQTL
#   approaches of the draft (Methods 4.6) side by side, per editing site and
#   cell type:
#     1. bulk tensorQTL                       ("standard edQTL pipeline")
#     2. tensorQTL on caNRD-edit / TCA cell-type matrices
#     3. end-to-end caNRD-editQTL on the bulk ratios
#   and counts, per cell type, edSites found by a cell-type approach but not
#   by the bulk scan ("cell-type-specific edSites").
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   1. Read the combined tensorQTL table and the caNRD-editQTL lead table.
#   2. Long table: one row per (site, approach, cell type) with the call,
#      effect size, p-value and lead variant.
#   3. Bulk edSites; per approach x cell type: calls, calls also in bulk,
#      calls not in bulk.
#
# Inputs:
#   --tensorqtl-results   <prefix>.tensorqtl_cis.all_analyses.tsv.gz
#   --canrd-lead          <prefix>.canrd_editqtl.lead_site_celltype.tsv
#   --celltypes           celltypes.tsv (celltype, slug) from prepare_cohort.R
#   --prefix, --output-dir
# Outputs (in --output-dir):
#   <prefix>.edsite_calls_long.tsv.gz
#   <prefix>.edsite_overview.tsv
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
  make_option("--tensorqtl-results", type = "character", default = ""),
  make_option("--canrd-lead", type = "character", default = ""),
  make_option("--celltypes", type = "character"),
  make_option("--prefix", type = "character", default = "cohort"),
  make_option("--output-dir", type = "character", default = "outputs")
)
opt <- clean_option_names(parse_args(OptionParser(option_list = option_list)))
require_options(opt, "celltypes")
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(suffix) file.path(opt$output_dir, paste0(opt$prefix, suffix))
log_info("summarize_edqtl_results.R (caEditR_edQTL ", pipeline_version(), ")")
celltype_table <- fread(opt$celltypes, sep = "\t")

## ---- 1-2. long table -------------------------------------------------------------------
long_rows <- list()
if (nzchar(opt$tensorqtl_results) && file.exists(opt$tensorqtl_results)) {
  tq <- fread(cmd = stream_command(opt$tensorqtl_results), sep = "\t", colClasses = list(character = c("site_id", "variant_id")))
  # one tensorQTL source per set: nominal (default) is preferred when both were run
  if ("source" %in% names(tq) && "cis_nominal" %in% tq$source) tq <- tq[source == "cis_nominal"]
  tq[, celltype := celltype_table$celltype[match(celltype_slug, celltype_table$slug)]]
  tq[analysis == "bulk", celltype := "bulk"]
  long_rows[[length(long_rows) + 1]] <- tq[, .(
    site_id,
    approach = paste0(analysis, "_tensorqtl"),
    celltype,
    lead_variant_id = variant_id,
    effect = slope,
    p_value = p_for_q,
    q_value = q,
    significant = edsite
  )]
}
if (nzchar(opt$canrd_lead) && file.exists(opt$canrd_lead) && file.info(opt$canrd_lead)$size > 0) {
  cl <- fread(opt$canrd_lead, sep = "\t", colClasses = list(character = c("site_id", "variant_id")))
  # only cell types that were actually tested at the site (gated / aliased ones are not calls)
  cl <- cl[status == "tested"]
  if (nrow(cl) > 0) {
    long_rows[[length(long_rows) + 1]] <- cl[, .(
      site_id,
      approach = "canrd_editqtl",
      celltype,
      lead_variant_id = variant_id,
      effect = beta,
      p_value = p,
      q_value = q_site,
      significant = significant_celltype %in% TRUE
    )]
  }
}
calls <- rbindlist(long_rows, fill = TRUE)
write_tsv(calls, out(".edsite_calls_long.tsv.gz"))

## ---- 3. overview -------------------------------------------------------------------------
bulk_edsites <- if (nrow(calls)) calls[approach == "bulk_tensorqtl" & significant == TRUE, unique(site_id)] else character(0)
overview <- if (nrow(calls)) {
  calls[, .(
    n_tested = .N,
    n_edsites = sum(significant, na.rm = TRUE),
    n_edsites_also_bulk = sum(significant & site_id %in% bulk_edsites, na.rm = TRUE),
    n_edsites_not_in_bulk = sum(significant & !site_id %in% bulk_edsites, na.rm = TRUE)
  ), by = .(approach, celltype)]
} else {
  data.table(approach = character(0), celltype = character(0))
}
setorder(overview, approach, celltype)
write_tsv(overview, out(".edsite_overview.tsv"))
log_info("summarize_edqtl_results.R finished: ", length(bulk_edsites), " bulk edSites")
