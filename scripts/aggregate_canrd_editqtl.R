#!/usr/bin/env Rscript
################################################################################
# scripts/aggregate_canrd_editqtl.R
#
# Description:
#   Gather stage of the end-to-end caNRD-editQTL scan (WDL task
#   AggregateCanrdEditqtl): combines the per-chunk outputs and applies the
#   multiple-testing correction, producing the finalized summary-statistic
#   tables.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   1. Read every chunk's site summary, lead-variant table and per site x cell
#      type minimum-p table (file lists, one path per line).
#   2. Site level ("is there an edQTL at this site in any cell type?"):
#        p_site_lead  = exact LRT p of the lead variant (K-df joint test)
#        p_site_adj   = min(1, p_site_lead x n_variants_tested)  (Bonferroni
#                       over the site's cis variants, as eigenMT-style
#                       gene-level correction)
#        q_site       = Storey q-value of p_site_adj across sites
#        edSite       = q_site <= --fdr
#   3. Cell type at the lead variant (hierarchical): for edSites, a cell type
#      is called when its exact LRT p x (cell types tested at the site)
#      <= --celltype-alpha. Weakly identifiable effects (VIF > 10) are flagged.
#   4. Cell type scan level (independent per cell type): min scan p over the
#      site's variants x n variants tested -> Storey q per cell type.
#   5. Write the three tables and a summary of counts per cell type.
#
# Inputs:
#   --site-summary-list, --lead-list, --site-celltype-list  text files of paths
#   --fdr (0.05) --celltype-alpha (0.05) --qvalue-lambda (optional, e.g. 0.85)
#   --prefix, --output-dir
# Outputs (in --output-dir):
#   <prefix>.canrd_editqtl.site_level.tsv
#   <prefix>.canrd_editqtl.lead_site_celltype.tsv
#   <prefix>.canrd_editqtl.celltype_scan_level.tsv
#   <prefix>.canrd_editqtl.summary.tsv
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
  make_option("--site-summary-list", type = "character"),
  make_option("--lead-list", type = "character"),
  make_option("--site-celltype-list", type = "character"),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--celltype-alpha", type = "double", default = 0.05),
  make_option("--qvalue-lambda", type = "character", default = ""),
  make_option("--prefix", type = "character", default = "cohort"),
  make_option("--output-dir", type = "character", default = "outputs")
)
opt <- clean_option_names(parse_args(OptionParser(option_list = option_list)))
require_options(opt, c("site_summary_list", "lead_list", "site_celltype_list"))
qvalue_lambda <- if (nzchar(opt$qvalue_lambda)) as.numeric(opt$qvalue_lambda) else NULL
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(suffix) file.path(opt$output_dir, paste0(opt$prefix, ".canrd_editqtl", suffix))
log_info("aggregate_canrd_editqtl.R (caEditR_edQTL ", pipeline_version(), ")")

read_all <- function(list_file) {
  paths <- readLines(list_file, warn = FALSE)
  paths <- paths[nzchar(paths)]
  tables <- lapply(paths, function(p) {
    if (file.info(p)$size == 0) return(NULL)
    fread(p, sep = "\t", header = TRUE, colClasses = list(character = "site_id"))
  })
  rbindlist(tables, fill = TRUE)
}

## ---- 1. inputs -------------------------------------------------------------------
site_summary <- read_all(opt$site_summary_list)
lead <- read_all(opt$lead_list)
site_celltype <- read_all(opt$site_celltype_list)
log_info("Sites: ", nrow(site_summary), "; lead rows: ", nrow(lead),
         "; site x cell type rows: ", nrow(site_celltype))

## ---- 2. site level ---------------------------------------------------------------------
site_level <- copy(site_summary)
site_level[, tested := n_variants_tested > 0 & is.finite(p_site_lead)]
site_level[, p_site_adjusted := ifelse(tested, bonferroni_within(p_site_lead, n_variants_tested), NA_real_)]
site_level[, q_site := storey_qvalue(p_site_adjusted, qvalue_lambda)]
q_method_site <- attr(site_level$q_site, "method")
site_level[, edsite := tested & !is.na(q_site) & q_site <= opt$fdr]
setorder(site_level, q_site, na.last = TRUE)
write_tsv(site_level, out(".site_level.tsv"))
log_info("edSites (q <= ", opt$fdr, ", ", q_method_site, "): ", sum(site_level$edsite, na.rm = TRUE))

## ---- 3. cell type at the lead variant ---------------------------------------------------
lead_celltype <- copy(lead)
if (nrow(lead_celltype) > 0) {
  lead_celltype[, n_celltypes_tested_at_lead := sum(status == "tested"), by = site_id]
  lead_celltype <- merge(lead_celltype,
                         site_level[, .(site_id, n_variants_tested, p_site_adjusted, q_site, edsite)],
                         by = "site_id", all.x = TRUE)
  lead_celltype[, p_celltype_adjusted := ifelse(status == "tested",
                                                pmin(1, p * pmax(n_celltypes_tested_at_lead, 1)),
                                                NA_real_)]
  lead_celltype[, significant_celltype := edsite %in% TRUE & status == "tested" &
                  !is.na(p_celltype_adjusted) & p_celltype_adjusted <= opt$celltype_alpha]
  setorder(lead_celltype, q_site, site_id, celltype, na.last = TRUE)
}
write_tsv(lead_celltype, out(".lead_site_celltype.tsv"))

## ---- 4. cell type scan level -------------------------------------------------------------
scan_level <- copy(site_celltype)
q_method_scan <- "none"
if (nrow(scan_level) > 0) {
  scan_level[, p_celltype_adjusted := bonferroni_within(min_p_scan, n_variants_tested)]
  scan_level[, q_celltype := {
    q <- storey_qvalue(p_celltype_adjusted, qvalue_lambda)
    q_method_scan <<- attr(q, "method")
    q
  }, by = celltype]
  scan_level[, significant := !is.na(q_celltype) & q_celltype <= opt$fdr]
  setorder(scan_level, celltype, q_celltype)
}
write_tsv(scan_level, out(".celltype_scan_level.tsv"))

## ---- 5. summary ---------------------------------------------------------------------------
celltypes <- sort(unique(c(lead_celltype$celltype, scan_level$celltype)))
summary_rows <- lapply(celltypes, function(h) {
  data.table(
    celltype = h,
    n_sites_tested_at_lead = if (nrow(lead_celltype)) lead_celltype[celltype == h & status == "tested", .N] else 0L,
    n_significant_lead_hierarchical = if (nrow(lead_celltype)) lead_celltype[celltype == h & significant_celltype == TRUE, .N] else 0L,
    n_weakly_identifiable_calls = if (nrow(lead_celltype)) lead_celltype[celltype == h & significant_celltype == TRUE & weakly_identifiable %in% TRUE, .N] else 0L,
    n_sites_tested_scan = if (nrow(scan_level)) scan_level[celltype == h, .N] else 0L,
    n_significant_scan_q = if (nrow(scan_level)) scan_level[celltype == h & significant == TRUE, .N] else 0L
  )
})
summary_table <- rbindlist(c(
  list(data.table(celltype = "any (site level)",
                  n_sites_tested_at_lead = site_level[tested == TRUE, .N],
                  n_significant_lead_hierarchical = site_level[edsite == TRUE, .N],
                  n_weakly_identifiable_calls = NA_integer_,
                  n_sites_tested_scan = NA_integer_,
                  n_significant_scan_q = NA_integer_)),
  summary_rows
))
summary_table[, fdr := opt$fdr]
summary_table[, q_method_site := q_method_site]
summary_table[, q_method_scan := q_method_scan]
write_tsv(summary_table, out(".summary.tsv"))
log_info("aggregate_canrd_editqtl.R finished")
