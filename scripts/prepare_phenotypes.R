#!/usr/bin/env Rscript
################################################################################
# scripts/prepare_phenotypes.R
#
# Description:
#   Prepares ONE editing matrix for tensorQTL (WDL task PreparePhenotypes):
#   the bulk QC-passing ratios ("standard edQTL pipeline", as Li et al. 2022)
#   or one caNRD-edit / TCA cell-type matrix. Produces the phenotype BED
#   (rank-based inverse normal transformed) and the matching tensorQTL
#   covariate file (supplied covariates + editing PCs of this matrix).
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   1. Read the aligned covariates (sample_id + numeric columns); stream the
#      merged BED (sites x cohort samples) in blocks of --block-size sites, so
#      memory does not grow with the number of sites (100k x 10k at AoU).
#   2. Per site: fraction missing (NA = no reads in bulk, or caNRD excluded
#      the cell type / estimate), sd. Drop sites with missing fraction >
#      --max-missing-fraction or no variation; write the report.
#   3. Impute remaining NAs with the site mean; INT each site across samples
#      (Blom offset, = RNOmni::RankNorm used by prepare_QTL). Logit then INT
#      (Li et al.) is identical because INT only uses ranks.
#   4. Editing PCs: top --n-pcs principal components of the INT matrix, from
#      the samples x samples Gram matrix accumulated over the blocks.
#   5. Covariates = supplied covariates (constant columns dropped) + PCs, in
#      tensorQTL layout (ID + one column per sample, same order as the BED).
#   6. Write the INT BED (bgzip + tabix), PCs, variance explained,
#      covariates and a summary.
#
# Inputs:
#   --input-bed             merged BED matrix (bgzipped)
#   --covariates            covariates_aligned.tsv (prepare_cohort.R)
#   --n-pcs                 editing PCs to compute (default 15, as the draft)
#   --max-missing-fraction  drop sites with more missing samples (default 0.5)
#   --label                 output label (default: input basename w/o .bed.gz)
#   --block-size            sites per streamed block (default 2000)
#   --output-dir
# Outputs (in --output-dir):
#   <label>.INT.bed.gz(.tbi)          tensorQTL phenotype BED
#   <label>.tensorqtl_covariates.tsv  tensorQTL covariates (ID + samples)
#   <label>.phenotype_PCs.tsv         sample_id + editing_PC1..k
#   <label>.pc_variance_explained.tsv
#   <label>.phenotype_filter_report.tsv
#   <label>.phenotype_summary.tsv
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
  make_option("--input-bed", type = "character"),
  make_option("--covariates", type = "character"),
  make_option("--n-pcs", type = "integer", default = 15L),
  make_option("--max-missing-fraction", type = "double", default = 0.5),
  make_option("--label", type = "character", default = ""),
  make_option("--block-size", type = "integer", default = 2000L),
  make_option("--output-dir", type = "character", default = "outputs")
)
opt <- clean_option_names(parse_args(OptionParser(option_list = option_list)))
require_options(opt, c("input_bed", "covariates"))
require_readable_file(opt$input_bed, "input_bed")
require_readable_file(opt$covariates, "covariates")
label <- if (nzchar(opt$label)) opt$label else sub("\\.bed\\.gz$", "", basename(opt$input_bed))
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(suffix) file.path(opt$output_dir, paste0(label, suffix))
log_info("prepare_phenotypes.R (caEditR_edQTL ", pipeline_version(), ") ", label)

## ---- 1. inputs -----------------------------------------------------------------
covariates_dt <- read_sample_table(opt$covariates, "sample_id")
covariates <- sample_table_to_matrix(covariates_dt, "sample_id")
samples <- rownames(covariates)
bed_columns <- read_bed_header(opt$input_bed)
check_bed_columns(bed_columns, "site_id", opt$input_bed)
missing_samples <- setdiff(samples, bed_columns[-(1:4)])
if (length(missing_samples) > 0) {
  stop_input(opt$input_bed, " lacks ", length(missing_samples), " covariate samples, e.g. ", missing_samples[1])
}

## ---- 2-3. stream the BED in blocks of sites --------------------------------------------
# Per block: site filter, mean imputation, INT, append to the output BED, and
# add the block to the samples x samples Gram matrix used for the PCs. Memory
# is one block + the Gram matrix, whatever the number of sites.
int_plain <- sub("\\.gz$", "", out(".INT.bed.gz"))
gram <- matrix(0, length(samples), length(samples), dimnames = list(samples, samples))
total_variance <- 0
report_blocks <- list()
n_input <- 0
n_kept <- 0
header_line <- paste(bed_columns, collapse = "\t")
connection <- pipe(stream_command(opt$input_bed), open = "r")
invisible(readLines(connection, n = 1))                       # skip the header
repeat {
  lines <- readLines(connection, n = opt$block_size)
  if (length(lines) == 0) break
  block <- fread(text = c(header_line, lines), sep = "\t", header = TRUE,
                 select = c(bed_columns[1:4], samples),
                 colClasses = list(character = c("#chr", "site_id"), integer = c("start", "end")),
                 na.strings = c("NA", "", "NaN"), check.names = FALSE, showProgress = FALSE)
  values <- as.matrix(block[, samples, with = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- block$site_id
  n_input <- n_input + nrow(values)

  # site filter
  report <- data.table(
    site_id = block$site_id,
    n_missing = rowSums(is.na(values)),
    fraction_missing = rowMeans(is.na(values)),
    mean_value = rowMeans(values, na.rm = TRUE),
    sd_value = apply(values, 1, stats::sd, na.rm = TRUE)
  )
  report[, reason := fcase(
    fraction_missing >= 1, "all_missing",
    fraction_missing > opt$max_missing_fraction, "too_many_missing",
    is.na(sd_value) | sd_value < 1e-8, "no_variation",
    default = "pass"
  )]
  report[, kept := reason == "pass"]
  report_blocks[[length(report_blocks) + 1]] <- report
  keep <- report$kept
  if (!any(keep)) next

  # impute + INT
  int_values <- inverse_normal_transform_rows(impute_row_means(values[keep, , drop = FALSE]))
  rm(values)

  # append to the output BED (header written with the first kept block)
  out_block <- cbind(block[keep, 1:4, with = FALSE],
                     as.data.table(signif(int_values, 6)))
  fwrite(out_block, int_plain, sep = "\t", na = "NA", quote = FALSE, append = n_kept > 0,
         col.names = n_kept == 0)
  n_kept <- n_kept + nrow(int_values)

  # Gram matrix of the site-centred INT values (samples x samples)
  centred <- int_values - rowMeans(int_values)
  gram <- gram + crossprod(centred)
  total_variance <- total_variance + sum(centred^2)
  rm(int_values, centred, out_block, block)
}
close(connection)
report <- rbindlist(report_blocks)
write_tsv(report, out(".phenotype_filter_report.tsv"))
log_info("Input: ", n_input, " sites x ", length(samples), " samples; kept ", n_kept, " sites")

# the WDL reads this number to skip tensorQTL for an empty phenotype set
writeLines(as.character(n_kept), out(".n_sites_kept.txt"))
if (n_kept == 0) {
  # header-only BED keeps the output contract
  writeLines(paste(c(bed_columns[1:4], samples), collapse = "\t"), int_plain)
}
system2("bgzip", c("-f", shQuote(int_plain)))
system2("tabix", c("-f", "-p", "bed", shQuote(out(".INT.bed.gz"))))

## ---- 4. editing PCs from the Gram matrix ---------------------------------------------------
pcs <- pcs_from_gram(gram, opt$n_pcs, total_variance, n_kept)
rm(gram)
log_info("Editing PCs: ", ncol(pcs$scores), " (", pcs$method, ")")
write_tsv(data.table(sample_id = rownames(pcs$scores), as.data.table(pcs$scores)),
          out(".phenotype_PCs.tsv"))
write_tsv(data.table(pc = colnames(pcs$scores), variance_explained = pcs$variance_explained),
          out(".pc_variance_explained.tsv"))

## ---- 5. covariates ---------------------------------------------------------------------
covariate_sd <- apply(covariates, 2, stats::sd)
constant_covariates <- names(covariate_sd)[!is.finite(covariate_sd) | covariate_sd == 0]
if (length(constant_covariates) > 0) {
  log_info("Dropping constant covariates: ", paste(constant_covariates, collapse = ", "))
}
all_covariates <- cbind(covariates[, setdiff(colnames(covariates), constant_covariates), drop = FALSE],
                        pcs$scores[rownames(covariates), , drop = FALSE])
write_tensorqtl_covariates(all_covariates, out(".tensorqtl_covariates.tsv"))

summary_table <- data.table(
  field = c("label", "pipeline_version", "n_sites_input", "n_sites_kept", "n_samples",
            "n_supplied_covariates", "n_constant_covariates_dropped", "n_editing_pcs", "pc_method",
            "max_missing_fraction", "block_size"),
  value = c(label, pipeline_version(), n_input, n_kept, length(samples),
            ncol(covariates), length(constant_covariates), ncol(pcs$scores), pcs$method,
            opt$max_missing_fraction, opt$block_size)
)
write_tsv(summary_table, out(".phenotype_summary.tsv"))
log_info("prepare_phenotypes.R finished ", label)
