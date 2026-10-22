#!/usr/bin/env Rscript
################################################################################
# scripts/convert_caeditr_csv_to_bed.R
#
# Description:
#   Adapter from the caEditR CSV matrix format (inst/cli/caedit.R:
#   sites x samples, first column = site id "chrom:pos:strand", first row =
#   sample ids) to the BED matrix format this pipeline reads
#   (#chr, start, end, site_id, samples...). Use it when editing levels /
#   coverage exist as caEditR CSVs instead of BEDs.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   1. Read the CSV (sample IDs kept as text).
#   2. Parse site ids -> chromosome, 1-based position; add the chromosome
#      prefix used by the genotypes (default "chr").
#   3. BED start = position - 1, end = position; sort; bgzip + tabix.
#
# Inputs:
#   --input-csv    caEditR CSV (ratios or coverage)
#   --output-bed   output .bed.gz
#   --chr-prefix   prefix added to the chromosome (default "chr")
# Outputs:
#   <output-bed> and <output-bed>.tbi
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

option_list <- list(
  make_option("--input-csv", type = "character"),
  make_option("--output-bed", type = "character"),
  make_option("--chr-prefix", type = "character", default = "chr")
)
opt <- clean_option_names(parse_args(OptionParser(option_list = option_list)))
require_options(opt, c("input_csv", "output_bed"))
require_readable_file(opt$input_csv, "input_csv")

## 1. read
csv <- fread(opt$input_csv, sep = ",", header = TRUE, check.names = FALSE,
             colClasses = list(character = 1))
setnames(csv, 1, "site_id")
values <- as.matrix(csv[, -1, with = FALSE])
storage.mode(values) <- "double"

## 2. coordinates
parsed <- parse_site_ids(csv$site_id)
meta <- data.table(`#chr` = paste0(opt$chr_prefix, parsed$chrom),
                   start = as.integer(parsed$pos - 1),
                   end = as.integer(parsed$pos),
                   site_id = csv$site_id)

## 3. sort and write
ord <- order(meta[["#chr"]], meta$start, meta$site_id)
write_bed_matrix(meta[ord], values[ord, , drop = FALSE], opt$output_bed, index = TRUE)
log_info("Wrote ", nrow(meta), " sites x ", ncol(values), " samples to ", opt$output_bed,
         " (run scripts/index_bed_matrix.sh or the WDL IndexBedMatrix task for natural chromosome order)")
