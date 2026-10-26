################################################################################
# R/samples.R
#
# Description:
#   Builds the analysis cohort: the samples present in EVERY input (editing
#   ratios, coverage, expression, proportions, genotypes, covariates and an
#   optional sample list), kept in the editing-BED column order.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   read_sample_list(): one ID per line; optional header research_id /
#     sample_id / ID (any case) is dropped; a tab-separated table (e.g. the
#     AoU eQTL sample lists) uses its research_id / sample_id column; blanks
#     and duplicates fail.
#   read_psam_sample_ids(): IID column of a PLINK2 .psam (#IID or #FID IID).
#   align_cohort_samples(): for each editing sample, record which inputs
#     contain it; keep samples present everywhere (and in the sample list).
#     Returns the kept IDs and a per-sample report table.
#
# Inputs:  sample ID vectors / file paths given by the calling script
# Outputs: none (returns R objects)
################################################################################

read_sample_list <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) {
    stop_input("sample list is empty: ", path)
  }
  id_headers <- c("research_id", "sample_id", "id", "iid", "#iid")
  if (any(grepl("\t", lines))) {
    # a table (e.g. research_id + ancestry): use the ID column, else the first column
    fields <- strsplit(lines, "\t", fixed = TRUE)
    header <- tolower(fields[[1]])
    column <- match(id_headers, header)
    column <- column[!is.na(column)]
    if (length(column) > 0) {
      ids <- vapply(fields[-1], `[`, character(1), column[1])
    } else {
      ids <- vapply(fields, `[`, character(1), 1)
    }
  } else {
    ids <- lines
    # drop an optional header line
    if (tolower(ids[1]) %in% id_headers) {
      ids <- ids[-1]
    }
  }
  ids <- trimws(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) {
    stop_input("sample list has a header but no sample IDs: ", path)
  }
  if (anyDuplicated(ids) > 0) {
    stop_input("sample list has duplicated IDs, e.g. ", ids[duplicated(ids)][1])
  }
  ids
}

read_psam_sample_ids <- function(path) {
  psam <- data.table::fread(path, sep = "\t", header = TRUE, colClasses = "character",
                            check.names = FALSE)
  id_column <- intersect(c("#IID", "IID"), names(psam))
  if (length(id_column) == 0) {
    stop_input(path, " has no #IID or IID column")
  }
  psam[[id_column[1]]]
}

# Intersect all inputs, keeping the order of the editing BED columns.
#   sample_sets: named list of character vectors (names = input labels);
#                the first element defines the order
#   sample_list: NULL or character vector restricting the cohort
align_cohort_samples <- function(sample_sets, sample_list = NULL) {
  reference_samples <- sample_sets[[1]]
  all_samples <- unique(c(unlist(sample_sets, use.names = FALSE), sample_list))

  report <- data.table::data.table(sample_id = all_samples)
  for (label in names(sample_sets)) {
    data.table::set(report, j = paste0("in_", label), value = all_samples %in% sample_sets[[label]])
  }
  if (!is.null(sample_list)) {
    report[, in_sample_list := sample_id %in% sample_list]
  }
  presence_columns <- grep("^in_", names(report), value = TRUE)
  report[, kept := Reduce(`&`, .SD), .SDcols = presence_columns]

  kept <- reference_samples[reference_samples %in% report[kept == TRUE, sample_id]]
  report[, kept := sample_id %in% kept]
  # order: kept samples first (in editing order), then the rest
  report[, order_key := match(sample_id, c(kept, setdiff(all_samples, kept)))]
  data.table::setorder(report, order_key)
  report[, order_key := NULL]

  list(samples = kept, report = report)
}
