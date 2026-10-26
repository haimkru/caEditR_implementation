################################################################################
# R/io_tables.R
#
# Description:
#   Reading and writing of the pipeline's tabular formats:
#     * BED matrices  (#chr, start, end, <id>, sample columns), plain or
#       bgzipped, optionally only one tabix region;
#     * sample-by-column TSVs (first column = sample_id);
#     * tensorQTL covariate files (first column = ID, covariates in rows).
#   All sample IDs are kept as CHARACTER (AoU research_ids look numeric).
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   read_bed_header()      -> column names of a BED (first line only)
#   read_bed_matrix()      -> list(meta = 4-column data.table,
#                                  values = numeric matrix rows x samples)
#                             reads either the whole file or one tabix region,
#                             and optionally only selected sample columns
#   write_bed_matrix()     -> sorted BED, bgzip, tabix index
#   read_sample_table()    -> data.table with sample_id first, as character
#                             (accepts sample_id / research_id / IID / ID anywhere)
#   write_tensorqtl_covariates() -> transposed covariates (ID + samples)
#
# Inputs:  file paths given by the calling script
# Outputs: files written by the calling script
################################################################################

# Command that streams a (possibly gzipped) text file to stdout.
stream_command <- function(path) {
  if (grepl("\\.gz$|\\.bgz$", path)) {
    paste("gzip -dc", shQuote(path))
  } else {
    paste("cat", shQuote(path))
  }
}

# Column names of a BED/TSV file (first line), without reading the data.
read_bed_header <- function(path) {
  header_line <- system(paste(stream_command(path), "| head -n 1"), intern = TRUE)
  if (length(header_line) == 0) {
    stop_input("file is empty: ", path)
  }
  strsplit(header_line, "\t", fixed = TRUE)[[1]]
}

# Validate the first four BED columns.
check_bed_columns <- function(columns, id_column, path) {
  expected <- c("#chr", "start", "end", id_column)
  if (length(columns) < 5 || !identical(columns[1:4], expected)) {
    stop_input(path, " must start with columns ", paste(expected, collapse = ", "),
               " followed by at least one sample column; found ",
               paste(utils::head(columns, 5), collapse = ", "))
  }
  samples <- columns[-(1:4)]
  if (anyDuplicated(samples) > 0) {
    stop_input(path, " has duplicated sample columns, e.g. ", samples[duplicated(samples)][1])
  }
  invisible(samples)
}

# Read a BED matrix.
#   path       : BED or BED.GZ (bgzipped + .tbi next to it when region is used)
#   id_column  : name of the 4th column ("site_id" or "gene_id")
#   region     : NULL (whole file) or a tabix region "chr1:100-200"
#   samples    : NULL (all) or the sample columns to keep, in this order
# Returns list(meta = data.table(#chr, start, end, id), values = matrix)
read_bed_matrix <- function(path, id_column = "site_id", region = NULL, samples = NULL) {
  columns <- read_bed_header(path)
  file_samples <- check_bed_columns(columns, id_column, path)

  if (is.null(samples)) {
    samples <- file_samples
  }
  missing_samples <- setdiff(samples, file_samples)
  if (length(missing_samples) > 0) {
    stop_input(path, " is missing ", length(missing_samples), " requested samples, e.g. ",
               missing_samples[1])
  }

  # read only the needed columns: 4 metadata columns + requested samples
  keep_columns <- c(columns[1:4], samples)
  col_classes <- list(character = c("#chr", id_column), integer = c("start", "end"))

  if (is.null(region)) {
    dt <- data.table::fread(
      cmd = stream_command(path),
      sep = "\t",
      header = TRUE,
      select = keep_columns,
      colClasses = col_classes,
      na.strings = c("NA", "", "NaN"),
      check.names = FALSE,
      showProgress = FALSE
    )
  } else {
    # tabix -h prints the header line plus the rows of the region
    dt <- data.table::fread(
      cmd = paste("tabix -h", shQuote(path), shQuote(region)),
      sep = "\t",
      header = TRUE,
      select = keep_columns,
      colClasses = col_classes,
      na.strings = c("NA", "", "NaN"),
      check.names = FALSE,
      showProgress = FALSE
    )
  }

  meta <- dt[, 1:4, with = FALSE]
  values <- as.matrix(dt[, samples, with = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- meta[[id_column]]
  colnames(values) <- samples
  list(meta = meta, values = values)
}

# Write a BED matrix, bgzip it and (optionally) index it with tabix.
#   meta   : data.table with #chr, start, end, id (same row order as values)
#   values : numeric matrix rows x samples
#   path   : output path ending in .bed.gz
write_bed_matrix <- function(meta, values, path, digits = 6, index = TRUE) {
  if (nrow(meta) != nrow(values)) {
    stop("write_bed_matrix(): meta and values have different row counts", call. = FALSE)
  }
  value_dt <- data.table::as.data.table(values)
  if (is.double(values) && !is.null(digits)) {
    value_dt <- value_dt[, lapply(.SD, function(x) signif(x, digits))]
  }
  out <- cbind(data.table::copy(meta), value_dt)
  plain_path <- sub("\\.gz$", "", path)
  data.table::fwrite(out, plain_path, sep = "\t", na = "NA", quote = FALSE)
  bgzip_status <- system2("bgzip", c("-f", shQuote(plain_path)))
  if (bgzip_status != 0) {
    stop("bgzip failed for ", plain_path, call. = FALSE)
  }
  if (index) {
    tabix_status <- system2("tabix", c("-f", "-p", "bed", shQuote(path)))
    if (tabix_status != 0) {
      stop("tabix failed for ", path, call. = FALSE)
    }
  }
  invisible(path)
}

# Read a TSV whose first column is a sample identifier; IDs stay character.
#   id_column: expected name of the first column (e.g. "sample_id")
read_sample_table <- function(path, id_column = "sample_id") {
  # the ID column may be called sample_id, research_id, IID or ID (AoU files use
  # research_id; prepare_QTL genetic-PC files use sample_id) and need not be first
  columns <- read_bed_header(path)
  id_candidates <- c(id_column, "sample_id", "research_id", "#IID", "IID", "ID")
  found <- columns[tolower(columns) %in% tolower(id_candidates)]
  if (length(found) == 0) {
    stop_input(path, " needs a sample ID column (", paste(unique(id_candidates), collapse = ", "),
               "); found columns: ", paste(utils::head(columns, 6), collapse = ", "))
  }
  file_id_column <- found[1]
  col_classes <- list(character = file_id_column)
  dt <- data.table::fread(
    cmd = stream_command(path),
    sep = "\t",
    header = TRUE,
    colClasses = col_classes,
    check.names = FALSE,
    showProgress = FALSE
  )
  # standardize: ID column named id_column and placed first
  data.table::setnames(dt, file_id_column, id_column)
  data.table::setcolorder(dt, id_column)
  if (anyDuplicated(dt[[id_column]]) > 0) {
    stop_input(path, " has duplicated ", id_column, " values, e.g. ",
               dt[[id_column]][duplicated(dt[[id_column]])][1])
  }
  dt
}

# Convert a sample table to a numeric matrix (rows = samples).
sample_table_to_matrix <- function(dt, id_column = "sample_id") {
  values <- as.matrix(dt[, setdiff(names(dt), id_column), with = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- dt[[id_column]]
  values
}

# Write a samples x covariates matrix in tensorQTL's layout:
#   ID  sample_1  sample_2 ...
#   cov1  ...
write_tensorqtl_covariates <- function(covariate_matrix, path) {
  transposed <- t(covariate_matrix)
  out <- data.table::data.table(ID = rownames(transposed))
  out <- cbind(out, data.table::as.data.table(transposed))
  data.table::setnames(out, c("ID", colnames(transposed)))
  data.table::fwrite(out, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

# Read a tensorQTL covariate file back into a samples x covariates matrix.
read_tensorqtl_covariates <- function(path) {
  columns <- read_bed_header(path)
  if (columns[1] != "ID") {
    stop_input(path, " must have 'ID' as its first column (tensorQTL layout)")
  }
  dt <- data.table::fread(path, sep = "\t", header = TRUE, check.names = FALSE,
                          colClasses = list(character = "ID"))
  values <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- dt$ID
  colnames(values) <- columns[-1]
  t(values)
}

# Write a data.table as TSV (optionally gzip by .gz suffix).
write_tsv <- function(dt, path) {
  data.table::fwrite(dt, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}
