################################################################################
# R/constants_and_logging.R
#
# Description:
#   Pipeline-wide constants (version, file-name conventions) and the small
#   helpers every script uses: timestamped logging, input errors, cell-type
#   slugs and command-line option clean-up.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   - pipeline_version(): one string recorded in every output manifest.
#   - log_info(): prints "[UTC time] message" to stderr (captured by WDL logs).
#   - stop_input(): stops with "Input error: ..." so Terra logs are searchable.
#   - slugify_celltype(): "Monocyte/myeloid" -> "monocyte_myeloid"
#     (identical rule to CellTypeDeconvolution's slugify_cell_group()).
#   - clean_option_names(): optparse keeps dashes ("min-coverage"); we use
#     underscores inside R ("min_coverage").
#
# Inputs:  none (functions only)
# Outputs: none (functions only)
################################################################################

# Version of the caEditR_implementation pipeline; bump when outputs change.
pipeline_version <- function() {
  "0.1.0"
}

# Names of the three phenotype "analyses" produced by the deconvolution stage.
#   bulk       : quality-controlled bulk editing ratios (standard edQTL pipeline)
#   canrd_edit : caEditR::caNRD_edit() per-cell-type matrices
#   tca        : caEditR::TCA_Like() (CRAN TCA) per-cell-type matrices
analysis_names <- function() {
  c("bulk", "canrd_edit", "tca")
}

# Timestamped message to stderr. Every script logs through this function.
log_info <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  message("[", timestamp, "] ", ...)
}

# Stop with a clearly labelled input error (no R call stack in the message).
stop_input <- function(...) {
  stop(paste0("Input error: ", ...), call. = FALSE)
}

# File-name-safe, lower-case slug of a cell-type label.
slugify_celltype <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_", "", x)
  x <- gsub("_$", "", x)
  x
}

# optparse returns list names with dashes; convert them to underscores.
clean_option_names <- function(opt) {
  names(opt) <- gsub("-", "_", names(opt))
  opt
}

# Stop when a required option is missing or empty.
require_options <- function(opt, required) {
  for (name in required) {
    value <- opt[[name]]
    if (is.null(value) || length(value) == 0 || (is.character(value) && !nzchar(value[1]))) {
      stop_input("missing required option --", gsub("_", "-", name))
    }
  }
  invisible(TRUE)
}

# Stop when a file path does not exist or cannot be read.
require_readable_file <- function(path, label) {
  if (is.null(path) || !nzchar(path)) {
    stop_input(label, " path is empty")
  }
  if (grepl("^(gs|s3|https?)://", path)) {
    stop_input(label, " is an unresolved cloud URI (", path, "); WDL must localize it")
  }
  if (!file.exists(path)) {
    stop_input(label, " does not exist: ", path)
  }
  invisible(TRUE)
}

# Parse a comma-separated option into a character vector (empty -> character(0)).
split_comma_option <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    return(character(0))
  }
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

# Logical options arrive from WDL as "true"/"false" strings.
as_flag <- function(x) {
  if (is.logical(x)) {
    return(isTRUE(x))
  }
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}
