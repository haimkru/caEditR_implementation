################################################################################
# R/chunk_data.R
#
# Description:
#   Loads the editing data of ONE chunk of sites (one Terra scatter shard):
#   the bulk ratios and the coverage of the chunk's sites for the cohort
#   samples, plus the matching theta rows. Shared by the deconvolution and
#   the caNRD-editQTL chunk scripts, so both see identical inputs and the
#   same site QC.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   1. Look up the chunk's tabix region and its site list.
#   2. tabix-read ratios and coverage of that region (cohort columns only).
#   3. Keep exactly the chunk's sites (a region can touch a neighbouring
#      chunk's site at the same position), in chunk order.
#   4. Validate: same sites in both files, ratios in [0, 1], coverage >= 0;
#      a ratio with zero coverage is set to NA.
#   5. Return meta, ratio, coverage, theta (sites x cell types).
#
# Inputs:  paths / IDs from the calling script
# Outputs: none (returns R objects)
################################################################################

load_chunk_editing <- function(chunk_id, chunks_path, site_chunks_path, ratios_path,
                               coverage_path, site_theta_path, samples, celltypes) {
  chunks <- data.table::fread(chunks_path, sep = "\t", colClasses = list(character = "chunk_id"))
  site_chunks <- data.table::fread(site_chunks_path, sep = "\t", colClasses = "character")
  if (!chunk_id %in% chunks$chunk_id) {
    stop_input("chunk ", chunk_id, " is not listed in ", chunks_path)
  }
  region <- chunks$tabix_region[chunks$chunk_id == chunk_id]
  chunk_sites <- site_chunks$site_id[site_chunks$chunk_id == chunk_id]
  log_info("Chunk ", chunk_id, ": region ", region, ", ", length(chunk_sites), " sites")

  ratio_data <- read_bed_matrix(ratios_path, "site_id", region = region, samples = samples)
  coverage_data <- read_bed_matrix(coverage_path, "site_id", region = region, samples = samples)

  keep <- match(chunk_sites, ratio_data$meta$site_id)
  keep_coverage <- match(chunk_sites, coverage_data$meta$site_id)
  if (anyNA(keep) || anyNA(keep_coverage)) {
    stop_input("chunk ", chunk_id, ": sites missing from the editing BEDs, e.g. ",
               chunk_sites[is.na(keep) | is.na(keep_coverage)][1])
  }
  meta <- ratio_data$meta[keep]
  ratio <- ratio_data$values[keep, , drop = FALSE]
  coverage <- coverage_data$values[keep_coverage, , drop = FALSE]

  if (any(ratio < 0 | ratio > 1, na.rm = TRUE)) {
    stop_input("editing ratios must be in [0, 1] (chunk ", chunk_id, ")")
  }
  if (any(coverage < 0, na.rm = TRUE)) {
    stop_input("coverage must be non-negative (chunk ", chunk_id, ")")
  }
  coverage[is.na(coverage)] <- 0
  ratio[coverage == 0] <- NA_real_

  theta_dt <- data.table::fread(site_theta_path, sep = "\t", colClasses = list(character = "site_id"),
                                check.names = FALSE)
  theta_dt <- theta_dt[match(chunk_sites, theta_dt$site_id)]
  if (anyNA(theta_dt$site_id)) {
    stop_input("site_theta is missing sites of chunk ", chunk_id)
  }
  theta <- as.matrix(theta_dt[, celltypes, with = FALSE])
  rownames(theta) <- chunk_sites

  list(meta = meta, ratio = ratio, coverage = coverage, theta = theta, region = region)
}

# Read the cohort, cell types and proportions written by prepare_cohort.R.
load_cohort_tables <- function(samples_path, proportions_path) {
  samples <- read_sample_list(samples_path)
  proportions_dt <- read_sample_table(proportions_path, "sample_id")
  proportions <- sample_table_to_matrix(proportions_dt, "sample_id")
  if (!identical(rownames(proportions), samples)) {
    proportions <- proportions[samples, , drop = FALSE]
  }
  list(samples = samples, proportions = proportions, celltypes = colnames(proportions))
}
