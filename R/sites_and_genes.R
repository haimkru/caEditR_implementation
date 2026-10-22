################################################################################
# R/sites_and_genes.R
#
# Description:
#   Everything that is "per editing site" but not per sample:
#     * parsing caEditR site IDs ("chrom:pos:strand", no chr prefix, 1-based);
#     * reading gene records from a GTF;
#     * mapping each site to its host gene (needed for theta, the per-cell-type
#       host-gene expression that turns proportions into RNA shares);
#     * per-site coverage quality control (Li et al. 2022 style);
#     * splitting the sites into genomic chunks for the Terra scatter.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   parse_site_ids(): split "1:12345:+" into chrom "1", pos 12345, strand "+".
#   read_gtf_genes(): keep feature == "gene"; parse gene_id, gene_name,
#     gene_type (fallback gene_biotype).
#   map_sites_to_host_genes(): overlap each site with gene bodies on the same
#     chromosome; prefer genes on the site's strand; among candidates pick the
#     gene with the highest mean expression; record the rule used.
#   site_coverage_qc(): per site count samples with >= min_coverage reads,
#     mean / sd of the ratio over covered samples; pass when enough samples
#     are covered and the ratio varies.
#   define_site_chunks(): consecutive sites on one chromosome, at most
#     sites_per_chunk per chunk, never splitting two sites at the same
#     position; one tabix region per chunk.
#
# Inputs:  R objects / file paths given by the calling script
# Outputs: none (returns R objects)
################################################################################

parse_site_ids <- function(site_ids) {
  parts <- strsplit(site_ids, ":", fixed = TRUE)
  n_parts <- lengths(parts)
  if (any(n_parts < 2 | n_parts > 3)) {
    stop_input("site ids must look like 'chrom:pos:strand' or 'chrom:pos'; bad example: ",
               site_ids[which(n_parts < 2 | n_parts > 3)[1]])
  }
  chrom <- vapply(parts, `[`, character(1), 1)
  pos <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 2)))
  strand <- vapply(parts, function(p) if (length(p) == 3) p[3] else "*", character(1))
  if (any(!is.finite(pos))) {
    stop_input("site id position is not numeric, e.g. ", site_ids[which(!is.finite(pos))[1]])
  }
  data.table::data.table(site_id = site_ids, chrom = chrom, pos = pos, strand = strand)
}

# Extract one attribute (e.g. gene_id) from GTF attribute strings.
gtf_attribute <- function(attributes, key) {
  pattern <- paste0('.*', key, ' "([^"]*)".*')
  value <- ifelse(grepl(paste0(key, ' "'), attributes), sub(pattern, "\\1", attributes), NA_character_)
  value
}

read_gtf_genes <- function(path) {
  # stream only gene lines (fast on a full GENCODE GTF)
  command <- paste(stream_command(path), "| awk -F'\\t' '$0 !~ /^#/ && $3 == \"gene\"'")
  gtf <- data.table::fread(
    cmd = command,
    sep = "\t",
    header = FALSE,
    quote = "",
    col.names = c("seqname", "source", "feature", "start", "end", "score",
                  "strand", "frame", "attribute"),
    colClasses = list(character = c(1, 2, 3, 6, 7, 8, 9), integer = c(4, 5)),
    showProgress = FALSE
  )
  if (nrow(gtf) == 0) {
    stop_input("GTF has no 'gene' records: ", path)
  }
  genes <- data.table::data.table(
    chr = gtf$seqname,
    start = gtf$start,
    end = gtf$end,
    strand = gtf$strand,
    gene_id = gtf_attribute(gtf$attribute, "gene_id"),
    gene_name = gtf_attribute(gtf$attribute, "gene_name"),
    gene_type = gtf_attribute(gtf$attribute, "gene_type")
  )
  missing_type <- is.na(genes$gene_type)
  if (any(missing_type)) {
    genes$gene_type[missing_type] <- gtf_attribute(gtf$attribute[missing_type], "gene_biotype")
  }
  if (anyNA(genes$gene_id)) {
    stop_input("GTF gene records without gene_id attribute: ", path)
  }
  genes
}

# Map sites to host genes.
#   site_meta       : data.table with #chr, start, end, site_id (BED coordinates)
#   genes           : read_gtf_genes() output
#   gene_mean_expr  : named numeric vector (gene_id -> mean CPM); genes absent
#                     from the expression BED cannot be host genes
map_sites_to_host_genes <- function(site_meta, genes, gene_mean_expr) {
  parsed <- parse_site_ids(site_meta$site_id)
  sites <- data.table::data.table(
    site_id = site_meta$site_id,
    chr = site_meta[["#chr"]],
    start = site_meta$end,           # 1-based site position == BED end
    end = site_meta$end,
    site_strand = parsed$strand
  )
  # only genes that are present in the expression matrix can provide theta
  genes <- genes[gene_id %in% names(gene_mean_expr)]
  genes[, mean_expression := gene_mean_expr[gene_id]]

  data.table::setkey(genes, chr, start, end)
  overlaps <- data.table::foverlaps(sites, genes, by.x = c("chr", "start", "end"),
                                    type = "within", nomatch = NA)
  # overlaps columns: gene columns + i.start / i.end (site) + site columns
  overlaps[, same_strand := !is.na(gene_id) & site_strand %in% c("+", "-") & strand == site_strand]

  choose_gene <- function(dt) {
    candidates <- dt[!is.na(gene_id)]
    if (nrow(candidates) == 0) {
      return(data.table::data.table(host_gene_id = NA_character_, host_gene_name = NA_character_,
                                    n_candidate_genes = 0L, mapping_rule = "no_overlapping_expressed_gene"))
    }
    rule <- "overlap_any_strand"
    if (any(candidates$same_strand)) {
      candidates <- candidates[same_strand == TRUE]
      rule <- "overlap_same_strand"
    }
    best <- candidates[which.max(mean_expression)]
    if (nrow(candidates) > 1) {
      rule <- paste0(rule, "_highest_mean_expression")
    }
    data.table::data.table(host_gene_id = best$gene_id, host_gene_name = best$gene_name,
                           n_candidate_genes = nrow(candidates), mapping_rule = rule)
  }
  mapping <- overlaps[, choose_gene(.SD), by = site_id]
  mapping <- mapping[match(site_meta$site_id, mapping$site_id)]
  mapping
}

# Per-site coverage QC.
#   ratio, coverage : sites x samples matrices (same dimnames)
# A sample is "covered" at a site when coverage >= min_coverage.
site_coverage_qc <- function(ratio, coverage, min_coverage = 20, min_covered_samples = 60,
                             min_covered_fraction = 0, min_sd = 1e-6) {
  covered <- !is.na(coverage) & coverage >= min_coverage & !is.na(ratio)
  any_reads <- !is.na(coverage) & coverage > 0 & !is.na(ratio)
  n_samples <- ncol(ratio)

  n_covered <- rowSums(covered)
  ratio_covered <- ratio
  ratio_covered[!covered] <- NA_real_
  mean_ratio <- rowMeans(ratio_covered, na.rm = TRUE)
  sd_ratio <- apply(ratio_covered, 1, stats::sd, na.rm = TRUE)
  mean_coverage <- rowMeans(coverage, na.rm = TRUE)

  qc <- data.table::data.table(
    site_id = rownames(ratio),
    n_samples = n_samples,
    n_samples_with_reads = rowSums(any_reads),
    n_samples_covered = n_covered,
    fraction_covered = n_covered / n_samples,
    mean_coverage = mean_coverage,
    mean_ratio_covered = mean_ratio,
    sd_ratio_covered = sd_ratio
  )
  qc[, reason := data.table::fcase(
    n_samples_covered < min_covered_samples, "too_few_covered_samples",
    fraction_covered < min_covered_fraction, "too_low_covered_fraction",
    is.na(sd_ratio_covered) | sd_ratio_covered < min_sd, "no_variation",
    default = "pass"
  )]
  qc[, pass := reason == "pass"]
  qc
}

# Split sites (already sorted by chromosome and position) into chunks.
define_site_chunks <- function(site_meta, sites_per_chunk) {
  n_sites <- nrow(site_meta)
  chunk_index <- integer(n_sites)
  current_chunk <- 1L
  in_chunk <- 0L
  for (i in seq_len(n_sites)) {
    new_chromosome <- i > 1 && site_meta[["#chr"]][i] != site_meta[["#chr"]][i - 1]
    same_position <- i > 1 && !new_chromosome && site_meta$end[i] == site_meta$end[i - 1]
    if (i > 1 && (new_chromosome || (in_chunk >= sites_per_chunk && !same_position))) {
      current_chunk <- current_chunk + 1L
      in_chunk <- 0L
    }
    chunk_index[i] <- current_chunk
    in_chunk <- in_chunk + 1L
  }
  chunk_ids <- sprintf("chunk_%05d", chunk_index)
  site_chunks <- data.table::data.table(site_id = site_meta$site_id, chunk_id = chunk_ids)

  chunks <- data.table::data.table(
    chunk_id = chunk_ids,
    chr = site_meta[["#chr"]],
    site_start = site_meta$start,
    site_end = site_meta$end
  )[, .(chr = chr[1],
        first_site_bed_start = min(site_start),
        last_site_bed_end = max(site_end),
        n_sites = .N), by = chunk_id]
  # tabix regions are 1-based and inclusive
  chunks[, tabix_region := paste0(chr, ":", first_site_bed_start + 1, "-", last_site_bed_end)]
  list(chunks = chunks, site_chunks = site_chunks)
}
