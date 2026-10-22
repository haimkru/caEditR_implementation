################################################################################
# R/genotypes.R
#
# Description:
#   Genotype access for the caNRD-editQTL scan: extract ALT-allele dosages of
#   the cis variants of one chunk of editing sites from PLINK2 files (with
#   plink2, so the full 1M-variant x 10k-sample matrix is never loaded), and
#   build the (site, variant) cis pairs.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   extract_cis_dosages():
#     1. write the cohort sample IDs to a keep file;
#     2. plink2 --pfile <prefix> --keep <file> --chr <chr>
#               --from-bp <a> --to-bp <b> --maf <maf> --export A-transpose;
#     3. read the .traw (variants x samples), rename sample columns to IIDs,
#        make sure the counted allele is ALT (flip 2 - x otherwise), so betas
#        have the same sign convention as tensorQTL.
#   build_cis_pairs(): every variant within +/- window of every site.
#
# Inputs:  PLINK2 prefix, sample IDs, region, thresholds (from the script)
# Outputs: none (returns R objects; temporary plink2 files in a work dir)
################################################################################

extract_cis_dosages <- function(pfile_prefix, samples, chr, from_bp, to_bp, maf, work_dir) {
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  keep_file <- file.path(work_dir, "keep_samples.txt")
  writeLines(samples, keep_file)
  out_prefix <- file.path(work_dir, "cis_dosages")

  args <- c(
    "--pfile", shQuote(pfile_prefix),
    "--keep", shQuote(keep_file),
    "--chr", shQuote(chr),
    "--from-bp", format(max(1, floor(from_bp)), scientific = FALSE),
    "--to-bp", format(ceiling(to_bp), scientific = FALSE),
    "--maf", maf,
    "--export", "A-transpose",
    "--out", shQuote(out_prefix)
  )
  plink_output <- suppressWarnings(system2("plink2", args, stdout = TRUE, stderr = TRUE))
  traw_file <- paste0(out_prefix, ".traw")
  if (!file.exists(traw_file)) {
    # plink2 stops with an error when no variant is left in the region
    if (any(grepl("No variants remaining|no variants|Empty", plink_output, ignore.case = TRUE))) {
      log_info("No variants in ", chr, ":", from_bp, "-", to_bp, " after filters")
      empty <- matrix(numeric(0), 0, length(samples), dimnames = list(NULL, samples))
      return(list(dosages = empty,
                  variants = data.table::data.table(variant_id = character(0), chr = character(0),
                                                    pos = numeric(0))))
    }
    stop("plink2 failed:\n", paste(plink_output, collapse = "\n"), call. = FALSE)
  }

  traw <- data.table::fread(traw_file, sep = "\t", header = TRUE, check.names = FALSE,
                            na.strings = c("NA", "."), showProgress = FALSE)
  info_columns <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")
  sample_columns <- setdiff(names(traw), info_columns)
  # plink2 names sample columns FID_IID (or IID when there is no FID column)
  sample_ids <- sample_columns
  strip_prefix <- !all(sample_ids %in% samples)
  if (strip_prefix) {
    sample_ids <- sub("^[^_]*_", "", sample_ids)
  }
  if (!all(samples %in% sample_ids)) {
    stop("could not match plink2 .traw sample columns to cohort sample IDs", call. = FALSE)
  }

  dosages <- as.matrix(traw[, sample_columns, with = FALSE])
  storage.mode(dosages) <- "double"
  colnames(dosages) <- sample_ids
  dosages <- dosages[, samples, drop = FALSE]

  # pvar columns: COUNTED is the counted allele, ALT the other one; we want ALT counts
  # (plink2 writes the counted allele in COUNTED; by default that is REF)
  pvar <- data.table::fread(paste0(pfile_prefix, ".pvar"), sep = "\t", header = TRUE,
                            select = c("#CHROM", "ID", "REF", "ALT"), skip = "#CHROM",
                            colClasses = "character", showProgress = FALSE)
  pvar_row <- match(traw$SNP, pvar$ID)
  ref_allele <- pvar$REF[pvar_row]
  counted_is_ref <- !is.na(ref_allele) & traw$COUNTED == ref_allele
  if (any(counted_is_ref)) {
    dosages[counted_is_ref, ] <- 2 - dosages[counted_is_ref, ]
  }
  rownames(dosages) <- traw$SNP

  # chromosome names exactly as in the .pvar (the .traw drops the "chr" prefix)
  variants <- data.table::data.table(variant_id = traw$SNP, chr = pvar[["#CHROM"]][pvar_row],
                                     pos = as.numeric(traw$POS))
  unlink(list.files(work_dir, pattern = "^cis_dosages", full.names = TRUE))
  list(dosages = dosages, variants = variants)
}

# All (site, variant) pairs with |variant_pos - site_pos| <= window.
#   sites   : data.table(site_id, chr, pos)
#   variants: data.table(variant_id, chr, pos)
build_cis_pairs <- function(sites, variants, window) {
  if (nrow(sites) == 0 || nrow(variants) == 0) {
    return(data.table::data.table(site_id = character(0), variant_id = character(0),
                                  distance = numeric(0)))
  }
  pairs <- lapply(seq_len(nrow(sites)), function(i) {
    in_window <- variants$chr == sites$chr[i] & abs(variants$pos - sites$pos[i]) <= window
    if (!any(in_window)) {
      return(NULL)
    }
    data.table::data.table(site_id = sites$site_id[i],
                           variant_id = variants$variant_id[in_window],
                           distance = variants$pos[in_window] - sites$pos[i])
  })
  data.table::rbindlist(pairs)
}
