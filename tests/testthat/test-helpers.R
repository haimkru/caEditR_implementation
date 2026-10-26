################################################################################
# tests/testthat/test-helpers.R
#
# Description: unit tests of the shared helpers in R/ (no external tools).
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
# Pseudo code: one test_that() block per helper, small hand-made inputs with
#              known answers.
# Inputs / outputs: none (run with Rscript tests/testthat.R)
################################################################################

test_that("slugify_celltype matches CellTypeDeconvolution's rule", {
  expect_equal(slugify_celltype(c("Monocyte/myeloid", "CD4 T cells", "Gamma-delta T cells")),
               c("monocyte_myeloid", "cd4_t_cells", "gamma_delta_t_cells"))
})

test_that("parse_site_ids reads caEditR site ids", {
  parsed <- parse_site_ids(c("1:12345:+", "X:10:-", "2:99"))
  expect_equal(parsed$chrom, c("1", "X", "2"))
  expect_equal(parsed$pos, c(12345, 10, 99))
  expect_equal(parsed$strand, c("+", "-", "*"))
  expect_error(parse_site_ids("chr1-100"), "Input error")
})

test_that("inverse_normal_transform is the Blom rank transform and keeps NA", {
  x <- c(3, 1, NA, 2)
  out <- inverse_normal_transform(x)
  expect_true(is.na(out[3]))
  expect_equal(out[c(2, 4, 1)], stats::qnorm((1:3 - 0.375) / (3 + 0.25)))
})

test_that("impute_row_means fills NA with the site mean", {
  m <- matrix(c(1, NA, 3, 4, 5, NA), nrow = 2, byrow = TRUE)
  out <- impute_row_means(m)
  expect_equal(out[1, 2], 2)
  expect_equal(out[2, 3], 4.5)
})

test_that("site_coverage_qc applies coverage and variation rules", {
  ratio <- rbind(a = c(0.1, 0.2, 0.3, 0.4), b = c(0.1, 0.1, 0.1, 0.1), c = c(0.1, NA, 0.2, 0.3))
  coverage <- rbind(a = c(30, 30, 30, 30), b = c(30, 30, 30, 30), c = c(30, 0, 5, 30))
  qc <- site_coverage_qc(ratio, coverage, min_coverage = 20, min_covered_samples = 3)
  expect_equal(qc$reason, c("pass", "no_variation", "too_few_covered_samples"))
})

test_that("define_site_chunks never splits a chromosome or a position", {
  meta <- data.table::data.table(`#chr` = c("chr1", "chr1", "chr1", "chr2"),
                                 start = c(9L, 9L, 19L, 4L), end = c(10L, 10L, 20L, 5L),
                                 site_id = c("1:10:+", "1:10:-", "1:20:+", "2:5:+"))
  chunks <- define_site_chunks(meta, sites_per_chunk = 1)
  expect_equal(chunks$site_chunks$chunk_id, c("chunk_00001", "chunk_00001", "chunk_00002", "chunk_00003"))
  expect_equal(chunks$chunks$tabix_region, c("chr1:10-10", "chr1:20-20", "chr2:5-5"))
})

test_that("compute_phenotype_pcs returns one score per sample", {
  set.seed(1)
  m <- matrix(stats::rnorm(200), nrow = 20, dimnames = list(NULL, paste0("s", 1:10)))
  pcs <- compute_phenotype_pcs(m, 3)
  expect_equal(dim(pcs$scores), c(10, 3))
  expect_equal(rownames(pcs$scores), paste0("s", 1:10))
  expect_true(all(diff(pcs$variance_explained) <= 1e-12))
})

test_that("storey_qvalue falls back to BH for few tests", {
  q <- storey_qvalue(c(0.01, 0.2, NA))
  expect_equal(attr(q, "method"), "benjamini_hochberg")
  expect_true(is.na(q[3]))
  expect_equal(bonferroni_within(c(0.01, 0.5), c(10, 5)), c(0.1, 1))
})

test_that("build_cis_pairs keeps variants inside the window on the same chromosome", {
  sites <- data.table::data.table(site_id = "s1", chr = "chr1", pos = 1000)
  variants <- data.table::data.table(variant_id = c("v1", "v2", "v3"), chr = c("chr1", "chr1", "chr2"),
                                     pos = c(1500, 5000, 1000))
  pairs <- build_cis_pairs(sites, variants, window = 1000)
  expect_equal(pairs$variant_id, "v1")
  expect_equal(pairs$distance, 500)
})

test_that("pcs_from_gram reproduces compute_phenotype_pcs on block-accumulated data", {
  set.seed(2)
  m <- matrix(stats::rnorm(400), nrow = 40, dimnames = list(NULL, paste0("s", 1:10)))
  centred <- m - rowMeans(m)
  gram <- crossprod(centred[1:20, ]) + crossprod(centred[21:40, ])
  from_gram <- pcs_from_gram(gram, 3, sum(centred^2), 40)
  direct <- compute_phenotype_pcs(m, 3)
  expect_equal(abs(unname(from_gram$scores)), abs(unname(direct$scores)), tolerance = 1e-8)
  expect_equal(from_gram$variance_explained, direct$variance_explained, tolerance = 1e-8)
})

test_that("sample tables accept research_id anywhere and sample lists can be tables", {
  table_file <- tempfile(fileext = ".tsv")
  writeLines(c("PC1\tresearch_id\tPC2", "0.1\t0001\t0.2", "0.3\t0002\t0.4"), table_file)
  dt <- read_sample_table(table_file, "sample_id")
  expect_equal(names(dt)[1], "sample_id")
  expect_equal(dt$sample_id, c("0001", "0002"))
  list_file <- tempfile(fileext = ".tsv")
  writeLines(c("research_id\tancestry", "1000001\tEUR", "1000002\tAFR"), list_file)
  expect_equal(read_sample_list(list_file), c("1000001", "1000002"))
  plain_file <- tempfile(fileext = ".txt")
  writeLines(c("research_id", "1000001", "1000002"), plain_file)
  expect_equal(read_sample_list(plain_file), c("1000001", "1000002"))
})
