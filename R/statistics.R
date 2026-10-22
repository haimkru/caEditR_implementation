################################################################################
# R/statistics.R
#
# Description:
#   Multiple-testing helpers shared by the aggregation scripts.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   storey_qvalue(p): Storey q-values with the qvalue package (lambda fixed
#     when given, e.g. 0.85 as Li et al. 2022); falls back to Benjamini-
#     Hochberg (reported in the "method" attribute) when qvalue is missing or
#     fails (e.g. too few tests to estimate pi0).
#   bonferroni_within(p, n): min(1, p * n), the per-site correction for the
#     number of cis variants tested at that site.
#
# Inputs:  numeric vectors from the calling script
# Outputs: none (returns numeric vectors)
################################################################################

storey_qvalue <- function(p, lambda = NULL) {
  q <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  method <- "none"
  if (sum(ok) == 0) {
    attr(q, "method") <- method
    return(q)
  }
  fitted <- NULL
  if (requireNamespace("qvalue", quietly = TRUE) && sum(ok) >= 10) {
    fitted <- tryCatch({
      if (is.null(lambda)) {
        qvalue::qvalue(p[ok])
      } else {
        qvalue::qvalue(p[ok], lambda = lambda)
      }
    }, error = function(e) NULL)
  }
  if (!is.null(fitted)) {
    q[ok] <- fitted$qvalues
    method <- "storey_qvalue"
  } else {
    q[ok] <- stats::p.adjust(p[ok], method = "BH")
    method <- "benjamini_hochberg"
  }
  attr(q, "method") <- method
  q
}

bonferroni_within <- function(p, n_tests) {
  pmin(1, p * pmax(n_tests, 1))
}
