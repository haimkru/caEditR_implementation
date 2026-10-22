################################################################################
# R/phenotypes.R
#
# Description:
#   Turns an editing matrix (bulk ratios, or a per-cell-type matrix from
#   caNRD-edit / TCA) into a tensorQTL-ready phenotype: missing-value filter,
#   mean imputation, rank-based inverse normal transform (INT) per site, and
#   phenotype principal components used as covariates.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   inverse_normal_transform(x): Blom INT, qnorm((rank - 3/8) / (n + 1/4)),
#     ties averaged, NA kept (same formula as RNOmni::RankNorm used by
#     prepare_QTL).
#   impute_row_means(m): replace NA by the row (site) mean.
#   compute_phenotype_pcs(m, k): PCA of samples on the INT matrix; uses
#     irlba (truncated SVD) when installed, otherwise the eigen-decomposition
#     of the samples x samples cross-product (exact, fine up to ~20k samples).
#
# Inputs:  matrices given by the calling script
# Outputs: none (returns R objects)
################################################################################

inverse_normal_transform <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  n <- sum(ok)
  if (n == 0) {
    return(out)
  }
  ranks <- rank(x[ok], ties.method = "average")
  out[ok] <- stats::qnorm((ranks - 0.375) / (n + 0.25))
  out
}

# Apply the INT to every row of a sites x samples matrix.
inverse_normal_transform_rows <- function(m) {
  out <- t(apply(m, 1, inverse_normal_transform))
  if (nrow(m) == 1) {
    out <- matrix(out, nrow = 1)
  }
  dimnames(out) <- dimnames(m)
  out
}

impute_row_means <- function(m) {
  row_means <- rowMeans(m, na.rm = TRUE)
  missing <- which(is.na(m), arr.ind = TRUE)
  if (nrow(missing) > 0) {
    m[missing] <- row_means[missing[, 1]]
  }
  m
}

# Principal components of the samples.
#   m : sites x samples matrix without NA (rows already centred/INT)
#   k : number of PCs requested
# Returns list(scores = samples x k matrix, variance_explained = numeric(k))
compute_phenotype_pcs <- function(m, k) {
  n_samples <- ncol(m)
  n_sites <- nrow(m)
  k <- min(k, n_samples - 1L, n_sites - 1L)
  if (k < 1) {
    return(list(scores = matrix(numeric(0), n_samples, 0, dimnames = list(colnames(m), NULL)),
                variance_explained = numeric(0), method = "none"))
  }
  # samples x sites, centred per site
  x <- t(m)
  x <- sweep(x, 2, colMeans(x), `-`)
  total_variance <- sum(x^2)

  use_irlba <- requireNamespace("irlba", quietly = TRUE) && k < min(dim(x)) / 2
  if (use_irlba) {
    sv <- irlba::irlba(x, nv = k, nu = k)
    scores <- sv$u %*% diag(sv$d, nrow = k)
    variance <- sv$d^2
    method <- "irlba_truncated_svd"
  } else {
    gram <- tcrossprod(x)                      # samples x samples
    eig <- eigen(gram, symmetric = TRUE)
    values <- pmax(eig$values[seq_len(k)], 0)
    scores <- eig$vectors[, seq_len(k), drop = FALSE] %*% diag(sqrt(values), nrow = k)
    variance <- values
    method <- "eigen_gram_matrix"
  }
  colnames(scores) <- paste0("editing_PC", seq_len(k))
  rownames(scores) <- colnames(m)
  list(scores = scores, variance_explained = variance / total_variance, method = method)
}

# Principal components from an accumulated samples x samples Gram matrix
# (sum over site blocks of crossprod(site-centred block)); identical scores to
# compute_phenotype_pcs() on the full matrix, without holding it.
#   gram           : samples x samples matrix
#   k              : number of PCs requested
#   total_variance : sum of squares of the centred values (for variance explained)
#   n_sites        : number of sites that went into gram
pcs_from_gram <- function(gram, k, total_variance, n_sites) {
  n_samples <- nrow(gram)
  k <- min(k, n_samples - 1L, n_sites - 1L)
  if (k < 1) {
    return(list(scores = matrix(numeric(0), n_samples, 0, dimnames = list(rownames(gram), NULL)),
                variance_explained = numeric(0), method = "none"))
  }
  if (requireNamespace("irlba", quietly = TRUE) && k < n_samples / 4) {
    eig <- irlba::partial_eigen(gram, n = k, symmetric = TRUE)
    values <- pmax(eig$values[seq_len(k)], 0)
    vectors <- eig$vectors[, seq_len(k), drop = FALSE]
    method <- "irlba_partial_eigen_gram"
  } else {
    eig <- eigen(gram, symmetric = TRUE)
    values <- pmax(eig$values[seq_len(k)], 0)
    vectors <- eig$vectors[, seq_len(k), drop = FALSE]
    method <- "eigen_gram"
  }
  scores <- vectors %*% diag(sqrt(values), nrow = k)
  dimnames(scores) <- list(rownames(gram), paste0("editing_PC", seq_len(k)))
  list(scores = scores, variance_explained = values / total_variance, method = method)
}
