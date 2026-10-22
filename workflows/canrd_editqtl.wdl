version 1.0

# ------------------------------------------------------------------------------
# workflows/canrd_editqtl.wdl
#
# Description: end-to-end cell-type edQTL discovery with caNRD-editQTL on the
#   bulk editing ratios (no per-sample deconvolution in the test).
#   Entry point 3 (Dockstore name canrd_editqtl). Takes the cohort outputs of
#   prepare_editing_phenotypes.wdl (same Terra data-table row).
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   scatter chunks: CanrdEditqtlChunk (cis dosages via plink2, scan + exact
#                   lead re-fit)
#   AggregateCanrdEditqtl: site-level and cell-type-level calls
# Inputs:  see the input block
# Outputs: see the output block
# ------------------------------------------------------------------------------

import "tasks/canrd_editqtl.wdl" as canrd_tasks

workflow CanrdEditqtl {
  input {
    Array[String] chunk_ids
    File chunks
    File site_chunks
    File editing_ratios
    File editing_ratios_index
    File editing_coverage
    File editing_coverage_index
    File site_theta
    File proportions
    File cohort_samples
    File covariates
    File? editing_pcs
    File plink_pgen
    File plink_pvar
    File plink_psam
    String output_prefix

    Int n_editing_pcs = 0
    Int cis_window = 1000000
    Float maf = 0.01
    Int min_minor_allele_samples = 10
    Int min_samples = 30
    Float min_mean_phi = 0.10
    Float theta_floor = 0.001
    Float min_coverage = 20
    Int min_covered_samples = 60
    Float min_covered_fraction = 0.0
    String refine = "lead"
    Boolean write_all_pairs = true
    Float fdr = 0.05
    Float celltype_alpha = 0.05
    Float? qvalue_lambda

    String docker_image = "ghcr.io/haimkru/caeditr_implementation:latest"
    Int preemptible_attempts = 2
    Int max_retries = 1
    Int chunk_cpu = 4
    String chunk_memory = "24 GB"
    Int chunk_disk_gb = 150
  }

  scatter (chunk_id in chunk_ids) {
    call canrd_tasks.CanrdEditqtlChunk {
      input:
        chunk_id = chunk_id,
        chunks = chunks,
        site_chunks = site_chunks,
        editing_ratios = editing_ratios,
        editing_ratios_index = editing_ratios_index,
        editing_coverage = editing_coverage,
        editing_coverage_index = editing_coverage_index,
        site_theta = site_theta,
        proportions = proportions,
        cohort_samples = cohort_samples,
        covariates = covariates,
        editing_pcs = editing_pcs,
        n_editing_pcs = n_editing_pcs,
        plink_pgen = plink_pgen,
        plink_pvar = plink_pvar,
        plink_psam = plink_psam,
        cis_window = cis_window,
        maf = maf,
        min_minor_allele_samples = min_minor_allele_samples,
        min_samples = min_samples,
        min_mean_phi = min_mean_phi,
        theta_floor = theta_floor,
        min_coverage = min_coverage,
        min_covered_samples = min_covered_samples,
        min_covered_fraction = min_covered_fraction,
        refine = refine,
        write_all_pairs = write_all_pairs,
        docker_image = docker_image,
        cpu = chunk_cpu,
        memory = chunk_memory,
        disk_gb = chunk_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  call canrd_tasks.AggregateCanrdEditqtl {
    input:
      site_summaries = CanrdEditqtlChunk.site_summary,
      lead_variants = CanrdEditqtlChunk.lead_variants,
      site_celltype_min_p = CanrdEditqtlChunk.site_celltype_min_p,
      output_prefix = output_prefix,
      fdr = fdr,
      celltype_alpha = celltype_alpha,
      qvalue_lambda = qvalue_lambda,
      docker_image = docker_image,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  output {
    Array[File] all_pairs = CanrdEditqtlChunk.all_pairs
    Array[File] lead_variants = CanrdEditqtlChunk.lead_variants
    Array[File] scan_info = CanrdEditqtlChunk.scan_info
    File site_level = AggregateCanrdEditqtl.site_level
    File lead_site_celltype = AggregateCanrdEditqtl.lead_site_celltype
    File celltype_scan_level = AggregateCanrdEditqtl.celltype_scan_level
    File summary = AggregateCanrdEditqtl.summary
  }
}
