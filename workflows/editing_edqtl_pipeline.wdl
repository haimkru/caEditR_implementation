version 1.0

# ------------------------------------------------------------------------------
# workflows/editing_edqtl_pipeline.wdl
#
# Description: END-TO-END cell-type-specific RNA-editing QTL pipeline for All
#   of Us (Dockstore name editing_edqtl_pipeline). Runs the three approaches
#   of the caEditR draft (Methods 4.6) on one cohort:
#     1. bulk edQTLs with tensorQTL (standard edQTL pipeline);
#     2. tensorQTL on caNRD-edit and TCA cell-type editing matrices;
#     3. end-to-end caNRD-editQTL on the bulk ratios;
#   and summarizes them side by side.
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   PrepareEditingPhenotypes  (index, cohort, deconvolution, phenotype sets)
#   TensorqtlEdqtl            (if run_tensorqtl)    approaches 1 + 2
#   CanrdEditqtl              (if run_canrd_editqtl) approach 3
#   SummarizeEdqtlResults     comparison table + counts
# Inputs:  see the input block; formats in docs/data-dictionary.md
# Outputs: see the output block; finalized tables are the *_summary,
#          *_site_level, *_lead_site_celltype and edsite_* files
# ------------------------------------------------------------------------------

import "prepare_editing_phenotypes.wdl" as prepare_workflow
import "tensorqtl_edqtl.wdl" as tensorqtl_workflow
import "canrd_editqtl.wdl" as canrd_workflow
import "tasks/summary.wdl" as summary_tasks

workflow EditingEdqtlPipeline {
  input {
    # ---- data -------------------------------------------------------------
    File editing_ratios
    File editing_coverage
    File expression
    File gtf
    File proportions
    File covariates
    File plink_pgen
    File plink_pvar
    File plink_psam
    File? sample_list
    String output_prefix

    # ---- what to run ------------------------------------------------------------
    Boolean run_tensorqtl = true
    Boolean run_tensorqtl_nominal = true
    Boolean run_tensorqtl_permutations = false
    Boolean run_canrd_editqtl = true

    # ---- shared analysis settings -----------------------------------------------
    Int sites_per_chunk = 500
    Float theta_floor = 0.001
    Float min_coverage = 20
    Int min_covered_samples = 60
    Float min_covered_fraction = 0.0
    Float canrd_min_coverage = 10
    Array[String] deconvolution_methods = ["canrd_edit", "tca"]
    Int n_editing_pcs = 15
    Float max_missing_fraction = 0.5
    Int cis_window = 1000000
    Float fdr = 0.05
    Float? qvalue_lambda
    Int? tensorqtl_seed

    # ---- caNRD-editQTL settings ---------------------------------------------------
    Int canrd_n_editing_pcs = 0
    Float canrd_maf = 0.01
    Int canrd_min_minor_allele_samples = 10
    Float canrd_min_mean_phi = 0.10
    String canrd_refine = "lead"
    Boolean canrd_write_all_pairs = true
    Float canrd_celltype_alpha = 0.05

    # ---- runtime ------------------------------------------------------------------
    String docker_image = "ghcr.io/haimkru/caeditr_implementation:latest"
    Int preemptible_attempts = 2
    Int max_retries = 1
    Int deconvolve_cpu = 1
    String deconvolve_memory = "6 GB"
    Int deconvolve_disk_gb = 100
    Int phenotype_cpu = 2
    String phenotype_memory = "8 GB"
    Int phenotype_disk_gb = 100
    Int canrd_chunk_cpu = 4
    String canrd_chunk_memory = "24 GB"
    Int canrd_chunk_disk_gb = 150
    Int tensorqtl_memory = 64
    Int tensorqtl_disk_space = 200
    Int tensorqtl_num_threads = 16
    Int tensorqtl_num_gpus = 1
    Int tensorqtl_num_preempt = 1
    String? tensorqtl_machine_type = "g2-standard-16"
    String tensorqtl_gpu_type = "nvidia-l4"
    Array[String] tensorqtl_zones = ["us-central1-c"]
    String tensorqtl_docker_image = "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"
  }

  call prepare_workflow.PrepareEditingPhenotypes as Prepare {
    input:
      editing_ratios = editing_ratios,
      editing_coverage = editing_coverage,
      expression = expression,
      gtf = gtf,
      proportions = proportions,
      covariates = covariates,
      plink_psam = plink_psam,
      sample_list = sample_list,
      output_prefix = output_prefix,
      sites_per_chunk = sites_per_chunk,
      theta_floor = theta_floor,
      min_coverage = min_coverage,
      min_covered_samples = min_covered_samples,
      min_covered_fraction = min_covered_fraction,
      canrd_min_coverage = canrd_min_coverage,
      deconvolution_methods = deconvolution_methods,
      n_editing_pcs = n_editing_pcs,
      max_missing_fraction = max_missing_fraction,
      docker_image = docker_image,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries,
      deconvolve_cpu = deconvolve_cpu,
      deconvolve_memory = deconvolve_memory,
      deconvolve_disk_gb = deconvolve_disk_gb,
      phenotype_cpu = phenotype_cpu,
      phenotype_memory = phenotype_memory,
      phenotype_disk_gb = phenotype_disk_gb
  }

  if (run_tensorqtl) {
    call tensorqtl_workflow.TensorqtlEdqtl as Tensorqtl {
      input:
        phenotype_labels = Prepare.phenotype_labels,
        phenotype_n_sites = Prepare.phenotype_n_sites,
        phenotype_int_beds = Prepare.phenotype_int_beds,
        phenotype_covariates = Prepare.phenotype_covariates,
        plink_pgen = plink_pgen,
        plink_pvar = plink_pvar,
        plink_psam = plink_psam,
        output_prefix = output_prefix,
        run_permutations = run_tensorqtl_permutations,
        run_nominal = run_tensorqtl_nominal,
        cis_window = cis_window,
        fdr = fdr,
        qvalue_lambda = qvalue_lambda,
        seed = tensorqtl_seed,
        tensorqtl_memory = tensorqtl_memory,
        tensorqtl_disk_space = tensorqtl_disk_space,
        tensorqtl_num_threads = tensorqtl_num_threads,
        tensorqtl_num_gpus = tensorqtl_num_gpus,
        tensorqtl_num_preempt = tensorqtl_num_preempt,
        tensorqtl_machine_type = tensorqtl_machine_type,
        tensorqtl_gpu_type = tensorqtl_gpu_type,
        tensorqtl_zones = tensorqtl_zones,
        tensorqtl_docker_image = tensorqtl_docker_image,
        docker_image = docker_image,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  if (run_canrd_editqtl) {
    call canrd_workflow.CanrdEditqtl as Canrd {
      input:
        chunk_ids = Prepare.chunk_ids,
        chunks = Prepare.chunks,
        site_chunks = Prepare.site_chunks,
        editing_ratios = Prepare.indexed_editing_ratios,
        editing_ratios_index = Prepare.indexed_editing_ratios_index,
        editing_coverage = Prepare.indexed_editing_coverage,
        editing_coverage_index = Prepare.indexed_editing_coverage_index,
        site_theta = Prepare.site_theta,
        proportions = Prepare.proportions_aligned,
        cohort_samples = Prepare.cohort_samples,
        covariates = Prepare.covariates_aligned,
        editing_pcs = Prepare.bulk_editing_pcs,
        n_editing_pcs = canrd_n_editing_pcs,
        plink_pgen = plink_pgen,
        plink_pvar = plink_pvar,
        plink_psam = plink_psam,
        output_prefix = output_prefix,
        cis_window = cis_window,
        maf = canrd_maf,
        min_minor_allele_samples = canrd_min_minor_allele_samples,
        min_mean_phi = canrd_min_mean_phi,
        theta_floor = theta_floor,
        min_coverage = min_coverage,
        min_covered_samples = min_covered_samples,
        min_covered_fraction = min_covered_fraction,
        refine = canrd_refine,
        write_all_pairs = canrd_write_all_pairs,
        fdr = fdr,
        celltype_alpha = canrd_celltype_alpha,
        qvalue_lambda = qvalue_lambda,
        docker_image = docker_image,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries,
        chunk_cpu = canrd_chunk_cpu,
        chunk_memory = canrd_chunk_memory,
        chunk_disk_gb = canrd_chunk_disk_gb
    }
  }

  call summary_tasks.SummarizeEdqtlResults {
    input:
      tensorqtl_results = Tensorqtl.tensorqtl_all_analyses,
      canrd_lead_site_celltype = Canrd.lead_site_celltype,
      celltypes = Prepare.celltypes,
      output_prefix = output_prefix,
      docker_image = docker_image,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  output {
    # cohort and QC
    File cohort_samples = Prepare.cohort_samples
    File sample_alignment_report = Prepare.sample_alignment_report
    File cohort_summary = Prepare.cohort_summary
    File site_to_gene = Prepare.site_to_gene
    File site_theta = Prepare.site_theta
    File site_qc = Prepare.site_qc
    File canrd_diagnostics = Prepare.canrd_diagnostics
    File deconvolution_summary = Prepare.deconvolution_summary

    # cell-type editing matrices and tensorQTL inputs
    File bulk_bed = Prepare.bulk_bed
    Array[File] celltype_beds = Prepare.celltype_beds
    File phenotype_inventory = Prepare.phenotype_inventory
    Array[String] phenotype_labels = Prepare.phenotype_labels
    Array[File] phenotype_int_beds = Prepare.phenotype_int_beds
    Array[File] phenotype_covariates = Prepare.phenotype_covariates

    # approaches 1 + 2: tensorQTL
    Array[File]? tensorqtl_cis_permutation_results = Tensorqtl.cis_permutation_results
    Array[Array[File]]? tensorqtl_cis_nominal_parquet = Tensorqtl.cis_nominal_parquet
    File? tensorqtl_all_analyses = Tensorqtl.tensorqtl_all_analyses
    File? tensorqtl_summary = Tensorqtl.tensorqtl_summary

    # approach 3: caNRD-editQTL
    Array[File]? canrd_editqtl_all_pairs = Canrd.all_pairs
    File? canrd_editqtl_site_level = Canrd.site_level
    File? canrd_editqtl_lead_site_celltype = Canrd.lead_site_celltype
    File? canrd_editqtl_celltype_scan_level = Canrd.celltype_scan_level
    File? canrd_editqtl_summary = Canrd.summary

    # comparison
    File edsite_calls_long = SummarizeEdqtlResults.edsite_calls_long
    File edsite_overview = SummarizeEdqtlResults.edsite_overview
  }
}
