version 1.0

# ------------------------------------------------------------------------------
# workflows/prepare_editing_phenotypes.wdl
#
# Description: from bulk RNA-editing ratios + coverage to tensorQTL-ready
#   phenotype sets: bulk, caNRD-edit x cell type and TCA x cell type.
#   Entry point 1 of caEditR_implementation (Dockstore name prepare_editing_phenotypes).
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   1. IndexBedMatrix x2   sort + bgzip + tabix ratios and coverage
#   2. PrepareCohort       cohort samples, site -> host gene, theta, chunks
#   3. scatter chunks:     DeconvolveChunk (site QC, bulk, caNRD-edit, TCA)
#   4. MergeChunkOutputs   genome-wide BED per (method, cell type)
#   5. PreparePhenotypes   bulk + every cell-type BED -> INT BED + covariates
# Inputs:  see the input block (file formats: docs/data-dictionary.md)
# Outputs: see the output block
# ------------------------------------------------------------------------------

import "tasks/cohort.wdl" as cohort_tasks
import "tasks/deconvolution.wdl" as deconvolution_tasks
import "tasks/phenotypes.wdl" as phenotype_tasks

workflow PrepareEditingPhenotypes {
  input {
    # ---- data -------------------------------------------------------------
    File editing_ratios
    File editing_coverage
    File expression
    File gtf
    File proportions
    File covariates
    File plink_psam
    File? sample_list
    String output_prefix

    # ---- analysis settings ----------------------------------------------------
    Int sites_per_chunk = 500
    Float theta_floor = 0.001
    Float min_coverage = 20
    Int min_covered_samples = 60
    Float min_covered_fraction = 0.0
    Float canrd_min_coverage = 10
    Array[String] deconvolution_methods = ["canrd_edit", "tca"]
    Int n_editing_pcs = 15
    Float max_missing_fraction = 0.5

    # ---- runtime ------------------------------------------------------------------
    String docker_image = "ghcr.io/haimkru/caeditr_implementation:latest"
    Int preemptible_attempts = 2
    Int max_retries = 1
    String index_memory = "16 GB"
    Int index_disk_gb = 200
    String prepare_cohort_memory = "32 GB"
    Int prepare_cohort_disk_gb = 200
    Int deconvolve_cpu = 1
    String deconvolve_memory = "6 GB"
    Int deconvolve_disk_gb = 100
    Int merge_disk_gb = 250
    Int phenotype_cpu = 2
    String phenotype_memory = "8 GB"
    Int phenotype_disk_gb = 100
  }

  call cohort_tasks.IndexBedMatrix as IndexEditingRatios {
    input:
      bed = editing_ratios,
      output_name = output_prefix + ".editing_ratios.bed.gz",
      docker_image = docker_image,
      memory = index_memory,
      disk_gb = index_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  call cohort_tasks.IndexBedMatrix as IndexEditingCoverage {
    input:
      bed = editing_coverage,
      output_name = output_prefix + ".editing_coverage.bed.gz",
      docker_image = docker_image,
      memory = index_memory,
      disk_gb = index_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  call cohort_tasks.PrepareCohort {
    input:
      editing_ratios = IndexEditingRatios.indexed_bed,
      editing_coverage = IndexEditingCoverage.indexed_bed,
      expression = expression,
      gtf = gtf,
      proportions = proportions,
      plink_psam = plink_psam,
      covariates = covariates,
      sample_list = sample_list,
      sites_per_chunk = sites_per_chunk,
      theta_floor = theta_floor,
      docker_image = docker_image,
      memory = prepare_cohort_memory,
      disk_gb = prepare_cohort_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  scatter (chunk_id in PrepareCohort.chunk_ids) {
    call deconvolution_tasks.DeconvolveChunk {
      input:
        chunk_id = chunk_id,
        chunks = PrepareCohort.chunks,
        site_chunks = PrepareCohort.site_chunks,
        editing_ratios = IndexEditingRatios.indexed_bed,
        editing_ratios_index = IndexEditingRatios.indexed_bed_index,
        editing_coverage = IndexEditingCoverage.indexed_bed,
        editing_coverage_index = IndexEditingCoverage.indexed_bed_index,
        site_theta = PrepareCohort.site_theta,
        proportions = PrepareCohort.proportions_aligned,
        cohort_samples = PrepareCohort.cohort_samples,
        min_coverage = min_coverage,
        min_covered_samples = min_covered_samples,
        min_covered_fraction = min_covered_fraction,
        theta_floor = theta_floor,
        canrd_min_coverage = canrd_min_coverage,
        methods = deconvolution_methods,
        docker_image = docker_image,
        cpu = deconvolve_cpu,
        memory = deconvolve_memory,
        disk_gb = deconvolve_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  call deconvolution_tasks.MergeChunkOutputs {
    input:
      chunk_beds = flatten(DeconvolveChunk.chunk_beds),
      chunk_reports = flatten(DeconvolveChunk.chunk_reports),
      output_prefix = output_prefix,
      docker_image = docker_image,
      disk_gb = merge_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  call phenotype_tasks.PreparePhenotypes as PrepareBulkPhenotypes {
    input:
      input_bed = MergeChunkOutputs.bulk_bed,
      covariates = PrepareCohort.covariates_aligned,
      n_pcs = n_editing_pcs,
      max_missing_fraction = max_missing_fraction,
      docker_image = docker_image,
      cpu = phenotype_cpu,
      memory = phenotype_memory,
      disk_gb = phenotype_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  scatter (celltype_bed in MergeChunkOutputs.celltype_beds) {
    call phenotype_tasks.PreparePhenotypes as PrepareCelltypePhenotypes {
      input:
        input_bed = celltype_bed,
        covariates = PrepareCohort.covariates_aligned,
        n_pcs = n_editing_pcs,
        max_missing_fraction = max_missing_fraction,
        docker_image = docker_image,
        cpu = phenotype_cpu,
        memory = phenotype_memory,
        disk_gb = phenotype_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  output {
    # cohort definition (inputs of the caNRD-editQTL entry point)
    File cohort_samples = PrepareCohort.cohort_samples
    File sample_alignment_report = PrepareCohort.sample_alignment_report
    File celltypes = PrepareCohort.celltypes
    File proportions_aligned = PrepareCohort.proportions_aligned
    File covariates_aligned = PrepareCohort.covariates_aligned
    File site_to_gene = PrepareCohort.site_to_gene
    File gene_theta = PrepareCohort.gene_theta
    File site_theta = PrepareCohort.site_theta
    File chunks = PrepareCohort.chunks
    File site_chunks = PrepareCohort.site_chunks
    Array[String] chunk_ids = PrepareCohort.chunk_ids
    File cohort_summary = PrepareCohort.cohort_summary
    File indexed_editing_ratios = IndexEditingRatios.indexed_bed
    File indexed_editing_ratios_index = IndexEditingRatios.indexed_bed_index
    File indexed_editing_coverage = IndexEditingCoverage.indexed_bed
    File indexed_editing_coverage_index = IndexEditingCoverage.indexed_bed_index

    # genome-wide matrices (linear editing-ratio scale)
    File bulk_bed = MergeChunkOutputs.bulk_bed
    Array[File] celltype_beds = MergeChunkOutputs.celltype_beds
    File phenotype_inventory = MergeChunkOutputs.phenotype_inventory
    File site_qc = MergeChunkOutputs.site_qc
    File canrd_diagnostics = MergeChunkOutputs.canrd_diagnostics
    File deconvolution_summary = MergeChunkOutputs.deconvolution_summary

    # tensorQTL inputs: element 1 = bulk, then one per method x cell type
    Array[String] phenotype_labels = flatten([[PrepareBulkPhenotypes.phenotype_label], PrepareCelltypePhenotypes.phenotype_label])
    Array[Int] phenotype_n_sites = flatten([[PrepareBulkPhenotypes.n_sites_kept], PrepareCelltypePhenotypes.n_sites_kept])
    Array[File] phenotype_int_beds = flatten([[PrepareBulkPhenotypes.int_bed], PrepareCelltypePhenotypes.int_bed])
    Array[File] phenotype_covariates = flatten([[PrepareBulkPhenotypes.tensorqtl_covariates], PrepareCelltypePhenotypes.tensorqtl_covariates])
    Array[File] phenotype_filter_reports = flatten([[PrepareBulkPhenotypes.phenotype_filter_report], PrepareCelltypePhenotypes.phenotype_filter_report])
    Array[File] phenotype_summaries = flatten([[PrepareBulkPhenotypes.phenotype_summary], PrepareCelltypePhenotypes.phenotype_summary])
    File bulk_editing_pcs = PrepareBulkPhenotypes.phenotype_pcs
    File bulk_pc_variance_explained = PrepareBulkPhenotypes.pc_variance_explained
  }
}
