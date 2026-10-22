version 1.0

# ------------------------------------------------------------------------------
# workflows/tensorqtl_edqtl.wdl
#
# Description: "standard" edQTL mapping with tensorQTL on every phenotype set
#   from prepare_editing_phenotypes.wdl (bulk; caNRD-edit and TCA cell-type
#   matrices). Entry point 2 (Dockstore name tensorqtl_edqtl).
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   scatter phenotype sets (skip sets with 0 sites):
#     tensorqtl_cis_nominal      (default; all pairs, parquet per chromosome)
#     tensorqtl_cis_permutations (optional; lead variant per site, beta-approx p)
#   AggregateTensorqtl: lead variant per site and set, raw p-values, counts
# Inputs:  arrays from PrepareEditingPhenotypes (same order), PLINK2 genotypes
# Outputs: see the output block
# ------------------------------------------------------------------------------

import "tasks/tensorqtl.wdl" as tensorqtl_tasks

workflow TensorqtlEdqtl {
  input {
    Array[String] phenotype_labels
    Array[Int] phenotype_n_sites
    Array[File] phenotype_int_beds
    Array[File] phenotype_covariates
    File plink_pgen
    File plink_pvar
    File plink_psam
    String output_prefix

    # nominal = every site x variant pair with raw p-values (default; the user
    # does the multiple-testing correction). Permutations only add tensorQTL's
    # empirical per-site p-values and cost much more GPU time.
    Boolean run_nominal = true
    Boolean run_permutations = false
    Int cis_window = 1000000
    Float fdr = 0.05
    Float? qvalue_lambda
    # lambda for tensorQTL's OWN q-value step only: 0 (pi0 = 1, BH-like) never fails, whereas
    # qvalue's default smoothing errors on small phenotype sets (e.g. a rare cell type with
    # few sites). The reported q-values are recomputed in AggregateTensorqtl (qvalue_lambda).
    Float tensorqtl_internal_qvalue_lambda = 0.0
    Int? seed

    # tensorQTL runtime (same defaults as AoU-Multiomics-Analysis/tensorQTL)
    Int tensorqtl_memory = 64
    Int tensorqtl_disk_space = 200
    Int tensorqtl_num_threads = 16
    Int tensorqtl_num_gpus = 1
    Int tensorqtl_num_preempt = 1
    String? tensorqtl_machine_type = "g2-standard-16"
    String tensorqtl_gpu_type = "nvidia-l4"
    Array[String] tensorqtl_zones = ["us-central1-c"]
    String tensorqtl_docker_image = "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"

    String docker_image = "ghcr.io/haimkru/caeditr_implementation:latest"
    Int preemptible_attempts = 2
    Int max_retries = 1
  }

  scatter (i in range(length(phenotype_int_beds))) {
    Boolean has_sites = phenotype_n_sites[i] > 0

    if (run_permutations && has_sites) {
      call tensorqtl_tasks.tensorqtl_cis_permutations {
        input:
          plink_pgen = plink_pgen,
          plink_pvar = plink_pvar,
          plink_psam = plink_psam,
          phenotype_bed = phenotype_int_beds[i],
          covariates = phenotype_covariates[i],
          prefix = phenotype_labels[i],
          fdr = fdr,
          qvalue_lambda = tensorqtl_internal_qvalue_lambda,
          seed = seed,
          cis_window = cis_window,
          memory = tensorqtl_memory,
          disk_space = tensorqtl_disk_space,
          num_threads = tensorqtl_num_threads,
          num_gpus = tensorqtl_num_gpus,
          num_preempt = tensorqtl_num_preempt,
          machine_type = tensorqtl_machine_type,
          gpu_type = tensorqtl_gpu_type,
          zones = tensorqtl_zones,
          docker_image = tensorqtl_docker_image
      }
    }

    if (run_nominal && has_sites) {
      call tensorqtl_tasks.tensorqtl_cis_nominal {
        input:
          plink_pgen = plink_pgen,
          plink_pvar = plink_pvar,
          plink_psam = plink_psam,
          phenotype_bed = phenotype_int_beds[i],
          covariates = phenotype_covariates[i],
          prefix = phenotype_labels[i],
          cis_window = cis_window,
          memory = tensorqtl_memory,
          disk_space = tensorqtl_disk_space,
          num_threads = tensorqtl_num_threads,
          num_gpus = tensorqtl_num_gpus,
          num_preempt = tensorqtl_num_preempt,
          machine_type = tensorqtl_machine_type,
          gpu_type = tensorqtl_gpu_type,
          zones = tensorqtl_zones,
          docker_image = tensorqtl_docker_image
      }
    }
  }

  Array[File] cis_qtl_results = select_all(tensorqtl_cis_permutations.cis_qtl)
  Array[File] cis_nominal_files = flatten(select_all(tensorqtl_cis_nominal.chr_parquet))

  if (run_permutations || run_nominal) {
    call tensorqtl_tasks.AggregateTensorqtl {
      input:
        cis_nominal_parquet = cis_nominal_files,
        cis_qtl_results = cis_qtl_results,
        output_prefix = output_prefix,
        fdr = fdr,
        qvalue_lambda = qvalue_lambda,
        docker_image = docker_image,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  output {
    Array[File] cis_permutation_results = cis_qtl_results
    Array[File] cis_permutation_logs = select_all(tensorqtl_cis_permutations.log)
    Array[Array[File]] cis_nominal_parquet = select_all(tensorqtl_cis_nominal.chr_parquet)
    File? tensorqtl_all_analyses = AggregateTensorqtl.all_analyses
    File? tensorqtl_summary = AggregateTensorqtl.summary
  }
}
