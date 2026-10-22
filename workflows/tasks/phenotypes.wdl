version 1.0

# ------------------------------------------------------------------------------
# workflows/tasks/phenotypes.wdl
#
# Description: turns one merged editing matrix (bulk or one method x cell
#              type) into tensorQTL inputs: INT phenotype BED + covariates.
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   PreparePhenotypes: site missingness filter, mean imputation, INT,
#                      editing PCs, covariates (supplied + PCs) in tensorQTL
#                      layout; reports the number of kept sites so the
#                      workflow can skip empty sets.
# Inputs / outputs: see the input and output blocks.
# ------------------------------------------------------------------------------

task PreparePhenotypes {
  input {
    File input_bed
    File covariates
    Int n_pcs = 15
    Float max_missing_fraction = 0.5
    String docker_image
    Int cpu = 2
    String memory = "8 GB"
    Int disk_gb = 100
    Int preemptible_attempts = 2
    Int max_retries = 1
  }

  String label = basename(input_bed, ".bed.gz")

  command <<<
    set -euo pipefail
    stage="prepare_phenotypes"
    log="$stage.log"
    printf 'stage=%s start_time=%s label=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" '~{label}' | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/caeditr_edqtl/scripts/prepare_phenotypes.R \
      --input-bed '~{input_bed}' \
      --covariates '~{covariates}' \
      --n-pcs ~{n_pcs} \
      --max-missing-fraction ~{max_missing_fraction} \
      --label '~{label}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "sites_kept:$(cat 'outputs/~{label}.n_sites_kept.txt')" \
      "int_bed,tensorqtl_covariates,phenotype_pcs,filter_report" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    String phenotype_label = label
    Int n_sites_kept = read_int("outputs/" + label + ".n_sites_kept.txt")
    File int_bed = "outputs/" + label + ".INT.bed.gz"
    File int_bed_index = "outputs/" + label + ".INT.bed.gz.tbi"
    File tensorqtl_covariates = "outputs/" + label + ".tensorqtl_covariates.tsv"
    File phenotype_pcs = "outputs/" + label + ".phenotype_PCs.tsv"
    File pc_variance_explained = "outputs/" + label + ".pc_variance_explained.tsv"
    File phenotype_filter_report = "outputs/" + label + ".phenotype_filter_report.tsv"
    File phenotype_summary = "outputs/" + label + ".phenotype_summary.tsv"
    File log = "prepare_phenotypes.log"
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: memory
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}
