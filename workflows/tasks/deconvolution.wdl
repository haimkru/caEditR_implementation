version 1.0

# ------------------------------------------------------------------------------
# workflows/tasks/deconvolution.wdl
#
# Description: per-chunk cell-type deconvolution (caNRD-edit, TCA) and the
#              genome-wide merge of the chunk outputs.
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   DeconvolveChunk:   one scatter shard: tabix-read the chunk's sites, site
#                      QC, bulk ratios, caNRD_edit, TCA_Like (chunk BEDs).
#   MergeChunkOutputs: concatenate chunk BEDs per (method, cell type), index,
#                      concatenate reports.
# Inputs / outputs: see each task's input and output blocks.
# ------------------------------------------------------------------------------

task DeconvolveChunk {
  input {
    String chunk_id
    File chunks
    File site_chunks
    File editing_ratios
    File editing_ratios_index
    File editing_coverage
    File editing_coverage_index
    File site_theta
    File proportions
    File cohort_samples
    Float min_coverage = 20
    Int min_covered_samples = 60
    Float min_covered_fraction = 0.0
    Float theta_floor = 0.001
    Float canrd_min_coverage = 10
    Array[String] methods = ["canrd_edit", "tca"]
    String docker_image
    Int cpu = 1
    String memory = "6 GB"
    Int disk_gb = 100
    Int preemptible_attempts = 2
    Int max_retries = 1
  }

  command <<<
    set -euo pipefail
    stage="deconvolve_chunk"
    log="$stage.log"
    printf 'stage=%s start_time=%s chunk=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" '~{chunk_id}' | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    # Cromwell may localize a BED and its .tbi into different directories:
    # link both under one directory so tabix finds the index.
    mkdir -p localized
    ln -s '~{editing_ratios}' localized/editing_ratios.bed.gz
    ln -s '~{editing_ratios_index}' localized/editing_ratios.bed.gz.tbi
    ln -s '~{editing_coverage}' localized/editing_coverage.bed.gz
    ln -s '~{editing_coverage_index}' localized/editing_coverage.bed.gz.tbi
    Rscript /opt/caeditr_implementation/scripts/deconvolve_chunk.R \
      --chunk-id '~{chunk_id}' \
      --chunks '~{chunks}' \
      --site-chunks '~{site_chunks}' \
      --editing-ratios localized/editing_ratios.bed.gz \
      --editing-coverage localized/editing_coverage.bed.gz \
      --site-theta '~{site_theta}' \
      --proportions '~{proportions}' \
      --samples '~{cohort_samples}' \
      --min-coverage ~{min_coverage} \
      --min-covered-samples ~{min_covered_samples} \
      --min-covered-fraction ~{min_covered_fraction} \
      --theta-floor ~{theta_floor} \
      --canrd-min-coverage ~{canrd_min_coverage} \
      --methods '~{sep="," methods}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "chunk_beds:$(ls outputs/beds | wc -l | tr -d ' ')" \
      "chunk_beds,chunk_reports" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    Array[File] chunk_beds = glob("outputs/beds/*.bed.gz")
    Array[File] chunk_reports = glob("outputs/reports/*.tsv")
    File log = "deconvolve_chunk.log"
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

task MergeChunkOutputs {
  input {
    Array[File] chunk_beds
    Array[File] chunk_reports
    String output_prefix
    String docker_image
    Int cpu = 4
    String memory = "16 GB"
    Int disk_gb = 250
    Int preemptible_attempts = 2
    Int max_retries = 1
  }

  command <<<
    set -euo pipefail
    stage="merge_chunk_outputs"
    log="$stage.log"
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    # file lists are generated at command time and stay File-typed (Terra path rules)
    cp '~{write_lines(chunk_beds)}' bed_list.txt
    cp '~{write_lines(chunk_reports)}' report_list.txt
    bash /opt/caeditr_implementation/scripts/merge_chunk_outputs.sh \
      --bed-list bed_list.txt \
      --report-list report_list.txt \
      --prefix '~{output_prefix}' \
      --threads ~{cpu} \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "phenotype_sets:$(tail -n +2 outputs/phenotype_inventory.tsv | wc -l | tr -d ' ')" \
      "bulk_bed,celltype_beds,reports,phenotype_inventory" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File bulk_bed = "outputs/merged/bulk/" + output_prefix + ".bulk.all.bed.gz"
    File bulk_bed_index = "outputs/merged/bulk/" + output_prefix + ".bulk.all.bed.gz.tbi"
    Array[File] celltype_beds = glob("outputs/merged/celltype/*.bed.gz")
    Array[File] celltype_bed_indexes = glob("outputs/merged/celltype/*.bed.gz.tbi")
    File phenotype_inventory = "outputs/phenotype_inventory.tsv"
    File site_qc = "outputs/reports/" + output_prefix + ".site_qc.tsv.gz"
    File canrd_diagnostics = "outputs/reports/" + output_prefix + ".canrd_diagnostics.tsv.gz"
    File deconvolution_summary = "outputs/reports/" + output_prefix + ".deconvolution_summary.tsv"
    File log = "merge_chunk_outputs.log"
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
