version 1.0

# ------------------------------------------------------------------------------
# workflows/tasks/cohort.wdl
#
# Description: input indexing and cohort preparation tasks of caEditR_implementation.
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   IndexBedMatrix: sort + bgzip + tabix one BED matrix (ratios or coverage),
#                   so later tasks read one genomic chunk at a time.
#   PrepareCohort:  shared cohort definition (samples, site -> host gene,
#                   theta, chunks); reads only coordinates of the big BEDs.
# Inputs / outputs: see each task's input and output blocks.
# ------------------------------------------------------------------------------

task IndexBedMatrix {
  input {
    File bed
    String output_name
    String docker_image
    Int cpu = 4
    String memory = "16 GB"
    Int disk_gb = 200
    Int preemptible_attempts = 2
    Int max_retries = 1
  }

  command <<<
    set -euo pipefail
    stage="index_bed_matrix"
    log="$stage.log"
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    output_name='~{output_name}'
    [[ "$output_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*\.bed\.gz$ ]] || { echo "Input error: output_name must be a file name ending in .bed.gz" | tee -a "$log"; exit 1; }
    mkdir -p outputs
    bash /opt/caeditr_implementation/scripts/index_bed_matrix.sh \
      --input '~{bed}' \
      --output "outputs/$output_name" \
      --threads ~{cpu} 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "rows:$(gzip -dc "outputs/$output_name" | tail -n +2 | wc -l | tr -d ' ')" \
      "indexed_bed,indexed_bed_index" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File indexed_bed = "outputs/" + output_name
    File indexed_bed_index = "outputs/" + output_name + ".tbi"
    File log = "index_bed_matrix.log"
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

task PrepareCohort {
  input {
    File editing_ratios
    File editing_coverage
    File expression
    File gtf
    File proportions
    File plink_psam
    File covariates
    File? sample_list
    Int sites_per_chunk = 500
    Float theta_floor = 0.001
    String docker_image
    Int cpu = 4
    String memory = "32 GB"
    Int disk_gb = 200
    Int preemptible_attempts = 2
    Int max_retries = 1
  }

  command <<<
    set -euo pipefail
    stage="prepare_cohort"
    log="$stage.log"
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    sample_list_path='~{default="" sample_list}'
    optional_args=()
    if [[ -n "$sample_list_path" ]]; then optional_args+=(--sample-list "$sample_list_path"); fi
    Rscript /opt/caeditr_implementation/scripts/prepare_cohort.R \
      --editing-ratios '~{editing_ratios}' \
      --editing-coverage '~{editing_coverage}' \
      --expression '~{expression}' \
      --gtf '~{gtf}' \
      --proportions '~{proportions}' \
      --psam '~{plink_psam}' \
      --covariates '~{covariates}' \
      --sites-per-chunk ~{sites_per_chunk} \
      --theta-floor ~{theta_floor} \
      ${optional_args[@]+"${optional_args[@]}"} \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "samples:$(wc -l < outputs/cohort_samples.txt | tr -d ' '),chunks:$(wc -l < outputs/chunk_ids.txt | tr -d ' ')" \
      "cohort_samples,proportions,covariates,site_theta,chunks" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File cohort_samples = "outputs/cohort_samples.txt"
    File sample_alignment_report = "outputs/sample_alignment_report.tsv"
    File celltypes = "outputs/celltypes.tsv"
    File proportions_aligned = "outputs/proportions_aligned.tsv"
    File covariates_aligned = "outputs/covariates_aligned.tsv"
    File site_to_gene = "outputs/site_to_gene.tsv"
    File gene_theta = "outputs/gene_theta.tsv"
    File site_theta = "outputs/site_theta.tsv"
    File chunks = "outputs/chunks.tsv"
    File site_chunks = "outputs/site_chunks.tsv"
    Array[String] chunk_ids = read_lines("outputs/chunk_ids.txt")
    File cohort_summary = "outputs/cohort_summary.tsv"
    File log = "prepare_cohort.log"
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
