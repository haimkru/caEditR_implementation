version 1.0

# ------------------------------------------------------------------------------
# workflows/tasks/summary.wdl
#
# Description: side-by-side comparison of the three edQTL approaches
#              (bulk tensorQTL, tensorQTL on cell-type matrices,
#              end-to-end caNRD-editQTL).
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   SummarizeEdqtlResults: long table of calls per (site, approach, cell
#                          type) and counts incl. cell-type edSites not found
#                          by the bulk scan.
# Inputs / outputs: see the input and output blocks.
# ------------------------------------------------------------------------------

task SummarizeEdqtlResults {
  input {
    File? tensorqtl_results
    File? canrd_lead_site_celltype
    File celltypes
    String output_prefix
    String docker_image
    Int cpu = 1
    String memory = "8 GB"
    Int disk_gb = 20
    Int preemptible_attempts = 2
    Int max_retries = 1
  }

  command <<<
    set -euo pipefail
    stage="summarize_edqtl_results"
    log="$stage.log"
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/caeditr_edqtl/scripts/summarize_edqtl_results.R \
      --tensorqtl-results '~{default="" tensorqtl_results}' \
      --canrd-lead '~{default="" canrd_lead_site_celltype}' \
      --celltypes '~{celltypes}' \
      --prefix '~{output_prefix}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "overview_rows:$(tail -n +2 'outputs/~{output_prefix}.edsite_overview.tsv' | wc -l | tr -d ' ')" \
      "edsite_calls_long,edsite_overview" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File edsite_calls_long = "outputs/" + output_prefix + ".edsite_calls_long.tsv.gz"
    File edsite_overview = "outputs/" + output_prefix + ".edsite_overview.tsv"
    File log = "summarize_edqtl_results.log"
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
