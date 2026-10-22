version 1.0

# ------------------------------------------------------------------------------
# workflows/tasks/canrd_editqtl.wdl
#
# Description: end-to-end cell-type edQTL scan with caEditR::caNRD_editQTL
#              (scan engine + exact lead re-fit) and its aggregation.
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   CanrdEditqtlChunk:     one chunk of sites; plink2 extracts the chunk's cis
#                          dosages; scan all (site, variant) pairs.
#   AggregateCanrdEditqtl: Bonferroni over each site's variants, Storey q
#                          across sites, cell-type calls; finalized tables.
# Inputs / outputs: see each task's input and output blocks.
# ------------------------------------------------------------------------------

task CanrdEditqtlChunk {
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
    File covariates
    File? editing_pcs
    Int n_editing_pcs = 0
    File plink_pgen
    File plink_pvar
    File plink_psam
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
    String docker_image
    Int cpu = 4
    String memory = "24 GB"
    Int disk_gb = 150
    Int preemptible_attempts = 2
    Int max_retries = 1
  }

  command <<<
    set -euo pipefail
    stage="canrd_editqtl_chunk"
    log="$stage.log"
    printf 'stage=%s start_time=%s chunk=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" '~{chunk_id}' | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    # one directory per index pair / PLINK trio (Cromwell may localize them apart)
    mkdir -p localized
    ln -s '~{editing_ratios}' localized/editing_ratios.bed.gz
    ln -s '~{editing_ratios_index}' localized/editing_ratios.bed.gz.tbi
    ln -s '~{editing_coverage}' localized/editing_coverage.bed.gz
    ln -s '~{editing_coverage_index}' localized/editing_coverage.bed.gz.tbi
    ln -s '~{plink_pgen}' localized/genotypes.pgen
    ln -s '~{plink_pvar}' localized/genotypes.pvar
    ln -s '~{plink_psam}' localized/genotypes.psam
    editing_pcs_path='~{default="" editing_pcs}'
    optional_args=()
    if [[ -n "$editing_pcs_path" ]]; then optional_args+=(--editing-pcs "$editing_pcs_path"); fi
    Rscript /opt/caeditr_implementation/scripts/canrd_editqtl_chunk.R \
      --chunk-id '~{chunk_id}' \
      --chunks '~{chunks}' \
      --site-chunks '~{site_chunks}' \
      --editing-ratios localized/editing_ratios.bed.gz \
      --editing-coverage localized/editing_coverage.bed.gz \
      --site-theta '~{site_theta}' \
      --proportions '~{proportions}' \
      --samples '~{cohort_samples}' \
      --covariates '~{covariates}' \
      --n-editing-pcs ~{n_editing_pcs} \
      --pfile-prefix localized/genotypes \
      --cis-window ~{cis_window} \
      --maf ~{maf} \
      --min-minor-allele-samples ~{min_minor_allele_samples} \
      --min-samples ~{min_samples} \
      --min-mean-phi ~{min_mean_phi} \
      --theta-floor ~{theta_floor} \
      --min-coverage ~{min_coverage} \
      --min-covered-samples ~{min_covered_samples} \
      --min-covered-fraction ~{min_covered_fraction} \
      --refine '~{refine}' \
      --write-all-pairs ~{write_all_pairs} \
      ${optional_args[@]+"${optional_args[@]}"} \
      --output-dir outputs 2>&1 | tee -a "$log"
    # an empty all-pairs file keeps the output contract when write_all_pairs = false
    all_pairs='outputs/~{chunk_id}.canrd_editqtl.all_pairs.tsv.gz'
    [[ -e "$all_pairs" ]] || printf '' | gzip -c > "$all_pairs"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "lead_rows:$(tail -n +2 'outputs/~{chunk_id}.canrd_editqtl.lead_variants.tsv' | wc -l | tr -d ' ')" \
      "all_pairs,lead_variants,site_summary,site_celltype_min_p,scan_info,site_qc" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File all_pairs = "outputs/" + chunk_id + ".canrd_editqtl.all_pairs.tsv.gz"
    File lead_variants = "outputs/" + chunk_id + ".canrd_editqtl.lead_variants.tsv"
    File site_summary = "outputs/" + chunk_id + ".canrd_editqtl.site_summary.tsv"
    File site_celltype_min_p = "outputs/" + chunk_id + ".canrd_editqtl.site_celltype_min_p.tsv"
    File scan_info = "outputs/" + chunk_id + ".canrd_editqtl.scan_info.tsv"
    File site_qc = "outputs/" + chunk_id + ".canrd_editqtl.site_qc.tsv"
    File log = "canrd_editqtl_chunk.log"
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

task AggregateCanrdEditqtl {
  input {
    Array[File] site_summaries
    Array[File] lead_variants
    Array[File] site_celltype_min_p
    String output_prefix
    Float fdr = 0.05
    Float celltype_alpha = 0.05
    Float? qvalue_lambda
    String docker_image
    Int cpu = 2
    String memory = "16 GB"
    Int disk_gb = 50
    Int preemptible_attempts = 2
    Int max_retries = 1
  }

  command <<<
    set -euo pipefail
    stage="aggregate_canrd_editqtl"
    log="$stage.log"
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    cp '~{write_lines(site_summaries)}' site_summary_list.txt
    cp '~{write_lines(lead_variants)}' lead_list.txt
    cp '~{write_lines(site_celltype_min_p)}' site_celltype_list.txt
    Rscript /opt/caeditr_implementation/scripts/aggregate_canrd_editqtl.R \
      --site-summary-list site_summary_list.txt \
      --lead-list lead_list.txt \
      --site-celltype-list site_celltype_list.txt \
      --fdr ~{fdr} \
      --celltype-alpha ~{celltype_alpha} \
      --qvalue-lambda '~{default="" qvalue_lambda}' \
      --prefix '~{output_prefix}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "sites:$(tail -n +2 'outputs/~{output_prefix}.canrd_editqtl.site_level.tsv' | wc -l | tr -d ' ')" \
      "site_level,lead_site_celltype,celltype_scan_level,summary" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File site_level = "outputs/" + output_prefix + ".canrd_editqtl.site_level.tsv"
    File lead_site_celltype = "outputs/" + output_prefix + ".canrd_editqtl.lead_site_celltype.tsv"
    File celltype_scan_level = "outputs/" + output_prefix + ".canrd_editqtl.celltype_scan_level.tsv"
    File summary = "outputs/" + output_prefix + ".canrd_editqtl.summary.tsv"
    File log = "aggregate_canrd_editqtl.log"
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
