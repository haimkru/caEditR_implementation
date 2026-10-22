version 1.0

# ------------------------------------------------------------------------------
# workflows/tasks/tensorqtl.wdl
#
# Description: tensorQTL cis mapping tasks + aggregation of their results.
#   tensorqtl_cis_permutations and tensorqtl_cis_nominal are COPIED from
#   https://github.com/AoU-Multiomics-Analysis/tensorQTL
#     cis_permutations/workflows/tensorqtl_cis_permutations.wdl
#     cis_nominal/tensorqtl_cis_nominal.wdl
#   (author Francois Aguet), unchanged except: (1) a docker_image input with
#   the same default image, (2) the stage/start_time/completion_time/
#   dimensions/outputs log lines required by tools/check_wdl_logging.py,
#   (3) a cis_window input for the nominal task (the source task always used
#   tensorQTL's built-in 1 Mb), so every scan uses the same window.
#   Keep them in sync with the source repository.
#   AggregateTensorqtl is new (caEditR_edQTL).
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude (copied tasks: Francois Aguet)
#
# Pseudo code:
#   tensorqtl_cis_permutations: python3 -m tensorqtl --mode cis (one phenotype set)
#   tensorqtl_cis_nominal:      python3 -m tensorqtl --mode cis_nominal
#   AggregateTensorqtl:         lead variant per site and phenotype set from the
#                               nominal parquet (default) and/or permutation
#                               tables; raw p-values kept; counts.
# Inputs / outputs: see each task's input and output blocks.
# ------------------------------------------------------------------------------

task tensorqtl_cis_permutations {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        String prefix
        File? phenotype_groups
        Float? fdr
        Float? qvalue_lambda
        Float? pval_threshold
        Int? seed
        String? flags
        Int cis_window = 1000000
        Int memory
        Int disk_space
        Int num_threads
        Int num_gpus
        Int num_preempt
        String? machine_type = "g2-standard-16"
        String gpu_type = "nvidia-l4"
        Array[String] zones = ["us-central1-c"]
        # added in caEditR_edQTL: image override (same default as the source task)
        String docker_image = "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"
    }

    command <<<
        set -euo pipefail
        # added in caEditR_edQTL: standard AoU logging fields
        printf 'stage=%s start_time=%s prefix=%s\n' "tensorqtl_cis_permutations" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" '~{prefix}' >&2
        log() { echo "[$(date -u +%FT%TZ)] stage=tensorqtl_cis_permutations $*" >&2; }
        log "Validate localized inputs."
        plink_pgen='~{sub(plink_pgen, "'", "'\"'\"'")}'
        plink_pvar='~{sub(plink_pvar, "'", "'\"'\"'")}'
        plink_psam='~{sub(plink_psam, "'", "'\"'\"'")}'
        phenotype_bed='~{sub(phenotype_bed, "'", "'\"'\"'")}'
        covariates='~{sub(covariates, "'", "'\"'\"'")}'
        phenotype_groups='~{sub(select_first([phenotype_groups, ""]), "'", "'\"'\"'")}'
        prefix='~{sub(prefix, "'", "'\"'\"'")}'
        for path in "$plink_pgen" "$plink_pvar" "$plink_psam" "$phenotype_bed" "$covariates"; do
            case "$path" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        done
        if [[ -n "$phenotype_groups" ]]; then
            path="$phenotype_groups"
            case "$path" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        fi
        [[ "$prefix" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { log "Input error: prefix must be a filename prefix."; exit 1; }
        # Cromwell can localize the three PLINK files into different directories.
        mkdir -p localized_genotypes
        ln -s "$plink_pgen" localized_genotypes/input.pgen
        ln -s "$plink_pvar" localized_genotypes/input.pvar
        ln -s "$plink_psam" localized_genotypes/input.psam
        args=(localized_genotypes/input "$phenotype_bed" "$prefix" --mode cis
              --covariates "$covariates" --window ~{cis_window})
        if [[ -n "$phenotype_groups" ]]; then args+=(--phenotype_groups "$phenotype_groups"); fi
        fdr='~{default="" fdr}'
        qvalue_lambda='~{default="" qvalue_lambda}'
        pval_threshold='~{default="" pval_threshold}'
        seed='~{default="" seed}'
        if [[ -n "$fdr" ]]; then args+=(--fdr "$fdr"); fi
        if [[ -n "$qvalue_lambda" ]]; then args+=(--qvalue_lambda "$qvalue_lambda"); fi
        if [[ -n "$pval_threshold" ]]; then args+=(--pval_threshold "$pval_threshold"); fi
        if [[ -n "$seed" ]]; then args+=(--seed "$seed"); fi
        flags='~{sub(select_first([flags, ""]), "'", "'\"'\"'")}'
        if [[ -n "$flags" ]]; then
            # Split shell-style arguments without executing shell expressions.
            python3 -c 'import shlex, sys; sys.stdout.buffer.write(b"".join(arg.encode() + b"\0" for arg in shlex.split(sys.argv[1])))' "$flags" > extra_flags.list
            while IFS= read -r -d '' flag; do args+=("$flag"); done < extra_flags.list
        fi
        log "Start cis permutation mapping."
        python3 -m tensorqtl "${args[@]}"
        log "Mapping completed."
        printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' "tensorqtl_cis_permutations" "phenotype_bed:$(basename "$phenotype_bed")" "cis_qtl,log" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >&2
    >>>

    runtime {
        docker: docker_image
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        bootDiskSizeGb: 25
        cpu: num_threads
        preemptible: num_preempt
        # Batch requires String. An empty string lets Batch select the machine.
        predefinedMachineType: select_first([machine_type, ""])
        gpuType: gpu_type
        gpuCount: num_gpus
        zones: zones
    }

    output {
        File cis_qtl = prefix + ".cis_qtl.txt.gz"
        File log = prefix + ".tensorQTL.cis.log"
    }

    parameter_meta {
        machine_type: "GCP machine type. Default: g2-standard-16 with L4. Set null or an empty string to let Batch select the machine from cpu, memory, gpu_type, and num_gpus."
        gpu_type: "GCP GPU type. Default: nvidia-l4. The GPU type and count must be compatible with the selected machine."
        zones: "GCP zones that support the selected machine and GPU. Default: us-central1-c."
    }

    meta {
        author: "Francois Aguet"
    }
}

task tensorqtl_cis_nominal {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        String prefix
        File? interaction
        File? phenotype_groups
        Int memory
        Int disk_space
        Int num_threads
        Int num_gpus = 1
        Int num_preempt
        # added in caEditR_edQTL: same cis window as the permutation task and caNRD-editQTL
        Int cis_window = 1000000
        String? machine_type = "g2-standard-16"
        String gpu_type = "nvidia-l4"
        Array[String] zones = ["us-central1-c"]
        # added in caEditR_edQTL: image override (same default as the source task)
        String docker_image = "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"
    }

    command <<<
        set -euo pipefail
        # added in caEditR_edQTL: standard AoU logging fields
        printf 'stage=%s start_time=%s prefix=%s\n' "tensorqtl_cis_nominal" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" '~{prefix}' >&2
        log() { echo "[$(date -u +%FT%TZ)] stage=tensorqtl_cis_nominal $*" >&2; }
        log "Validate localized inputs."
        plink_pgen='~{sub(plink_pgen, "'", "'\"'\"'")}'
        plink_pvar='~{sub(plink_pvar, "'", "'\"'\"'")}'
        plink_psam='~{sub(plink_psam, "'", "'\"'\"'")}'
        phenotype_bed='~{sub(phenotype_bed, "'", "'\"'\"'")}'
        covariates='~{sub(covariates, "'", "'\"'\"'")}'
        interaction='~{sub(select_first([interaction, ""]), "'", "'\"'\"'")}'
        phenotype_groups='~{sub(select_first([phenotype_groups, ""]), "'", "'\"'\"'")}'
        prefix='~{sub(prefix, "'", "'\"'\"'")}'
        for path in "$plink_pgen" "$plink_pvar" "$plink_psam" "$phenotype_bed" "$covariates"; do
            case "$path" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        done
        for path in "$interaction" "$phenotype_groups"; do
            [[ -n "$path" ]] || continue
            case "$path" in gs://*|s3://*|https://*|http://*) log "Input localization error: unresolved URI $path"; exit 1;; esac
            [[ -r "$path" ]] || { log "Input localization error: unreadable file $path"; exit 1; }
        done
        [[ "$prefix" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { log "Input error: prefix must be a filename prefix."; exit 1; }
        # Cromwell can localize the three PLINK files into different directories.
        mkdir -p localized_genotypes
        ln -s "$plink_pgen" localized_genotypes/input.pgen
        ln -s "$plink_pvar" localized_genotypes/input.pvar
        ln -s "$plink_psam" localized_genotypes/input.psam
        args=(localized_genotypes/input "$phenotype_bed" "$prefix" --mode cis_nominal
              --covariates "$covariates" --window ~{cis_window})
        if [[ -n "$interaction" ]]; then args+=(--interaction "$interaction"); fi
        if [[ -n "$phenotype_groups" ]]; then args+=(--phenotype_groups "$phenotype_groups"); fi
        log "Start nominal cis mapping; interaction enabled: ~{defined(interaction)}."
        python3 -m tensorqtl "${args[@]}"
        log "Mapping completed."
        printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' "tensorqtl_cis_nominal" "phenotype_bed:$(basename "$phenotype_bed")" "chr_parquet,log" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >&2
    >>>

    runtime {
        docker: docker_image
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        bootDiskSizeGb: 25
        cpu: num_threads
        preemptible: num_preempt
        # Batch requires String. An empty string lets Batch select the machine.
        predefinedMachineType: select_first([machine_type, ""])
        gpuType: gpu_type
        gpuCount: num_gpus
        zones: zones
    }


    output {
        Array[File] chr_parquet = glob(prefix + "*.parquet")
        File log = glob(prefix + "*.log")[0]
    }
    parameter_meta {
        machine_type: "GCP machine type. Default: g2-standard-16 with L4. Set null or an empty string to let Batch select the machine from cpu, memory, gpu_type, and num_gpus."
        gpu_type: "GCP GPU type. Default: nvidia-l4. The GPU type and count must be compatible with the selected machine."
        zones: "GCP zones that support the selected machine and GPU. Default: us-central1-c."
    }

    meta {
        author: "Francois Aguet"
    }
}

task AggregateTensorqtl {
  input {
    Array[File] cis_nominal_parquet = []
    Array[File] cis_qtl_results = []
    String output_prefix
    Float fdr = 0.05
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
    stage="aggregate_tensorqtl"
    log="$stage.log"
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    cp '~{write_lines(cis_nominal_parquet)}' cis_nominal_list.txt
    cp '~{write_lines(cis_qtl_results)}' cis_permutation_list.txt
    Rscript /opt/caeditr_edqtl/scripts/aggregate_tensorqtl.R \
      --cis-nominal-list cis_nominal_list.txt \
      --cis-permutation-list cis_permutation_list.txt \
      --fdr ~{fdr} \
      --qvalue-lambda '~{default="" qvalue_lambda}' \
      --prefix '~{output_prefix}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "nominal_files:$(wc -l < cis_nominal_list.txt | tr -d ' '),permutation_files:$(wc -l < cis_permutation_list.txt | tr -d ' ')" \
      "all_analyses,summary" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File all_analyses = "outputs/" + output_prefix + ".tensorqtl_cis.all_analyses.tsv.gz"
    File summary = "outputs/" + output_prefix + ".tensorqtl_cis.summary.tsv"
    File log = "aggregate_tensorqtl.log"
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
