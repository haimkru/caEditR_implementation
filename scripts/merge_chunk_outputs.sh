#!/usr/bin/env bash
################################################################################
# scripts/merge_chunk_outputs.sh
#
# Description:
#   Gather stage after the DeconvolveChunk scatter (WDL task
#   MergeChunkOutputs). Concatenates the per-chunk BED matrices into one
#   genome-wide, bgzipped, tabix-indexed BED per (method, cell type), and the
#   per-chunk reports into one table per report type. Streams text (no R), so
#   memory stays flat at 100k sites x 10k samples.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   1. Read the list of chunk BED files (one path per line). File names are
#      <method>.<celltype_slug>.<chunk_id>.bed.gz.
#   2. For every <method>.<celltype_slug> group: sort its files by chunk_id
#      (chunks are in genomic order), write the header of the first file and
#      the body of every file, bgzip, tabix.
#      Bulk goes to merged/bulk/, cell types to merged/celltype/.
#   3. Write phenotype_inventory.tsv (method, slug, file, sites, samples).
#   4. Concatenate reports of the same type (<chunk>.<report>.tsv) keeping one
#      header: site_qc, canrd_diagnostics, deconvolution_summary.
#
# Inputs:
#   --bed-list     text file, one chunk BED path per line
#   --report-list  text file, one chunk report path per line
#   --prefix       output file prefix (e.g. "aou_wgs_v8"); letters, digits, _ -
#   --output-dir   output directory
#   --threads      bgzip threads (default 1)
# Outputs:
#   merged/bulk/<prefix>.bulk.all.bed.gz(.tbi)
#   merged/celltype/<prefix>.<method>.<slug>.bed.gz(.tbi)
#   reports/<prefix>.site_qc.tsv.gz, <prefix>.canrd_diagnostics.tsv.gz,
#   reports/<prefix>.deconvolution_summary.tsv, phenotype_inventory.tsv
################################################################################
set -euo pipefail

bed_list=""
report_list=""
prefix=""
output_dir="outputs"
threads=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bed-list) bed_list="$2"; shift 2 ;;
    --report-list) report_list="$2"; shift 2 ;;
    --prefix) prefix="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    *) echo "Input error: unknown argument $1" >&2; exit 1 ;;
  esac
done
[[ -r "$bed_list" ]] || { echo "Input error: unreadable --bed-list $bed_list" >&2; exit 1; }
[[ -r "$report_list" ]] || { echo "Input error: unreadable --report-list $report_list" >&2; exit 1; }
[[ "$prefix" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || {
  echo "Input error: --prefix must start with a letter/digit and contain only letters, digits, _ or -" >&2
  exit 1
}

mkdir -p "$output_dir/merged/bulk" "$output_dir/merged/celltype" "$output_dir/reports"
inventory="$output_dir/phenotype_inventory.tsv"
printf 'method\tcelltype_slug\tfile\tn_sites\tn_samples\tn_chunks\n' > "$inventory"

# ---- 1-2. BED groups -----------------------------------------------------------
# table: group <tab> chunk_id <tab> path   (sorted by group, then chunk_id)
group_table="$output_dir/bed_groups.tsv"
while IFS= read -r path; do
  [[ -n "$path" ]] || continue
  [[ -r "$path" ]] || { echo "Input error: unreadable chunk BED $path" >&2; exit 1; }
  name="$(basename "$path" .bed.gz)"          # method.slug.chunk_id
  chunk_id="${name##*.}"
  group="${name%.*}"
  printf '%s\t%s\t%s\n' "$group" "$chunk_id" "$path"
done < "$bed_list" | LC_ALL=C sort -t "$(printf '\t')" -k1,1 -k2,2 > "$group_table"

for group in $(cut -f1 "$group_table" | uniq); do
  method="${group%%.*}"
  slug="${group#*.}"
  if [[ "$method" == "bulk" ]]; then
    merged="$output_dir/merged/bulk/${prefix}.${group}.bed.gz"
  else
    merged="$output_dir/merged/celltype/${prefix}.${group}.bed.gz"
  fi
  files=()
  while IFS= read -r f; do files+=("$f"); done < <(awk -F'\t' -v g="$group" '$1 == g { print $3 }' "$group_table")
  {
    # header of the first chunk; head closes the pipe early, so ignore SIGPIPE here
    gzip -dc "${files[0]}" | head -n 1 || true
    for f in "${files[@]}"; do
      gzip -dc "$f" | tail -n +2
    done
  } | bgzip -@ "$threads" -c > "$merged"
  tabix -f -p bed "$merged"
  n_sites="$(gzip -dc "$merged" | tail -n +2 | wc -l | tr -d ' ')"
  n_samples="$(gzip -dc "$merged" | head -n 1 | awk -F'\t' '{ print NF - 4 }' || true)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$method" "$slug" "$(basename "$merged")" \
    "$n_sites" "$n_samples" "${#files[@]}" >> "$inventory"
  echo "Merged ${#files[@]} chunks into $merged ($n_sites sites)" >&2
done
rm -f "$group_table"

# ---- 4. reports -----------------------------------------------------------------
merge_reports() {
  local report_type="$1"
  local output="$2"
  reports=()
  while IFS= read -r r; do reports+=("$r"); done < <(grep -E "\.${report_type}\.tsv$" "$report_list" | LC_ALL=C sort || true)
  if [[ ${#reports[@]} -eq 0 ]]; then
    : > "$output"
    return
  fi
  {
    head -n 1 "${reports[0]}"
    for r in "${reports[@]}"; do
      tail -n +2 "$r"
    done
  } > "${output%.gz}"
  if [[ "$output" == *.gz ]]; then
    gzip -f "${output%.gz}"
  fi
}
merge_reports "site_qc" "$output_dir/reports/${prefix}.site_qc.tsv.gz"
merge_reports "canrd_diagnostics" "$output_dir/reports/${prefix}.canrd_diagnostics.tsv.gz"
merge_reports "deconvolution_summary" "$output_dir/reports/${prefix}.deconvolution_summary.tsv"
echo "Merged reports into $output_dir/reports" >&2
