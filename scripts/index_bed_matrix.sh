#!/usr/bin/env bash
################################################################################
# scripts/index_bed_matrix.sh
#
# Description:
#   Sorts a BED matrix (#chr, start, end, id, sample columns ...) by genomic
#   position, compresses it with bgzip and indexes it with tabix, so later
#   tasks can read one genomic chunk at a time (tabix region queries) instead
#   of loading 100k sites x 10k samples into memory.
#   Accepts plain or gzipped input. Used for the editing-ratio and the
#   editing-coverage BEDs.
#
# Version:      0.1.0
# Date written: 2026-10-22
# Author:       Haim Krupkin+Claude
#
# Pseudo code:
#   1. Stream the input (gzip -dc for .gz, cat otherwise).
#   2. Keep the header line; sort the body by chromosome (natural order:
#      chr1, chr2, ..., chr10), start, end and id (id breaks ties, so two
#      files with the same sites always get the same row order).
#   3. bgzip the result, tabix -p bed it, check the index exists.
#
# Inputs:
#   --input   BED or BED.GZ matrix
#   --output  output path ending in .bed.gz
#   --threads bgzip threads (default 1)
# Outputs:
#   <output> and <output>.tbi
################################################################################
set -euo pipefail

input=""
output=""
threads=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) input="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    *) echo "Input error: unknown argument $1" >&2; exit 1 ;;
  esac
done

[[ -r "$input" ]] || { echo "Input error: unreadable input BED: $input" >&2; exit 1; }
[[ "$output" == *.bed.gz ]] || { echo "Input error: --output must end in .bed.gz" >&2; exit 1; }

# 1. a command that streams the input as text
if [[ "$input" == *.gz || "$input" == *.bgz ]]; then
  reader=(gzip -dc "$input")
else
  reader=(cat "$input")
fi

# 2. header + sorted body; a local temp dir keeps big sorts on the task disk
sort_tmp="$(dirname "$output")/sort_tmp"
mkdir -p "$sort_tmp"
tab="$(printf '\t')"
# head closes the pipe early (SIGPIPE for the reader): tolerate it for the header only
header_line="$("${reader[@]}" 2>/dev/null | head -n 1 || true)"
[[ "$header_line" == "#chr${tab}start${tab}end${tab}"* ]] || {
  echo "Input error: $input must start with the header #chr, start, end, <id>, samples..." >&2
  exit 1
}
{
  printf '%s\n' "$header_line"
  "${reader[@]}" | tail -n +2 \
    | LC_ALL=C sort -t "$tab" -k1,1V -k2,2n -k3,3n -k4,4 -S 40% -T "$sort_tmp"
} | bgzip -@ "$threads" -c > "$output"
rm -rf "$sort_tmp"

# 3. tabix index (BED preset: 0-based start, end in column 3)
tabix -f -p bed "$output"
[[ -s "${output}.tbi" ]] || { echo "tabix index was not created for $output" >&2; exit 1; }

echo "Indexed $(gzip -dc "$output" | tail -n +2 | wc -l | tr -d ' ') rows into $output" >&2
