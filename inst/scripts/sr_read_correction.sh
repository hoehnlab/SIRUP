#!/bin/bash

# Summarize CB/CR correction stats from Cell Ranger BAM
# Usage: ./sr_read_correction.sh input.bam output_prefix

set -euo pipefail

BAM_FILE=${1:-}
OUT_PREFIX=${2:-}
THREADS=${THREADS:-4}

if [ -z "$BAM_FILE" ] || [ -z "$OUT_PREFIX" ]; then
  echo "Usage: $0 <input.bam> <output_prefix>" >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "Error: samtools not found on PATH." >&2
  exit 1
fi
if ! command -v csvtk >/dev/null 2>&1; then
  echo "Error: csvtk not found on PATH." >&2
  exit 1
fi

TSV_RAW="${OUT_PREFIX}_raw.tsv"
SUMMARY_CSV="${OUT_PREFIX}_barcode_summary.csv"
PAIR_CSV="${OUT_PREFIX}_correction_pairs.csv"

samtools view -@ "$THREADS" "$BAM_FILE" | awk '
BEGIN { OFS="\t" }
{
    cb = cr = ""
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^CB:Z:/) cb = substr($i, 6)
        else if ($i ~ /^CR:Z:/) cr = substr($i, 6)
    }
    if (cb != "" && cr != "")
        print cb, cr
}' > "$TSV_RAW"

# Summary: for each CB, count raw==corrected (raw), raw!=corrected (corrected), total
awk '
{
    cb = $1
    cr = $2
    total[cb]++
    if (cb == cr) raw[cb]++
    else corrected[cb]++
}
END {
    print "barcode,raw_ct,corrected_ct,total_ct"
    for (k in total) {
        r = (k in raw ? raw[k] : 0)
        c = (k in corrected ? corrected[k] : 0)
        print k "," r "," c "," total[k]
    }
}' "$TSV_RAW" > "$SUMMARY_CSV"

# Pair counts (CB, CR)
csvtk -t freq -f 1,2 "$TSV_RAW" | csvtk -t rename -f 1:corrected_barcode,2:original_barcode,3:count > "$PAIR_CSV"

echo "Done. Output files:"
echo "  Summary table:        $SUMMARY_CSV"
echo "  Correction pairs:     $PAIR_CSV"
