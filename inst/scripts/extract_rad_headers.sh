#!/bin/bash
# Extract RAD tags (CB/CR/UB/TS) from FASTQ headers to CSV.
# Optional: include read length column when --length is provided.
# Usage: ./extract_rad_headers.sh input.fastq.gz output.csv [--length]

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_fastq> <output_csv> [--length]" >&2
    exit 1
fi

if ! command -v seqkit >/dev/null 2>&1; then
    echo "Error: seqkit not found on PATH." >&2
    exit 1
fi

input_fastq=$1
output_csv=$2
with_length=${3:-""}

if [ "$with_length" == "--length" ]; then
  seqkit fx2tab -n -l "$input_fastq" \
  | awk -F'\t' '
BEGIN { 
    OFS=","; 
    print "id","final_barcode","original_barcode","strand","umi","read_length" 
}
{
    read_length = $NF
    n = split($1, fields, " ")
    raw_id = fields[1]
    gsub(/-[FR](-CT)?$/, "", raw_id)
    id = raw_id
    cb = ""
    cr = ""
    ub = ""
    ts = ""

    for (i = 2; i < NF; i++) {
        if ($i ~ /^CB:Z:/) cb = substr($i, 6)
        else if ($i ~ /^CR:Z:/) cr = substr($i, 6)
        else if ($i ~ /^UB:Z:/) ub = substr($i, 6)
        else if ($i ~ /^TS:A:/) ts = substr($i, 6)
    }

    print id, cb, cr, ts, ub, read_length
}' > "$output_csv"
else
  seqkit seq -n "$input_fastq" \
  | sed 's/^@//' \
  | awk -F'\t' '
BEGIN { 
    OFS=","; 
    print "id","final_barcode","original_barcode","strand","umi" 
}
{
    raw_id = $1
    gsub(/-[FR](-CT)?$/, "", raw_id)
    id = raw_id
    cb = ""
    cr = ""
    ub = ""
    ts = ""

    for (i = 2; i <= NF; i++) {
        if ($i ~ /^CB:Z:/) cb = substr($i, 6)
        else if ($i ~ /^CR:Z:/) cr = substr($i, 6)
        else if ($i ~ /^UB:Z:/) ub = substr($i, 6)
        else if ($i ~ /^TS:A:/) ts = substr($i, 6)
    }

    print id, cb, cr, ts, ub
}' > "$output_csv"
fi

echo "Headers extracted to $output_csv"
