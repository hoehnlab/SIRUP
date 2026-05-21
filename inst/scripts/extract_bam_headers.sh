#!/bin/bash
# Usage: ./extract_bam_headers.sh input.bam output.csv

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.bam> <output.csv>" >&2
    exit 1
fi

input_bam=$1
output_csv=$2
threads=${THREADS:-4}

# Check if samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "samtools not found. Please install it to proceed." >&2
    exit 1
fi

samtools view -@ "$threads" "$input_bam" \
| awk -F'\t' '
BEGIN {
    OFS = ",";
    print "id", "final_barcode", "original_barcode", "strand", "umi";
}
{
    id = $1;
    cb = cr = ub = ts = "";

    for (i = 12; i <= NF; i++) {
        if ($i ~ /^CB:Z:/) cb = substr($i, 6);
        else if ($i ~ /^CR:Z:/) cr = substr($i, 6);
        else if ($i ~ /^UB:Z:/) ub = substr($i, 6);
        else if ($i ~ /^TS:A:/) ts = substr($i, 6);
    }

    print id, cb, cr, ts, ub;
}' > "$output_csv"

echo "Extracted barcode data saved to $output_csv"
