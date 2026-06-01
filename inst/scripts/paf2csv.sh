#!/bin/bash

# Convert PAF to CSV and split embedded UMI/barcode from read IDs
# Usage: paf2csv.sh <input_file.paf> <output_directory>

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_file.paf> <output_directory>" >&2
    exit 1
fi

INPUT_PAF="$1"
OUTPUT_DIR="$2"

if [ ! -f "$INPUT_PAF" ]; then
    echo "Error: Input file '$INPUT_PAF' does not exist." >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

CSV_OUTPUT="$OUTPUT_DIR/output.csv"
PROCESSED_OUTPUT="$OUTPUT_DIR/output_processed.csv"

# PAF columns:
# 1 qname, 2 qlen, 3 qstart, 4 qend, 5 strand, 6 tname, 7 tlen,
# 8 tstart, 9 tend, 10 nmatch, 11 alen, 12 mapq
awk -F '\t' 'BEGIN {OFS=","; print "qname","qlen","qstart","qend","strand","tname","tlen","tstart","tend","nmatch","alen","mapq","rest"}
{
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$0
}' "$INPUT_PAF" > "$CSV_OUTPUT"

# Process CSV to split id into umi, barcode, and id columns
awk -F ',' 'BEGIN {OFS=","}
NR==1 {print "id","umi","barcode", $0; next}
NR>1 {
    split($1, parts, /\|/)
    id_part = parts[1]
    umi_part = ""
    barcode_part = ""
    for (i=2; i<=length(parts); i++) {
        if (parts[i] ~ /^umi:/) umi_part = substr(parts[i], 5)
        else if (parts[i] ~ /^barcode:/) barcode_part = substr(parts[i], 9)
    }
    print id_part, umi_part, barcode_part, $0
}' "$CSV_OUTPUT" > "$PROCESSED_OUTPUT"

echo "Processing complete. Outputs:"
echo "  $CSV_OUTPUT"
echo "  $PROCESSED_OUTPUT"
