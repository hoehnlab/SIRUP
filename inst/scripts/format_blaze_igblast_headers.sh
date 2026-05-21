#!/bin/bash

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <igblast_output.tsv>" >&2
    exit 1
fi

INPUT_TSV="$1"
OUTPUT_TSV="${INPUT_TSV%.tsv}_processed.tsv"

echo "Processing $INPUT_TSV to extract id, umi, and barcode from sequence_id..."

awk -F '\t' 'BEGIN {OFS="\t"}
NR==1 {
    print "id", "umi", "barcode", $0
    next
}
{
    # Expected format: barcode_umi#readid_strand
    split($1, prepost, "#")
    split(prepost[1], bu, "_")
    split(prepost[2], rs, "_")

    barcode = bu[1]
    umi = bu[2]
    id = rs[1]

    print id, umi, barcode, $0
}' "$INPUT_TSV" > "$OUTPUT_TSV"

echo "Finished processing. Output saved to $OUTPUT_TSV."
