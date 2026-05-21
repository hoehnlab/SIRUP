#!/bin/bash

# Usage: ./extract_blaze_headers.sh input.fastq.gz output.csv

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fastq> <output_csv>" >&2
    exit 1
fi

if ! command -v seqkit >/dev/null 2>&1; then
    echo "Error: seqkit not found on PATH." >&2
    exit 1
fi

input_fastq=$1
output_csv=$2

seqkit seq -n "$input_fastq" \
  | sed 's/^@//' \
  | awk -F'[_#]' '
      BEGIN { OFS=","; print "id","barcode","umi" }
      { print $3, $1, $2 }
    ' \
  > "$output_csv"

echo "BLAZE headers extracted to $output_csv"
