#!/bin/bash

# Reformat FASTQ headers to embed barcode in read name
# Usage: ./reformat_headers_for_minimap.sh input.fastq.gz output.fastq.gz

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <input.fastq.gz> <output.fastq.gz>" >&2
    echo "Reformats headers from:" >&2
    echo "  @read_id CB:Z:BARCODE ..." >&2
    echo "To:" >&2
    echo "  @read_id-BARCODE" >&2
    exit 1
fi

if ! command -v seqkit >/dev/null 2>&1; then
    echo "Error: seqkit not found on PATH." >&2
    exit 1
fi

INPUT=$1
OUTPUT=$2

seqkit replace -p "^(\S+)\s+CB:Z:(\w+).*" -r '${1}-${2}' "$INPUT" -o "$OUTPUT"
