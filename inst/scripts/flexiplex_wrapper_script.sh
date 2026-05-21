#!/bin/bash

# flex_runner.sh
# Wrapper for flexiplex and flexiplex-filter with default settings

set -euo pipefail

# Known barcode whitelist
WHITELIST="/Users/cmv/Desktop/rad_paper/cellranger_whitelists/3M-february-2018-3v3.txt"

# Check arguments
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input.fastq[.gz]> [output_prefix]"
    echo "Example: $0 reads.fastq.gz /path/to/output/sample1"
    exit 1
fi

INPUT_FASTQ="$1"
OUTPUT_PREFIX="${2:-barcodes}"
# Normalize trailing slash so suffixes append correctly
OUTPUT_PREFIX="${OUTPUT_PREFIX%/}"
# Fallback if the user passed only a slash
if [[ -z "$OUTPUT_PREFIX" ]]; then
    OUTPUT_PREFIX="barcodes"
fi

# Extract output directory from prefix
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
PREFIX_NAME=$(basename "$OUTPUT_PREFIX")

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "=== Flexiplex Wrapper ==="
echo "Input FASTQ: $INPUT_FASTQ"
echo "Output directory: $OUTPUT_DIR"
echo "Output prefix: $PREFIX_NAME"
echo "Whitelist: $WHITELIST"
echo ""

# Run flexiplex with -k (known barcodes) and -f 0 (no filtering, discovery mode)
echo "[1/2] Running flexiplex (discovery mode)..."
if [[ "$INPUT_FASTQ" == *.gz ]]; then
    echo "[Info] Detected gzipped FASTQ input."
    # Use gzcat on macOS (BSD gzip), otherwise gzip -cd
    if [[ "$(uname -s)" == "Darwin" ]]; then
        gunzip -c "$INPUT_FASTQ" | flexiplex \
        -d 10x3v3 \
        -p 4 \
        -f 0 \
        -n "${OUTPUT_PREFIX}_flexiplex" \
        2>&1 | tee "${OUTPUT_PREFIX}_flexiplex.log"
        fi
else
    flexiplex \
        -f 0 \
        -p 4 \
        -d 10x3v3 \
        -n "${OUTPUT_PREFIX}_flexiplex" \
        "$INPUT_FASTQ" \
        2>&1 | tee "${OUTPUT_PREFIX}_flexiplex.log"
fi

echo ""

# Run flexiplex-filter on the output
# modified to use -u 0 for the filtering step because defaults just don't work
echo "[2/2] Running flexiplex-filter..."
flexiplex-filter -w "$WHITELIST" \
    --verbose \
    -u 0 \
    -l 1 \
    "${OUTPUT_PREFIX}_flexiplex_barcodes_counts.txt" \
    > "${OUTPUT_PREFIX}_whitelist.txt" \
    2>&1 | tee "${OUTPUT_PREFIX}_filter.log"

echo ""
echo "=== Done ==="
echo "Outputs:"
echo "  Barcode counts: ${OUTPUT_PREFIX}_flexiplex_barcodes_counts.txt"
echo "  Filtered whitelist: ${OUTPUT_PREFIX}_whitelist.txt"
echo "  Logs: ${OUTPUT_PREFIX}_flexiplex.log, ${OUTPUT_PREFIX}_filter.log"
