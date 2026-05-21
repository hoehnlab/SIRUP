#!/bin/bash

# Split BAM by barcode (extracted from read names) and generate per‑base depth
# Usage: ./generate_per_base_mappings.sh input.bam output_directory

set -euo pipefail

INPUT_BAM=${1:-}
OUTPUT_DIR=${2:-}
THREADS=${THREADS:-4}

if [ -z "$INPUT_BAM" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <input.bam> <output_directory>" >&2
    echo "Splits BAM by barcode suffix in read names (format: readid-BARCODE)" >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "Error: samtools not found on PATH." >&2
    exit 1
fi

if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM not found: $INPUT_BAM" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Extract unique barcodes
samtools view -@ "$THREADS" "$INPUT_BAM" | \
awk '{
    split($1, parts, "-")
    barcode = parts[length(parts)]
    print barcode
}' | sort | uniq > "$OUTPUT_DIR/barcodes.txt"

BARCODE_COUNT=$(wc -l < "$OUTPUT_DIR/barcodes.txt")

# Split BAM by barcode
counter=0
while read -r barcode; do
    counter=$((counter + 1))
    echo "Processing barcode $counter/$BARCODE_COUNT: $barcode"

    samtools view -@ "$THREADS" -h "$INPUT_BAM" | \
    awk -v bc="$barcode" '
    /^@/ {print; next}
    {
        split($1, parts, "-")
        read_barcode = parts[length(parts)]
        if (read_barcode == bc) print
    }' | \
    samtools view -@ "$THREADS" -bS - > "$OUTPUT_DIR/${barcode}.bam"

    samtools index "$OUTPUT_DIR/${barcode}.bam"

done < "$OUTPUT_DIR/barcodes.txt"

# Generate depth files
counter=0
while read -r barcode; do
    counter=$((counter + 1))
    echo "Generating depth $counter/$BARCODE_COUNT: $barcode"
    samtools depth -aa "$OUTPUT_DIR/${barcode}.bam" | gzip > "$OUTPUT_DIR/${barcode}_depth.txt.gz"

done < "$OUTPUT_DIR/barcodes.txt"

echo "BAM splitting and depth analysis complete. Output files in $OUTPUT_DIR"
