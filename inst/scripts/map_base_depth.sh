#!/bin/bash

# Usage: ./map_base_depth.sh ref.fa reads.fq.gz output.csv

set -euo pipefail

REF=${1:-}
READS=${2:-}
OUTPUT=${3:-depth.csv}
THREADS=${THREADS:-4}

if [ -z "$REF" ] || [ -z "$READS" ]; then
  echo "Usage: $0 <ref.fa> <reads.fq.gz> [output.csv]" >&2
  exit 1
fi

if ! command -v minimap2 >/dev/null 2>&1; then
  echo "Error: minimap2 not found on PATH." >&2
  exit 1
fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "Error: samtools not found on PATH." >&2
  exit 1
fi

BAM_OUT=$(mktemp -t map_base_depth.XXXXXX.bam)

minimap2 -ax splice -C0 -06,24 --splice-flank=no -t "$THREADS" "$REF" "$READS" \
  | samtools view -@ "$THREADS" -bS - \
  | samtools sort -@ "$THREADS" -o "$BAM_OUT"

samtools index "$BAM_OUT"

samtools depth -a "$BAM_OUT" \
  | awk 'BEGIN {OFS=","; print "ref,pos,count"} {print $1, $2, $3}' > "$OUTPUT"

rm -f "$BAM_OUT" "$BAM_OUT.bai"

echo "Coverage written to: $OUTPUT"
