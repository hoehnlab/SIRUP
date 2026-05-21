#!/bin/bash
set -euo pipefail

# Usage:
#   generate_ig_mappings.sh [ig|filter] <rad|blaze|raw> <input_fastq.gz> [output_csv]
#
# Examples:
#   generate_ig_mappings.sh rad reads.fastq.gz
#   generate_ig_mappings.sh filter rad reads.fastq.gz

MODE="ig"
if [[ "${1:-}" == "ig" || "${1:-}" == "filter" ]]; then
  MODE="$1"
  shift
fi

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 [ig|filter] <rad|blaze|raw> <input_fastq.gz> [output_csv]" >&2
  exit 1
fi

TYPE=$1
INPUT_FASTQ=$2
OUTPUT_CSV=${3:-""}

if ! command -v minimap2 >/dev/null 2>&1; then
  echo "Error: minimap2 not found on PATH." >&2
  exit 1
fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "Error: samtools not found on PATH." >&2
  exit 1
fi

# Choose decompression command
if [[ "$INPUT_FASTQ" == *.gz ]]; then
  if command -v gzcat >/dev/null 2>&1; then
    CAT_CMD="gzcat"
  else
    CAT_CMD="zcat"
  fi
else
  CAT_CMD="cat"
fi

# UUID extraction logic based on type
if [ "$TYPE" == "rad" ]; then
  UUID_EXTRACT='if ($F[0] =~ /^([a-f0-9\-]+)[\-][FR]/) { $uuid = $1; } else { $uuid = "NA"; }'
elif [ "$TYPE" == "blaze" ]; then
  UUID_EXTRACT='if ($F[0] =~ /\#([a-f0-9\-]+)_/) { $uuid = $1; } else { $uuid = "NA"; }'
elif [ "$TYPE" == "raw" ]; then
  UUID_EXTRACT='if ($F[0] =~ /^([a-f0-9\-]+)/) { $uuid = $1; } else { $uuid = "NA"; }'
else
  echo "Error: Unrecognized type '$TYPE'. Must be 'rad', 'blaze', or 'raw'." >&2
  exit 1
fi

if [ "$MODE" == "filter" ]; then
  REF_FASTA=~/Desktop/filter_seqs.fa
  OUTPUT_CSV=${OUTPUT_CSV:-"$(dirname "$INPUT_FASTQ")/filter_mapping_results.csv"}
  STATUS_LOGIC='    if ($flag & 4) {
        $status = "unmapped";
    } elsif ($ref =~ /_RR/) {
        $status = "ribosomal";
    } elsif ($ref =~ /MT_DNA/) {
        $status = "mitochondrial";
    } else {
        $status = "unmapped";
    }'
  DEDUP='    next if $seen{$uuid}++;'
else
  REF_FASTA=~/Desktop/hg38_bcrtcr.fa
  OUTPUT_CSV=${OUTPUT_CSV:-"$(dirname "$INPUT_FASTQ")/mapping_results.csv"}
  STATUS_LOGIC='    if ($flag & 4) {
        $status = "unmapped";
    } elsif ($ref =~ /^TR[ABDG]/) {
        $status = "mapped_tcr";
    } elsif ($ref =~ /^IG[HKL]/) {
        $status = "mapped_bcr";
    } else {
        $status = "unmapped";
    }'
  DEDUP=''
fi

if [ ! -f "$REF_FASTA" ]; then
  echo "Error: Reference file '$REF_FASTA' does not exist." >&2
  exit 1
fi

# Escape $ for perl so shell doesn't expand
UUID_EXTRACT_ESC=${UUID_EXTRACT//$/\\$}
STATUS_LOGIC_ESC=${STATUS_LOGIC//$/\\$}
DEDUP_ESC=${DEDUP//$/\\$}

$CAT_CMD "$INPUT_FASTQ" \
| minimap2 -ax map-ont "$REF_FASTA" - \
| samtools view -h - \
| perl -lane "BEGIN { print \"uuid,status\"; %seen = (); }
    next if /^@/;
    $UUID_EXTRACT_ESC
$DEDUP_ESC
    \$ref = \$F[2];
    \$flag = \$F[1];
$STATUS_LOGIC_ESC
    print \"\\$uuid,\\$status\";" \
> "$OUTPUT_CSV"

echo "Mapping summary written to: $OUTPUT_CSV"
