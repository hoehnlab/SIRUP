#!/usr/bin/env bash
# extract_flexiplex_headers
# Usage: ./extract_flexiplex_headers.sh <input.fastq[.gz]> <output.csv>

set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_fastq> <output_csv>" >&2
  exit 1
fi

if ! command -v seqkit >/dev/null 2>&1; then
  echo "Error: seqkit not found on PATH." >&2
  exit 1
fi

input_fastq="$1"
output_csv="$2"

seqkit seq -n "$input_fastq" \
| sed 's/^@//' \
| awk -F'\t' '
BEGIN {
  OFS=","; print "id","barcode","umi"
}
{
  # Field 1 looks like: <...>#<ID>_<...>
  hdr = $1
  id  = hdr
  n = split(hdr, a, /#/) 
  if (n >= 2) {
    rest = a[2]
    m = split(rest, b, /_/) 
    if (m >= 1) id = b[1]
  }

  cb = ""; ub = ""
  for (i = 2; i <= NF; i++) {
    if ($i ~ /^CB:Z:/) cb = substr($i, 6)
    else if ($i ~ /^UB:Z:/) ub = substr($i, 6)
  }

  print id, cb, ub
}' > "$output_csv"

echo "Headers extracted to $output_csv"
