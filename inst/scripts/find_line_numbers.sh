#!/usr/bin/env bash
# Usage: ./find_line_numbers_simple.sh file.gz ids.txt

gzfile="$1"
idsfile="$2"

while IFS= read -r id; do
  [[ -z "$id" ]] && continue
  zgrep -n "$id" "$gzfile" | awk -F: -v id="$id" '{print id "," $1}'
done < "$idsfile"
