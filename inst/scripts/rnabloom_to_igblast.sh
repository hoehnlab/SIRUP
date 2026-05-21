#!/usr/bin/env bash
set -euo pipefail

# Usage
if [ "$#" -lt 2 ]; then
  echo "Usage: $(basename "$0") <input.fa.gz|input.fq.gz> <output_dir> [nthreads] [dry_run]"
  exit 1
fi

INPUT_GZ="$1"
OUTPUT_DIR="$2"
NUM_THREADS="${3:-96}"
DRY_RUN="${4:-FALSE}"

# Prefer local IgBLAST installs; fall back to shared copy
# Set IGBLAST_HOME in the env to force a specific install.
if [ -z "${IGBLAST_HOME:-}" ]; then
  IGBLAST_CANDIDATES=(
    "/Users/cmv/Desktop/igblast"
    "/dartfs-hpc/rc/home/3/f005f43/igblast"
    "/dartfs/rc/home/3/f005f43/igblast"
  )

  for candidate in "${IGBLAST_CANDIDATES[@]}"; do
    if [ -x "$candidate/bin/igblastn" ] \
      && [ -d "$candidate/database/Homo_sapiens_imgt" ] \
      && [ -f "$candidate/optional_file/human_gl.aux" ]; then
      IGBLAST_HOME="$candidate"
      break
    fi
  done
fi

if [ -z "${IGBLAST_HOME:-}" ]; then
  echo "Error: igblast installation not found. Checked: ${IGBLAST_CANDIDATES[*]}" >&2
  exit 1
fi

echo "Using IgBLAST from: $IGBLAST_HOME"

IGBLAST_BIN="$IGBLAST_HOME/bin/igblastn"

# Prefer local seqkit/csvtk; fall back to shared copies if missing
SEQKIT_BIN="/Users/cmv/miniconda3/bin/seqkit"
CSVTK_BIN="/Users/cmv/miniconda3/bin/csvtk"
if [ ! -x "$SEQKIT_BIN" ]; then
  SEQKIT_BIN="$(command -v seqkit || true)"
fi
if [ ! -x "$CSVTK_BIN" ]; then
  CSVTK_BIN="$(command -v csvtk || true)"
fi
if [ -z "$SEQKIT_BIN" ] || [ ! -x "$SEQKIT_BIN" ]; then
  echo "Error: seqkit not found on PATH or at /Users/cmv/miniconda3/bin/seqkit" >&2
  exit 1
fi
if [ -z "$CSVTK_BIN" ] || [ ! -x "$CSVTK_BIN" ]; then
  echo "Error: csvtk not found on PATH or at /Users/cmv/miniconda3/bin/csvtk" >&2
  exit 1
fi

# Ig/TCR DB paths
BCR_V_DB="$IGBLAST_HOME/database/Homo_sapiens_imgt/IG_V"
BCR_D_DB="$IGBLAST_HOME/database/Homo_sapiens_imgt/IG_D"
BCR_J_DB="$IGBLAST_HOME/database/Homo_sapiens_imgt/IG_J"
BCR_C_DB="$IGBLAST_HOME/database/Homo_sapiens_imgt/IG_C"

TCR_V_DB="$IGBLAST_HOME/database/Homo_sapiens_imgt/TR_V"
TCR_D_DB="$IGBLAST_HOME/database/Homo_sapiens_imgt/TR_D"
TCR_J_DB="$IGBLAST_HOME/database/Homo_sapiens_imgt/TR_J"

AUX_DATA="$IGBLAST_HOME/optional_file/human_gl.aux"

for required in "$IGBLAST_BIN" "$AUX_DATA" "$BCR_V_DB" "$BCR_D_DB" "$BCR_J_DB" "$BCR_C_DB" "$TCR_V_DB" "$TCR_D_DB" "$TCR_J_DB"; do
  if [ ! -e "$required" ]; then
    echo "Error: expected IgBLAST resource missing: $required" >&2
    exit 1
  fi
done

# Outputs
mkdir -p "$OUTPUT_DIR"
BCR_OUTPUT="$OUTPUT_DIR/igblast_bcr.tsv"
TCR_OUTPUT="$OUTPUT_DIR/igblast_tcr.tsv"
BCR_CLONOTYPE="$OUTPUT_DIR/igblast_bcr_clonotypes.csv"
TCR_CLONOTYPE="$OUTPUT_DIR/igblast_tcr_clonotypes.csv"
BCR_FILTERED="$OUTPUT_DIR/filtered_igblast_bcr_complete.tsv"
TCR_FILTERED="$OUTPUT_DIR/filtered_igblast_tcr_complete.tsv"

# ---- Build the streaming command with seqkit ----
# - If FASTQ, convert to FASTA with fq2fa
# - If FASTA, normalize wrapping to 0 (one-line seqs)
# - Add dry-run sampling with seqkit head if requested
ext="${INPUT_GZ##*.}"
# naive sniff: if name contains "fq" or "fastq" → FASTQ, else treat as FASTA
if [[ "$INPUT_GZ" =~ \.f(ast)?q(\.gz)?$ ]]; then
  DECOMP_CMD="$SEQKIT_BIN fq2fa -j $NUM_THREADS \"$INPUT_GZ\""
else
  DECOMP_CMD="$SEQKIT_BIN seq -j $NUM_THREADS -w 0 \"$INPUT_GZ\""
fi

if [[ "$DRY_RUN" == "TRUE" ]]; then
  SAMPLE_CMD="$SEQKIT_BIN head -n 1000"
else
  SAMPLE_CMD="cat"  # no-op passthrough
fi

PIPE_CMD="$DECOMP_CMD | $SAMPLE_CMD"

# ---- IgBLAST BCR ----
echo "Running IgBLAST (BCR) with $NUM_THREADS threads..."
bash -lc "cd \"$IGBLAST_HOME\" && { $PIPE_CMD; } | \"$IGBLAST_BIN\" \
    -ig_seqtype Ig \
    -organism human \
    -evalue 1e-4 \
    -germline_db_V \"$BCR_V_DB\" \
    -germline_db_D \"$BCR_D_DB\" \
    -germline_db_J \"$BCR_J_DB\" \
    -c_region_db  \"$BCR_C_DB\" \
    -auxiliary_data \"$AUX_DATA\" \
    -outfmt 19 \
    -out \"$BCR_OUTPUT\" \
    -clonotype_out \"$BCR_CLONOTYPE\" \
    -num_threads \"$NUM_THREADS\""

# ---- IgBLAST TCR ----
echo "Running IgBLAST (TCR) with $NUM_THREADS threads..."
bash -lc "cd \"$IGBLAST_HOME\" && { $PIPE_CMD; } | \"$IGBLAST_BIN\" \
    -ig_seqtype TCR \
    -organism human \
    -evalue 1e-4 \
    -germline_db_V \"$TCR_V_DB\" \
    -germline_db_D \"$TCR_D_DB\" \
    -germline_db_J \"$TCR_J_DB\" \
    -auxiliary_data \"$AUX_DATA\" \
    -outfmt 19 \
    -out \"$TCR_OUTPUT\" \
    -clonotype_out \"$TCR_CLONOTYPE\" \
    -num_threads \"$NUM_THREADS\""

# ---- Filter complete VDJ with csvtk ----
echo "Filtering for complete VDJ..."
"$CSVTK_BIN" filter2 -j "$NUM_THREADS" -E -T -t -f '$complete_vdj == "T"' "$BCR_OUTPUT" > "$BCR_FILTERED"
"$CSVTK_BIN" filter2 -j "$NUM_THREADS" -E -T -t -f '$complete_vdj == "T"' "$TCR_OUTPUT" > "$TCR_FILTERED"

# ---- Split id into id/umi/barcode in filtered outputs ----
# Expected id format: id|umi:XXXX|barcode:YYYY
if [ -f "$BCR_FILTERED" ]; then
  awk -F '\t' 'BEGIN {OFS="\t"} NR==1 {print "id","umi","barcode", $0; next}
  {
    split($1, parts, /\|/)
    id_part = parts[1]
    umi_part = ""
    barcode_part = ""
    for (i=2; i<=length(parts); i++) {
      if (parts[i] ~ /^umi:/) umi_part = substr(parts[i], 5)
      else if (parts[i] ~ /^barcode:/) barcode_part = substr(parts[i], 9)
    }
    print id_part, umi_part, barcode_part, $0
  }' "$BCR_FILTERED" > "${BCR_FILTERED}_processed"
  mv "${BCR_FILTERED}_processed" "$BCR_FILTERED"
fi

if [ -f "$TCR_FILTERED" ]; then
  awk -F '\t' 'BEGIN {OFS="\t"} NR==1 {print "id","umi","barcode", $0; next}
  {
    split($1, parts, /\|/)
    id_part = parts[1]
    umi_part = ""
    barcode_part = ""
    for (i=2; i<=length(parts); i++) {
      if (parts[i] ~ /^umi:/) umi_part = substr(parts[i], 5)
      else if (parts[i] ~ /^barcode:/) barcode_part = substr(parts[i], 9)
    }
    print id_part, umi_part, barcode_part, $0
  }' "$TCR_FILTERED" > "${TCR_FILTERED}_processed"
  mv "${TCR_FILTERED}_processed" "$TCR_FILTERED"
fi

echo "Done. Outputs:"
echo "  BCR: $BCR_OUTPUT"
echo "  TCR: $TCR_OUTPUT"
echo "  BCR (filtered): $BCR_FILTERED"
echo "  TCR (filtered): $TCR_FILTERED"
