# Real 500-read equivalence check

This opt-in check compares the reviewed `processSpatialSample()` implementation
against the manual workflow plus the old six-command pRESTO pipeline from SIRAP
commit `bb1b912ebc97f9daadd0c338d0d30031b1d1e1b3`. It calls real IgBLAST,
pRESTO, Change-O, vsearch and MUSCLE. It does not download tools, references or
reads, and it is separate from the mocked tests run during ordinary package checks.

Use FASTA exported by rrad with `to_fasta = TRUE, header_format = "presto"`.
Headers must retain `|BARCODE=...` and may contain `|UMI=...`, before any whitespace.
Run the following from the SIRAP checkout, replacing paths with your installed
tools, reference databases and prepared reads. The baseline commit must exist
in the checkout's Git history. `--sirap-repo .` uses `pkgload` to load the code
under review, including uncommitted fixes; omit it to test installed SIRAP.

```sh
Rscript inst/scripts/check_rrad_equivalence.R \
  --sirap-repo . --baseline-repo . \
  --input /data/rrad_prepared.fasta.gz \
  --output-dir /data/validation/rrad_ig_500 \
  --igblast-base /software/igblast \
  --igblast-db-dir /references/raw_igblast \
  --aux-file /software/igblast/optional_file/human_gl.aux \
  --germline-glob '/references/imgt/IG*.fasta' \
  --muscle /software/muscle3 \
  --receptor ig --organism human --threads 1
```

The output directory must not already exist. IgBLAST and MUSCLE paths are
explicit; `ClusterSets.py`, `ParseHeaders.py`, `AlignSets.py`,
`BuildConsensus.py`, `AssignGenes.py`, `MakeDb.py` and `vsearch` must be on `PATH`.
The raw databases default to `IG_V`, `IG_D`, `IG_J`, `IG_C` for `ig`, or
`TR_V`, `TR_D`, `TR_J` for `tr`. Per-database path overrides and an IgBLAST
executable override are available through `--help`. `igblast_base` must also
contain the databases and internal files expected by `AssignGenes.py` for the
second IgBLAST pass. Use the matching IMGT FASTA references for `MakeDb.py`.

The runner selects the first 500 complete FASTA records in input order, without
altering the input or collapsing duplicate sequences. Fewer than 500 records or
duplicate bare read identifiers cause failure. To test a particular 500-read
cohort, provide that cohort as the input. Single-threaded processing is the
default for a reproducible correctness check; `--threads` can exercise the same
comparison with parallel processing.

Both arms use the same selected input and reference databases. The legacy arm
calls raw `igblastn -outfmt 19` directly, independently parses the canonical
barcode/UMI tags, applies `complete_vdj` filtering and the old `filterAirr()`
with count filtering disabled, and writes `sequence_alignment` to spatial
FASTA. It then calls the old `processSpatialSample(input_fasta = ...)` from the
pinned commit. The current arm calls the reviewed
`processSpatialSample(input = ..., filter_args = list(min_consensus_count = NULL))`
from the same starting FASTA. This checks both IgBLAST passes and all six
consensus-processing commands, including the explicit `NULL` override.

Four checkpoints must match: raw AIRR, filtered AIRR, spatial FASTA and final
`db_pass`. The comparison preserves every column, value and row multiplicity;
it ignores row order and FASTA line wrapping. No scientific fields, barcode
fields, read identifiers or generated cluster identifiers are removed to make
a comparison pass. Empty results fail. External-tool failures and mismatches
produce a nonzero exit code rather than a skipped or successful check.

The result directory contains both arms' outputs, `equivalence.tsv`, the
selected read IDs, `reads_500.fasta`, and a manifest with Git references,
working-tree status, package/tool versions, parameters, original input path,
size and modification time, the selected 500-read checksum, and IMGT reference
checksums. The original input is not scanned in full just to compute a digest;
the saved subset and checksum identify the exact reads used by both arms.
Failures after setup are recorded in the manifest.
Preserve the directory as review evidence. Check that any exported data and
identifiers are suitable for sharing before attaching them to a public issue.

Repeat with a real TCR cohort using `--receptor tr` and the corresponding TR
germline FASTA glob. A successful BCR run does not establish TCR equivalence.
If the cohort contains mixed BCR/TCR reads, run each receptor separately and
report the selected and retained counts; each final result must be nonempty.
The runner validates the FASTA handoff and downstream processing; it does not
establish demultiplexing accuracy or general equivalence across datasets.
