# 500-record equivalence results: 2026-09-07

Both real-tool comparisons passed. Each used 500 distinct original contig IDs,
prepared by rrad 1.2.1, and compared the old manual workflow with the reviewed
`processSpatialSample()` workflow. No external tools were mocked.

- Legacy processing source: `bb1b912ebc97f9daadd0c338d0d30031b1d1e1b3`.
- Tested review source: `5e2297147c431b184bd8f55abcd17cc2960d2041`, clean worktree.
- Both arms explicitly disabled the minimum-consensus-count filter. The current
  arm used `filter_args = list(min_consensus_count = NULL)`, exercising the fix
  that preserves explicit `NULL` overrides.

## Results

Counts below are identical between the legacy and current arms. PASS means all
columns, values and duplicate-row multiplicities matched, not merely the counts.
Only row/column ordering and FASTA line wrapping were ignored.

| Checkpoint | BCR records, each arm | BCR | TCR records, each arm | TCR |
| --- | ---: | --- | ---: | --- |
| Raw IgBLAST AIRR | 500 | PASS | 500 | PASS |
| Filtered AIRR | 461 | PASS | 396 | PASS |
| Spatial FASTA | 461 | PASS | 396 | PASS |
| Final `db_pass` | 458 | PASS | 396 | PASS |

Final BCR annotations were 458 IGH records. Final TCR annotations were 190 TRA,
203 TRB and 3 TRD records. Both runs completed with exit status zero.

## Inputs and preparation

The BCR input was the first 500 IGH rows of `BCR_data.tsv` from the
[Immcantation 10x tutorial](https://github.com/immcantation/immcantation/blob/master/training/10x_tutorial.md)
and its [data deposit](https://doi.org/10.5281/zenodo.10028129).
We used the original `sequence`, not `sequence_alignment`. Read IDs were
`sample_id__sequence_id`; barcodes were `sample_id__cell_id`, preventing
cross-sample collisions. The separate 477-record tutorial FASTA and replicated
benchmark fixtures were not used to inflate the sample to 500.

The TCR input was the first 500 records of the pre-existing
`tcr_filtered_contig.fasta` benchmark file, which contains 1,198 records.
Barcodes were recovered from original contig identifiers by removing the final
`_contig_N`. Its original download URL was not recorded, so its provenance is
limited to the retained local file and checksum below.

For both cohorts, the original sequences and identifiers were retained. Barcode
metadata was supplied as `CB:Z:` FASTA comments, then converted by the actual
`rrad::rad_reformat(..., to_fasta = TRUE, header_format = "presto")` call.
The resulting `|BARCODE=...` headers, sequences and identifiers were checked
against the input. No UMI values were invented; UMI tags were absent.

SHA-256 identifiers for the exact input files:

```text
BCR source BCR_data.tsv:
5705ffd2929434ffced8fba6c1640385b9597e0cca3ec1879fcf87f6f6bed9f1
BCR prepared 500-record FASTA:
05ba6fbf458b355a504fcc18a5a0fb7915e2d5a76c17e03e5e72a2180c776940
TCR source tcr_filtered_contig.fasta:
be675d9fba7aef9ee63329a1c4943880d20d5c4cd8667d1344046cc92dec1c37
TCR prepared 500-record FASTA:
b0a44d77d4d280ddddf8de526be38e98f737cb65bc3b0c50c2c72f495cf65350
```

## Tools and references

The runs used R 4.5.2, SIRAP 0.0.0.9001, rrad 1.2.1, IgBLAST 1.22.0,
pRESTO 0.7.9, Change-O 1.3.4, MUSCLE 3.8.1551 and vsearch 2.31.0.
The runner's `--threads` setting was 1; MakeDb retained its own default worker
count. This is a correctness test, not a performance benchmark.

Both arms shared the same human IMGT references, prepared using the
[Change-O reference preparation procedure](https://changeo.readthedocs.io/en/1.3.1/examples/igblast.html).
The source was the [IMGT gapped nucleotide aggregate](https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP),
retrieved on the test date (SHA-256
`e71cd5a35e25d523ff8c3ec81e63ab2a2b80780e601260a2349023990f8f5429`).
V/D/J references retained the appropriate regions. Constant-reference resources
used authentic first constant exons (IGH CH1, TR EX1) or light-chain C-REGION,
not concatenated artificial constant transcripts. The raw BCR pass received
the IG constant database; the raw TCR pass did not receive a constant database.

## Reproduction and retained evidence

Use the [runner instructions](rrad-equivalence.md) with the prepared FASTA and
matching references. They document the two independent arms and all four
comparisons. The runner rejects missing records, empty checkpoints, failed tools
and mismatches with a nonzero exit status.

Local evidence from these runs is under
`/tmp/sirap-equivalence-input.7zVFnH/`: `bcr_results/` and `tcr_results/` contain
the selected FASTA, read IDs, both arms' outputs, `equivalence.tsv` and manifests.
`prepare_inputs.R`, `run_equivalence.R`, and `input_provenance.rds` record the
exact preparation and invocation. Per-reference checksums and database-building
commands are in
`/tmp/sirap-equivalence-resources.Ykgjfz/reference_manifest.json`.
These temporary artifacts are local, not included in Git, and may be removed
by operating-system cleanup. Preserve them separately before sharing or
attempting an exact replay; the committed runner and this report persist.

## Limits

These were alternative real contig cohorts, not a replay of the original
Discovery engblom 500-read sample, which was not accessible in this environment.
Every final consensus group had `consensus_count = 1`. Consequently these runs
exercise both IgBLAST passes, the real pRESTO/Change-O commands, barcode handoff,
filtering and output equivalence, but do not validate nontrivial multi-read
consensus construction, UMI-bearing real data or raw-read demultiplexing
accuracy. A multi-read raw cohort, ideally the original requested sample,
remains necessary for those claims. Success here does not imply equivalence on
all datasets or validate fresh-install dependency setup.
