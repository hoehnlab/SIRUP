# rrad FASTA handoff fixture

`canonical.fasta` was generated from `tagged.fastq` with public rrad v1.2.1:

```r
rrad::rad_reformat(
  "tagged.fastq",
  output_fastq = "canonical.fasta",
  to_fasta = TRUE,
  header_format = "presto"
)
```

The second read deliberately has a barcode but no UMI. SIRAP's automated
tests consume this fixed output without installing rrad. External IgBLAST
and pRESTO calls are mocked in the handoff tests; these tests verify metadata
preservation and pipeline ordering, not biological annotation accuracy.
