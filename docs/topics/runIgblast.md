# runIgblast

## Annotate raw reads with IgBLAST

## Description

`runIgblast` runs `airrblast::assignGenes()` on raw reads. It writes an AIRR table.

## Usage

``` R
runIgblast(
  input,
  output_tsv,
  receptor = c("ig", "tr"),
  organism = "human",
  igblast_base = "~/software/igblast",
  igblast_refs = NULL,
  num_threads = 8
)
```

## Arguments

`input`

:   Path of a FASTA or FASTQ file, optionally gzipped. A plain FASTA file must end in `.fasta` or `.fa`.

`output_tsv`

:   Path of the AIRR TSV file. It must end in `.tsv`.

`receptor`

:   `"ig"` for BCR or `"tr"` for TCR.

`organism`

:   `"human"`, `"mouse"`, or `"rhesus_monkey"`.

`igblast_base`

:   IgBLAST installation folder. It must hold `bin/igblastn`, `internal_data`, and `optional_file`.

`igblast_refs`

:   Folder with the IgBLAST `⁠database/⁠` and `⁠fasta/⁠` folders (the `refs` argument of `airrblast::assignGenes()`). If `NULL`, the function uses `igblast_base`.

`num_threads`

:   Number of threads, passed to `assignGenes()` as `nproc`.

## Details

The input can be FASTA or FASTQ, optionally gzipped. The function first writes gzipped and FASTQ inputs as a plain FASTA next to `output_tsv`. It drops the FASTQ quality strings, because `igblastn` accepts only FASTA. Each call runs one receptor group: `ig` or `tr`.

`assignGenes()` does not search the C-region database, so the tables of raw reads have no `c_call` from this step.

## Value

`output_tsv`, invisibly.

## See Also

[`processSpatialSample()`](processSpatialSample.md), [`parseRadSequenceIds()`](parseRadSequenceIds.md)

## Examples

``` R
## Not run: 
# This example needs IgBLAST and the IgBLAST reference files.
runIgblast(
  "sample1_reads.fasta.gz",
  output_tsv = "sample1_igblast.tsv",
  receptor = "ig",
  organism = "human",
  igblast_base = "~/software/igblast"
)

## End(Not run)
```
