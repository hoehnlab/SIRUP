# writeFastqas

## Write sequences to a FASTQ or FASTA file

## Description

`writeFastqas` writes a table of sequences to a FASTQ or FASTA file.

## Usage

``` R
writeFastqas(
  df,
  fn,
  type = "fq",
  header_col = "id",
  append = FALSE,
  nthreads = 1
)
```

## Arguments

`df`

:   Data frame or data.table with a `seq` column. For FASTQ, it also needs a `qual` column. The `header_col` column holds the record headers.

`fn`

:   Path of the output file.

`type`

:   File type: `"fq"` (default), `"fastq"`, `"fa"`, or `"fasta"`.

`header_col`

:   Name of the column that holds the record headers (default `"id"`). Use `"full_id"` to write the full header.

`append`

:   Add the records to an existing file (default `FALSE`)?

`nthreads`

:   Number of threads for `data.table::fwrite()` (default 1).

## Value

`fn`, invisibly.

## See Also

[`processSpatialSample()`](processSpatialSample.md)

## Examples

``` R
reads <- data.frame(id = c("read1", "read2"), seq = c("ACGTACGT", "TTGACCGA"))
path <- tempfile(fileext = ".fasta")
writeFastqas(reads, path, type = "fasta")
readLines(path)
```
