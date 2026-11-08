# processSpatialSamples

## Process several spatial samples

## Description

`processSpatialSamples` runs [`processSpatialSample()`](processSpatialSample.md) on each input in turn.

## Usage

``` R
processSpatialSamples(inputs, output_dirs = NULL, sample_names = NULL, ...)
```

## Arguments

`inputs`

:   List or character vector of inputs, one for each sample. Each input goes to [`processSpatialSample()`](processSpatialSample.md). It can be a reads file, the path of an AIRR TSV, or an AIRR data frame.

`output_dirs`

:   Vector of output directories, one for each input. If `NULL`, each sample uses the default of [`processSpatialSample()`](processSpatialSample.md).

`sample_names`

:   Vector of sample names, one for each input. If `NULL`, the function uses the file name of each input, or `⁠sample_<i>⁠` for a data frame.

`...`

:   Other arguments for [`processSpatialSample()`](processSpatialSample.md).

## Value

A list with one list of output paths for each sample, invisibly.

## See Also

[`processSpatialSample()`](processSpatialSample.md)

## Examples

``` R
## Not run: 
# This example needs IgBLAST, pRESTO, vsearch, and MUSCLE.
results <- processSpatialSamples(
  c("sample1_reads.fasta.gz", "sample2_reads.fasta.gz"),
  output_dirs = c("sample1_output", "sample2_output"),
  receptor = "ig",
  imgt_dir = "~/share/germlines/imgt"
)

## End(Not run)
```
