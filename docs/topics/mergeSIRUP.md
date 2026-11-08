# mergeSIRUP

## Merge SIRUP objects into one multi-sample object

## Description

`mergeSIRUP` combines several SIRUP objects into one multi-sample SIRUP object.

## Usage

``` R
mergeSIRUP(
  ...,
  samples,
  sample_field = "sample_id",
  keep_analyses = TRUE,
  validate = TRUE
)
```

## Arguments

`...`

:   Two or more SIRUP objects to merge.

`samples`

:   Character vector of sample names. Give one unique, non-empty name for each SIRUP object.

`sample_field`

:   Name of the receptor-table column that holds the sample name (default `"sample_id"`). The function creates the column or checks it.

`keep_analyses`

:   Keep the stored analyses (default `TRUE`)? Each analysis becomes a named list with one entry for each name in `samples`.

`validate`

:   Validate the merged object (default `TRUE`)?

## Details

Each input object becomes one spatial sample, named by `samples`. The function does these steps:

- It sets `sample_field` in each receptor table, then binds the rows.

- It stores the Seurat objects in a named list.

- It adds the sample name to each stored polygon (if any) and keeps the polygons.

- It stores each analysis as a list with one entry per sample, if `keep_analyses` is `TRUE`.

## Value

A multi-sample SIRUP object.

## See Also

[`createSIRUP()`](createSIRUP.md), [`subsetSIRUP()`](subsetSIRUP.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
merged <- mergeSIRUP(x, x, samples = c("a", "b"))
spatialSamples(merged)
```
