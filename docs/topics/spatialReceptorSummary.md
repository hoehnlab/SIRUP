# spatialReceptorSummary

## Count the sequences and spots of a receptor table

## Description

`spatialReceptorSummary` returns the basic counts of a receptor table.

## Usage

``` R
spatialReceptorSummary(airr_df, spatial_id = "spatial_id_2um", slot = NULL)
```

## Arguments

`airr_df`

:   AIRR-style data frame, or a SIRUP object.

`spatial_id`

:   Name of the spatial ID column to count (default `"spatial_id_2um"`). For a SIRUP object, if you omit it, the function uses the stored spatial ID field. If the object has a sample field, the function counts each spatial ID once for each sample.

`slot`

:   Name of the receptor slot to use if `airr_df` is a SIRUP object. Leave it `NULL` if the object has one receptor table. You must give it if the object has several slots.

## Value

A one-row tibble with the columns `n_sequences`, `n_spatial_ids`, `n_loci`, and `n_productive`. `n_loci` and `n_productive` are `NA` if the table has no `locus` or `productive` column.

## See Also

[`filterAirr()`](filterAirr.md), [`spatialCloneDominance()`](spatialCloneDominance.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
spatialReceptorSummary(x)
```
