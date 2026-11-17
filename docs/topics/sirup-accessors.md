# sirup-accessors

## Extract data from a SIRUP object

## Description

`getSpatial` returns the Seurat object(s). `spatialSamples` returns the spatial sample names. `getReceptors` returns one receptor table. `receptorSlots` returns the names of the receptor slots.

## Usage

``` R
getSpatial(x, sample_id = NULL)

spatialSamples(x)

getReceptors(x, sample_id = NULL, slot = NULL)

receptorSlots(x)
```

## Arguments

`x`

:   A SIRUP object.

`sample_id`

:   Name of one spatial sample. Use it only if the object holds a named list of Seurat objects.

`slot`

:   Name of the receptor slot. If the object has one receptor table, leave it `NULL` (the slot is called `"only_slot"`). If the object has several slots, you must give one.

## Value

`getSpatial` returns the Seurat object. For a multi-sample object, it returns the named list of Seurat objects, or one Seurat object if you give `sample_id`. `spatialSamples` returns a character vector of sample names. `getReceptors` returns the receptor tibble of one slot. If you give `sample_id`, it returns only the rows for that sample. `receptorSlots` returns a character vector of slot names.

## See Also

[`getPolygons()`](sirup-polygons.md), [`getAnalysis()`](sirup-analyses.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
getSpatial(x)
head(getReceptors(x))
spatialSamples(x)
receptorSlots(x)
```
