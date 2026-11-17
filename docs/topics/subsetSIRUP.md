# subsetSIRUP

## Subset a SIRUP object

## Description

`subsetSIRUP` keeps selected samples, spatial IDs, and receptor rows of a SIRUP object. `filterSample` keeps only the samples you name.

## Usage

``` R
subsetSIRUP(
  x,
  sample_id = NULL,
  spatial_id = NULL,
  receptor_rows = NULL,
  slot = NULL,
  keep_nonspatial_receptors = FALSE,
  keep_analyses = FALSE,
  drop_empty_spatial = TRUE,
  validate = TRUE
)

filterSample(x, sample_id, keep_analyses = FALSE, validate = TRUE)
```

## Arguments

`x`

:   A SIRUP object.

`sample_id`

:   Name or names of the spatial samples to keep. Use it only if the object holds a named list of Seurat objects.

`spatial_id`

:   Spatial IDs (barcodes) to keep. In a multi-sample object, the function matches the IDs in each selected sample.

`receptor_rows`

:   Logical vector with length `nrow(getReceptors(x, slot = slot))`, or a function. The function takes the receptor tibble of `slot` and returns such a vector. The function keeps the rows where the vector is `TRUE`. It applies this filter to `slot` only. The other filters apply to every slot.

`slot`

:   Name of the receptor slot that `receptor_rows` applies to. Leave it `NULL` if the object has one receptor table. You must give it if the object has several slots and you use `receptor_rows`.

`keep_nonspatial_receptors`

:   Keep receptor rows that have no spatial ID after the other filters (default `FALSE`)? The function uses this argument only if you give `spatial_id`.

`keep_analyses`

:   Keep the stored analyses (default `FALSE`)? Most analyses depend on the original spatial and receptor data.

`drop_empty_spatial`

:   Drop samples that have no cells left after the `spatial_id` filter (default `TRUE`)? The function uses this argument only for multi-sample objects.

`validate`

:   Validate the subset object (default `TRUE`)?

## Details

The function subsets the Seurat data, the receptor rows, the stored polygons, and the analyses. It applies the stored polygons again to the subset, so the polygon columns stay consistent.

## Value

A SIRUP object that contains the subset.

## See Also

[`mergeSIRUP()`](mergeSIRUP.md), [`createSIRUP()`](createSIRUP.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Keep three spots.
subsetSIRUP(x, spatial_id = c("spot01", "spot02", "spot03"))

# Keep the receptor rows of clone 2.
subsetSIRUP(x, receptor_rows = function(receptors) receptors$clone_id == "2")

# Keep one sample of a multi-sample object.
merged <- mergeSIRUP(x, x, samples = c("a", "b"))
filterSample(merged, "a")
```
