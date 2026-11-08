# spatialCloneSpread

## Summarize the spatial spread of clones

## Description

`spatialCloneSpread` gives the distances between the spots of each clone in each sample.

## Usage

``` R
spatialCloneSpread(
  object,
  receptors = NULL,
  clone = NULL,
  spatial_id = NULL,
  clone_val = NULL,
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  min_spatial_ids = 2L,
  collapse_within_spot = FALSE
)
```

## Arguments

`object`

:   SIRUP object, or Seurat object with spatial images. For a Seurat object, give `receptors`.

`receptors`

:   AIRR-style receptor table. Omit it if `object` is a SIRUP object.

`clone`

:   Name of the column that holds the clone IDs. If `NULL`, the function uses the stored clone field of a SIRUP object, if there is one.

`spatial_id`

:   Name of the column that holds the spatial IDs. If `NULL`, the function uses the stored spatial ID field of a SIRUP object.

`clone_val`

:   Clone IDs to keep.

`images`

:   Names of the images or samples to include.

`image_scale`

:   Image scale for the Seurat image coordinates.

`swap_coords`

:   Swap the x and y coordinates before the function calculates distances?

`min_spatial_ids`

:   Smallest number of observations that a clone needs in a sample (default 2). If `collapse_within_spot` is `TRUE`, the function counts distinct spots. Otherwise, it counts sequences.

`collapse_within_spot`

:   Merge the sequences that share a spatial ID into one point before the function calculates distances (default `FALSE`)? If `FALSE`, the function counts every sequence, so sequences in the same spot give pairs with distance zero. Use `TRUE` to measure the spread over distinct spots only.

## Details

For each sample and clone, the function calculates the distance between every pair of points. It reports the mean and the maximum distance.

## Value

A tibble with one row for each sample and clone. The columns are `sample`, `clone_id`, `n_spatial_ids`, `mean_pairwise_distance`, and `max_pairwise_distance`.

## See Also

[`spatialCloneJoinCount()`](spatialCloneJoinCount.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
spatialCloneSpread(x)
```
