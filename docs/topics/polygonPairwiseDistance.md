# polygonPairwiseDistance

## Calculate distances between polygon centroids

## Description

`polygonPairwiseDistance` calculates the distance between the centroids of each pair of polygons in a sample.

## Usage

``` R
polygonPairwiseDistance(
  object,
  polygons = NULL,
  normalize_by_diagonal = TRUE,
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL
)
```

## Arguments

`object`

:   SIRUP object with stored polygons.

`polygons`

:   Table of polygon vertices.

`normalize_by_diagonal`

:   Divide each distance by the diagonal of the coordinate range of its sample (default `TRUE`)?

`images`

:   Names of the images or spatial samples to use for the coordinate range.

`image_scale`

:   Image scale for the Seurat image coordinates.

`swap_coords`

:   Swap the x and y coordinates before the function calculates the coordinate range?

## Value

A tibble with one row for each pair of polygons in a sample. The columns are `sample`, `polygon_id_1`, `polygon_id_2`, `polygon_label_1`, `polygon_label_2`, `distance`, `normalized_distance`, and `section_diagonal`.

## See Also

[`polygonCentroids()`](polygonCentroids.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
polygonPairwiseDistance(x)
```
