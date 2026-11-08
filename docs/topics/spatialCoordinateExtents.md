# spatialCoordinateExtents

## Summarize spatial coordinate ranges

## Description

`spatialCoordinateExtents` gives the x and y range of each sample or image.

## Usage

``` R
spatialCoordinateExtents(
  object,
  by = c("sample", "image"),
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL
)
```

## Arguments

`object`

:   SIRUP object, or Seurat object with spatial images.

`by`

:   Group the ranges by `sample` or by `image`.

`images`

:   Names of the images or samples to include.

`image_scale`

:   Image scale for the Seurat image coordinates.

`swap_coords`

:   Swap the x and y coordinates?

## Value

A tibble with the ranges (`x_min`, `x_max`, `y_min`, `y_max`), the `width`, `height`, and `diagonal`, and the number of spots (`n_spatial_ids`).

## See Also

[`spatialCoordinateTable()`](spatialCoordinateTable.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
spatialCoordinateExtents(x)
```
