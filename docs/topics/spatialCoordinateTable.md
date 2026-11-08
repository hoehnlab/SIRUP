# spatialCoordinateTable

## Get spatial plotting coordinates

## Description

`spatialCoordinateTable` returns the coordinate table that SIRUP uses for spatial statistics and plots.

## Usage

``` R
spatialCoordinateTable(
  object,
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL
)
```

## Arguments

`object`

:   Seurat object.

`images`

:   Names of the images to include.

`image_scale`

:   Image scale for the Seurat image coordinates.

`swap_coords`

:   Swap the x and y coordinates?

## Details

The function takes the coordinates from the Seurat images. If there are none, it takes them from the coordinate columns in the Seurat metadata. The `.sirup_y` column holds the y-coordinate that SIRUP uses for plots.

## Value

A data frame with at least these columns: `cell`, `x`, `y`, `.sirup_y`, and `image`.

## See Also

[`spatialCoordinateExtents()`](spatialCoordinateExtents.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
head(spatialCoordinateTable(getSpatial(x), image_scale = "hires", swap_coords = FALSE))
```
