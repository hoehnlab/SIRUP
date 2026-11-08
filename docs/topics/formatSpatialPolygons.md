# formatSpatialPolygons

## Format point sets as polygon vertices

## Description

`formatSpatialPolygons` converts a table of points into the vertex table that [`selectSpatialPolygons()`](selectSpatialPolygons.md) returns.

## Usage

``` R
formatSpatialPolygons(
  object,
  polygon_points,
  spatial_id = "spatial_id",
  polygon = "polygon",
  x = "x",
  y = "y",
  image = "image",
  point_type = c("interior", "boundary"),
  buffer_radius = NULL,
  order = NULL,
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  close_polygons = TRUE
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`polygon_points`

:   Data frame of points. With `point_type = "interior"`, it holds all points inside each polygon. With `point_type = "boundary"`, it holds the boundary points.

`spatial_id`

:   Name of the column in `polygon_points` that holds the spatial IDs (default `"spatial_id"`). You need it if `polygon_points` has no `x` and `y` columns.

`polygon`

:   Name of the column in `polygon_points` that identifies the polygons (default `"polygon"`).

`x`

:   Name of the x-coordinate column in `polygon_points`, if it exists (default `"x"`).

`y`

:   Name of the y-coordinate column in `polygon_points`, if it exists (default `"y"`).

`image`

:   Name of the image column in `polygon_points`, if it exists (default `"image"`).

`point_type`

:   Type of the points in `polygon_points`: `interior` (all points inside) or `boundary`.

`buffer_radius`

:   Radius for the buffer around the points. If `NULL`, the function estimates it from the spacing between nearest neighbors in the full coordinate grid.

`order`

:   Name of a column that orders the points in each polygon. It is most useful with `point_type = "boundary"`.

`images`

:   Names of the Seurat images to use.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object.

`close_polygons`

:   Close each polygon by repeating the first vertex, if needed (default `TRUE`)?

## Details

By default, `polygon_points` holds all points inside each polygon. The function merges the points with `sf`. Then it maps the outer ring back to the nearest selected spatial ID coordinates. The boundary keeps a valid order for concave polygons, and the vertices are spatial ID coordinates. If the input already holds boundary points, use `point_type = "boundary"`.

## Value

A data frame with the columns `spatial_id`, `x`, `y`, `polygon`, and `image`. The function drops polygons with fewer than three points and gives a warning.

## See Also

[`selectSpatialPolygons()`](selectSpatialPolygons.md), [`addPolygons()`](sirup-polygons.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Four spots that lie inside one polygon.
points <- data.frame(
  spatial_id = c("spot04", "spot05", "spot09", "spot10"),
  polygon = "block"
)
formatSpatialPolygons(x, points, image_scale = "hires", swap_coords = FALSE)
```
