# seuratSubsetPolygon

## Subset a Seurat object to the spots inside polygons

## Description

`seuratSubsetPolygon` keeps only the spatial IDs inside the polygons. It also records the polygon of each spot in the metadata.

## Usage

``` R
seuratSubsetPolygon(
  object,
  polygon = NULL,
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  polygon_val = NULL,
  x = "x",
  y = "y",
  polygon_id = "polygon",
  metadata = "polygon",
  include_boundary = TRUE
)
```

## Arguments

`object`

:   Seurat object, or a named list of Seurat objects.

`polygon`

:   Data frame of polygon vertices in order. It needs the columns set by `x` and `y`. If the column set by `polygon_id` exists, the function treats each polygon on its own.

`images`

:   Names of the Seurat images to use.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`.

`swap_coords`

:   Swap the two Seurat coordinate axes? You must give a value.

`polygon_val`

:   IDs of the polygons to use. If `NULL`, the function uses all polygons.

`x`

:   Name of the x-coordinate column in `polygon` (default `"x"`).

`y`

:   Name of the y-coordinate column in `polygon` (default `"y"`).

`polygon_id`

:   Name of the column in `polygon` that groups the vertices by polygon (default `"polygon"`).

`metadata`

:   Name of the metadata column that holds the primary polygon of each spot (default `"polygon"`). A spot in overlapping polygons gets the smallest polygon by area.

`include_boundary`

:   Count spots on a polygon edge as inside the polygon (default `TRUE`)?

## Value

A Seurat object with only the spatial IDs inside the selected polygons. For a list of Seurat objects, a list of the subset objects. The function drops an object with no spot inside a polygon. It adds the `metadata` column to `⁠@meta.data⁠`.

## See Also

[`assignSpatialPolygons()`](assignSpatialPolygons.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
flat <- data.frame(
  polygon = "block",
  x = c(3.5, 5.5, 5.5, 3.5),
  y = c(0.5, 0.5, 2.5, 2.5)
)

# Keep the four spots inside the polygon.
subset <- seuratSubsetPolygon(
  getSpatial(x),
  polygon = flat,
  image_scale = "hires",
  swap_coords = FALSE
)
ncol(subset)
```
