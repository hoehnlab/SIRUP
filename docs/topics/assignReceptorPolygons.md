# assignReceptorPolygons

## Add polygon IDs to a receptor table

## Description

`assignReceptorPolygons` adds polygon columns to a receptor table. It keeps all rows.

## Usage

``` R
assignReceptorPolygons(
  receptors,
  object = NULL,
  polygon = NULL,
  column = "polygon_id",
  spatial_id = NULL,
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  polygon_val = NULL,
  x = "x",
  y = "y",
  receptor_x = "x",
  receptor_y = "y",
  polygon_id = "polygon",
  include_boundary = TRUE,
  slot = NULL
)
```

## Arguments

`receptors`

:   AIRR-style receptor table. You can give a SIRUP object here if you omit `object`.

`object`

:   Seurat object or SIRUP object.

`polygon`

:   Data frame of polygon vertices in order. It needs the columns set by `x` and `y`. If the column set by `polygon_id` exists, the function treats each polygon on its own. If `object` is a SIRUP object and `polygon` is `NULL`, the function uses the stored polygons.

`column`

:   Name of the column to add to the receptor table (default `"polygon_id"`).

`spatial_id`

:   Name of the `receptors` column that holds the spatial IDs. If `NULL`, the function uses the stored field of a SIRUP object, or looks for a standard spatial ID column.

`images`

:   Names of the Seurat images to use.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`polygon_val`

:   IDs of the polygons to use. If `NULL`, the function uses all polygons.

`x`

:   Name of the x-coordinate column in `polygon` (default `"x"`).

`y`

:   Name of the y-coordinate column in `polygon` (default `"y"`).

`receptor_x`

:   Name of the x-coordinate column in `receptors`, for `object = NULL` (default `"x"`).

`receptor_y`

:   Name of the y-coordinate column in `receptors`, for `object = NULL` (default `"y"`).

`polygon_id`

:   Name of the column in `polygon` that groups the vertices by polygon (default `"polygon"`).

`include_boundary`

:   Count spots on a polygon edge as inside the polygon (default `TRUE`)?

`slot`

:   Name of the receptor slot to use. Leave it `NULL` if the SIRUP object has one receptor table. You must give it if the object has several slots. The function ignores `slot` if you give a plain receptor table.

## Details

If you give `object`, a receptor row gets `NA` if its spatial ID is outside the selected polygons or does not match the spatial object. If `object` is `NULL`, the function finds the polygons from the coordinate columns in `receptors`.

A row in overlapping polygons gets the smallest polygon by area in the primary column. The function also adds one logical column `⁠polygon_<ID>⁠` for each polygon. These columns show membership in every polygon.

If `receptors` is a SIRUP object and you omit `object` and `polygon`, the function uses the stored spatial object and polygons.

## Value

A tibble with the `column` added.

## See Also

[`assignSpatialPolygons()`](assignSpatialPolygons.md), [`withPolygons()`](withPolygons.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Use the spatial object and the polygons stored in the SIRUP object.
receptors <- assignReceptorPolygons(x)
table(receptors$polygon_id, useNA = "ifany")
```
