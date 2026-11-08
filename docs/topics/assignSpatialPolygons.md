# assignSpatialPolygons

## Add polygon IDs to the spatial metadata

## Description

`assignSpatialPolygons` adds polygon columns to the Seurat metadata. It keeps all spots.

## Usage

``` R
assignSpatialPolygons(
  object,
  polygon = NULL,
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  polygon_val = NULL,
  x = "x",
  y = "y",
  polygon_id = "polygon",
  metadata = "polygon_id",
  include_boundary = TRUE
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`polygon`

:   Data frame of polygon vertices in order. It needs the columns set by `x` and `y`. If the column set by `polygon_id` exists, the function treats each polygon on its own. If `object` is a SIRUP object and `polygon` is `NULL`, the function uses the stored polygons.

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

`polygon_id`

:   Name of the column in `polygon` that groups the vertices by polygon (default `"polygon"`).

`metadata`

:   Name of the metadata column that holds the primary polygon of each spot (default `"polygon_id"`).

`include_boundary`

:   Count spots on a polygon edge as inside the polygon (default `TRUE`)?

## Details

A spot outside the selected polygons gets `NA` in the primary column. A spot in overlapping polygons gets the smallest polygon by area. The function also adds one logical column `⁠polygon_<ID>⁠` for each polygon. These columns show membership in every polygon.

## Value

A Seurat object with the `metadata` column added to `⁠@meta.data⁠`.

## See Also

[`seuratSubsetPolygon()`](seuratSubsetPolygon.md), [`assignReceptorPolygons()`](assignReceptorPolygons.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
flat <- data.frame(
  polygon = "block",
  x = c(3.5, 5.5, 5.5, 3.5),
  y = c(0.5, 0.5, 2.5, 2.5)
)
spatial <- assignSpatialPolygons(
  getSpatial(x),
  polygon = flat,
  image_scale = "hires",
  swap_coords = FALSE
)
table(spatial$polygon_id, useNA = "ifany")
```
