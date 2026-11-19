# sirup-polygons

## Add or retrieve polygon geometry from a SIRUP object

## Description

`addPolygons` stores a polygon table in a SIRUP object. It replaces any stored polygons. `getPolygons` returns the stored table.

## Usage

``` R
getPolygons(x, sample_id = NULL)

addPolygons(
  x,
  polygons = NULL,
  polygons_path = NULL,
  sample_id = NULL,
  include_boundary = TRUE,
  validate = TRUE
)
```

## Arguments

`x`

:   A SIRUP object.

`sample_id`

:   Sample name or names to add to the polygons from `polygons_path`. Use it only if the paths have no names.

`polygons`

:   Polygon table.

`polygons_path`

:   Path to a flat CSV file of polygon vertices. For a multi-sample object, give a named vector or list of paths, with one path for each sample.

`include_boundary`

:   Count spots on a polygon edge as members of the polygon?

`validate`

:   Validate the updated object?

## Details

The input to `addPolygons` can be a nested tibble with `polygon_id` and `vertices` columns. It can also be a flat vertex table, such as the CSV file from [`selectSpatialPolygons()`](selectSpatialPolygons.md).

The function adds a `polygon_id` column to both the Seurat metadata and the receptor table. If a barcode lies in overlapping or nested polygons, it gets the ID of the smallest polygon by area. The function also adds one logical `⁠polygon_<ID>⁠` column for each polygon. These columns show membership in every polygon as booleans.

## Value

`getPolygons` returns a polygon tibble. `addPolygons` returns the updated SIRUP object.

## See Also

[`selectSpatialPolygons()`](selectSpatialPolygons.md), [`createSIRUP()`](createSIRUP.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
getPolygons(x)

# Replace the stored polygons with one polygon in the middle of the data.
center <- tibble::tibble(
  polygon_id = "center",
  vertices = list(data.frame(x = c(2.5, 3.5, 3.5, 2.5), y = c(2.5, 2.5, 3.5, 3.5)))
)
x2 <- addPolygons(x, center)
getPolygons(x2)$polygon_id
```
