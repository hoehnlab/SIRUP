# addQuPathPolygons

## Add QuPath annotations to a SIRUP object as polygons

## Description

`addQuPathPolygons` reads a GeoJSON file that QuPath exported. It converts the annotations to the coordinates of a spatial sample in a SIRUP object. It adds a buffer to each polygon, so that every spot whose circle touches the polygon becomes a member of it. Then it stores the polygons in the object.

## Usage

``` R
addQuPathPolygons(
  x,
  path,
  sample_id = NULL,
  class = NULL,
  keep_unclassified = TRUE,
  image_scale = NULL,
  swap_coords = NULL,
  image = NULL,
  buffer = NULL,
  replace = FALSE,
  include_boundary = TRUE,
  validate = TRUE
)
```

## Arguments

`x`

:   A SIRUP object.

`path`

:   Path to the QuPath GeoJSON file.

`sample_id`

:   Name of the spatial sample that the annotations belong to. A SIRUP object with one spatial sample does not need it. A multi-sample object does.

`class`

:   Class names to keep. If `NULL`, the function keeps all classes.

`keep_unclassified`

:   Keep annotations that have no class (default `TRUE`)? The function ignores this if you set `class`.

`image_scale`

:   Image scale of the polygon coordinates, for example `"hires"`. If `NULL`, the function uses the value stored in `x`. If `x` has no value, you must give one, and the function stores it in `x`. A value that differs from the stored one gives an error.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in `x`. The rules are the same as for `image_scale`.

`image`

:   Name of the Seurat image to use. If `NULL`, the function uses the only image of the sample.

`buffer`

:   Distance to grow each polygon, in the coordinates of the object. If `NULL` (default), the function uses the spot radius. Use `0` for no buffer.

`replace`

:   Replace the stored polygons of this sample (default `FALSE`)? If `FALSE`, the function adds the new polygons to the stored ones, and the polygon IDs must not repeat.

`include_boundary`

:   Count spots on a polygon edge as members of the polygon?

`validate`

:   Validate the updated object?

## Details

QuPath writes coordinates in pixels of the full-resolution image. The function finds the factor that converts them to the coordinates of the object. It compares `SeuratObject::GetTissueCoordinates()` at the `image_scale` of the object with the unscaled coordinates. This works for `"hires"` and `"lowres"`. The image that you annotated in QuPath must be the full-resolution image that the spots were registered to.

`image_scale` and `swap_coords` come from the SIRUP object. If the object has no `image_scale`, you must give it. The function swaps the polygon axes in the same way as the spots when `swap_coords` is `TRUE`.

A spot is a circle. By default, the function grows each polygon by the spot radius, which it reads from the scale factors of the Seurat image. If the sample has no Seurat image, the function uses the QuPath coordinates as they are and sets the default buffer to 0. Then [`assignSpatialPolygons()`](assignSpatialPolygons.md) counts a spot as inside when its centre is inside the grown polygon. That is the same as counting every spot whose circle touches the original polygon. Set `buffer = 0` to count only spots with the centre inside the QuPath polygon. For Visium HD, check that the stored spot size suits your bins, or give `buffer`.

A SIRUP polygon cannot have holes. The function drops holes and gives a warning. See [`readQuPathPolygons()`](readQuPathPolygons.md) for the rules about shapes, classes, and IDs.

## Value

The updated SIRUP object. Spatial metadata and receptor records get the new `polygon_id` values.

## See Also

[`readQuPathPolygons()`](readQuPathPolygons.md), [`addPolygons()`](sirup-polygons.md), [`plotSpatialPolygons()`](plotSpatialPolygons.md)

## Examples

``` R
# Write a small GeoJSON file with one annotation. QuPath makes files like it.
geojson <- tempfile(fileext = ".geojson")
writeLines(paste0(
  '{"type":"FeatureCollection","features":[{"type":"Feature",',
  '"geometry":{"type":"Polygon","coordinates":',
  '[[[3.5,0.5],[5.5,0.5],[5.5,2.5],[3.5,2.5],[3.5,0.5]]]},',
  '"properties":{"objectType":"annotation","name":"GC1"}}]}'
), geojson)

x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- addQuPathPolygons(x, geojson, replace = TRUE)
getPolygons(x)
table(getSpatial(x)$polygon_id, useNA = "ifany")
```
