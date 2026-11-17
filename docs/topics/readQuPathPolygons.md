# readQuPathPolygons

## Read QuPath annotations as SIRUP polygons

## Description

`readQuPathPolygons` reads a GeoJSON file exported from QuPath. It returns a polygon table. Pass the table to [`addPolygons()`](sirup-polygons.md).

## Usage

``` R
readQuPathPolygons(
  path,
  sample_id = NULL,
  class = NULL,
  keep_unclassified = TRUE,
  scale = 1
)
```

## Arguments

`path`

:   Path to the GeoJSON file.

`sample_id`

:   Sample name to add to all polygons. A multi-sample SIRUP object needs it. Use `NULL` for a single-sample object.

`class`

:   Class names to keep, for example `"GC"`. If `NULL`, the function keeps all classes.

`keep_unclassified`

:   Keep annotations that have no class (default `TRUE`)? The function ignores this if you set `class`.

`scale`

:   Number that multiplies all x and y values. See Details.

## Details

In QuPath, use `⁠File > Object data... > Export as GeoJSON⁠`. Select the annotations first, or export all objects. A script can also call `exportObjectsToGeoJson(getAnnotationObjects(), path, "FEATURE_COLLECTION")`.

QuPath writes coordinates in pixels of the full-resolution image. The origin is the top-left corner, and y increases downward. SIRUP polygons must use the same coordinates as `SeuratObject::GetTissueCoordinates()` at the `image_scale` of the object. If the object uses the `hires` or `lowres` image, set `scale` to the matching scale factor from `scalefactors_json.json`. Use `scale = 1` for full-resolution coordinates.

The function accepts a FeatureCollection, a single Feature, and the bare list of Features that older QuPath versions write. It keeps features with `objectType` equal to `"annotation"`, or with no `objectType`. It reads only `Polygon` and `MultiPolygon` shapes. A `MultiPolygon` becomes one SIRUP polygon for each part. The IDs of the parts end in `⁠_p1⁠`, `⁠_p2⁠`, and so on.

A SIRUP polygon has one ring and cannot have holes. If an annotation has holes, the function drops them and gives a warning. Spots in a hole then count as inside the polygon.

Polygon IDs come from the annotation name. If an annotation has no name, the ID comes from the class name and a number, for example `ENDO_01`. An annotation with neither gets the label `annotation`.

## Value

A tibble with the columns `polygon_id`, `class`, `vertices`, and `sample` (only if you set `sample_id`). `vertices` is a list-column of tables with the columns `x` and `y`.

## See Also

[`addPolygons()`](sirup-polygons.md), [`plotSpatialPolygons()`](plotSpatialPolygons.md)

## Examples

``` R
geojson <- tempfile(fileext = ".geojson")
writeLines(paste0(
  '{"type":"FeatureCollection","features":[{"type":"Feature",',
  '"geometry":{"type":"Polygon","coordinates":',
  '[[[2.5,2.5],[3.5,2.5],[3.5,3.5],[2.5,3.5],[2.5,2.5]]]},',
  '"properties":{"objectType":"annotation",',
  '"classification":{"name":"GC","color":[255,0,0]}}}]}'
), geojson)
polygons <- readQuPathPolygons(geojson)
polygons

x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- addPolygons(x, polygons)
getPolygons(x)$polygon_id
```
