# withPolygons

## Apply a function to the receptors of each polygon

## Description

`withPolygons` splits the receptor table by polygon and applies a function to the rows of each polygon. Or it applies the function once to all receptors in any polygon.

## Usage

``` R
withPolygons(
  receptors,
  .f,
  ...,
  polygons = NULL,
  polygon_id = "polygon",
  mode = c("by_polygon", "all_polygons"),
  pass_polygon_id = FALSE
)
```

## Arguments

`receptors`

:   SIRUP object, or an AIRR-style receptor table.

`.f`

:   Function to apply.

`...`

:   Other arguments for `.f`.

`polygons`

:   Polygon table that sets the order of the output polygons. If `NULL`, the function takes the polygon IDs from the SIRUP object, or from the `polygon_id` values in `receptors`.

`polygon_id`

:   Name of the column that holds the primary polygon of each receptor (default `"polygon"`).

`mode`

:   `"by_polygon"` runs `.f` once for each polygon. `"all_polygons"` runs `.f` once on all receptors in any polygon.

`pass_polygon_id`

:   Pass the current polygon ID to `.f` as the argument `polygon_id`, if `.f` accepts it (default `FALSE`)?

## Value

For `mode = "by_polygon"`, a named list with one result for each polygon. Otherwise, the result of `.f`.

## See Also

[`assignReceptorPolygons()`](assignReceptorPolygons.md), [`polygonAbundance()`](polygonAbundance.md), [`polygonDiversity()`](polygonDiversity.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Count the receptors in each polygon.
withPolygons(x, nrow)

# Count the receptors in all polygons together.
withPolygons(x, nrow, mode = "all_polygons")
```
