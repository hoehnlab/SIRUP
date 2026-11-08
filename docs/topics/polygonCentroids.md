# polygonCentroids

## Calculate polygon centroids

## Description

`polygonCentroids` returns the centroid of each polygon.

## Usage

``` R
polygonCentroids(
  object = NULL,
  polygons = NULL,
  polygon_id = "polygon",
  x = "x",
  y = "y"
)
```

## Arguments

`object`

:   SIRUP object with stored polygons, or a table of polygon vertices.

`polygons`

:   Table of polygon vertices. If you give it, you can omit `object`, or use `object` only for context.

`polygon_id`

:   Name of the polygon ID column, for flat polygon input (default `"polygon"`).

`x`

:   Name of the x-coordinate column, for flat polygon input (default `"x"`).

`y`

:   Name of the y-coordinate column, for flat polygon input (default `"y"`).

## Value

A tibble with one row for each polygon. The columns are `sample`, `polygon_id`, `polygon_label`, `x`, `y`, `.sirup_y`, and `n_vertices`.

## See Also

[`polygonPairwiseDistance()`](polygonPairwiseDistance.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
polygonCentroids(x)
```
