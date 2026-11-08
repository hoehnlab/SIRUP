# expandPolygons

## Expand polygons outward from their centroids

## Description

`expandPolygons` makes each polygon larger by moving its vertices away from the polygon centroid.

## Usage

``` R
expandPolygons(
  polygons,
  delta = NULL,
  x = "x",
  y = "y",
  polygon_id = "polygon"
)
```

## Arguments

`polygons`

:   Polygons to expand. Give a flat data frame of vertices (with `x`, `y`, and, if you want, a polygon ID column). You can also pass a SIRUP object. The function then expands the stored polygons.

`delta`

:   Distance that each vertex moves outward, in coordinate units. Required. A negative value moves the vertices inward, so the polygon becomes smaller.

`x`

:   Name of the x-coordinate column, if `polygons` is a data frame (default `"x"`).

`y`

:   Name of the y-coordinate column, if `polygons` is a data frame (default `"y"`).

`polygon_id`

:   Name of the column that separates the polygons, if `polygons` is a data frame (default `"polygon"`). If the column does not exist, the function uses a `polygon` or `polygon_id` column. If neither exists, all rows form one polygon.

## Details

Each vertex moves along the line from the centroid through the vertex. It moves by the distance `delta`. The polygon keeps its shape and its number of vertices. The polygons must be convex.

## Value

The expanded polygons in the form of the input: a data frame for a data frame, or the updated SIRUP object for a SIRUP object.

## See Also

[`addPolygons()`](sirup-polygons.md), [`plotSpatialPolygons()`](plotSpatialPolygons.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Make the stored polygons larger by 0.25 coordinate units.
x2 <- expandPolygons(x, delta = 0.25)

# Use a data frame of vertices.
flat <- data.frame(
  polygon = "block",
  x = c(3.5, 5.5, 5.5, 3.5),
  y = c(0.5, 0.5, 2.5, 2.5)
)
expandPolygons(flat, delta = 0.25)
```
