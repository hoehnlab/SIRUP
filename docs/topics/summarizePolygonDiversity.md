# summarizePolygonDiversity

## Summarize polygon diversity at one value of q

## Description

`summarizePolygonDiversity` takes the diversity estimate for one diversity order q from a stored analysis. It makes one table row for each polygon and stores the table as a new analysis.

## Usage

``` R
summarizePolygonDiversity(
  object,
  q,
  diversity_analysis = "polygon_diversity",
  analysis_name = NULL,
  tolerance = sqrt(.Machine$double.eps),
  validate = TRUE
)
```

## Arguments

`object`

:   SIRUP object.

`q`

:   Diversity order to take, for example `2` for Simpson diversity.

`diversity_analysis`

:   Name of the stored analysis from [`polygonDiversity()`](polygonDiversity.md) (default `"polygon_diversity"`).

`analysis_name`

:   Name for the stored table. If `NULL`, the function uses `polygon_diversity_q` followed by the value of `q`.

`tolerance`

:   Tolerance for the match of `q` (default `sqrt(.Machine$double.eps)`).

`validate`

:   Validate the updated SIRUP object (default `TRUE`)?

## Details

Polygons that [`polygonDiversity()`](polygonDiversity.md) left out, because they had fewer than `min_n` receptors, get a row with `NA` values. [`plotSpatialPolygons()`](plotSpatialPolygons.md) draws these polygons in grey. The function reads both kinds of result: one `DiversityCurve` for all polygons, and the list that `polygonDiversity(common_size = FALSE)` makes.

## Value

The updated SIRUP object, with the diversity summary as an analysis.

## See Also

[`polygonDiversity()`](polygonDiversity.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- polygonDiversity(x, min_n = 2, max_q = 2, step_q = 1)
x <- summarizePolygonDiversity(x, q = 1)
listAnalyses(x)
```
