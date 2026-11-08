# polygonDiversity

## Estimate clonal diversity inside polygons

## Description

`polygonDiversity` runs [`alakazam::alphaDiversity()`](https://alakazam.readthedocs.io/en/stable/topics/alphaDiversity/) on the receptors of each polygon.

## Usage

``` R
polygonDiversity(
  receptors,
  polygons = NULL,
  polygon_id = "polygon",
  mode = c("by_polygon", "all_polygons"),
  analysis_name = "polygon_diversity",
  clone = "clone_id",
  copy = NULL,
  group = NULL,
  min_n = 30,
  min_q = 0,
  max_q = 4,
  step_q = 0.1,
  ci = 0.95,
  ...,
  validate = TRUE
)
```

## Arguments

`receptors`

:   SIRUP object, or an AIRR-style receptor table.

`polygons`

:   Polygon table that sets the order of the output polygons. If `NULL`, the function takes the polygon IDs from the SIRUP object, or from the `polygon_id` values in `receptors`.

`polygon_id`

:   Name of the column that holds the primary polygon of each receptor (default `"polygon"`).

`mode`

:   `"by_polygon"` runs the calculation once for each polygon. `"all_polygons"` runs it once on all receptors in any polygon.

`analysis_name`

:   Name for the stored result in a SIRUP object (default `"polygon_diversity"`).

`clone`

:   Name of the column that holds the clone IDs (default `"clone_id"`).

`copy`

:   Name of the column that holds the copy number of each sequence (default `NULL`). If `NULL`, the function measures clone abundance as the number of sequences. Otherwise, it uses the sum of the copy numbers in each clone.

`group`

:   Name of the column that holds the group IDs (default `NULL`). If `NULL`, the function does not group, and the `group` column of the result is `NA`.

`min_n`

:   Smallest number of observations in a group (default 30). The function drops a group with fewer observations.

`min_q`

:   Smallest diversity order q (default 0).

`max_q`

:   Largest diversity order q (default 4).

`step_q`

:   Step between the values of q (default 0.1).

`ci`

:   Confidence interval, between 0 and 1 (default 0.95).

`...`

:   Other arguments for [`alakazam::alphaDiversity()`](https://alakazam.readthedocs.io/en/stable/topics/alphaDiversity/).

`validate`

:   Validate the updated SIRUP object (default `TRUE`)?

## Details

With `mode = "by_polygon"`, the function runs once for each polygon. With `mode = "all_polygons"`, it runs once on all receptors in any polygon. A polygon with fewer than `min_n` receptor rows returns `NULL`. For a SIRUP object, the function stores the result as an analysis and returns the updated object.

The arguments from `clone` to `...` go to [`alakazam::alphaDiversity()`](https://alakazam.readthedocs.io/en/stable/topics/alphaDiversity/).

## Value

For a SIRUP object, the updated SIRUP object. Otherwise, the diversity results.

## See Also

[`summarizePolygonDiversity()`](summarizePolygonDiversity.md), [`polygonAbundance()`](polygonAbundance.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- polygonDiversity(x, min_n = 2, max_q = 2, step_q = 1)
listAnalyses(x)
```
