# polygonAbundance

## Estimate clonal abundance inside polygons

## Description

`polygonAbundance` runs [`alakazam::estimateAbundance()`](https://alakazam.readthedocs.io/en/stable/topics/estimateAbundance/) on the receptors of each polygon.

## Usage

``` R
polygonAbundance(
  receptors,
  polygons = NULL,
  polygon_id = "polygon",
  mode = c("by_polygon", "all_polygons"),
  analysis_name = "polygon_abundance",
  clone = "clone_id",
  copy = NULL,
  group = NULL,
  min_n = 5,
  max_n = NULL,
  uniform = TRUE,
  ci = 0.95,
  nboot = 2000,
  cell_id = NULL,
  progress = FALSE,
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

:   Name for the stored result in a SIRUP object (default `"polygon_abundance"`).

`clone`

:   Name of the column that holds the clone IDs (default `"clone_id"`).

`copy`

:   Name of the column that holds the copy number of each sequence (default `NULL`). If `NULL`, the function measures clone abundance as the number of sequences. Otherwise, it uses the sum of the copy numbers in each clone.

`group`

:   Name of the column that holds the group IDs (default `NULL`). If `NULL`, the function does not group, and the `group` column of the result is `NA`.

`min_n`

:   Smallest number of observations in a group (default 5). The function drops a group with fewer observations.

`max_n`

:   Largest number of observations to sample (default `NULL`). If `NULL`, there is no maximum.

`uniform`

:   Resample each group to the same number of observations (default `TRUE`)? If `FALSE`, the function resamples each group to its original size, or to `max_n` if you give it.

`ci`

:   Confidence interval, between 0 and 1 (default 0.95).

`nboot`

:   Number of bootstrap samples (default 2000).

`cell_id`

:   Name of the column that holds the cell IDs (default `NULL`). If `NULL`, the function assumes bulk data.

`progress`

:   Show a progress bar (default `FALSE`)?

`validate`

:   Validate the updated SIRUP object (default `TRUE`)?

## Details

With `mode = "by_polygon"`, the function runs once for each polygon. With `mode = "all_polygons"`, it runs once on all receptors in any polygon. A polygon with fewer than `min_n` receptor rows returns `NULL`. For a SIRUP object, the function stores the result as an analysis and returns the updated object.

The arguments from `clone` to `progress` go to [`alakazam::estimateAbundance()`](https://alakazam.readthedocs.io/en/stable/topics/estimateAbundance/).

## Value

For a SIRUP object, the updated SIRUP object. Otherwise, the abundance results.

## See Also

[`polygonDiversity()`](polygonDiversity.md), [`withPolygons()`](withPolygons.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Use all polygons together. The example polygons have few receptors.
x <- polygonAbundance(x, mode = "all_polygons", min_n = 2, nboot = 50)
listAnalyses(x)
```
