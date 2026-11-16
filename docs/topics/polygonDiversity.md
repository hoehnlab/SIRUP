# polygonDiversity

## Estimate clonal diversity inside polygons

## Description

`polygonDiversity` runs [`alakazam::alphaDiversity()`](https://alakazam.readthedocs.io/en/stable/topics/alphaDiversity/) once on the receptors of all polygons, with the polygon as the group. alakazam then resamples every polygon to the same number of receptors, so the diversity of different polygons is comparable. With `common_size = FALSE`, the function analyses each polygon alone instead.

## Usage

``` R
polygonDiversity(
  receptors,
  polygons = NULL,
  polygon_id = "polygon",
  analysis_name = "polygon_diversity",
  clone = "clone_id",
  copy = NULL,
  common_size = TRUE,
  min_n = 10,
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

:   Polygon table that sets the polygons to use. If `NULL`, the function takes the polygon IDs from the SIRUP object, or from the `polygon_id` values in `receptors`.

`polygon_id`

:   Name of the column that holds the primary polygon of each receptor (default `"polygon"`). For a SIRUP object, the default is `"polygon_id"`.

`analysis_name`

:   Name for the stored result in a SIRUP object (default `"polygon_diversity"`).

`clone`

:   Name of the column that holds the clone IDs (default `"clone_id"`).

`copy`

:   Name of the column that holds the copy number of each sequence (default `NULL`). If `NULL`, the function measures clone abundance as the number of sequences. Otherwise, it uses the sum of the copy numbers in each clone.

`common_size`

:   Resample all polygons to the same size (default `TRUE`)? If `FALSE`, the function analyses each polygon alone, as in earlier versions of SIRUP, and returns a list with one result for each polygon.

`min_n`

:   Smallest number of receptors that a polygon needs (default 10). The function drops a polygon with fewer receptors. With `common_size = TRUE`, the smallest polygon that stays sets the size of the resampling, so a larger `min_n` gives less noisy estimates and keeps fewer polygons.

`min_q`

:   Smallest diversity order q (default 0).

`max_q`

:   Largest diversity order q (default 4).

`step_q`

:   Step between the values of q (default 0.1).

`ci`

:   Confidence interval, between 0 and 1 (default 0.95).

`...`

:   Other arguments for [`alakazam::alphaDiversity()`](https://alakazam.readthedocs.io/en/stable/topics/alphaDiversity/), such as `nboot`, `max_n` and `uniform`. You cannot give `group`, because the function uses the polygon as the group.

`validate`

:   Validate the updated SIRUP object (default `TRUE`)?

## Details

The function puts the receptors of all polygons in one table and passes the polygon column to [`alakazam::alphaDiversity()`](https://alakazam.readthedocs.io/en/stable/topics/alphaDiversity/) as `group`. alakazam resamples each polygon to the size of the smallest polygon that has at least `min_n` receptors. Use `max_n` (through `...`) to set a smaller size. The function leaves out polygons with fewer than `min_n` receptors. If no polygon has enough receptors, the function gives an error.

The result is one `DiversityCurve` object. Its `diversity` table has one column of polygon IDs. [`alakazam::plotDiversityCurve()`](https://alakazam.readthedocs.io/en/stable/topics/plotDiversityCurve/) plots all polygons together. For a SIRUP object, the function stores the result as an analysis and returns the updated object.

With `common_size = FALSE`, the function runs [`alakazam::alphaDiversity()`](https://alakazam.readthedocs.io/en/stable/topics/alphaDiversity/) once for each polygon. Each polygon is then resampled to its own size, so polygons with more receptors can look more diverse. A polygon with fewer than `min_n` receptors gives `NULL`. The result is a named list with one `DiversityCurve` for each polygon.

The arguments from `clone` to `...` go to [`alakazam::alphaDiversity()`](https://alakazam.readthedocs.io/en/stable/topics/alphaDiversity/).

## Value

For a SIRUP object, the updated SIRUP object. Otherwise, the diversity results: one `DiversityCurve` object, or a list with one object for each polygon if `common_size = FALSE`.

## See Also

[`summarizePolygonDiversity()`](summarizePolygonDiversity.md), [`polygonAbundance()`](polygonAbundance.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- polygonDiversity(x, min_n = 2, max_q = 2, step_q = 1)
listAnalyses(x)
```
