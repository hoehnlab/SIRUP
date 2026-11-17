# polygonAbundance

## Estimate clonal abundance inside polygons

## Description

`polygonAbundance` runs [`alakazam::estimateAbundance()`](https://alakazam.readthedocs.io/en/stable/topics/estimateAbundance/) once on the receptors of all polygons, with the polygon as the group. alakazam then resamples every polygon to the same number of receptors, so the curves of different polygons are comparable. With `common_size = FALSE`, the function analyses each polygon alone instead.

## Usage

``` R
polygonAbundance(
  receptors,
  polygons = NULL,
  polygon_id = "polygon",
  analysis_name = NULL,
  clone = "clone_id",
  copy = NULL,
  common_size = TRUE,
  slot = NULL,
  min_n = 10,
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

:   Polygon table that sets the polygons to use. If `NULL`, the function takes the polygon IDs from the SIRUP object, or from the `polygon_id` values in `receptors`.

`polygon_id`

:   Name of the column that holds the primary polygon of each receptor (default `"polygon"`). For a SIRUP object, the default is `"polygon_id"`.

`analysis_name`

:   Name for the stored result in a SIRUP object (default `"polygon_abundance"`). If the object has several receptor slots, the default also ends with `⁠_⁠` and the slot name.

`clone`

:   Name of the column that holds the clone IDs (default `"clone_id"`).

`copy`

:   Name of the column that holds the copy number of each sequence (default `NULL`). If `NULL`, the function measures clone abundance as the number of sequences. Otherwise, it uses the sum of the copy numbers in each clone.

`common_size`

:   Resample all polygons to the same size (default `TRUE`)? If `FALSE`, the function analyses each polygon alone, as in earlier versions of SIRUP, and returns a list with one result for each polygon.

`slot`

:   Name of the receptor slot to use. Leave it `NULL` if the SIRUP object has one receptor table. You must give it if the object has several slots. The function ignores `slot` if you give a plain receptor table.

`min_n`

:   Smallest number of receptors that a polygon needs (default 10). The function drops a polygon with fewer receptors. With `common_size = TRUE`, the smallest polygon that stays sets the size of the resampling, so a larger `min_n` gives less noisy curves and keeps fewer polygons.

`max_n`

:   Largest number of receptors to sample from each polygon (default `NULL`). If `NULL`, the size is the size of the smallest polygon that stays.

`uniform`

:   Resample each polygon to the same number of receptors (default `TRUE`)? Passed to [`alakazam::estimateAbundance()`](https://alakazam.readthedocs.io/en/stable/topics/estimateAbundance/).

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

The function puts the receptors of all polygons in one table and passes the polygon column to [`alakazam::estimateAbundance()`](https://alakazam.readthedocs.io/en/stable/topics/estimateAbundance/) as `group`. With `uniform = TRUE` (the default), alakazam resamples each polygon to the size of the smallest polygon that has at least `min_n` receptors. Use `max_n` to set a smaller size. The function leaves out polygons with fewer than `min_n` receptors. If no polygon has enough receptors, the function gives an error.

The result is one `AbundanceCurve` object. Its `abundance` table has one column of polygon IDs. [`alakazam::plotAbundanceCurve()`](https://alakazam.readthedocs.io/en/stable/topics/plotAbundanceCurve/) plots all polygons together. For a SIRUP object, the function stores the result as an analysis and returns the updated object.

With `common_size = FALSE`, the function runs [`alakazam::estimateAbundance()`](https://alakazam.readthedocs.io/en/stable/topics/estimateAbundance/) once for each polygon. Each polygon is then resampled to its own size, so polygons with more receptors can look more diverse. A polygon with fewer than `min_n` receptors gives `NULL`. The result is a named list with one `AbundanceCurve` for each polygon.

The arguments from `clone` to `progress` go to [`alakazam::estimateAbundance()`](https://alakazam.readthedocs.io/en/stable/topics/estimateAbundance/).

## Value

For a SIRUP object, the updated SIRUP object. Otherwise, the abundance results: one `AbundanceCurve` object, or a list with one object for each polygon if `common_size = FALSE`.

## See Also

[`polygonDiversity()`](polygonDiversity.md), [`withPolygons()`](withPolygons.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# The example polygons have few receptors, so lower min_n.
x <- polygonAbundance(x, min_n = 2, nboot = 50)
listAnalyses(x)
```
