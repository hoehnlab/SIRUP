# polygonReceptorSummary

## Summarize columns in receptors data inside polygons

## Description

`polygonReceptorSummary` makes one summary row for each polygon.

## Usage

``` R
polygonReceptorSummary(
  object,
  polygons = NULL,
  spatial_id = NULL,
  polygon_id = "polygon",
  clone = "clone_id",
  weight = NULL,
  summary_columns = NULL,
  summary_operations = NULL,
  include_total_spots = FALSE,
  analysis_name = NULL,
  slot = NULL,
  validate = TRUE
)
```

## Arguments

`object`

:   SIRUP object, or an AIRR-style receptor table.

`polygons`

:   Polygon table that sets the order of the output polygons. If `NULL`, the function takes the polygon IDs from the SIRUP object, or from the `polygon_id` values in `object`.

`spatial_id`

:   Name of the column that holds the spatial IDs. If `NULL`, the function uses the stored field of a SIRUP object, or looks for a standard spatial ID column in the receptor table.

`polygon_id`

:   Name of the column that holds the primary polygon of each receptor (default `"polygon"`).

`clone`

:   Name of the column that holds the clone IDs (default `"clone_id"`). The function uses it for `n_clones` and `top_clone_fraction`. If the column does not exist, the function leaves out these two summaries.

`weight`

:   Name of a numeric column that holds the receptor weights. If you give it, `n_receptors`, `top_clone_fraction`, and the numeric custom summaries use the weights.

`summary_columns`

:   Names of the receptor columns to summarize.

`summary_operations`

:   Operations for `summary_columns`. Give one operation, or one for each column. The operations are `count`, `sum`, `mean`, `median`, `variance`, `var`, `sd`, `min`, `max`, `any`, `all`, and `cat`.

`include_total_spots`

:   Add `n_total_spatial_ids` (default `FALSE`)? This column holds the number of all spots in each polygon, including spots with no receptor rows. It works only for SIRUP objects.

`analysis_name`

:   Name for the stored summary in a SIRUP object (default `"polygon_receptor_summary"`). If the object has several receptor slots, the default also ends with `⁠_⁠` and the slot name.

`slot`

:   Name of the receptor slot to use. Leave it `NULL` if the SIRUP object has one receptor table. You must give it if the object has several slots. The function ignores `slot` if you give a plain receptor table.

`validate`

:   Validate the updated SIRUP object (default `TRUE`)?

## Details

The input can be a SIRUP object or a receptor table that already has polygon columns. For a SIRUP object, the function stores the summary as an analysis and returns the updated object. For a receptor table, it returns the summary tibble.

## Value

For a SIRUP object, the updated SIRUP object. Otherwise, a summary tibble.

## See Also

[`assignReceptorPolygons()`](assignReceptorPolygons.md), [`polygonDiversity()`](polygonDiversity.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- polygonReceptorSummary(x)
getAnalysis(x, "polygon_receptor_summary")
```
