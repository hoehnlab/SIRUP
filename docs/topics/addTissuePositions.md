# addTissuePositions

## Add tissue coordinates to an AIRR table

## Description

`addTissuePositions` joins the tissue positions of the spots to an AIRR table. It adds the columns `x` and `y`.

## Usage

``` R
addTissuePositions(
  airr_df,
  positions,
  spatial_id = "spatial_id_8um",
  position_spatial_id = "barcode",
  x = "pxl_col_in_fullres",
  y = "pxl_row_in_fullres",
  keep_position_cols = TRUE
)
```

## Arguments

`airr_df`

:   AIRR-style data frame with a spatial ID column.

`positions`

:   Data frame of tissue positions, or the path of a CSV or parquet file.

`spatial_id`

:   Name of the column in `airr_df` to match (default `"spatial_id_8um"`).

`position_spatial_id`

:   Name of the spatial ID column in `positions` (default `"barcode"`).

`x`

:   Name of the column in `positions` to copy to `x` (default `"pxl_col_in_fullres"`).

`y`

:   Name of the column in `positions` to copy to `y` (default `"pxl_row_in_fullres"`).

`keep_position_cols`

:   Keep all columns of `positions` (default `TRUE`)? If `FALSE`, the result has the columns of `airr_df`, `x`, and `y`.

## Value

A data frame with the columns `x` and `y` added.

## See Also

[`readSpatialIdMapping()`](readSpatialIdMapping.md), [`extractSpatialIds()`](extractSpatialIds.md)

## Examples

``` R
reads <- data.frame(sequence_id = c("a", "b"), spatial_id_8um = c("s1", "s2"))
positions <- data.frame(
  barcode = c("s1", "s2"),
  pxl_col_in_fullres = c(100, 200),
  pxl_row_in_fullres = c(50, 60)
)
addTissuePositions(reads, positions)
```
