# assignTipColors

## Assign colors to the tips of clones

## Description

`assignTipColors` gives each tip of each clone a color ID and a color. Tree plots and the matching spatial plots use the same colors.

## Usage

``` R
assignTipColors(
  clones,
  clone_val = NULL,
  clone = "clone_id",
  id = "sequence_id",
  shared_color = NULL,
  palette = NULL
)
```

## Arguments

`clones`

:   Data frame of AIRR data with clone IDs.

`clone_val`

:   Clone IDs to keep before the function assigns the colors. If `NULL`, the function keeps all clones.

`clone`

:   Name of the column that holds the clone IDs (default `"clone_id"`).

`id`

:   Name of the column that holds the sequence or tip IDs (default `"sequence_id"`).

`shared_color`

:   Name of a column. Tips that share a value in this column get the same color within a clone, for example a spatial ID column.

`palette`

:   Vector of colors. If `NULL`, the function uses `pals::glasbey()` if the package is installed. Otherwise, it uses the base R colors from `grDevices::hcl.colors()`.

## Details

The function numbers the tips of each clone in the order of `id`. If you give `shared_color`, tips with the same value of that column get the same number within a clone. For example, tips in the same spot get the same color.

## Value

A data frame with the columns `tip_color_id` and `tip_color` added.

## See Also

[`plotSpatialCloneTreePairs()`](plotSpatialCloneTreePairs.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
colors <- assignTipColors(getReceptors(x), shared_color = "spatial_id")
head(colors[, c("clone_id", "spatial_id", "tip_color_id", "tip_color")])
```
