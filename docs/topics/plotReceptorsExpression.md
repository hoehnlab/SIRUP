# plotReceptorsExpression

## Plot receptor groups over gene expression

## Description

`plotReceptorsExpression` plots one feature together with one or more groups of receptors on a tissue image. For example, it shows the BCR-positive and TCR-positive spots on one expression map.

## Usage

``` R
plotReceptorsExpression(
  object,
  receptors = NULL,
  group = NULL,
  group_val = NULL,
  spatial_id = NULL,
  feature = NULL,
  gene = NULL,
  colors = NULL,
  overlap = "average",
  value = NULL,
  value_colors = c("#FFFFFF", "#000080"),
  value_limits = c(0, NA),
  value_combine = c("mean", "median", "mode", "max", "min", "sum"),
  receptor_size = 1,
  receptor_stroke = 0.7,
  receptor_alpha = 1,
  feature_colors = c("#939393", "#FFFF00"),
  feature_size = 1,
  feature_alpha = 1,
  zero_tol = 0,
  title = NULL,
  order_points = TRUE,
  assay = NULL,
  layer = "counts",
  image_scale = NULL,
  swap_coords = NULL,
  images = NULL,
  sample_id = NULL,
  image_alpha = 1,
  zoom = FALSE,
  receptor_legend_title = NULL,
  show_receptor_legend = TRUE,
  combine = TRUE,
  combine_sample_plots = FALSE
)
```

## Arguments

`object`

:   Seurat object, or a SIRUP object with spatial data and receptors.

`receptors`

:   AIRR-style receptor table with a spatial ID column and a group column. Omit it if `object` is a SIRUP object.

`group`

:   Name of the `receptors` column that defines the receptor groups, for example `"receptor_type"` with the values `"BCR"` and `"TCR"`. You can omit it if you give `value`. If you give both, `group` and `group_val` only select the receptors to draw.

`group_val`

:   Group values to plot, in the order of the plot. If `NULL`, the function uses all non-missing values of `group`.

`spatial_id`

:   Name of the `receptors` column that holds the Seurat spot names. If `NULL`, the function uses the stored field of a SIRUP object, or looks for a standard spatial ID column.

`feature`

:   Feature to plot under the spots with receptors.

`gene`

:   Alias for `feature`.

`colors`

:   Colors of the receptor groups. Named colors match the group values. Unnamed colors follow the order of the values.

`overlap`

:   Color of a spot that holds more than one receptor group. The function draws such a spot once, as a combined group with the member names joined by `" + "`. With `"average"` (default), the spot gets the mean color of its members. With one color, all combined spots get that color.

`value`

:   Color the receptor spots on a continuous gradient and not by group. Give the name of a numeric `receptors` column, or `"count"` for the number of receptors in the spot.

`value_colors`

:   Colors of the `value` gradient (default white to navy). Two colors give a low-to-high gradient. More colors give a gradient with several stops.

`value_limits`

:   Limits of the `value` scale (default `c(0, NA)`).

`value_combine`

:   How to combine `value` if a spot holds more than one receptor: `"mean"` (default), `"median"`, `"mode"`, `"max"`, `"min"`, or `"sum"`. The function ignores it if `value = "count"`.

`receptor_size`

:   Size of the receptor circles (default 1).

`receptor_stroke`

:   Line width of the receptor circles (default 0.7).

`receptor_alpha`

:   Opacity of the receptor circles (default 1).

`feature_colors`

:   Colors that make the expression scale (default `c("#939393", "#FFFF00")`).

`feature_size`

:   Size of the expression points (default 1).

`feature_alpha`

:   Opacity of the expression points (default 1).

`zero_tol`

:   Values at or below this threshold count as zero (default 0).

`title`

:   Plot title.

`order_points`

:   Draw the low expression values first, so that the high values show on top (default `TRUE`)? The function draws the receptor groups in the order of `group_val`.

`assay`

:   Assay to use.

`layer`

:   Seurat expression layer to use (default `"counts"`).

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`.

`swap_coords`

:   Swap the two Seurat coordinate axes?

`images`

:   Names of the Seurat images to use.

`sample_id`

:   Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`image_alpha`

:   Opacity of the tissue image.

`zoom`

:   Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot shows the full tissue area that the spatial coordinates cover.

`receptor_legend_title`

:   Title of the receptor legend. If `NULL`, the function uses `group`.

`show_receptor_legend`

:   Show the receptor legend (default `TRUE`)?

`combine`

:   Combine the image panels with patchwork (default `TRUE`)?

`combine_sample_plots`

:   Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

## Details

The function draws the spots with expression above zero with a continuous fill scale. It draws the spots with receptors on top as hollow circles.

If the object has multiple images or samples, the function draws one plot for each image or sample.

## Value

A ggplot object or a patchwork object. If `combine = FALSE`, a list.

## See Also

[`plotReceptors()`](plotReceptors.md), [`spatialExpressionPlot()`](spatialExpressionPlot.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
plotReceptorsExpression(x, group = "clone_id", feature = "gene1")
```
