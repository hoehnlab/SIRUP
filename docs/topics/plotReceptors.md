# plotReceptors

## Plot receptor groups on a tissue image

## Description

`plotReceptors` plots the spots that hold receptors. It colors the spots by receptor group, or by a continuous value.

## Usage

``` R
plotReceptors(
  object,
  receptors = NULL,
  group = NULL,
  group_val = NULL,
  spatial_id = NULL,
  colors = NULL,
  overlap = "average",
  value = NULL,
  value_colors = c("#FFFFFF", "#000080"),
  value_limits = c(0, NA),
  value_combine = c("mean", "median", "mode", "max", "min", "sum"),
  receptor_size = 1,
  receptor_alpha = 1,
  title = NULL,
  order_points = TRUE,
  image_scale = NULL,
  swap_coords = NULL,
  images = NULL,
  sample_id = NULL,
  image_alpha = 1,
  zoom = FALSE,
  receptor_legend_title = NULL,
  show_receptor_legend = TRUE,
  polygons = FALSE,
  polygon_x = "x",
  polygon_y = "y",
  polygon_id = "polygon",
  polygon_val = NULL,
  polygon_alpha = 0.2,
  polygon_fill = "white",
  polygon_color = "white",
  polygon_linewidth = 0.6,
  polygon_label = FALSE,
  polygon_label_by = NULL,
  polygon_label_size = 3,
  polygon_label_color = "white",
  combine = TRUE,
  combine_sample_plots = FALSE,
  slot = NULL
)
```

## Arguments

`object`

:   Seurat object, or a SIRUP object with spatial data and receptors.

`receptors`

:   AIRR-style receptor table with a spatial ID column and a group column. Omit it if `object` is a SIRUP object.

`group`

:   Name of the `receptors` column that selects and colors the receptor groups, for example `"clone_id"`.

`group_val`

:   Values of `group` to plot. If `NULL`, the function uses all non-missing values.

`spatial_id`

:   Name of the `receptors` column that holds the Seurat spot names. If `NULL`, the function uses the stored field of a SIRUP object, or uses the first standard spatial ID column that it finds.

`colors`

:   Colors of the receptor groups. Named colors match the group values. Unnamed colors follow the order of the values.

`overlap`

:   Color of a spot that holds more than one receptor group. The function draws such a spot once, as a combined group with the member names joined by `" + "`. With `"average"` (default), the spot gets the mean color of its members. With one color, all combined spots get that color.

`value`

:   Color each spot on a continuous gradient and not by group. Give the name of a numeric `receptors` column, or `"count"` for the number of receptors in the spot. You can omit `group` in this mode. If you give it (with `group_val`), the function uses only the matching receptors.

`value_colors`

:   Colors of the `value` gradient (default white to navy). Two colors give a low-to-high gradient. More colors give a gradient with several stops, so you can give a palette such as `viridisLite::viridis(9)`.

`value_limits`

:   Limits of the `value` scale (default `c(0, NA)`, which starts the scale at zero).

`value_combine`

:   How to combine `value` if a spot holds more than one receptor: `"mean"` (default), `"median"`, `"mode"`, `"max"`, `"min"`, or `"sum"`. The function ignores it if `value = "count"`.

`receptor_size`

:   Size of the receptor points (default 1).

`receptor_alpha`

:   Opacity of the receptor points (default 1).

`title`

:   Plot title. Give one title for each plot, in the order of the plots, to title every plot of a multi-sample object.

`order_points`

:   Draw the early receptor groups first, so that the later groups show on top (default `TRUE`)?

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

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

`polygons`

:   Draw polygons under the points on each spatial panel. With `FALSE` (default), the function draws no polygons. With `TRUE`, it uses the polygons stored in the SIRUP object. You can also give a data frame of polygon vertices, or the path of a CSV file of polygon vertices. A vertex table needs the columns set by `polygon_x` and `polygon_y`. It needs the column set by `polygon_id` if it holds more than one polygon.

`polygon_x`

:   Name of the x-coordinate column in the polygon table (default `"x"`).

`polygon_y`

:   Name of the y-coordinate column in the polygon table (default `"y"`).

`polygon_id`

:   Name of the column that separates the polygons in the polygon table (default `"polygon"`). If the column does not exist, the function uses a `polygon` or `polygon_id` column. If neither exists, all rows form one polygon.

`polygon_val`

:   IDs of the polygons to draw. If `NULL`, the function draws all polygons.

`polygon_alpha`

:   Opacity of the polygon fill (default 0.2).

`polygon_fill`

:   Fill color of the polygons (default `"white"`).

`polygon_color`

:   Color of the polygon outlines (default `"white"`).

`polygon_linewidth`

:   Width of the polygon outlines (default 0.6).

`polygon_label`

:   Write a label at the center of each polygon (default `FALSE`)?

`polygon_label_by`

:   Name of the column that holds the labels. If `NULL`, the function uses the polygon ID.

`polygon_label_size`

:   Text size of the polygon labels (default 3).

`polygon_label_color`

:   Text color of the polygon labels (default `"white"`).

`combine`

:   Combine the image panels with patchwork (default `TRUE`)?

`combine_sample_plots`

:   Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

`slot`

:   Name of the receptor slot to use. Leave it `NULL` if the SIRUP object has one receptor table. You must give it if the object has several slots. The function ignores `slot` if you give a plain receptor table.

## Details

#\' If the object has multiple images or samples, the function draws one plot for each image or sample.

## Value

A ggplot object or a patchwork object. If `combine = FALSE`, a list.

## See Also

[`plotReceptorsExpression()`](plotReceptorsExpression.md), [`plotSpatialColorTips()`](plotSpatialColorTips.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
plotReceptors(x, group = "clone_id")

# Color the spots by the number of receptors in each spot.
plotReceptors(x, value = "count")
```
