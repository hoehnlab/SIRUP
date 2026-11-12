# plotSpatialColorTips

## Plot tip colors on a tissue image

## Description

`plotSpatialColorTips` plots the spots of tree tips. Each spot gets the color of its tip.

## Usage

``` R
plotSpatialColorTips(
  object,
  tip_data,
  spatial_id = "spatial_id_8um",
  sample_col = NULL,
  color = "tip_color",
  label = NULL,
  legend_title = "Tip color",
  legend_colors = NULL,
  title = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  point_size = 1.6,
  point_alpha = 1,
  order_points = TRUE,
  images = NULL,
  sample_id = NULL,
  image_alpha = 1,
  zoom = FALSE,
  show_legend = FALSE,
  polygons = FALSE,
  polygon_x = "x",
  polygon_y = "y",
  polygon_id = "polygon",
  polygon_val = NULL,
  polygon_alpha = 0.2,
  polygon_fill = "#0072B2",
  polygon_color = "#0072B2",
  polygon_linewidth = 0.6,
  polygon_label = FALSE,
  polygon_label_by = NULL,
  polygon_label_size = 3,
  polygon_label_color = "black",
  combine = TRUE,
  combine_sample_plots = FALSE
)
```

## Arguments

`object`

: Seurat object or SIRUP object.

`tip_data`

: Data frame with a spatial ID column and a color column. For example, use the result of [`assignTipColors()`](assignTipColors.md).

`spatial_id`

: Name of the column in `tip_data` that holds the Seurat spot names (default `"spatial_id_8um"`).

`sample_col`

: Name of the column in `tip_data` that matches the names of a list of spatial objects.

`color`

: Name of the column that holds the colors (default `"tip_color"`).

`label`

: Name of the column that holds the legend labels of the colors.

`legend_title`

: Title of the color legend (default `"Tip color"`).

`legend_colors`

: Named vector of colors for the legend.

`title`

: Plot title.

`image_scale`

: Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

: Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`point_size`

: Size of the spots. The default is 1.6.

`point_alpha`

: Opacity of the spots.

`order_points`

: Draw the early color groups first, so that the later groups show on top (default `TRUE`)?

`images`

: Names of the Seurat images to use.

`sample_id`

: Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`image_alpha`

: Opacity of the tissue image.

`zoom`

: Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot shows the full tissue area that the spatial coordinates cover.

`show_legend`

: Show the color legend (default `FALSE`)?

`polygons`

: Draw polygons under the points on each spatial panel. With `FALSE` (default), the function draws no polygons. With `TRUE`, it uses the polygons stored in the SIRUP object. You can also give a data frame of polygon vertices, or the path of a CSV file of polygon vertices. A vertex table needs the columns set by `polygon_x` and `polygon_y`. It needs the column set by `polygon_id` if it holds more than one polygon.

`polygon_x`

: Name of the x-coordinate column in the polygon table (default `"x"`).

`polygon_y`

: Name of the y-coordinate column in the polygon table (default `"y"`).

`polygon_id`

: Name of the column that separates the polygons in the polygon table (default `"polygon"`). If the column does not exist, the function uses a `polygon` or `polygon_id` column. If neither exists, all rows form one polygon.

`polygon_val`

: IDs of the polygons to draw. If `NULL`, the function draws all polygons.

`polygon_alpha`

: Opacity of the polygon fill (default 0.2).

`polygon_fill`

: Fill color of the polygons (default `"#0072B2"`).

`polygon_color`

: Color of the polygon outlines (default `"#0072B2"`).

`polygon_linewidth`

: Width of the polygon outlines (default 0.6).

`polygon_label`

: Write a label at the center of each polygon (default `FALSE`)?

`polygon_label_by`

: Name of the column that holds the labels. If `NULL`, the function uses the polygon ID.

`polygon_label_size`

: Text size of the polygon labels (default 3).

`polygon_label_color`

: Text color of the polygon labels (default `"black"`).

`combine`

: Combine the image panels with patchwork (default `TRUE`)?

`combine_sample_plots`

: Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

## Details

If the object has multiple images or samples, the function draws one plot for each image or sample.

## Value

A ggplot object or a patchwork object.

## See Also

[`plotSpatialCloneColorTips()`](plotSpatialCloneColorTips.md), [`plotSpatialCloneTreePairs()`](plotSpatialCloneTreePairs.md), [`assignTipColors()`](assignTipColors.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
colors <- assignTipColors(getReceptors(x), shared_color = "spatial_id")
plotSpatialColorTips(x, tip_data = colors, spatial_id = "spatial_id")
```
