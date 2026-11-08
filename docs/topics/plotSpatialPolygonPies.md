# plotSpatialPolygonPies

## Plot clonal dominance results within polygons as pies or bars

## Description

`plotSpatialPolygonPies` summarizes one receptor column for each polygon. It plots the result as pie charts or stacked bars.

## Usage

``` R
plotSpatialPolygonPies(
  dominance,
  object,
  fill = "clone_id",
  mode = c("side_by_side", "spatial_clipped", "dominance_only", "stacked_bar"),
  top_n = 10,
  label_top_n = 0,
  other_label = "Other",
  other_color = "#D3D3D3",
  plot_other = TRUE,
  images = NULL,
  sample_id = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  polygon_val = NULL,
  polygon_alpha = 0.2,
  image_alpha = 1,
  pie_alpha = 0.9,
  label_size = 3,
  pies_per_page = 6L,
  pie_order = c("receptors", "polygon"),
  spatial_legend = FALSE,
  zoom = FALSE,
  combine_sample_plots = FALSE
)
```

## Arguments

`dominance`

:   Result of [`spatialCloneDominance()`](spatialCloneDominance.md).

`object`

:   Seurat object or SIRUP object.

`fill`

:   Name of the column in `dominance$region_receptors` to summarize. If you omit it, the function uses the clone column that [`spatialCloneDominance()`](spatialCloneDominance.md) recorded (default `"clone_id"`).

`mode`

:   Plot style: `side_by_side`, `spatial_clipped`, `dominance_only`, or `stacked_bar`. See Details.

`top_n`

:   Number of the most frequent `fill` levels to find in each polygon (default 10). The function shows the union of these levels on their own in all polygons. It groups the other levels under `other_label`.

`label_top_n`

:   Number of levels to label in each polygon (default 0).

`other_label`

:   Label for the grouped levels (default `"Other"`).

`other_color`

:   Color of the grouped levels (default `"#D3D3D3"`).

`plot_other`

:   Plot the grouped levels (default `TRUE`)?

`images`

:   Names of the Seurat images to use.

`sample_id`

:   Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`polygon_val`

:   IDs of the polygons to select from `dominance`.

`polygon_alpha`

:   Opacity of the polygon fill in the numbered map (default 0.2).

`image_alpha`

:   Opacity of the tissue image.

`pie_alpha`

:   Opacity of the pies (default 0.9).

`label_size`

:   Text size of the polygon labels and pie labels (default 3).

`pies_per_page`

:   Largest number of pie charts on one page, if `mode = "side_by_side"` (default 6). If there are more polygons, the function returns a named list with one plot for each page.

`pie_order`

:   Order of the pie charts, if `mode = "side_by_side"`: `"receptors"` (default) puts the polygon with the most receptors first, then the next, and so on. `"polygon"` keeps the order of the numbered map. Each pie keeps its map label. The function ignores this argument if `mode = "spatial_clipped"`.

`spatial_legend`

:   Show one shared fill legend for the `spatial_clipped` plots (default `FALSE`)?

`zoom`

:   Zoom to the extent of the plotted polygons (default `FALSE`)? If `FALSE`, the plot shows the full tissue image.

`combine_sample_plots`

:   Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

## Details

The function uses the receptor rows in the result of [`spatialCloneDominance()`](spatialCloneDominance.md). The `mode` argument sets the plot style:

- `side_by_side`: pie charts next to a numbered map of the polygons.

- `spatial_clipped`: pie charts inside the polygon boundaries on the tissue. It needs the `sf` package.

- `dominance_only`: the pie charts only, with no map.

- `stacked_bar`: one stacked bar for each polygon in one panel, with no map. Each bar stacks the `fill` levels in dominance order.

If there are multiple samples or images, the function creates separate plots for each.

## Value

A ggplot object, a patchwork object, or a named list of ggplot objects.

## See Also

[`spatialCloneDominance()`](spatialCloneDominance.md), [`plotSpatialRegionReceptors()`](plotSpatialRegionReceptors.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
dominance <- spatialCloneDominance(x)
plotSpatialPolygonPies(dominance, x)
plotSpatialPolygonPies(dominance, x, mode = "stacked_bar")
```
