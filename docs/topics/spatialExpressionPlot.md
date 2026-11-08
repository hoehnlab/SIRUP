# spatialExpressionPlot

## Plot gene expression on a tissue image

## Description

`spatialExpressionPlot` plots the expression of one or more features. By default, it draws only the spots with expression above zero.

## Usage

``` R
spatialExpressionPlot(
  object,
  features,
  assay = NULL,
  layer = "counts",
  images = NULL,
  sample_id = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  feature_colors = NULL,
  feature_color_palette = "viridis",
  zero_tol = 0,
  point_size = 0.6,
  point_alpha = 1,
  order_points = TRUE,
  title = NULL,
  image_alpha = 1,
  zoom = FALSE,
  combine = TRUE,
  combine_sample_plots = FALSE
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`features`

:   Names of the features to plot.

`assay`

:   Assay to use.

`layer`

:   Seurat expression layer to use (default `"counts"`).

`images`

:   Names of the Seurat images to use.

`sample_id`

:   Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`feature_colors`

:   Colors for the continuous expression scale, from the low value to the high value. If you give them, they replace `feature_color_palette`.

`feature_color_palette`

:   Color palette for the expression scale, if `feature_colors` is `NULL`: `"viridis"` (default), `"inferno"`, `"plasma"`, `"cividis"`, `"rocket"`, or `"mako"`.

`zero_tol`

:   Values at or below this threshold count as zero (default 0). Set to -Inf to plot all spots.

`point_size`

:   Size of the spots. The default is 0.6.

`point_alpha`

:   Opacity of the spots.

`order_points`

:   Draw the low values first, so that the high values show on top (default `TRUE`)?

`title`

:   Plot title. For one panel, it replaces the feature name. For several panels that the function combines, it is the title of the whole patchwork, and each panel keeps its feature name.

`image_alpha`

:   Opacity of the tissue image.

`zoom`

:   Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot shows the full tissue area that the spatial coordinates cover.

`combine`

:   Combine the panels into one figure with patchwork, if there are several images or features (default `TRUE`)?

`combine_sample_plots`

:   Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

## Details

If the object has multiple images or samples, the function draws one plot for each image or sample. If given several features, the function draws one panel for each feature.

## Value

A ggplot object, a patchwork object, or a list of plots. The function returns `NULL` if all spots have zero expression.

## See Also

[`spatialExpressionBlendPlot()`](spatialExpressionBlendPlot.md), [`spatialExpressionMultiPlot()`](spatialExpressionMultiPlot.md), [`plotReceptorsExpression()`](plotReceptorsExpression.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
spatialExpressionPlot(x, "gene1")

# Plot two features. The function combines the panels with patchwork.
spatialExpressionPlot(x, c("gene1", "gene2"))
```
