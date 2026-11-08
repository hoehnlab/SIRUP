# spatialExpressionMultiPlot

## Plot several features as ordered color layers

## Description

`spatialExpressionMultiPlot` draws each feature in its own color on one tissue image. The features are layers. A feature that the function draws later shows on top.

## Usage

``` R
spatialExpressionMultiPlot(
  object,
  features,
  feat_order = NULL,
  combine = TRUE,
  combine_sample_plots = FALSE,
  assay = NULL,
  layer = "counts",
  colors = NULL,
  breadth_metric = c("n_expressing"),
  images = NULL,
  sample_id = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  point_size = 0.6,
  point_alpha = 1,
  order_points = TRUE,
  zero_tol = 0,
  image_alpha = 1,
  zoom = FALSE
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`features`

:   Names of the features to plot.

`feat_order`

:   Order of the layers. Put broad features (low priority) first and rare features (high priority) last. If `NULL`, the function orders the features with `breadth_metric`.

`combine`

:   Combine the spatial panels and the legends with patchwork (default `TRUE`)?

`combine_sample_plots`

:   Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

`assay`

:   Assay to use.

`layer`

:   Seurat expression layer to use (default `"counts"`).

`colors`

:   Colors, one for each feature.

`breadth_metric`

:   Measure that sets the default order of the features (default `"n_expressing"`, the number of spots that express the feature).

`images`

:   Names of the Seurat images to use.

`sample_id`

:   Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`point_size`

:   Size of the spots. The default is 0.6.

`point_alpha`

:   Opacity of the spots.

`order_points`

:   Draw the low-priority features first, so that the high-priority features show on top (default `TRUE`)?

`zero_tol`

:   Values at or below this threshold count as zero (default 0). Set to -Inf to plot all spots.

`image_alpha`

:   Opacity of the tissue image.

`zoom`

:   Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot shows the full tissue area that the spatial coordinates cover.

## Details

If the object has multiple images or samples, the function draws one plot for each image or sample.

## Value

A ggplot object or a patchwork object. If `combine = FALSE`, a list.

## See Also

[`spatialExpressionPlot()`](spatialExpressionPlot.md), [`spatialExpressionBlendPlot()`](spatialExpressionBlendPlot.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
spatialExpressionMultiPlot(x, c("gene1", "gene2", "gene3"))
```
