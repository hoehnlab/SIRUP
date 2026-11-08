# spatialExpressionBlendPlot

## Blend the expression of two features on a tissue image

## Description

`spatialExpressionBlendPlot` colors each spot by the expression of two features at once. By default, it leaves out the spots where both features are zero.

## Usage

``` R
spatialExpressionBlendPlot(
  object,
  features,
  combine = TRUE,
  combine_sample_plots = FALSE,
  feature_1_alt_name = NULL,
  feature_2_alt_name = NULL,
  assay = NULL,
  layer = "counts",
  only_blended = FALSE,
  bottom_left = "#FFFFFF",
  bottom_right = "#0000FF",
  top_left = "#FFFF00",
  top_right = "#00FF00",
  zero_tol = 0,
  images = NULL,
  sample_id = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  point_size = 0.6,
  point_alpha = 1,
  order_points = TRUE,
  image_alpha = 1,
  zoom = FALSE
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`features`

:   Names of two features.

`combine`

:   Combine the panels and the legend with patchwork (default `TRUE`)?

`combine_sample_plots`

:   Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

`feature_1_alt_name`

:   Label of the first feature in the blend legend. If `NULL`, the function uses the feature name.

`feature_2_alt_name`

:   Label of the second feature in the blend legend. If `NULL`, the function uses the feature name.

`assay`

:   Assay to use.

`layer`

:   Seurat expression layer to use (default `"counts"`).

`only_blended`

:   Return only the blended panel and the legend (default `FALSE`)?

`bottom_left`

:   Color of the bottom-left corner of the blend grid (low in both features). The default is `"#FFFFFF"`.

`bottom_right`

:   Color of the bottom-right corner of the blend grid (high in the first feature). The default is `"#0000FF"`.

`top_left`

:   Color of the top-left corner of the blend grid (high in the second feature). The default is `"#FFFF00"`.

`top_right`

:   Color of the top-right corner of the blend grid (high in both features). The default is `"#00FF00"`.

`zero_tol`

:   Values at or below this threshold count as zero (default 0). Set to -Inf to plot all spots.

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

:   Draw the spots with low blended expression first, so that the spots with high expression show on top (default `TRUE`)?

`image_alpha`

:   Opacity of the tissue image.

`zoom`

:   Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot shows the full tissue area that the spatial coordinates cover.

## Details

#\' If the object has multiple images or samples, the function draws one plot for each image or sample.

## Value

A ggplot object or a patchwork object. If `combine = FALSE`, a list.

## See Also

[`spatialExpressionPlot()`](spatialExpressionPlot.md), [`spatialExpressionMultiPlot()`](spatialExpressionMultiPlot.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
spatialExpressionBlendPlot(x, c("gene1", "gene2"))
```
