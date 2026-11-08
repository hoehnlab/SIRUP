# spatialDiscretePlot

## Plot a metadata column on a tissue image

## Description

`spatialDiscretePlot` colors each spot by the value of a metadata column. It leaves out the spots that have the \"zero_value\" value.

## Usage

``` R
spatialDiscretePlot(
  object,
  meta_col,
  zero_value = 0,
  images = NULL,
  sample_id = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  cols = NULL,
  point_size = 1.6,
  point_alpha = 1,
  order_points = TRUE,
  image_alpha = 1,
  zoom = FALSE,
  combine = TRUE,
  combine_sample_plots = FALSE
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`meta_col`

:   Name of the metadata column to plot.

`zero_value`

:   Value that means \"not present\" (default 0). The function does not draw these spots.

`images`

:   Names of the Seurat images to use.

`sample_id`

:   Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`cols`

:   Named vector of colors, one for each value.

`point_size`

:   Size of the spots. The default is 1.6.

`point_alpha`

:   Opacity of the spots.

`order_points`

:   Draw the low numeric values, or the early factor levels, first, so that the high or late values show on top (default `TRUE`)?

`image_alpha`

:   Opacity of the tissue image.

`zoom`

:   Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot shows the full tissue area that the spatial coordinates cover.

`combine`

:   Combine the panels with patchwork, if possible (default `TRUE`)?

`combine_sample_plots`

:   Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

## Details

#\' If the object has multiple images or samples, the function draws one plot for each image or sample.

## Value

A ggplot object, a patchwork object, or a list of plots. The function returns `NULL` if all spots have the blank value.

## See Also

[`spatialExpressionPlot()`](spatialExpressionPlot.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
spatialDiscretePlot(x, "cell_type")
```
