# plotSpatialLocalMoran

## Plot local Moran\'s I on a tissue image

## Description

`plotSpatialLocalMoran` colors each spot by its local Moran\'s I result.

## Usage

``` R
plotSpatialLocalMoran(
  object,
  local_moran,
  value_col = "z_i",
  images = NULL,
  sample_id = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  colors = c("#2c7bb6", "#f7f7f7", "#d7191c"),
  midpoint = 0,
  point_size = 0.7,
  point_alpha = 1,
  order_points = TRUE,
  image_alpha = 1,
  zoom = FALSE,
  title = NULL,
  combine = TRUE,
  combine_sample_plots = FALSE
)
```

## Arguments

`object`

:   SIRUP object, or Seurat object with spatial images.

`local_moran`

:   Result of [`spatialLocalMoran()`](spatialLocalMoran.md).

`value_col`

:   Name of the column to color by (default `"z_i"`).

`images`

:   Names of the Seurat images to use.

`sample_id`

:   Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`colors`

:   Three colors for the low, middle, and high values.

`midpoint`

:   Value that gets the middle color.

`point_size`

:   Size of the spots.

`point_alpha`

:   Opacity of the spots.

`order_points`

:   Draw low values first, so that high values show on top (default `TRUE`)?

`image_alpha`

:   Opacity of the tissue image.

`zoom`

:   Zoom to the plotted spots (default `FALSE`)?

`title`

:   Plot title.

`combine`

:   Combine the image panels with patchwork (default `TRUE`)?

`combine_sample_plots`

:   Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

## Details

If the data come from more than one image, or a sirup object with multiple samples, then the function builds a list of plots, one for each sample or image.

## Value

A ggplot or patchwork object.

## See Also

[`spatialLocalMoran()`](spatialLocalMoran.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
local_moran <- spatialLocalMoran(x, "gene1", k = 4)
plotSpatialLocalMoran(x, local_moran)
```
