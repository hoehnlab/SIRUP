# plotSpatialRegionReceptors

## Plot receptors inside polygons using the results of [`spatialCloneDominance()`](spatialCloneDominance.md)

## Description

Plot receptors inside polygons using the results of [`spatialCloneDominance()`](spatialCloneDominance.md)

## Usage

``` R
plotSpatialRegionReceptors(
  dominance,
  object,
  color = "clone_id",
  images = NULL,
  sample_id = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  polygon_val = NULL,
  point_size = 1,
  point_alpha = 0.5,
  polygon_alpha = 0.12,
  polygon_fill = "#FFFFFF",
  polygon_color = "black",
  image_alpha = 1,
  top_n = NULL,
  other_label = "Other",
  other_color = "#D3D3D3",
  plot_other = TRUE,
  zoom = FALSE,
  combine = TRUE,
  combine_sample_plots = FALSE
)
```

## Arguments

`dominance`

:   Result of [`spatialCloneDominance()`](spatialCloneDominance.md).

`object`

:   Seurat object or SIRUP object.

`color`

:   Name of the column in `dominance$region_receptors` that sets the color of the points. If you omit it, the function uses the clone column that [`spatialCloneDominance()`](spatialCloneDominance.md) recorded (default `"clone_id"`).

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

`point_size`

:   Size of the receptor points (default 1).

`point_alpha`

:   Opacity of the receptor points (default 0.5).

`polygon_alpha`

:   Opacity of the polygon fill (default 0.12).

`polygon_fill`

:   Fill color of the polygons (default `"#FFFFFF"`).

`polygon_color`

:   Color of the polygon outlines (default `"black"`).

`image_alpha`

:   Opacity of the tissue image.

`top_n`

:   Number of the most frequent `color` levels to show. The function groups the other levels under `other_label`. If `NULL`, it shows all levels.

`other_label`

:   Label for the grouped levels (default `"Other"`).

`other_color`

:   Color of the grouped levels (default `"#D3D3D3"`).

`plot_other`

:   Plot the grouped levels (default `TRUE`)?

`zoom`

:   Zoom to the extent of the selected polygons (default `FALSE`)? If `FALSE`, the plot shows the full tissue image.

`combine`

:   Combine the plots with patchwork, if possible (default `TRUE`)?

`combine_sample_plots`

:   Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

## Details

`plotSpatialRegionReceptors` plots the receptors inside the polygons on the tissue image. The result of [`spatialCloneDominance()`](spatialCloneDominance.md) supplies the receptors. If there are multiple samples or images, the function creates separate plots for each.

## Value

A ggplot or patchwork object. If `combine = FALSE`, a named list of ggplot objects.

## See Also

[`spatialCloneDominance()`](spatialCloneDominance.md), [`plotSpatialPolygonPies()`](plotSpatialPolygonPies.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
dominance <- spatialCloneDominance(x)
plotSpatialRegionReceptors(dominance, x)
```
