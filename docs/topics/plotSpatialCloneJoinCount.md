# plotSpatialCloneJoinCount

## Plot a clone join-count result

## Description

`plotSpatialCloneJoinCount` plots the tissue map of one tested clone. It also shows the observed and expected number of joins between spots of the same clone.

## Usage

``` R
plotSpatialCloneJoinCount(
  join_count,
  object,
  clone_val,
  images = NULL,
  sample_id = NULL,
  polygon_val = NULL,
  polygon = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  show_neighbors = TRUE,
  clone_color = "#D55E00",
  other_color = "grey75",
  point_size = 1.8,
  point_alpha = 0.9,
  edge_width = 0.6,
  edge_alpha = 0.75,
  polygon_color = "white",
  polygon_linewidth = 0.6,
  image_alpha = 1,
  zoom = FALSE,
  title = NULL,
  combine_sample_plots = FALSE
)
```

## Arguments

`join_count`

:   Result of [`spatialCloneJoinCount()`](spatialCloneJoinCount.md). Run that function with `return_data = TRUE`.

`object`

:   Seurat object with spatial images, or a SIRUP object.

`clone_val`

:   Clone ID to plot.

`images`

:   Names of the images to plot.

`sample_id`

:   Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`polygon_val`

:   IDs of the polygons to include.

`polygon`

:   Data frame of polygon vertices to draw on the tissue.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object.

`show_neighbors`

:   Draw a line between each pair of neighboring spots that both contain the clone (default `TRUE`)?

`clone_color`

:   Color of the spots and lines for the clone.

`other_color`

:   Color of the other spots that have receptors.

`point_size`

:   Size of the spots.

`point_alpha`

:   Opacity of the spots.

`edge_width`

:   Width of the neighbor lines.

`edge_alpha`

:   Opacity of the neighbor lines.

`polygon_color`

:   Color of the polygon outlines (default `"white"`).

`polygon_linewidth`

:   Width of the polygon outlines.

`image_alpha`

:   Opacity of the tissue image.

`zoom`

:   Zoom to the spots in the test (default `FALSE`)?

`title`

:   Title of the combined plot.

`combine_sample_plots`

:   Combine the plots of all samples into one figure (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

## Details

The map uses the same spots and the same neighbor graph as the call to [`spatialCloneJoinCount()`](spatialCloneJoinCount.md).

If the data come from more than one image, or a sirup object with multiple samples, then the function builds a list of plots, one for each sample or image.

## Value

A ggplot object. For a multi-sample object and `combine_sample_plots = FALSE`, a named list of plots.

## See Also

[`spatialCloneJoinCount()`](spatialCloneJoinCount.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
join_count <- spatialCloneJoinCount(
  x,
  use_stored_polygons = FALSE,
  k = 4,
  return_data = TRUE
)
plotSpatialCloneJoinCount(join_count, x, clone_val = "2")
```
