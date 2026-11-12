# plotSpatialPolygons

## Plot polygons on a tissue image

## Description

Plot polygons on a tissue image

## Usage

``` R
plotSpatialPolygons(
  object,
  polygon = NULL,
  images = NULL,
  sample_id = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  polygon_val = NULL,
  x = "x",
  y = "y",
  polygon_id = "polygon",
  polygon_alpha = 0.2,
  polygon_fill = "#0072B2",
  fill_analysis = NULL,
  fill_column = NULL,
  analysis_fill = c("#D3D3D3", "#FFFF00"),
  polygon_color = "#0072B2",
  linewidth = 0.6,
  image_alpha = 1,
  label = TRUE,
  label_by = "polygon",
  label_size = 3,
  label_color = "black",
  label_outside = FALSE,
  label_offset = 0.08,
  label_leader = TRUE,
  label_leader_color = "grey40",
  zoom = FALSE,
  combine = TRUE,
  combine_sample_plots = FALSE
)
```

## Arguments

`object`

: Seurat object or SIRUP object.

`polygon`

: Data frame of polygon vertices in order. It needs the columns set by `x` and `y`. If the column set by `polygon_id` exists, the function treats each polygon on its own. If `object` is a SIRUP object and `polygon` is `NULL`, the function uses the stored polygons.

`images`

: Names of the Seurat images to use.

`sample_id`

: Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`image_scale`

: Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

: Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`polygon_val`

: IDs of the polygons to use. If `NULL`, the function uses all polygons.

`x`

: Name of the x-coordinate column in `polygon` (default `"x"`).

`y`

: Name of the y-coordinate column in `polygon` (default `"y"`).

`polygon_id`

: Name of the column in `polygon` that groups the vertices by polygon (default `"polygon"`).

`polygon_alpha`

: Opacity of the polygon fill (default 0.2).

`polygon_fill`

: Fill color of the polygons, if no analysis is mapped (default `"#0072B2"`).

`fill_analysis`

: Name of a stored SIRUP analysis that has one row for each polygon. The function colors the polygons with it.

`fill_column`

: Name of the column that sets the fill color. If `fill_analysis` is `NULL`, the function reads the column from the polygon data. Otherwise, it reads the column from the analysis. The analysis must have a `polygon_id` or `polygon` column as the polygon key.

`analysis_fill`

: Two colors for the low and high end of the gradient for numeric analysis values (default `c("#D3D3D3", "#FFFF00")`). The function uses the Glasbey palette for other values.

`polygon_color`

: Color of the polygon outlines (default `"#0072B2"`).

`linewidth`

: Width of the polygon outlines (default 0.6).

`image_alpha`

: Opacity of the tissue image.

`label`

: Write a label at the center of each polygon (default `TRUE`)?

`label_by`

: Name of the column in the polygon data that holds the labels (default `"polygon"`).

`label_size`

: Text size of the labels (default 3).

`label_color`

: Text color of the labels (default `"black"`).

`label_outside`

: Put each label outside its polygon (default `FALSE`)? The function pushes each label away from the centroid of all plotted polygons, so that labels spread outward from a group of small polygons. The label is `label_offset` away from the polygon.

`label_offset`

: Gap between a polygon edge and its label, if `label_outside` is `TRUE` (default 0.08). The value is a fraction of the plotted coordinate range.

`label_leader`

: Draw a line from each polygon to its label, if `label_outside` is `TRUE` (default `TRUE`)?

`label_leader_color`

: Color of the leader lines (default `"grey40"`).

`zoom`

: Zoom to the extent of the polygons, not the full tissue (default `FALSE`)?

`combine`

: Combine the plots with patchwork, if possible (default `TRUE`)?

`combine_sample_plots`

: Combine the plots of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns a named list with one plot for each sample.

## Details

`plotSpatialPolygons` plots polygons over the Seurat tissue image. If there are multiple samples or images, the function creates separate plots for each.

## Value

A ggplot or patchwork object. If `combine = FALSE`, a named list of ggplot objects.

## See Also

[`selectSpatialPolygons()`](selectSpatialPolygons.md), [`spatialCloneDominance()`](spatialCloneDominance.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
plotSpatialPolygons(x)
plotSpatialPolygons(x, label_outside = TRUE)
```
