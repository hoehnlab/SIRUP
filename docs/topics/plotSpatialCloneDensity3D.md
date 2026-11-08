# plotSpatialCloneDensity3D

## Plot 3D clone-density surfaces on a tissue image

## Description

`plotSpatialCloneDensity3D` draws a smoothed density surface for each clone above the tissue image in an interactive rgl scene.

## Usage

``` R
plotSpatialCloneDensity3D(
  dominance,
  object,
  clone = "clone_id",
  weight = NULL,
  images = NULL,
  sample_id = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  polygon_val = NULL,
  top_n = 10,
  other_label = "Other",
  other_color = "#D3D3D3",
  plot_other = TRUE,
  bandwidth = NULL,
  grid_n = 80,
  plot_type = c("density", "histogram"),
  hist_bins = 30,
  zoom = FALSE,
  height_scale = 0.22,
  z_offset = 0.04,
  density_alpha = 0.55,
  density_floor = 0.03,
  image_alpha = 1,
  show_points = FALSE,
  point_size = 5,
  point_alpha = 0.75,
  show_polygon = TRUE,
  polygon_color = "white",
  polygon_lwd = 3,
  label_density_peaks = FALSE,
  label_top_n = 5,
  label_cex = 0.8,
  show_axes = TRUE,
  texture_dir = NULL,
  as_widget = TRUE,
  return_data = FALSE
)
```

## Arguments

`dominance`

:   Result of [`spatialCloneDominance()`](spatialCloneDominance.md).

`object`

:   Seurat object or SIRUP object.

`clone`

:   Name of the column in `dominance$region_receptors` that defines the clone groups. If you omit it, the function uses the clone column that [`spatialCloneDominance()`](spatialCloneDominance.md) recorded (default `"clone_id"`).

`weight`

:   Name of a numeric column that gives the weight of each receptor. If `NULL`, each receptor row counts as one.

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

`top_n`

:   Number of the most abundant clone groups to show on their own (default 10). The function groups the other clone groups under `other_label`.

`other_label`

:   Label for the grouped clone groups (default `"Other"`).

`other_color`

:   Color of the grouped clone groups (default `"#D3D3D3"`).

`plot_other`

:   Plot the grouped clone groups (default `TRUE`)?

`bandwidth`

:   Bandwidth of the Gaussian kernel, in coordinate units. Give one value for equal smoothing in x and y, or two values (x, y). If `NULL`, the function estimates it from the coordinates in the regions.

`grid_n`

:   Number of grid points for the density surfaces (default 80). Give one value for both axes, or two values (x, y).

`plot_type`

:   Type of plot: `density` (smoothed surfaces) or `histogram` (binned 3D bars).

`hist_bins`

:   Number of x and y bins, if `plot_type = "histogram"` (default 30).

`zoom`

:   Zoom the image plane and the 3D view to the selected region (default `FALSE`)? If `FALSE`, the plot shows the full image.

`height_scale`

:   Maximum height of the density, as a fraction of the plotted range (default 0.22).

`z_offset`

:   Vertical gap between the density surfaces and the tissue image, as a fraction of the plotted range (default 0.04).

`density_alpha`

:   Opacity of the density surfaces (default 0.55).

`density_floor`

:   Relative density below which the function hides surface cells (default 0.03). This removes the low Gaussian tails, which can look like a colored plane.

`image_alpha`

:   Opacity of the tissue image plane.

`show_points`

:   Show the receptor locations above the tissue plane (default `FALSE`)?

`point_size`

:   Size of the receptor points (default 5).

`point_alpha`

:   Opacity of the receptor points (default 0.75).

`show_polygon`

:   Show the polygon boundaries above the tissue plane (default `TRUE`)?

`polygon_color`

:   Color of the polygon boundary lines (default `"white"`).

`polygon_lwd`

:   Width of the polygon boundary lines (default 3).

`label_density_peaks`

:   Label the highest point of the first density surfaces (default `FALSE`)?

`label_top_n`

:   Number of surfaces to label, if `label_density_peaks` is `TRUE` (default 5).

`label_cex`

:   Text scale of the peak labels (default 0.8).

`show_axes`

:   Show the 3D bounding box and the axis labels (default `TRUE`)?

`texture_dir`

:   Directory to write the tissue texture PNG file to. Use a permanent directory if you save an HTML widget.

`as_widget`

:   Return an rgl htmlwidget (default `TRUE`)? If `FALSE`, the function draws the scene and returns the ID of the current rgl device.

`return_data`

:   Return the density data and do not draw the scene (default `FALSE`)? This does not need `rgl`.

## Details

The function uses the receptor rows in the result of [`spatialCloneDominance()`](spatialCloneDominance.md). It needs the `rgl` package, except if you set `return_data = TRUE`.

## Value

An rgl widget or an rgl device ID. If `return_data = TRUE`, a list of the plot data.

## See Also

[`spatialCloneDominance()`](spatialCloneDominance.md), [`plotSpatialRegionReceptors()`](plotSpatialRegionReceptors.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
dominance <- spatialCloneDominance(x)

# Return the density data. This needs no rgl window.
density <- plotSpatialCloneDensity3D(dominance, x, return_data = TRUE)
names(density)
```
