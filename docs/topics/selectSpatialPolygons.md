# selectSpatialPolygons

## Select polygon vertices in an interactive app

## Description

`selectSpatialPolygons` opens a Shiny app. You click polygon vertices on a plot of the spatial coordinates.

## Usage

``` R
selectSpatialPolygons(
  object,
  feature = NULL,
  assay = "Spatial",
  layer = "counts",
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  save_path = NULL,
  plot_image = FALSE,
  downsample_n = NULL,
  downsample_seed = 1L,
  point_size = 0.8,
  point_alpha = 0.9,
  show_zeros = TRUE,
  launch_browser = TRUE,
  host = getOption("shiny.host", "127.0.0.1"),
  port = getOption("shiny.port", NULL),
  quiet = TRUE
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`feature`

:   Feature to show as a color scale under the vertices. You can change it in the app. The function looks in the assay first. If the feature is not there, it looks in the object metadata. This lets you show module scores (from `AddModuleScore()`, `UCell`, or any numeric column) without adding them to an assay.

`assay`

:   Assay to set before the function gets the feature values (default `"Spatial"`).

`layer`

:   Expression layer to use for `feature` (default `"counts"`).

`images`

:   Names of the Seurat images to use.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`save_path`

:   Path of a CSV file to save the vertices to.

`plot_image`

:   Show the Seurat tissue image under the points (default `FALSE`)?

`downsample_n`

:   Number of points to show in the app. If `NULL`, the app shows all points. A click still snaps to the nearest spatial ID in the full coordinate table.

`downsample_seed`

:   Random seed for the downsampling (default 1).

`point_size`

:   Size of the points (default 0.8).

`point_alpha`

:   Starting opacity of the points (default 0.9).

`show_zeros`

:   Show the points where `feature` is zero or `NA` (default `TRUE`)? On sparse Visium HD data, the zeros can hide the signal. Use `FALSE`, or clear \"Show zero points\" in the app, to show only the positive points. This changes the display only. A click can still select a hidden point as a vertex.

`launch_browser`

:   Passed to `shiny::runApp()` (default `TRUE`).

`host`

:   Host for `shiny::runApp()`. Set it when you run the function in an interactive HPC session through an SSH tunnel.

`port`

:   Port for `shiny::runApp()`. Set it when you run the function in an interactive HPC session through an SSH tunnel.

`quiet`

:   Hide warnings while the app runs (default `TRUE`)?

## Details

Each click snaps to the nearest spatial ID. Click \"Close polygon\" or \"New polygon\" before you finish. This repeats the first vertex at the end of each polygon.

The vertices can describe overlapping or nested polygons. The function that receives the vertices decides the membership: [`assignSpatialPolygons()`](assignSpatialPolygons.md), [`assignReceptorPolygons()`](assignReceptorPolygons.md), [`createSIRUP()`](createSIRUP.md), or [`addPolygons()`](sirup-polygons.md). The primary `polygon_id` gets the ID of the smallest matching polygon by area. The logical `⁠polygon_<ID>⁠` columns show membership in each polygon.

## Value

A data frame with the columns `spatial_id`, `x`, `y`, `num`, `polygon`, and `image`.

## See Also

[`formatSpatialPolygons()`](formatSpatialPolygons.md), [`addPolygons()`](sirup-polygons.md), [`plotSpatialPolygons()`](plotSpatialPolygons.md)

## Examples

``` R

x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Click the vertices in the app. Use "Close polygon" when you finish.
polygons <- selectSpatialPolygons(
  getSpatial(x),
  image_scale = "hires",
  swap_coords = FALSE
)
```
