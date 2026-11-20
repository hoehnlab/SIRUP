# plotSpatialCloneTreePairs

## Plot each clone tree next to its tissue map

## Description

`plotSpatialCloneTreePairs` makes two plots for each clone: the lineage tree, and a tissue map of the tips. Each tip has one color in both plots, so you can find a tip of the tree on the tissue.

## Usage

``` R
plotSpatialCloneTreePairs(
  trees,
  object,
  spatial_id = NULL,
  clone_val = NULL,
  id = "sequence_id",
  sequence_id = id,
  clone = "clone_id",
  color_by = "identity",
  palette = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  images = NULL,
  sample_id = NULL,
  zoom = FALSE,
  tipsize = NULL,
  scale = 0.01,
  layout = "rectangular",
  title = TRUE,
  shape_by = ".sirup_spatial_status",
  shape_values = c(spatial = 21, not_spatial = 22),
  tree_layout = NULL,
  tree_tip_size = NULL,
  tree_linewidth = 0.5,
  tip_stroke = 0.3,
  treescale_fontsize = 3.88,
  treescale = NULL,
  treescale_position = c("bottomleft", "topleft", "midleft"),
  spatial_tip_size = 1.6,
  spatial_tip_alpha = 1,
  order_points = FALSE,
  show_spatial_legend = FALSE,
  legend_theme = NULL,
  combine_sample_plots = FALSE,
  polygons = FALSE,
  polygon_x = "x",
  polygon_y = "y",
  polygon_id = "polygon",
  polygon_val = NULL,
  polygon_alpha = 0.2,
  polygon_fill = "white",
  polygon_color = "white",
  polygon_linewidth = 0.6,
  polygon_label = FALSE,
  polygon_label_by = NULL,
  polygon_label_size = 3,
  polygon_label_color = "white",
  nproc = 1L,
  slot = NULL,
  ...
)
```

## Arguments

`trees`

:   Result of [`dowser::getTrees()`](https://dowser.readthedocs.io/en/latest/topics/getTrees/), or a list with an element `trees` that holds this result.

`object`

:   Seurat object or SIRUP object.

`spatial_id`

:   Name of the column in the dowser clone data that holds the Seurat spot names. Omit it if `object` is a SIRUP object.

`clone_val`

:   IDs of the clones to plot. If `NULL`, the function plots all clones.

`id`

:   Name of the column in the dowser clone data that matches the tree tip labels (default `"sequence_id"`).

`sequence_id`

:   Name of the column in the receptor table of `object` that decides which tree tips belong to the current spatial object (default: the value of `id`). The function uses it only if `object` is a SIRUP object. A tip that is not in this column stays in the tree plot, but the map does not show it.

`clone`

:   Name of the column that holds the clone IDs (default `"clone_id"`).

`color_by`

:   Column that sets the colors of the tips in both plots. With `"identity"` (default), each tip gets its own color. With another column, such as `"tree_cluster"` or a spatial ID column, tips that share a value get one color.

`palette`

:   Colors to use. If you give fewer colors than the function needs, it makes more colors by averaging the colors that you give.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`images`

:   Names of the Seurat images to use.

`sample_id`

:   Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`zoom`

:   Zoom the maps to the clone (default `FALSE`)? If `FALSE`, the maps show the full tissue area that the spatial coordinates cover.

`tipsize`

:   Size of the tip points of the tree.

`scale`

:   Distance that the scale bar of the tree shows, in mutations for each site (default 0.01).

`layout`

:   Layout of the tree, for `ggtree::ggtree()` (default `"rectangular"`).

`title`

:   Title of each tree plot. With `TRUE` (default), the title is the clone ID. Give a character string to use it as the title.

`shape_by`

:   Name of the column in the dowser clone data that sets the shape of the tree tips. The default marks a tip as `"spatial"` if it has a spatial ID, and as `"not_spatial"` if it has none.

`shape_values`

:   Named vector of shapes, for [`ggplot2::scale_shape_manual()`](https://ggplot2.tidyverse.org/reference/scale_shape_manual.html) (default `c(spatial = 21, not_spatial = 22)`). Shapes 21 to 25 are filled shapes: the tip color fills the shape, and the outline is black. Other shapes get the tip color on the whole mark.

`tree_layout`

:   Alias for `layout`.

`tree_tip_size`

:   Alias for `tipsize`.

`tree_linewidth`

:   Line width of the tree branches, for `ggtree::ggtree()` (default 0.5).

`tip_stroke`

:   Outline width of the filled tree tips (default 0.3).

`treescale_fontsize`

:   Font size of the label of the scale bar (default 3.88).

`treescale`

:   Alias for `scale`.

`treescale_position`

:   Position of the scale bar: `"bottomleft"` (default), `"topleft"` (near the root, at the top of the tips), or `"midleft"` (near the root, halfway down the tips).

`spatial_tip_size`

:   Size of the tip points on the map (default 1.6).

`spatial_tip_alpha`

:   Opacity of the tip points on the map (default 1).

`order_points`

:   Draw the early color groups first on the map, so that the later groups show on top (default `FALSE`)?

`show_spatial_legend`

:   Show a color legend below the map, based on the tip colors of the tree (default `FALSE`)?

`legend_theme`

:   A ggplot2 theme, made with `ggplot2::theme()`, that changes the look of the shared legend (default `NULL`). Use it, for example, to make the legend text and keys smaller, so that a long legend fits on one line. The function uses it only if `show_spatial_legend` is `TRUE`.

`combine_sample_plots`

:   Combine the maps of all samples into one patchwork (default `FALSE`)? The function uses it only for multi-sample objects. If `FALSE`, it returns one tree and map pair for each sample.

`polygons`

:   Draw polygons under the points on each spatial panel. With `FALSE` (default), the function draws no polygons. With `TRUE`, it uses the polygons stored in the SIRUP object. You can also give a data frame of polygon vertices, or the path of a CSV file of polygon vertices. A vertex table needs the columns set by `polygon_x` and `polygon_y`. It needs the column set by `polygon_id` if it holds more than one polygon.

`polygon_x`

:   Name of the x-coordinate column in the polygon table (default `"x"`).

`polygon_y`

:   Name of the y-coordinate column in the polygon table (default `"y"`).

`polygon_id`

:   Name of the column that separates the polygons in the polygon table (default `"polygon"`). If the column does not exist, the function uses a `polygon` or `polygon_id` column. If neither exists, all rows form one polygon.

`polygon_val`

:   IDs of the polygons to draw. If `NULL`, the function draws all polygons.

`polygon_alpha`

:   Opacity of the polygon fill (default 0.2).

`polygon_fill`

:   Fill color of the polygons (default `"white"`).

`polygon_color`

:   Color of the polygon outlines (default `"white"`).

`polygon_linewidth`

:   Width of the polygon outlines (default 0.6).

`polygon_label`

:   Write a label at the center of each polygon (default `FALSE`)?

`polygon_label_by`

:   Name of the column that holds the labels. If `NULL`, the function uses the polygon ID.

`polygon_label_size`

:   Text size of the polygon labels (default 3).

`polygon_label_color`

:   Text color of the polygon labels (default `"white"`).

`nproc`

:   Number of processes for building the plot pairs (default 1).

`slot`

:   Name of the receptor slot to use. Leave it `NULL` if the SIRUP object has one receptor table. You must give it if the object has several slots. The function ignores `slot` if you give a plain receptor table.

`...`

:   Other arguments for [`plotSpatialColorTips()`](plotSpatialColorTips.md).

## Details

A tip with no spatial ID, or with a spatial ID that is not in `object`, shows in the tree plot only. The map shows only the tips whose spatial ID matches a spot in `object`. In the tree, by default tips with a spot are circles, and tips without a spot are squares.

If the object has multiple images or samples, the function draws one paired plot for each image or sample.

To match the tips to receptors, see [`formatSpatialClones()`](formatSpatialClones.md).

## Value

A named list with one element for each clone. Each element is a list with the items `tree` and `spatial_clone` (the two plots), and `legend` and `legend_rel_height` (for the optional legend). Print an element to see the tree and the map side by side, with the legend beneath them. Use the `tree` and `spatial_clone` items to make your own layout.

## See Also

[`spatialCloneTreePairsToPDF()`](spatialCloneTreePairsToPDF.md), [`clusterTreesTips()`](clusterTreesTips.md), [`getSpatialTrees()`](getSpatialTrees.md), [`plotSpatialColorTips()`](plotSpatialColorTips.md), [`plotSpatialCloneColorTips()`](plotSpatialCloneColorTips.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- formatSpatialClones(
  x,
  id = "sequence_id",
  germ = "germline_alignment_d_mask",
  traits = "spatial_id"
)
x <- getSpatialTrees(x)
trees <- getAnalysis(x, "trees")

pairs <- plotSpatialCloneTreePairs(trees, x, clone_val = "4")
names(pairs)

# Show the tree and the map side by side.
pairs[["4"]]
```
