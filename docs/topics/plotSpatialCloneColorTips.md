# plotSpatialCloneColorTips

## Plot the tip colors of one clone on a tissue image

## Description

`plotSpatialCloneColorTips` plots the spots of one clone. Each spot gets the color of its tip. The function calls [`plotSpatialColorTips()`](plotSpatialColorTips.md).

## Usage

``` R
plotSpatialCloneColorTips(
  colored_clones,
  clone_val,
  tissue_positions = NULL,
  clone = "clone_id",
  spatial_id = "spatial_id_8um",
  image_scale = NULL,
  swap_coords = NULL,
  sample_id = NULL,
  ...
)
```

## Arguments

`colored_clones`

:   Result of [`assignTipColors()`](assignTipColors.md), or a SIRUP object whose receptor table has the tip colors.

`clone_val`

:   ID of the clone to plot.

`tissue_positions`

:   Seurat object or SIRUP object. Omit it if `colored_clones` is a SIRUP object.

`clone`

:   Name of the column that holds the clone IDs (default `"clone_id"`). If `colored_clones` is a SIRUP object and you omit this argument, the function uses the stored clone field.

`spatial_id`

:   Name of the column in `colored_clones` that holds the spatial IDs (default `"spatial_id_8um"`). If `colored_clones` is a SIRUP object and you omit this argument, the function uses the stored spatial ID field.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object.

`sample_id`

:   Names of the samples to plot. Use it for multi-sample SIRUP objects. If `NULL`, the function plots all samples.

`...`

:   Other arguments for [`plotSpatialColorTips()`](plotSpatialColorTips.md).

## Value

A ggplot object.

## See Also

[`plotSpatialColorTips()`](plotSpatialColorTips.md), [`assignTipColors()`](assignTipColors.md), [`plotSpatialCloneTreePairs()`](plotSpatialCloneTreePairs.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- updateReceptors(x, assignTipColors, shared_color = "spatial_id")
plotSpatialCloneColorTips(x, clone_val = "1")
```
