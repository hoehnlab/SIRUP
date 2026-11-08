# spatialRipleyCrossK

## Calculate Ripley\'s cross K

## Description

`spatialRipleyCrossK` measures how spots of mark `i` cluster around spots of mark `j` at different distances. Similar to Lee\'s L but better for factor data.

## Usage

``` R
spatialRipleyCrossK(
  object,
  mark_col,
  i,
  j,
  nsim = 0,
  max_r = NULL,
  r = NULL,
  image_scale = NULL,
  swap_coords = NULL
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`mark_col`

:   Name of the metadata column that holds the point types.

`i`

:   First mark level.

`j`

:   Second mark level.

`nsim`

:   Number of random-label simulations for the envelope (default 0). Use `0` for the observed K only.

`max_r`

:   Largest radius to calculate. The function ignores it if you give `r`.

`r`

:   Vector of distances for `spatstat.explore::Kcross()`.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

## Details

The function calls `spatstat.explore::Kcross()` with the border correction. If `nsim` is greater than 0, it also builds an envelope from random label simulations.

If the data come from more than one image, or a sirup object with multiple samples, then the function builds one graph and runs one test for each sample or image.

## Value

A list with these items: `observed` (the observed K), `envelope` (the simulation envelope, if any), `plot_data`, `points` (the point data), `i`, `j`, `correction`, `nsim`, `max_r`, and `r`.

## See Also

[`plotRipleyCrossK()`](plotRipleyCrossK.md)

## Examples

``` R

x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Compare spots of type A with spots of type B. Use 19 random-label
# simulations for the envelope.
set.seed(1)
kcross <- spatialRipleyCrossK(
  x,
  mark_col = "cell_type",
  i = "A",
  j = "B",
  nsim = 19
)
head(kcross$plot_data)
```
