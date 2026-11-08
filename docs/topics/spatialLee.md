# spatialLee

## Calculate global Lee\'s L

## Description

`spatialLee` tests the spatial association between two continuous variables.

## Usage

``` R
spatialLee(
  object,
  x,
  y,
  x_type = c("feature", "metadata"),
  y_type = c("feature", "metadata"),
  assay = NULL,
  layer = "counts",
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  radius = NULL,
  k = 6,
  alternative = "greater",
  method = c("analytical", "mc"),
  nsim = 999,
  progress = FALSE
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`x`

:   Name of the first variable. It is a feature or a numeric metadata column.

`y`

:   Name of the second variable. It is a feature or a numeric metadata column.

`x_type`

:   Type of `x`: `feature` or `metadata`.

`y_type`

:   Type of `y`: `feature` or `metadata`.

`assay`

:   Assay to use if `x_type` or `y_type` is `"feature"`.

`layer`

:   Seurat expression layer to use for features (default `"counts"`).

`images`

:   Names of the Seurat images to use.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`radius`

:   Distance that defines a neighborhood. If `NULL`, the function uses `k`.

`k`

:   Number of nearest neighbors, used if `radius` is `NULL` (default 6). Use `k = 8` for Visium HD, where the bins form a square grid. Use `k = 6` for standard Visium, where the spots form a hexagonal grid.

`alternative`

:   Alternative hypothesis for the significance test (default `"greater"`).

`method`

:   Significance test: `analytical` or `mc` (permutation). See Details.

`nsim`

:   Number of permutations for `method = "mc"` (default 999).

`progress`

:   Show a progress bar for the permutations (default `FALSE`)? The function uses it only if `method = "mc"`. Each permutation re-lags both variables over every spot. For Visium HD, a large `nsim` takes minutes for each pair.

## Details

The function builds one neighbor graph from all spots, with `radius` or `k` as in [`spatialMoran()`](spatialMoran.md). For example, give two gene names to test if the expression of one gene relates to the expression of the second gene in neighboring spots.

`method = "analytical"` calls `spdep::lee.test()`. This function builds a dense n by n weights matrix. The memory use grows with the square of the number of spots, so be careful with larger datasets, like with Visium HD. For large datasets, use `method = "mc"`. It calls `spdep::lee.mc()`, which permutes the labels. It keeps the weights sparse and needs only a few GB.

If the data come from more than one image, or a sirup object with multiple samples, then the function builds one graph and runs one test for each sample or image.

## Value

A list with these items: `test` (the test result), `L` (the observed Lee\'s L), `method`, `data` (the values tested), `neighbors` and `weights` (from `spdep`), `radius`, `k`, `x`, `y`, `x_type`, and `y_type`.

## See Also

[`spatialMoran()`](spatialMoran.md), [`spatialLocalMoran()`](spatialLocalMoran.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Test if gene1 and gene2 are associated.
spatialLee(x, "gene1", "gene2", k = 4)$L

# Use permutations. This method also works for large Visium HD data.
set.seed(1)
spatialLee(x, "gene1", "gene2", k = 4, method = "mc", nsim = 99)$test
```
