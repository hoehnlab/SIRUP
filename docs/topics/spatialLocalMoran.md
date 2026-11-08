# spatialLocalMoran

## Calculate local Moran\'s I (LISA)

## Description

`spatialLocalMoran` calculates local Moran\'s I for each spot. It shows where high or low values cluster.

## Usage

``` R
spatialLocalMoran(
  object,
  variable,
  variable_type = c("feature", "metadata", "receptor", "clone"),
  value = NULL,
  receptors = NULL,
  spatial_id = NULL,
  clone = NULL,
  assay = NULL,
  layer = "counts",
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  radius = NULL,
  k = 6
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`variable`

:   Name of the feature or the metadata column to test.

`variable_type`

:   Type of `variable`: `feature` (gene expression), `metadata` (Seurat metadata column), `receptor` (receptor table column), or `clone` (clone ID).

`value`

:   Metadata value to convert to a 0/1 indicator (default `NULL`).

`receptors`

:   AIRR-style receptor table. You must give it if `variable_type` is `"receptor"` or `"clone"` and `object` is not a SIRUP object.

`spatial_id`

:   Name of the `receptors` column that holds the spatial IDs. If `object` is a SIRUP object and this argument is `NULL`, the function uses the stored spatial ID field. Give a value if you use separate Seurat and receptor objects.

`clone`

:   Name of the `receptors` column that holds the clone IDs. The function uses it only if `variable_type` is `"clone"`. If `object` is a SIRUP object and this argument is `NULL`, the function uses the stored clone field, if there is one.

`assay`

:   Assay to use if `variable_type` is `"feature"`.

`layer`

:   Seurat expression layer to use for features (default `"counts"`).

`images`

:   Names of the Seurat images to use.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`radius`

:   Distance that defines a neighborhood. If `NULL`, the function uses the `k` nearest neighbors.

`k`

:   Number of nearest neighbors, used if `radius` is `NULL` (default 6).

## Details

The function builds one neighbor graph from all spots, with `radius` or `k` as in [`spatialMoran()`](spatialMoran.md). Set `variable_type = "clone"` to use a 0/1 value for each spot: 1 if the spot contains the clone, 0 if not.

If the data come from more than one image, or a sirup object with multiple samples, then the function builds one graph and runs one test for each sample or image.

## Value

A data frame with the coordinates, the values, local Moran\'s I, the z-scores, and the p-values.

## See Also

[`plotSpatialLocalMoran()`](plotSpatialLocalMoran.md), [`spatialMoran()`](spatialMoran.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
local_moran <- spatialLocalMoran(x, "gene1", k = 4)
head(local_moran[, c("cell", "local_i", "z_i", "p_value")])
```
