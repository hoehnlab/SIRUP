# spatialMoran

## Calculate global Moran\'s I

## Description

`spatialMoran` tests one variable for spatial autocorrelation across spots using Moran\'s I.

## Usage

``` R
spatialMoran(
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
  k = 6,
  alternative = "greater",
  slot = NULL
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

:   Distance that defines a neighborhood. If `NULL`, the function uses `k`.

`k`

:   Number of nearest neighbors, used if `radius` is `NULL` (default 6).

`alternative`

:   Alternative hypothesis for `spdep::moran.test()` (default `"greater"`).

`slot`

:   Name of the receptor slot to use. Leave it `NULL` if the SIRUP object has one receptor table. You must give it if the object has several slots and the function reads receptors. The function ignores `slot` if you give a plain receptor table.

## Details

The function builds a neighbor graph from the spot coordinates. It uses all spots within `radius`. If `radius` is `NULL`, it uses the `k` nearest spots (or fewer, if the image has fewer spots). Then it calls `spdep::moran.test()` on the variable values. The function drops spots with non-finite values.

If the data come from more than one image, or a sirup object with multiple samples, then the function builds one graph and runs one test for each sample or image.

Set `variable_type = "clone"` to test where one clone occurs. The function makes a 0/1 value for each spot: 1 if sthe spot contains the clone, 0 if not.

## Value

For one image, a list with these items: `test` (the `spdep::moran.test()` result), `data` (the values tested), `neighbors` and `weights` (from `spdep`), `radius`, `k`, and `variable`. For several images, the list has `tests`, `neighbors`, and `weights` as named lists with one entry for each image. It also has `summary`, a data frame with one row for each image (`image`, `statistic`, `p_value`, `moran_i`, `expectation`, `variance`), and the items `data`, `radius`, `k`, and `variable`.

## See Also

[`spatialLocalMoran()`](spatialLocalMoran.md), [`spatialLee()`](spatialLee.md), [`spatialCloneJoinCount()`](spatialCloneJoinCount.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Test gene1 for spatial autocorrelation.
result <- spatialMoran(x, "gene1", k = 4)
result$test
```
