# spatialCloneJoinCount

## Test clone localization with join counts

## Description

`spatialCloneJoinCount` tests if each clone clusters in space.

## Usage

``` R
spatialCloneJoinCount(
  object,
  receptors = NULL,
  clone = NULL,
  spatial_id = NULL,
  clone_val = NULL,
  polygon = NULL,
  polygon_val = NULL,
  images = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  x = "x",
  y = "y",
  polygon_id = "polygon",
  include_boundary = TRUE,
  min_spatial_ids = 2,
  radius = NULL,
  k = 6,
  alternative = "greater",
  other_label = "Other",
  use_stored_polygons = NULL,
  return_data = FALSE
)
```

## Arguments

`object`

:   Seurat object or SIRUP object.

`receptors`

:   AIRR-style receptor table with spatial IDs and clone IDs. Omit if `object` is a SIRUP object.

`clone`

:   Name of the `receptors` column that holds the clone IDs. If `object` is a SIRUP object and this argument is `NULL`, the function uses the stored clone field, if there is one.

`spatial_id`

:   Name of the `receptors` column that holds the spatial IDs. If `object` is a SIRUP object and this argument is `NULL`, the function uses the stored spatial ID field. Give a value if you use separate Seurat and receptor objects.

`clone_val`

:   Clone IDs to test. If `NULL`, the function tests all clones that occur in at least `min_spatial_ids` spatial IDs.

`polygon`

:   Data frame of polygon vertices. The function tests only the spots inside these polygons (default `NULL`).

`polygon_val`

:   Polygon IDs to select from `polygon`.

`images`

:   Names of the Seurat images to use.

`image_scale`

:   Image scale for `SeuratObject::GetTissueCoordinates()`. If `NULL`, the function uses the value stored in the SIRUP object.

`swap_coords`

:   Swap the two Seurat coordinate axes? If `NULL`, the function uses the value stored in the SIRUP object. Give a value if `object` is a Seurat object.

`x`

:   Name of the x-coordinate column in `polygon` (default `"x"`).

`y`

:   Name of the y-coordinate column in `polygon` (default `"y"`).

`polygon_id`

:   Name of the column in `polygon` that groups the vertices by polygon (default `"polygon"`).

`include_boundary`

:   Count spots on a polygon edge as inside the polygon (default `TRUE`)? The function ignores this argument if `object` is a SIRUP object with stored polygons and `polygon` is `NULL`. In that case, it uses the stored polygon annotations.

`min_spatial_ids`

:   Smallest number of spatial IDs that a clone needs to be tested (default 2).

`radius`

:   Distance that defines a neighborhood. If `NULL`, the function uses `k`.

`k`

:   Number of nearest neighbors, used if `radius` is `NULL` (default 6).

`alternative`

:   Alternative hypothesis for `spdep::joincount.test()` (default `"greater"`).

`other_label`

:   Label for the spatial IDs that do not contain the tested clone (default `"Other"`).

`use_stored_polygons`

:   Run one test for each stored polygon of a SIRUP object? Use `TRUE` to require stored polygons. Use `FALSE` to build one neighbor graph for each sample and ignore stored polygons. If `NULL` (default), the function uses `TRUE` when a SIRUP object has stored polygons and you give no `polygon`. Otherwise, it uses `FALSE`.

`return_data`

:   Return the data used for each test (default `FALSE`)?

## Details

The function calls `spdep::joincount.test()` for each clone. The test asks if neighboring spots with receptors share the clone more often than random labels would give.

The function skips a clone if it occurs in fewer than `min_spatial_ids` spatial IDs, or in every spatial ID of the region.

If `object` is a SIRUP object with stored polygons and you give no `polygon`, the function runs one test for each stored polygon. Set `use_stored_polygons` to FALSE to change this.

If the data come from more than one image, or a sirup object with multiple samples, then the function builds one graph and runs one test for each sample or image.

## Value

A list with these items: `summary` (a table with one row for each tested clone and region, ordered by p-value), `tests` (the `spdep` test results), `data` (the data for each test, or `NULL` if `return_data` is `FALSE`), `neighbors`, `weights`, `radius`, `k`, `alternative`, and `min_spatial_ids`.

## See Also

[`plotSpatialCloneJoinCount()`](plotSpatialCloneJoinCount.md), [`spatialMoran()`](spatialMoran.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Test each clone on the whole tissue.
join_count <- spatialCloneJoinCount(
  x,
  use_stored_polygons = FALSE,
  k = 4,
  return_data = TRUE
)
join_count$summary[, c("clone_id", "observed_join_count", "expected_join_count", "p_value")]
```
