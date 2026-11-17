# createSIRUP

## Create a SIRUP object

## Description

`createSIRUP` joins a Seurat spatial object and an AIRR receptor table into one SIRUP object.

## Usage

``` R
createSIRUP(
  spatial,
  receptors,
  spatial_id = NULL,
  image = NULL,
  sample_id = NULL,
  clone = NULL,
  locus = "locus",
  count = NULL,
  image_scale = NULL,
  swap_coords = NULL,
  metadata = list(),
  polygons = NULL,
  polygons_path = NULL,
  include_boundary = TRUE,
  analyses = list(),
  validate = TRUE
)
```

## Arguments

`spatial`

:   Seurat object with spatial expression data, or a named list of Seurat objects (one per spatial sample).

`receptors`

:   AIRR-style receptor table (a data frame or similar), or a named list of such tables. Each table in a list is one slot, for example `list(bcr = bcr_table, tcr = tcr_table)`.

`spatial_id`

:   Name of the `receptors` column that holds the spatial cell or spot IDs. If `NULL`, the function uses the first column found from this list: `spatial_id`, `spatial_id_8um`, `spatial_id_2um`, `spatial_id_16um`, `barcode_8um`, `barcode_2um`, `cell`, `cell_barcode`. It stores the name in the metadata and gives an error if none exists.

`image`

:   Name of the `receptors` column that identifies the Seurat image (default `NULL`).

`sample_id`

:   Name of the `receptors` column that holds the sample name. If `spatial` is a named list, this argument is required. If a Seurat object has no such metadata column, the function fills it with the list name.

`clone`

:   Name of the `receptors` column that holds the clone ID (default `NULL`). The function gives an error if the column does not exist.

`locus`

:   Name of the `receptors` column that holds the receptor locus (default `"locus"`). If you use the default and the column does not exist, the function stores `NULL`. If you give another name that does not exist, the function gives an error.

`count`

:   Name of the `receptors` column that holds the observation count or weight (default `NULL`).

`image_scale`

:   Default Seurat image scale (default `NULL`). Plotting and spatial statistics functions use it when their own `image_scale` is `NULL`. For multi-sample objects, give a named vector or list with one value for each sample, or one unnamed value that applies to all samples.

`swap_coords`

:   Swap the two Seurat coordinate axes by default (default `NULL`)? Plotting and spatial statistics functions use this value when their own `swap_coords` is `NULL`. For multi-sample objects, give a named vector or list with one value for each sample, or one unnamed value that applies to all samples.

`metadata`

:   Named list of dataset-level metadata.

`polygons`

:   Polygon table (default `NULL`). The table has one row per polygon, with a `polygon_id` column and a list column `vertices`. Each `vertices` table needs columns `x` and `y`. It can also have `spatial_id` and `image`. The function adds a `polygon_id` column and one logical `⁠polygon_<ID>⁠` column per polygon to both the Seurat metadata and the receptor table.

`polygons_path`

:   Path to a flat CSV file of polygon vertices (default `NULL`). The file needs `x` and `y` columns and a `polygon_id` or `polygon` column. The function annotates these polygons in the same way as `polygons`.

`include_boundary`

:   Count spots on a polygon edge as members of the polygon (default `TRUE`)?

`analyses`

:   Named list of derived analysis results.

`validate`

:   Validate the new object (default `TRUE`)?

## Details

`spatial` can be one Seurat object or a named list of Seurat objects. If you give a list, the column of `receptors` that `sample_id` names must name the corresponding list entry for each row. The function keeps rows with no matching entry as nonspatial records. For example, these rows can hold bulk or single-cell data.

The function stores key receptor column names in `metadata$receptor_fields` and converts `receptors` to a [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html). The optional `analyses` list keeps derived results.

`receptors` can also be a named list of tables, for example one table of BCRs and one of TCRs. Each table is a *slot*. All slots share the same spatial data, polygons and field mappings. The spatial ID and sample columns must be in every table. Each other mapped column, such as `clone_id`, must be in at least one table. Clone IDs only need to be unique within one slot. Functions that use receptors analyse one slot at a time, chosen with their `slot` argument. A SIRUP object with one table has the slot name `"only_slot"`.

The examples in this package use a small synthetic SIRUP object. It has 25 spots, 10 genes, 5 clones, and 2 polygons. The spots and gene counts are simulated. The receptor rows are real BCR sequences from a human thymus section, assigned to synthetic spots. Load the object with `readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))`. The script `data-raw/make_synthetic_example.R` builds it.

## Value

A SIRUP object (S3 class `sirup`).

## See Also

[`validateSIRUP()`](validateSIRUP.md), [`mergeSIRUP()`](mergeSIRUP.md), [`saveSIRUP()`](sirup-io.md)

## Examples

``` R
# Use the small synthetic object that comes with the package.
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x

# Build a SIRUP object from a Seurat object and a receptor table.
sirup <- createSIRUP(
  spatial = getSpatial(x),
  receptors = getReceptors(x),
  spatial_id = "spatial_id",
  clone = "clone_id",
  image_scale = "hires",
  swap_coords = FALSE
)
sirup
```
