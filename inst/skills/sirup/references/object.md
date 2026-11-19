# The SIRUP object

A SIRUP object keeps the spatial data, the receptor table, the tissue regions, and your results in one place. SIRUP checks that these parts agree each time you change the object.

## Get a small test object

Use this object to test code before you run it on the user's data. It has 25 spots, 10 genes, 5 clones, and 2 polygons. The spots and gene counts come from a simulation. The receptor rows are real BCR sequences from a human thymus section.

```r
library(SIRUP)
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x
summary(x)
```

## The five parts

| Part | Holds | Read it with |
|---|---|---|
| `spatial` | The Seurat object, with expression and the H&E image | `getSpatial(x)` |
| `receptors` | The AIRR receptor table, as a tibble | `getReceptors(x)` |
| `metadata` | Column names and image settings | `x$metadata` |
| `polygons` | Tissue regions, one row for each polygon, with a list column `vertices` | `getPolygons(x)` |
| `analyses` | Any R object that you store under a name | `getAnalysis(x, name)` |

The accessors return ordinary R objects. Pass them to any function from another package. In a multi-sample object, `getSpatial(x)` returns a named list of Seurat objects. Give `sample_id = "name"` to get one sample.

## Field names

SIRUP stores which receptor column plays which role in `x$metadata$receptor_fields`. There are six fields.

| Field | The receptor column that holds... |
|---|---|
| `spatial_id` | the name of the spot that holds the receptor |
| `image` | the name of the Seurat image, when the object has several images |
| `sample` | the sample name, in a multi-sample object |
| `clone` | the clone ID |
| `locus` | the receptor locus, such as `IGH` or `TRB` |
| `count` | the number of reads or UMIs of the sequence |

You set the fields in `createSIRUP()`. If the clone column has another name than `clone_id`, give `clone = "your_column"`. To change a field later, give `receptor_fields = list(clone = "your_column")` to `updateReceptors()`. A value of `NULL` clears an optional field.

`createSIRUP()` picks the spatial ID column from this list when you give none: `spatial_id`, `spatial_id_8um`, `spatial_id_2um`, `spatial_id_16um`, `barcode_8um`, `barcode_2um`, `cell`, `cell_barcode`. Give `spatial_id = "column"` when your column has another name.

## Rules that SIRUP checks

`validateSIRUP(x, verbose = TRUE)` gives an error when one of these rules fails:

1. The object has class `sirup`.
2. At least one receptor spatial ID matches a spot name in the Seurat object.
3. Each column that `receptor_fields` names exists in the receptor table.
4. If you set the `sample` field, the column exists in both the receptor table and the Seurat metadata.

`createSIRUP()`, `updateReceptors()`, `updateSpatial()`, `updateAnalysis()`, `setAnalysis()`, and `addPolygons()` run this check before they return. A change that breaks a rule fails and leaves your object as it was.

The error "No receptor IDs match the spatial Seurat object." means that no receptor spatial ID equals a name in `colnames(getSpatial(x))`. Compare a few values from each side. On 10x Visium, spot names end in `-1`. Add it with the `suffix` argument of `extractSpatialIds()`.

SIRUP keeps receptor rows that match no spot. Bulk and blood sequences are examples. The rules need only one match.

## Image settings

Two settings tell SIRUP how to put spots on the H&E image. Both live in `x$metadata`.

- `image_scale` is the Seurat image scale, `"lowres"` or `"hires"`. Use the image that the Seurat object stores. `Seurat::Load10X_Spatial()` stores the low-resolution image by default.
- `swap_coords` says whether to swap the x and y coordinates. Some older Seurat objects store the coordinates in the opposite order to the image. To test it, plot `SeuratObject::GetTissueCoordinates(spatial)` and compare the shape with the tissue image. Set `swap_coords = TRUE` when the spots sit rotated or mirrored over the image.

Plots and spatial statistics read these settings when their own `image_scale` and `swap_coords` arguments are `NULL`.

## Change the object

Three functions run any function on one part of the object and store the result.

- `updateReceptors(x, fun, ...)` calls `fun(getReceptors(x), ...)` and stores the table that `fun` returns.
- `updateSpatial(x, fun, ...)` calls `fun(getSpatial(x), ...)` and stores the Seurat object that `fun` returns.
- `updateAnalysis(x, name, fun, ...)` does the same for a stored analysis.

```r
# Remove the sequences of one clone.
x2 <- updateReceptors(x, dplyr::filter, clone_id != "3")

# Add a column to the Seurat metadata.
x3 <- updateSpatial(x, function(spatial) {
  spatial$gene_rich <- spatial$nFeature_Spatial >= 8
  spatial
})

# Add a column to the receptor table.
x4 <- updateReceptors(x, function(receptors) {
  dplyr::mutate(receptors, junction_nchar = nchar(junction))
})

# Call clones again. scoper stops when the table already has a clone_id column, so remove it first.
x5 <- updateReceptors(x, function(receptors) {
  receptors$clone_id <- NULL
  scoper::hierarchicalClones(receptors, threshold = 0.1, summarize_clones = FALSE, nproc = 1)
})
```

The same pattern runs germline reconstruction (`dowser::createGermlines()`) and mutation counting (`shazam::observedMutations()`). The function must return a data frame (for `updateReceptors`) or a Seurat object (for `updateSpatial`). If it returns fewer rows, the object keeps only those rows.

## Regions

A region is a polygon on the tissue image. `addPolygons(x, polygons = data_frame)` or `addPolygons(x, polygons_path = "file.csv")` stores the polygons. Then SIRUP marks the polygon of each spot and each receptor. It adds a `polygon_id` column and one logical `polygon_<ID>` column for each polygon, to both the Seurat metadata and the receptor table. A spot on the edge of a polygon counts as a member when `include_boundary = TRUE` (the default).

The polygon CSV file needs the columns `x`, `y`, and `polygon` (or `polygon_id`). Each row is one vertex, in order. `selectSpatialPolygons()` writes this file from an interactive app. An agent cannot use the app. Ask the user to draw the polygons.

## Several samples

`createSIRUP()` accepts one Seurat object or a named list of Seurat objects. With a list, the `sample` column of the receptor table names the list entry of each row. Rows with no matching entry stay in the table as non-spatial records. Bulk and blood sequences use this feature. Run `spatialSamples(x)` to list the samples. `mergeSIRUP()` joins several objects, and `subsetSIRUP()` keeps chosen samples, spots, or receptor rows.

## Keep and save results

```r
x <- setAnalysis(x, "n_clones", 5)     # store any R object under a name
listAnalyses(x)
getAnalysis(x, "n_clones")
x <- removeAnalysis(x, "n_clones")

path <- saveSIRUP(x, tempfile(fileext = ".sirup.rds"))   # write one RDS file
x_again <- readSIRUP(path)             # read it back and validate it
```

`saveSIRUP()` records the versions of SIRUP and SeuratObject. `readSIRUP()` does four things. It converts an object saved under the old package name (SIRAP). It upgrades an older object to the current layout. It warns you if the package versions differ. It validates the object.
