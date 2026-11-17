# sirup-transform

## Update data stored in a SIRUP object

## Description

`updateReceptors`, `updateSpatial`, and `updateAnalysis` apply a function to one part of a SIRUP object. They return the updated SIRUP object.

## Usage

``` R
updateReceptors(
  x,
  fun,
  ...,
  sample_id = NULL,
  slot = NULL,
  receptor_fields = NULL,
  validate = TRUE
)

updateSpatial(x, fun, ..., sample_id = NULL, validate = TRUE)

updateAnalysis(x, name, fun, ..., sample_id = NULL, validate = TRUE)
```

## Arguments

`x`

:   A SIRUP object.

`fun`

:   Function, or the name of a function. The function receives the selected part of the SIRUP object as its first argument.

`...`

:   Other arguments for `fun`.

`sample_id`

:   Name of the spatial sample to update. Use it only for multi-sample objects.

`slot`

:   Name of the receptor slot that `updateReceptors` changes. Leave it `NULL` if the object has one receptor table. You must give it if the object has several slots.

`receptor_fields`

:   Named list of receptor field mappings to change. The names can be `spatial_id`, `image`, `sample`, `clone`, `locus`, or `count`. A `NULL` value clears an optional mapping.

`validate`

:   Validate the updated SIRUP object?

`name`

:   Name of the stored analysis. Give one non-empty character string.

## Details

`updateReceptors` needs a function that returns a data frame. The function stores the result as a tibble and keeps the existing field mappings. If the object has no clone field and the result has a `clone_id` column, the function records `clone_id` as the clone field. Use `receptor_fields` to set or clear a mapping yourself. If you give `sample_id`, the function applies `fun` to the rows of that sample only. It then puts the result back in the full table. If the object has several receptor slots, `slot` chooses the table that `fun` changes. The other slots stay as they are. The field mappings are shared by all slots. The `spatial_id` and `sample` columns must stay in every slot. The other mapped columns, such as `clone_id`, must stay in at least one slot.

`updateSpatial` needs a function that returns a Seurat object. In a multi-sample object, it applies `fun` to every Seurat object. If you give `sample_id`, it updates only that object.

`updateAnalysis` applies `fun` to the stored result called `name`. To store a new result, use [`setAnalysis()`](sirup-analyses.md).

## Value

The updated SIRUP object.

## See Also

[`getReceptors()`](sirup-accessors.md), [`getSpatial()`](sirup-accessors.md), [`setAnalysis()`](sirup-analyses.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Remove the sequences in clone 3.
x2 <- updateReceptors(x, dplyr::filter, clone_id != "3")
nrow(getReceptors(x2))

# Add a Seurat metadata column.
x3 <- updateSpatial(x, function(spatial) {
  spatial$gene_rich <- spatial$nFeature_Spatial >= 8
  spatial
})
table(getSpatial(x3)$gene_rich)

# Change a stored result.
x4 <- setAnalysis(x, "n_clones", 5)
x4 <- updateAnalysis(x4, "n_clones", function(value) value + 1)
getAnalysis(x4, "n_clones")
```
