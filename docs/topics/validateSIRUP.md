# validateSIRUP

## Validate a SIRUP object

## Description

`validateSIRUP` checks that the parts of a SIRUP object agree with each other.

## Usage

``` R
validateSIRUP(x, verbose = FALSE)
```

## Arguments

`x`

:   Object to validate.

`verbose`

:   Print a short summary of the result (default `FALSE`)?

## Details

The function gives an error if one of these checks fails:

- The object does not have class `sirup`.

- No receptor spatial ID matches a cell or spot name in the Seurat object. The function checks every receptor slot.

- A column named in `metadata$receptor_fields` does not exist in the receptor table. The `spatial_id` column must be in every receptor slot. Each other named column, such as `clone_id`, must be in at least one slot.

- The `sample` field is set, but the column is missing from a receptor table or from the Seurat metadata.

An object with several receptor slots can have a slot with no rows. For example, one section can have BCRs and no TCRs. The function skips the spatial ID checks for that slot.

## Value

`x`, invisibly, if the object is valid. If not, the function gives an error.

## See Also

[`createSIRUP()`](createSIRUP.md), [`isSIRUP()`](isSIRUP.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
validateSIRUP(x)
validateSIRUP(x, verbose = TRUE)
```
