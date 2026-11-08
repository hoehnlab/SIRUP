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

- No receptor spatial ID matches a cell or spot name in the Seurat object.

- A column named in `metadata$receptor_fields` does not exist in the receptor table.

- The `sample` field is set, but the column is missing from the receptor table or from the Seurat metadata.

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
