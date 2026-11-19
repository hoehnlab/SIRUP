# sirup-analyses

## Store and retrieve analyses in a SIRUP object

## Description

`getAnalysis` returns one stored result. `setAnalysis` adds or replaces a result. `listAnalyses` returns the names of the stored results. `removeAnalysis` deletes a result.

## Usage

``` R
getAnalysis(x, name, sample_id = NULL)

setAnalysis(x, name, value, sample_id = NULL, validate = TRUE)

listAnalyses(x)

removeAnalysis(x, name, validate = TRUE)
```

## Arguments

`x`

:   A SIRUP object.

`name`

:   Name of the analysis. Give one non-empty character string.

`sample_id`

:   Name of a spatial sample. Use it to store or get a result for one sample.

`value`

:   R object to store. It can be `NULL`.

`validate`

:   Validate the updated SIRUP object?

## Details

The `analyses` list of a SIRUP object holds any R object. Use it to keep derived results with the data that produced them. Examples are Dowser clones and trees, spatial statistics, and polygon summaries.

## Value

`getAnalysis` returns the stored object. `listAnalyses` returns a character vector. `setAnalysis` and `removeAnalysis` return the updated SIRUP object.

## See Also

[`updateAnalysis()`](sirup-transform.md), [`mergeSIRUP()`](mergeSIRUP.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- setAnalysis(x, "n_clones", 5)
listAnalyses(x)
getAnalysis(x, "n_clones")
x <- removeAnalysis(x, "n_clones")
listAnalyses(x)
```
