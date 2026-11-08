# isSIRUP

## Test whether an object is a SIRUP object

## Description

`isSIRUP` checks if an object has class `sirup`. It does not validate the object. Use [`validateSIRUP()`](validateSIRUP.md) for that.

## Usage

``` R
isSIRUP(x)
```

## Arguments

`x`

:   Object to test.

## Value

`TRUE` if `x` has class `sirup`. Otherwise, `FALSE`.

## See Also

[`validateSIRUP()`](validateSIRUP.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
isSIRUP(x)
isSIRUP(getReceptors(x))
```
