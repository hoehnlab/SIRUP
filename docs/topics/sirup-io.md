# sirup-io

## Save and read SIRUP objects

## Description

`saveSIRUP` writes a SIRUP object to one RDS file. `readSIRUP` reads the file back.

## Usage

``` R
saveSIRUP(x, file, ...)

readSIRUP(file)
```

## Arguments

`x`

:   A SIRUP object from [`createSIRUP()`](createSIRUP.md).

`file`

:   Path of the RDS file. A `.sirup.rds` extension is a good choice, but any path works.

`...`

:   Other arguments for [`saveRDS()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html), for example `compress`.

## Details

`saveSIRUP` validates the object. It stores the versions of SIRUP and SeuratObject that write the file.

`readSIRUP` does these steps:

1.  It converts an object that was saved under the old package name, SIRAP (class `sirap`), to the current names.

2.  It upgrades an older object to the current schema.

3.  It gives a warning if the versions of SIRUP or SeuratObject that wrote the file differ from the installed versions.

4.  It validates the object, then returns it.

## Value

`saveSIRUP` returns `file`, invisibly. `readSIRUP` returns the validated SIRUP object.

## See Also

[`createSIRUP()`](createSIRUP.md), [`validateSIRUP()`](validateSIRUP.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))

# Save the object to a temporary file, then read it back.
path <- saveSIRUP(x, tempfile(fileext = ".sirup.rds"))
x2 <- readSIRUP(path)
x2
```
