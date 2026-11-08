# readSpatialIdMapping

## Read a spatial ID mapping file

## Description

`readSpatialIdMapping` reads a table that maps spatial IDs to other values, such as tissue positions. The file can be CSV, TSV, or parquet.

## Usage

``` R
readSpatialIdMapping(path)
```

## Arguments

`path`

:   Path of the mapping file.

## Details

The function chooses the format from the file extension. A parquet file needs the `arrow` package.

## Value

A data frame.

## See Also

[`addTissuePositions()`](addTissuePositions.md)

## Examples

``` R
positions <- data.frame(
  barcode = c("s1", "s2"),
  pxl_col_in_fullres = c(100, 200),
  pxl_row_in_fullres = c(50, 60)
)
path <- tempfile(fileext = ".csv")
write.csv(positions, path, row.names = FALSE)
readSpatialIdMapping(path)
```
