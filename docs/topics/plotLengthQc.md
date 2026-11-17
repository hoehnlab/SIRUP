# plotLengthQc

## Plot length histograms of V, J, and junction regions

## Description

`plotLengthQc` plots one histogram each for the V length, the J length, and the junction length of the receptors. Used to choose length thresholds for [`filterAirr()`](filterAirr.md).

## Usage

``` R
plotLengthQc(airr_df, bin_width = 5, title = NULL, slot = NULL)
```

## Arguments

`airr_df`

:   AIRR-style receptor table, or a SIRUP object. The table needs the columns `v_sequence_start`, `v_sequence_end`, `j_sequence_start`, `j_sequence_end`, and `junction_length`.

`bin_width`

:   Width of the histogram bins (default 5).

`title`

:   Plot title.

`slot`

:   Name of the receptor slot to use. Leave it `NULL` if the SIRUP object has one receptor table. You must give it if the object has several slots. The function ignores `slot` if you give a plain receptor table.

## Value

A ggplot object with one panel for each length.

## See Also

[`filterAirr()`](filterAirr.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
plotLengthQc(x)
```
