# plotRipleyCrossK

## Plot Ripley\'s cross K against distance

## Description

`plotRipleyCrossK` plots the observed cross K, with the simulation envelope if there is one.

## Usage

``` R
plotRipleyCrossK(
  kcross,
  title = NULL,
  observed_color = "black",
  simulated_color = "grey30",
  envelope_fill = "grey75",
  envelope_alpha = 0.45,
  facet = TRUE,
  log_y = FALSE
)
```

## Arguments

`kcross`

:   Result of [`spatialRipleyCrossK()`](spatialRipleyCrossK.md).

`title`

:   Plot title.

`observed_color`

:   Color of the observed curve.

`simulated_color`

:   Color of the dotted simulated curve.

`envelope_fill`

:   Fill color of the simulation envelope.

`envelope_alpha`

:   Opacity of the simulation envelope.

`facet`

:   Make one panel for each image if the result has several images (default `TRUE`)?

`log_y`

:   Use a log10 scale for the y-axis (default `FALSE`)?

## Value

A ggplot object.

## See Also

[`spatialRipleyCrossK()`](spatialRipleyCrossK.md)

## Examples

``` R

x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
set.seed(1)
kcross <- spatialRipleyCrossK(
  x,
  mark_col = "cell_type",
  i = "A",
  j = "B",
  nsim = 19
)
plotRipleyCrossK(kcross)
```
