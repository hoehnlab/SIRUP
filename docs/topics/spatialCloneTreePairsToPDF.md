# spatialCloneTreePairsToPDF

## Save the tree and map pairs to a PDF file

## Description

`spatialCloneTreePairsToPDF` writes one page for each clone to a PDF file. Each page shows the tree and the tissue map of the clone side by side.

## Usage

``` R
spatialCloneTreePairsToPDF(
  pairs,
  path,
  width = 12,
  height = 6,
  rel_widths = c(1.2, 1),
  nproc = 1L
)
```

## Arguments

`pairs`

:   Result of [`plotSpatialCloneTreePairs()`](plotSpatialCloneTreePairs.md).

`path`

:   Path of the PDF file.

`width`

:   Width of a page, in inches (default 12).

`height`

:   Height of a page, in inches (default 6).

`rel_widths`

:   Relative widths of the tree panel and the map panel (default `c(1.2, 1)`).

`nproc`

:   Number of processes for the pages (default 1). The function needs the `qpdf` package to join the pages that it makes in parallel. If `qpdf` is not installed, it makes the pages one after the other.

## Value

`path`, invisibly.

## See Also

[`plotSpatialCloneTreePairs()`](plotSpatialCloneTreePairs.md)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- formatSpatialClones(
  x,
  id = "sequence_id",
  germ = "germline_alignment_d_mask",
  traits = "spatial_id"
)
x <- getSpatialTrees(x)
trees <- getAnalysis(x, "trees")

pairs <- plotSpatialCloneTreePairs(trees, x, clone_val = c("1", "4"))
spatialCloneTreePairsToPDF(pairs, tempfile(fileext = ".pdf"))
```
