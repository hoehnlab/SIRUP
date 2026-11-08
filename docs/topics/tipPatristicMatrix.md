# tipPatristicMatrix

## Calculate patristic distances between tips in a tree

## Description

`tipPatristicMatrix` returns the matrix of patristic distances between the selected tips of a tree. The rows and columns follow the order of `tip_labels`.

## Usage

``` R
tipPatristicMatrix(tree, tip_labels)
```

## Arguments

`tree`

:   Phylogenetic tree with the fields `edge`, `edge.length`, and `tip.label`.

`tip_labels`

:   Character vector of the tip labels to use.

## Value

A numeric matrix of patristic distances.

## See Also

[`clusterTreeTips()`](clusterTreeTips.md)

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
tree <- trees$trees[[1]]

tipPatristicMatrix(tree, tree$tip.label[1:4])
```
