# chooseTreeTipK

## Choose the number of tip clusters for a tree

## Description

`chooseTreeTipK` picks a number of tip clusters for one tree.

## Usage

``` R
chooseTreeTipK(
  tree,
  k_range = 1:10,
  improvement_threshold = 0.1,
  ignore_tips = "Germline",
  method = c("cut_largest_branches", "split_largest_subtrees", "monophyletic_clades"),
  cut_terminal_edges = TRUE
)
```

## Arguments

`tree`

:   Phylogenetic tree with the fields `edge`, `edge.length`, and `tip.label`.

`k_range`

:   Values of `k` to test (default `1:10`).

`improvement_threshold`

:   Smallest relative improvement of the score that justifies a larger `k` (default 0.10).

`ignore_tips`

:   Tip labels to remove before the scoring (default `"Germline"`).

`method`

:   Clustering method for [`clusterTreeTips()`](clusterTreeTips.md). The default is `"cut_largest_branches"`.

`cut_terminal_edges`

:   Allow cuts of terminal branches (default `TRUE`)? The function passes the value to [`clusterTreeTips()`](clusterTreeTips.md).

## Details

For each value in `k_range`, the function runs [`clusterTreeTips()`](clusterTreeTips.md). It scores the result as the mean within cluster patristic distance. It finds the first `k` where the relative improvement of the score is below `improvement_threshold`, and picks the `k` before it. If the improvement never falls below the threshold, it picks the largest `k`.

## Value

A list with the chosen `k` and `scores`, a data frame with the columns `k`, `score`, and `relative_improvement`.

## See Also

[`clusterTreeTips()`](clusterTreeTips.md), [`clusterTreesTips()`](clusterTreesTips.md)

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

chooseTreeTipK(tree)
```
