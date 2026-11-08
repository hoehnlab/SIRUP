# clusterTreesTips

## Cluster the tips of all trees in a dowser result

## Description

`clusterTreesTips` runs [`clusterTreeTips()`](clusterTreeTips.md) on each tree of a dowser result. It adds the cluster column to the clone data of each tree.

## Usage

``` R
clusterTreesTips(
  trees,
  k = "auto",
  id = "sequence_id",
  cluster = "tree_cluster",
  prefix = "tree_cluster_",
  ignore_tips = "Germline",
  keep = NULL,
  keep_val = NULL,
  method = c("split_largest_subtrees", "cut_largest_branches", "monophyletic_clades"),
  cut_terminal_edges = FALSE,
  k_range = 1:8,
  improvement_threshold = 0.1
)
```

## Arguments

`trees`

:   Result of [`dowser::getTrees()`](https://dowser.readthedocs.io/en/latest/topics/getTrees/), or a list with an element `trees` that holds this result.

`k`

:   Number of subtrees to make in each tree. Use `"auto"` (default) to choose `k` for each tree with [`chooseTreeTipK()`](chooseTreeTipK.md).

`id`

:   Name of the column in each clone data frame that matches the tree tip labels (default `"sequence_id"`).

`cluster`

:   Name of the cluster column to add (default `"tree_cluster"`).

`prefix`

:   Prefix of the cluster labels (default `"tree_cluster_"`).

`ignore_tips`

:   Tip labels to remove before the clustering (default `"Germline"`).

`keep`

:   Name of a column in each clone data frame (default `NULL`). The function runs the clustering first. Then, for each tip where `keep` differs from `keep_val`, it replaces the cluster label with the value of `keep`.

`keep_val`

:   Value of `keep` that marks a tip to keep.

`method`

:   Clustering method for [`clusterTreeTips()`](clusterTreeTips.md). The default is `"split_largest_subtrees"`.

`cut_terminal_edges`

:   Allow cuts of terminal branches (default `FALSE`)? With `FALSE`, the function prefers internal splits of clades. With `TRUE`, branch cuts can remove single terminal tips.

`k_range`

:   Values of `k` to test, if `k = "auto"` (default `1:8`). The function passes it to [`chooseTreeTipK()`](chooseTreeTipK.md).

`improvement_threshold`

:   Smallest relative improvement that justifies a larger `k`, if `k = "auto"` (default 0.10). The function passes it to [`chooseTreeTipK()`](chooseTreeTipK.md).

## Value

`trees`, with the `cluster` column added to the clone data of each tree.

## See Also

[`clusterTreeTips()`](clusterTreeTips.md), [`chooseTreeTipK()`](chooseTreeTipK.md), [`plotSpatialCloneTreePairs()`](plotSpatialCloneTreePairs.md)

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

# Use 3 clusters for every tree.
clustered <- clusterTreesTips(trees, k = 3)
table(clustered$data[[1]]@data$tree_cluster)

# Choose k for each tree.
clustered <- clusterTreesTips(trees, k = "auto")
```
