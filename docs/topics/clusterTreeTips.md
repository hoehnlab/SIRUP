# clusterTreeTips

## Cluster the tips of a tree

## Description

`clusterTreeTips` groups the tips of a phylogenetic tree into `k` connected clusters.

## Usage

``` R
clusterTreeTips(
  tree,
  k,
  id = "sequence_id",
  cluster = "tree_cluster",
  ignore_tips = "Germline",
  prefix = "tree_cluster_",
  method = c("split_largest_subtrees", "cut_largest_branches", "monophyletic_clades"),
  cut_terminal_edges = TRUE
)
```

## Arguments

`tree`

:   Phylogenetic tree with the fields `edge`, `edge.length`, and `tip.label`.

`k`

:   Number of subtrees to make, before the function removes the ignored tips. If `k` is larger than the number of tips, the function gives a warning and uses the number of tips.

`id`

:   Name of the tip ID column in the result (default `"sequence_id"`).

`cluster`

:   Name of the cluster column in the result (default `"tree_cluster"`).

`ignore_tips`

:   Tip labels to remove **before** the clustering (default `"Germline"`). These tips do not form a cluster, do not use up a cluster, and do not appear in the result. The function always removes a tip named `"Germline"` (in any case). It prunes the tree with `ape::drop.tip()`. If `ape` is not installed, the function removes the tips from the result table after the clustering.

`prefix`

:   Prefix of the cluster labels (default `"tree_cluster_"`).

`method`

:   Clustering method: `"split_largest_subtrees"` (default), `"cut_largest_branches"`, or `"monophyletic_clades"`. See Details.

`cut_terminal_edges`

:   Allow cuts of terminal branches (default `TRUE`)? If `FALSE`, the function cuts internal branches first, so that it splits larger clades before it cuts off single tips.

## Details

The function groups the tips of the tree into `k` connected clusters. This helps with visualization: tips in the same cluster get the same label, and you can give them the same color. Each cluster is a subtree. The function makes the `k` subtrees by cutting `k - 1` branches. The `method` argument sets how the function chooses the branches to cut:

- `split_largest_subtrees` (default): The function splits one subtree at each step. It picks the subtree with the largest genetic span, which is the largest patristic distance between two of its tips. If two subtrees tie, it picks the subtree with more tips. In that subtree, a candidate branch must leave two subtrees that each have at least one tip. The function cuts the branch that makes the smaller of the two subtrees as large as possible (by tip count). If branches tie, it cuts the longest branch. If `cut_terminal_edges` is `FALSE`, the function considers internal branches first. It considers terminal branches (the branches that end in a tip) only if no internal branch is valid. This stops the function from cutting off long terminal branches before it splits the larger clades. The function stops after `k - 1` cuts, or when no valid cut remains.

- `cut_largest_branches`: The function cuts the `k - 1` longest branches of the whole tree. If `cut_terminal_edges` is `FALSE`, it cuts the longest internal branches first. It cuts terminal branches only if the tree has fewer than `k - 1` internal branches. This is the heuristic for fixed-k minimum spanning tree clustering, applied to a phylogenetic tree.

- `monophyletic_clades`: Each cluster is one complete clade, which is a node and all of its descendants. No cluster is a paraphyletic remainder. Use this method if each cluster must match a real clade, for example if each tip color in a figure must mark one monophyletic group. The function starts with the whole tree as one clade. At each step it picks the clade with the largest genetic span (the clade needs at least two tips). It replaces that clade with all of its child clades, which are the subtrees that start at the children of its root node. The function stops when `k` clades exist, or when no clade can split. `cut_terminal_edges` has no effect. For a bifurcating tree, the result has exactly `k` clades. At a multifurcation, one step adds all child clades, so the result can have more than `k` clades. If an ignored tip (such as the germline) forms its own clade, the function removes that clade. The result can then have fewer than `k` clusters.

The function labels the clusters with `prefix` and a number, such as `tree_cluster_1`.

## Value

A data frame with one row for each kept tip. It has the tip ID column (named by `id`) and the cluster label column (named by `cluster`).

## See Also

[`chooseTreeTipK()`](chooseTreeTipK.md), [`clusterTreesTips()`](clusterTreesTips.md)

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

# Cut the tree into 3 clusters.
clusterTreeTips(tree, k = 3)

# Make each cluster one complete clade.
clusterTreeTips(tree, k = 3, method = "monophyletic_clades")
```
