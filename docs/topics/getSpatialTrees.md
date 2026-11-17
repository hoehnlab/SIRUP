# getSpatialTrees

## Build trees from the formatted clones of a SIRUP object

## Description

`getSpatialTrees` runs [`dowser::getTrees()`](https://dowser.readthedocs.io/en/latest/topics/getTrees/) on the stored formatted clones. It stores the trees as an analysis.

## Usage

``` R
getSpatialTrees(
  object,
  f_clones_analysis = NULL,
  analysis_name = NULL,
  trait = NULL,
  id = NULL,
  dir = NULL,
  modelfile = NULL,
  build = "pratchet",
  exec = NULL,
  igphyml = NULL,
  fixtrees = FALSE,
  nproc = 1,
  quiet = 0,
  rm_temp = TRUE,
  palette = NULL,
  seq = NULL,
  collapse = FALSE,
  check_divergence = TRUE,
  ...,
  slot = NULL,
  validate = TRUE
)
```

## Arguments

`object`

:   A SIRUP object.

`f_clones_analysis`

:   Name of the stored formatted clones. If `NULL`, the function uses `"f_clones"`, or `"f_clones_<slot>"` if the object has several receptor slots.

`analysis_name`

:   Name for the stored trees (default `"trees"`). If the object has several receptor slots, the default also ends with `⁠_⁠` and the slot name.

`trait`

:   Trait for the parsimony models. You need it if you use `igphyml`.

`id`

:   Unique name for this analysis. You need it if you use `igphyml` or `dnapars`.

`dir`

:   Directory for the temporary files.

`modelfile`

:   File that sets the parsimony model.

`build`

:   Program that builds the trees: `pratchet` (default), `pml`, `dnapars`, `dnaml`, `igphyml`, or `raxml`.

`exec`

:   Location of the phylogenetic program.

`igphyml`

:   Location of the `igphyml` program, for parsimony.

`fixtrees`

:   Use the tree topologies that you supply (default `FALSE`)?

`nproc`

:   Number of cores for the calculations (default 1).

`quiet`

:   Amount of text that the function prints to the console (default 0).

`rm_temp`

:   Delete the temporary files (default `TRUE`)?

`palette`

:   Deprecated.

`seq`

:   Name of the column that holds the sequences.

`collapse`

:   Collapse the internal nodes that have identical sequences (default `FALSE`)?

`check_divergence`

:   Run `checkDivergence` on the trees (default `TRUE`)?

`...`

:   Other arguments for the tree-building programs.

`slot`

:   Name of the receptor slot that the trees come from. The function uses it only to find the default analysis names. Leave it `NULL` if the object has one receptor table. You must give it if the object has several slots and you use the default names.

`validate`

:   Validate the updated SIRUP object (default `TRUE`)?

## Details

[`dowser::getTrees()`](https://dowser.readthedocs.io/en/latest/topics/getTrees/) takes a tibble of formatted clones and returns a tibble with a `trees` column. `getSpatialTrees` does not change how dowser builds the trees. It reads the formatted clones from the SIRUP object by name (`f_clones_analysis`) and stores the trees in the object as an analysis. The trees and the clone data of each tree then stay with the spatial data. [`plotSpatialCloneTreePairs()`](plotSpatialCloneTreePairs.md) and [`clusterTreesTips()`](clusterTreesTips.md) can use them without more input.

By default, the function builds maximum parsimony trees with `phangorn::pratchet()`. It estimates the branch lengths with `ape::acctran()`. If you use `build = "igphyml"`, the function predicts the values of the internal nodes with maximum parsimony. In this case, give `dir` (a directory for the temporary files) and `id` (a unique name for this analysis). Use a different `id` for each call that shares a `dir`, so that the calls do not overwrite each other\'s files. The function writes `modelfile` if you do not give it. The file has no constraints. The function deletes the temporary files unless you set `rm_temp = FALSE`.

The arguments from `trait` to `...` go to [`dowser::getTrees()`](https://dowser.readthedocs.io/en/latest/topics/getTrees/).

## Value

The updated SIRUP object, with the trees as an analysis.

## See Also

[`formatSpatialClones()`](formatSpatialClones.md), [`clusterTreesTips()`](clusterTreesTips.md), [`plotSpatialCloneTreePairs()`](plotSpatialCloneTreePairs.md)

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
getAnalysis(x, "trees")
```
