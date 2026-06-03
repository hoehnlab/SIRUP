This vignette shows how to call clones, visualize BCR/TCR clone
locations over tissue, build BCR lineage trees with `dowser`, and pair
each tree with the spatial locations of its colored tips.

For now, this vignette uses the full local thymus dataset from the
walkthrough. We should replace this with a small packaged toy dataset
later.

This read-the-docs version can be rendered in two modes. By default,
code is shown but not evaluated. To generate plots from the full
dataset, render on the cluster with `SIRAP_RENDER_PLOTS=true`.

## Setup

## Call BCR and TCR Clones

SIRAP does not wrap `scoper` clone calling. Use `scoper` directly, then
pass the resulting data frames to SIRAP plotting functions.

``` r
bcrClones <- scoper::hierarchicalClones(
  filter(bcrTopBarcode, productive),
  threshold = 0.1,
  cell_id = "sequence_id",
  summarize_clones = FALSE,
  nproc = nproc
)
#> Warning in prepare_db(db = db, junction = junction, v_call = v_call, j_call =
#> j_call, : Removed 189 sequences with non ATCG charachters.
#> Running defineClonesScoper in single cell mode
#> Warning in defineClonesScoper(db = db, threshold = threshold, model =
#> "hierarchical", : Single cell mode requested, but `db` doesn't contain light
#> chain data. Skipping.

tcrClones <- scoper::identicalClones(
  filter(tcrTopBarcode, productive),
  cell_id = "sequence_id",
  summarize_clones = FALSE,
  nproc = nproc
)
#> Warning in prepare_db(db = db, junction = junction, v_call = v_call, j_call =
#> j_call, : Removed 136 sequences with non ATCG charachters.
#> Running defineClonesScoper in single cell mode
#> Warning in defineClonesScoper(db = db, method = match.arg(method), model =
#> "identical", : Single cell mode requested, but `db` doesn't contain light chain
#> data. Skipping.

saveRDS(bcrClones, bcrClonesPath)
saveRDS(tcrClones, tcrClonesPath)
```

If clone calling has already been run, load the saved RDS files instead.

``` r
bcrClones <- readRDS(bcrClonesPath)
tcrClones <- readRDS(tcrClonesPath)
```

## Choose Clones to Plot

For a first pass, choose a couple of clones with enough tips to be
interesting but not so many that the tree becomes hard to read.

``` r
bcrCloneCounts <- bcrClones %>%
  filter(!is.na(clone_id)) %>%
  count(clone_id, sort = TRUE)

tcrCloneCounts <- tcrClones %>%
  filter(!is.na(clone_id)) %>%
  count(clone_id, sort = TRUE)

bcrCloneIds <- bcrCloneCounts %>%
  filter(n >= 30, n <= 50) %>%
  slice_head(n = 2) %>%
  pull(clone_id)

tcrCloneIds <- tcrCloneCounts %>%
  filter(n >= 30, n <= 50) %>%
  slice_head(n = 2) %>%
  pull(clone_id)

bcrCloneCounts %>% filter(clone_id %in% bcrCloneIds)
#> # A tibble: 2 × 2
#>   clone_id     n
#>   <chr>    <int>
#> 1 1764        50
#> 2 364         50
tcrCloneCounts %>% filter(clone_id %in% tcrCloneIds)
#> # A tibble: 2 × 2
#>   clone_id     n
#>   <chr>    <int>
#> 1 3373        45
#> 2 5528        35
```

## Plot Clone Locations Over Tissue

`plotReceptorExpression()` can plot any receptor table, including a
subset of one clone.

``` r
bcrClonePlots <- lapply(
  bcrCloneIds,
  function(cloneId) {
    plotReceptorExpression(
      object = thymus,
      receptors = filter(bcrClones, .data$clone_id == .env$cloneId),
      barcode_col = "barcode_8um",
      receptor_color = "#2c7fb8",
      receptor_size = 1.7,
      title = paste("BCR clone", cloneId)
    )
  }
)

tcrClonePlots <- lapply(
  tcrCloneIds,
  function(cloneId) {
    plotReceptorExpression(
      object = thymus,
      receptors = filter(tcrClones, .data$clone_id == .env$cloneId),
      barcode_col = "barcode_8um",
      receptor_color = "#d95f02",
      receptor_size = 1.7,
      title = paste("TCR clone", cloneId)
    )
  }
)

cowplot::plot_grid(plotlist = bcrClonePlots, ncol = 2)
```

![](/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP/docs/readthedocs/markdown/spatial-clones-and-trees_files/figure-gfm/clone-location-plots-1.png)<!-- -->

``` r
cowplot::plot_grid(plotlist = tcrClonePlots, ncol = 2)
```

![](/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP/docs/readthedocs/markdown/spatial-clones-and-trees_files/figure-gfm/clone-location-plots-2.png)<!-- -->

## Build BCR Trees with Dowser

SIRAP keeps tree building in `dowser`, then focuses on coloring and
pairing trees with spatial plots.

``` r
germlineDir <- path.expand("~/software/germlines/imgt/human/vdj")
refs <- dowser::readIMGT(germlineDir)
#> [1] "Read in 1198 from 17 fasta files"

bcrTreeInput <- bcrClones %>%
  filter(clone_id %in% bcrCloneIds, !is.na(clone_id))

bcrTreeInput <- dowser::createGermlines(bcrTreeInput, refs, id = "sequence_id")

bcrFormattedClones <- dowser::formatClones(
  bcrTreeInput,
  id = "sequence_id",
  germ = "germline_alignment",
  traits = c("barcode_8um"),
  columns = c("x", "y", "barcode_8um"),
  filterstop = FALSE,
  collapse = TRUE,
  add_count = "consensus_count" %in% colnames(bcrTreeInput),
  num_fields = if ("consensus_count" %in% colnames(bcrTreeInput)) "consensus_count" else NULL
)
#> Warning in dowser::formatClones(bcrTreeInput, id = "sequence_id", germ =
#> "germline_alignment", : columns x y barcode_8um contain multiple values per
#> clone, setting to character, flattening with comma

bcrTrees <- dowser::getTrees(
  bcrFormattedClones,
  nproc = nproc
)

saveRDS(bcrTrees, bcrTreesPath)
```

If trees have already been built, load them from disk.

``` r
bcrTrees <- readRDS(bcrTreesPath)
```

## Pair BCR Trees with Spatial Tip Locations

`plotSpatialCloneTreePairs()` assigns tip colors, draws the BCR tree
with `ggtree`, and creates a matching spatial plot over the tissue. Tips
in the same spatial location share a color.

``` r
pairedTreeSpatialPlots <- plotSpatialCloneTreePairs(
  bcrTrees,
  object = thymus,
  barcode_col = "barcode_8um",
  clone_ids = bcrCloneIds,
  id_col = "sequence_id",
  zoom = FALSE
)

firstCloneId <- names(pairedTreeSpatialPlots)[1]
cowplot::plot_grid(
  pairedTreeSpatialPlots[[firstCloneId]]$tree,
  pairedTreeSpatialPlots[[firstCloneId]]$spatial_clone,
  ncol = 2,
  rel_widths = c(1.2, 1)
)
```

![](/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP/docs/readthedocs/markdown/spatial-clones-and-trees_files/figure-gfm/tree-spatial-pairs-1.png)<!-- -->

Export all paired tree/spatial plots to a single PDF.

``` r
spatialCloneTreePairsToPDF(
  pairedTreeSpatialPlots,
  file.path(plotDir, sprintf("%s_bcr_tree_spatial_pairs.pdf", sampleLower)),
  width = 12,
  height = 6
)
```
