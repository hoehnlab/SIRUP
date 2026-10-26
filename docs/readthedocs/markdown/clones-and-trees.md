Clones and lineage trees on the tissue
================

A clone is a group of B or T cells that descend from one ancestor cell
and share its receptor rearrangement. In germinal centers, B cell clones
gain somatic hypermutations. A lineage tree of a B cell clone shows how
its sequences relate to each other. Spatial AIRR-seq lets you put each
tip of that tree on the tissue.

SIRUP does not reimplement clone calling or tree building. It runs the
Immcantation tools (scoper, dowser, shazam) on the receptor table inside
the SIRUP object. It then adds the spatial part: clone maps, trees
paired with tissue plots, and distances between tips in and between
regions.

This vignette uses MG thymus section 2 and its eGC polygons. You need
the Immcantation packages to run it:

``` r
install.packages(c("scoper", "dowser", "shazam", "alakazam", "cowplot"))
BiocManager::install("ggtree")
```

## Set up

``` r
library(SIRUP)
library(dplyr)
library(ggplot2)
```

``` r
extdata <- system.file("extdata", package = "SIRUP")

thymus <- readRDS(file.path(extdata, "mg_thymus_slice2_spatial.rds"))
bcrs <- airr::read_rearrangement(file.path(extdata, "mg_thymus_slice2_bcr.tsv.gz")) |>
  extractSpatialIds(
    id = "sequence_id",
    spatial_id = "barcode_10x",
    regex = ",([ACGT]{16})$",
    suffix = "-1"
  ) |>
  mutate(sequence_id = gsub(",", "_", sequence_id)) |>
  filterAirr(
    locus_val = "IGH",
    min_v_length = 250,
    min_j_length = 30,
    min_junction_length = 10,
    productive_only = TRUE,
    min_consensus_count = NULL,
    spatial_id = "barcode_10x"
  ) |>
  filter(!grepl("[^ACGT]", junction))
#> filterAirr kept 2195 of 3801 rows (dropped 1606).

sirup <- createSIRUP(
  spatial = thymus,
  receptors = bcrs,
  spatial_id = "barcode_10x",
  image_scale = "lowres",
  swap_coords = FALSE,
  polygons_path = file.path(extdata, "mg_thymus_slice2_gc_polygons.csv")
)
```

The setup also replaces the comma in each `sequence_id` with `_`. dowser
writes trees in Newick format, where a comma separates tips, so it
renames tips that contain one. With the comma gone, the tip names match
the receptor table, and SIRUP can find each tip on the tissue.

## Call clones inside the object

The setup made a SIRUP object without clones. To call clones, run
`scoper::hierarchicalClones()` through `updateReceptors()`.
`updateReceptors()` gives the receptor table to the function as its
first argument, with your other arguments. It then stores the returned
table and validates the object:

``` r
sirup <- updateReceptors(
  sirup,
  scoper::hierarchicalClones,
  threshold = 0.1,
  summarize_clones = FALSE,
  nproc = 1
)
#> Warning: package 'scoper' was built under R version 4.5.3

sirup
#> <SIRUP object>
#> Spatial data:    Seurat, 3,203 spatial IDs, 449 genes, 1 image
#> Receptor data:   tibble, 2,105 receptor sequences x 67 columns
#> Spatial overlap: 797/868 receptor IDs in spatial data (91.8%)
#> Spatial ID field: barcode_10x
#> Clone ID field:  clone_id
#> Analyses:        0
#> Polygons:        9
dplyr::n_distinct(getReceptors(sirup)$clone_id)
#> [1] 953
```

`hierarchicalClones()` groups sequences with the same V gene, J gene,
and junction length whose junctions differ by less than `threshold`. The
SIRUP paper used `threshold = 0.1` for BCRs. TCRs do not mutate, so call
TCR clones with `scoper::identicalClones()`.

The function writes a `clone_id` column, and `updateReceptors()` records
it as the clone field. If a function writes a different column, record
it yourself with `receptor_fields = list(clone = "my_clone_column")`.

## Germlines

A lineage tree starts at the germline sequence of its clone.
`dowser::createGermlines()` reconstructs one germline per clone from the
IMGT reference sequences. The SIRUP paper made the germlines and removed
likely PCR chimeras in the object.

SIRUP does not include the IMGT references. Download them (for example
with the Immcantation script `fetch_imgtdb.sh`) and set `SIRUP_IMGT_DIR`
to the folder that holds the human V, D, and J FASTA files, such as
`~/share/germlines/imgt/human/vdj`. The germline and tree chunks below
run only when that folder exists.

``` r
references <- dowser::readIMGT(imgt_dir)
#> [1] "Read in 1198 from 17 fasta files"

sirup <- updateReceptors(
  sirup,
  dowser::createGermlines,
  references = references,
  id = "sequence_id",
  nproc = 1
)

sirup <- updateReceptors(sirup, function(receptors) {
  chimera <- shazam::slideWindowDb(
    receptors,
    sequenceColumn = "sequence_alignment",
    germlineColumn = "germline_alignment_d_mask",
    mutThresh = 6,
    windowSize = 10
  )
  receptors[!chimera, ]
})
```

## Somatic hypermutation

`shazam::observedMutations()` counts the mutations of each sequence
against a germline. With clonal germlines, use
`germline_alignment_d_mask`. Here we use `germline_alignment`, the
germline of each sequence that IgBLAST wrote, so the chunk runs without
the IMGT references. The result goes in `mu_freq`, the fraction of
mutated positions:

``` r
sirup <- updateReceptors(
  sirup,
  shazam::observedMutations,
  sequenceColumn = "sequence_alignment",
  germlineColumn = "germline_alignment",
  frequency = TRUE,
  combine = TRUE
)
summary(getReceptors(sirup)$mu_freq)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00000 0.01462 0.03750 0.04184 0.06377 0.13467
```

`plotReceptors()` can color spots by a numeric receptor column in place
of groups. Give the column to `value`. When a spot has several
receptors, `value_combine` sets how SIRUP combines them (the mean, by
default):

``` r
plotReceptors(
  sirup,
  group = "locus",
  value = "mu_freq",
  value_colors = c("white", "darkred"),
  receptor_size = 1.3,
  polygons = TRUE,
  polygon_alpha = 0,
  polygon_color = "black",
  title = "Mean SHM frequency per spot"
)
```

<img src="/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRUP/docs/readthedocs/markdown/clones-and-trees_files/figure-gfm/shm-map-1.jpeg" alt="" width="100%" />

## Choose clones

Trees are most useful for expanded clones. Here we choose clones in the
eGCs of this section that have between 10 and 60 sequences:

``` r
clone_sizes <- getReceptors(sirup) |>
  group_by(clone_id) |>
  summarise(
    n_sequences = n(),
    n_in_gc = sum(!is.na(polygon_id)),
    n_gcs = n_distinct(na.omit(polygon_id)),
    gcs = paste(sort(unique(na.omit(polygon_id))), collapse = ", ")
  ) |>
  filter(n_sequences >= 10, n_sequences <= 60, n_in_gc > 0) |>
  arrange(desc(n_in_gc))

head(clone_sizes)
#> # A tibble: 6 × 5
#>   clone_id n_sequences n_in_gc n_gcs gcs  
#>   <chr>          <int>   <int> <int> <chr>
#> 1 739               15       6     2 6, 7 
#> 2 810               13       6     1 5    
#> 3 196               40       5     2 3, 5 
#> 4 604               12       3     2 5, 6 
#> 5 428               11       2     1 5    
#> 6 598               13       1     1 7
tree_clones <- head(clone_sizes$clone_id, 3)
```

Show their spots on the tissue:

``` r
plotReceptors(
  sirup,
  group = "clone_id",
  group_val = tree_clones,
  receptor_size = 1.6,
  polygons = TRUE,
  polygon_alpha = 0,
  polygon_color = "black",
  zoom = TRUE,
  title = "Selected clones"
)
```

<img src="/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRUP/docs/readthedocs/markdown/clones-and-trees_files/figure-gfm/clone-map-1.jpeg" alt="" width="100%" />

## Build trees

Two functions run dowser on the object:

  - `formatSpatialClones()` runs `dowser::formatClones()` and stores the
    result as an analysis.
  - `getSpatialTrees()` runs `dowser::getTrees()` on a stored
    formatted-clone analysis and stores the trees as another analysis.

Both take the dowser arguments. Tree building takes time for large
clones, so first make a smaller object with only the clones you want.
`updateReceptors()` with `dplyr::filter()` does this:

``` r
tree_sirup <- updateReceptors(sirup, dplyr::filter, clone_id %in% tree_clones)

tree_sirup <- formatSpatialClones(
  tree_sirup,
  analysis_name = "formatted_clones",
  id = "sequence_id",
  germ = "germline_alignment_d_mask",
  traits = c("barcode_10x", "polygon_id")
)

tree_sirup <- getSpatialTrees(
  tree_sirup,
  f_clones_analysis = "formatted_clones",
  analysis_name = "trees",
  build = "pratchet"
)

trees <- getAnalysis(tree_sirup, "trees")
trees
#> # A tibble: 3 × 5
#>   clone_id data       locus  seqs trees  
#>   <chr>    <list>     <chr> <int> <list> 
#> 1 196      <airrClon> IGH      38 <phylo>
#> 2 739      <airrClon> IGH      15 <phylo>
#> 3 810      <airrClon> IGH      13 <phylo>
```

`traits` keeps the columns that the tree plots need. Keep the spatial
ID, so SIRUP can put each tip on the tissue, and any column that you
want to color by.

`build = "pratchet"` builds maximum parsimony trees, which is fast. For
the paper we used maximum likelihood (`build = "pml"`). For B cell
specific models, use IgPhyML (`build = "igphyml"` with the `igphyml`
path).

Store the trees in the main object:

``` r
sirup <- setAnalysis(sirup, "trees", trees)
```

## Pair each tree with the tissue

`plotSpatialCloneTreePairs()` makes two plots for each clone: the
lineage tree and a tissue map. Each tip gets one color in both plots, so
you can find a tip of the tree on the tissue:

``` r
pairs <- plotSpatialCloneTreePairs(
  trees,
  object = sirup,
  clone_val = tree_clones,
  id = "sequence_id",
  zoom = TRUE,
  polygons = TRUE,
  polygon_alpha = 0,
  polygon_color = "black"
)
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the ggtree package.
#>   Please report the issue at <https://github.com/YuLab-SMU/ggtree/issues>.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
names(pairs)
#> [1] "196" "739" "810"

cowplot::plot_grid(
  pairs[[1]]$tree,
  pairs[[1]]$spatial_clone,
  rel_widths = c(1.2, 1)
)
```

<img src="/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRUP/docs/readthedocs/markdown/clones-and-trees_files/figure-gfm/tree-pairs-1.jpeg" alt="" width="100%" />

The result is a list with one element per clone. Each element holds a
`tree` plot and a `spatial_clone` plot. Tips with a spot in this section
are circles. Tips without one, such as bulk sequences, are squares and
appear in the tree only. Sequences in the same spot share a point on the
map.

### Color tips by group

When a clone has many tips, one color per tip is hard to read.
`clusterTreesTips()` cuts each tree into `k` connected groups of related
tips. Here `k = 4` for every tree:

``` r
clustered <- clusterTreesTips(
  trees,
  k = 4,
  id = "sequence_id",
  cluster = "tree_cluster"
)

cluster_pairs <- plotSpatialCloneTreePairs(
  clustered,
  object = sirup,
  clone_val = tree_clones,
  id = "sequence_id",
  color_by = "tree_cluster",
  zoom = TRUE,
  polygons = TRUE,
  polygon_alpha = 0,
  polygon_color = "black",
  show_spatial_legend = TRUE
)

cowplot::plot_grid(
  cluster_pairs[[1]]$tree,
  cluster_pairs[[1]]$spatial_clone,
  rel_widths = c(1.2, 1)
)
```

<img src="/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRUP/docs/readthedocs/markdown/clones-and-trees_files/figure-gfm/tree-clusters-1.jpeg" alt="" width="100%" />

With `k = "auto"`, the function chooses `k` for each tree from
`k_range`. It stops when one more cut reduces the mean within-group
patristic distance by less than `improvement_threshold` (10% by
default). A tree with little structure can get `k = 1`, which gives all
its tips one color.

`color_by` accepts any column of the dowser clone data. For example,
`color_by = "polygon_id"` colors each tip by its eGC, because we kept
`polygon_id` in `traits`.

### Save every pair

`spatialCloneTreePairsToPDF()` writes one page per clone:

``` r
spatialCloneTreePairsToPDF(cluster_pairs, "tree_pairs.pdf", width = 12, height = 6)
```

## Tree distance and physical distance

Do sequences in the same eGC have shorter tree distances than sequences
in different eGCs? The patristic distance between two tips is the sum of
the branch lengths between them. `tipPatristicMatrix()` gives these
distances for a set of tips. Here we take the selected clone in the most
eGCs, and compare its tip pairs in the same eGC with its pairs in
different eGCs:

``` r
multi_gc_clone <- clone_sizes |>
  filter(clone_id %in% tree_clones) |>
  slice_max(n_gcs, n = 1, with_ties = FALSE) |>
  pull(clone_id)
multi_gc_clone
#> [1] "739"

index <- which(trees$clone_id == multi_gc_clone)
clone_data <- trees$data[[index]]@data
tree <- trees$trees[[index]]

gc_tips <- clone_data |>
  filter(!is.na(polygon_id), sequence_id %in% tree$tip.label) |>
  select(sequence_id, polygon_id)

distances <- tipPatristicMatrix(tree, gc_tips$sequence_id)

pairs_table <- as.data.frame(as.table(distances)) |>
  rename(tip_a = Var1, tip_b = Var2, distance = Freq) |>
  filter(as.character(tip_a) < as.character(tip_b)) |>
  mutate(
    gc_a = gc_tips$polygon_id[match(tip_a, gc_tips$sequence_id)],
    gc_b = gc_tips$polygon_id[match(tip_b, gc_tips$sequence_id)],
    pair = if_else(gc_a == gc_b, "same eGC", "different eGCs")
  )

pairs_table |>
  group_by(pair) |>
  summarise(n_pairs = n(), median_distance = median(distance))
#> # A tibble: 2 × 3
#>   pair           n_pairs median_distance
#>   <chr>            <int>           <dbl>
#> 1 different eGCs       5          0.0273
#> 2 same eGC            10          0.0178
```

In this clone, tips in the same eGC are closer in the tree than tips in
different eGCs. One clone is not enough for a conclusion, because tip
pairs share branches. In the SIRUP paper we ran this comparison over all
GC clones, with a mixed model and a random effect for each clone.

With `collapse = TRUE` (the default), `formatClones()` merges identical
sequences into one tip, so not every sequence keeps its own tip.

## Session information

``` r
sessionInfo()
#> R version 4.5.1 (2025-06-13)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Red Hat Enterprise Linux 8.10 (Ootpa)
#> 
#> Matrix products: default
#> BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.15.so;  LAPACK version 3.9.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: America/New_York
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] scoper_1.5.0     ggplot2_4.0.3    dplyr_1.2.1      SIRUP_0.0.0.9001
#> [5] testthat_3.3.2  
#> 
#> loaded via a namespace (and not attached):
#>   [1] fs_2.1.0                    matrixStats_1.5.0          
#>   [3] spatstat.sparse_3.2-0       bitops_1.1-0               
#>   [5] devtools_2.5.2              httr_1.4.8                 
#>   [7] RColorBrewer_1.1-3          doParallel_1.0.17          
#>   [9] tools_4.5.1                 sctransform_0.4.3          
#>  [11] utf8_1.2.6                  R6_2.6.1                   
#>  [13] lazyeval_0.2.3              uwot_0.2.5                 
#>  [15] withr_3.0.3                 sp_2.2-3                   
#>  [17] prettyunits_1.2.0           gridExtra_2.3.1            
#>  [19] progressr_1.0.0             cli_3.6.6                  
#>  [21] Biobase_2.70.0              phylotate_1.3              
#>  [23] spatstat.explore_3.8-1      fastDummies_1.7.6          
#>  [25] labeling_0.4.3              diptest_0.77-2             
#>  [27] Seurat_5.5.1                S7_0.2.2                   
#>  [29] spatstat.data_3.1-9         readr_2.2.0                
#>  [31] ggridges_0.5.7              pbapply_1.7-4              
#>  [33] Rsamtools_2.26.0            systemfonts_1.3.2          
#>  [35] yulab.utils_0.2.5           dichromat_2.0-1            
#>  [37] kedd_1.0.4                  parallelly_1.48.0          
#>  [39] sessioninfo_1.2.4           maps_3.4.3                 
#>  [41] pals_1.10                   rstudioapi_0.19.0          
#>  [43] generics_0.1.4              gridGraphics_0.5-1         
#>  [45] ica_1.0-3                   spatstat.random_3.5-0      
#>  [47] vroom_1.7.1                 Matrix_1.7-3               
#>  [49] S4Vectors_0.48.1            abind_1.4-8                
#>  [51] lifecycle_1.0.5             yaml_2.3.12                
#>  [53] SummarizedExperiment_1.40.0 SparseArray_1.8.0          
#>  [55] Rtsne_0.17                  grid_4.5.1                 
#>  [57] alakazam_1.4.3              promises_1.5.0             
#>  [59] crayon_1.5.3                pwalign_1.6.0              
#>  [61] miniUI_0.1.2                lattice_0.22-7             
#>  [63] cowplot_1.2.0               cigarillo_1.0.0            
#>  [65] mapproj_1.2.12              pillar_1.11.1              
#>  [67] knitr_1.51                  GenomicRanges_1.62.1       
#>  [69] future.apply_1.20.2         codetools_0.2-20           
#>  [71] fastmatch_1.1-8             glue_1.8.1                 
#>  [73] ggiraph_0.9.6               ggfun_0.2.1                
#>  [75] spatstat.univar_3.2-0       fontLiberation_0.1.0       
#>  [77] data.table_1.18.4           vctrs_0.7.3                
#>  [79] png_0.1-9                   treeio_1.34.0              
#>  [81] spam_2.11-4                 gtable_0.3.6               
#>  [83] cachem_1.1.0                xfun_0.60                  
#>  [85] S4Arrays_1.10.1             mime_0.13                  
#>  [87] Seqinfo_1.0.0               coda_0.19-4.1              
#>  [89] survival_3.8-9              iterators_1.0.14           
#>  [91] ellipsis_0.3.3              fitdistrplus_1.2-6         
#>  [93] ROCR_1.0-12                 nlme_3.1-168               
#>  [95] ggtree_4.0.5                usethis_3.2.1              
#>  [97] bit64_4.8.2                 fontquiver_0.2.1           
#>  [99] progress_1.2.3              RcppAnnoy_0.0.23           
#> [101] rprojroot_2.1.1             irlba_2.3.7                
#> [103] airr_2.0.0                  KernSmooth_2.23-26         
#> [105] otel_0.2.0                  colorspace_2.1-3           
#> [107] BiocGenerics_0.56.0         ade4_1.7-24                
#> [109] phangorn_2.12.1             tidyselect_1.2.1           
#> [111] bit_4.6.0                   compiler_4.5.1             
#> [113] desc_1.4.3                  fontBitstreamVera_0.1.1    
#> [115] DelayedArray_0.36.1         plotly_4.12.1              
#> [117] scales_1.4.0                lmtest_0.9-40              
#> [119] quadprog_1.5-8              rappdirs_0.3.4             
#> [121] stringr_1.6.0               digest_0.6.39              
#> [123] goftest_1.2-3               spatstat.utils_3.2-4       
#> [125] dowser_2.5.1.999            rmarkdown_2.31             
#> [127] XVector_0.50.0              htmltools_0.5.9            
#> [129] pkgconfig_2.0.3             MatrixGenerics_1.22.0      
#> [131] fastmap_1.2.0               rlang_1.3.0                
#> [133] htmlwidgets_1.6.4           shiny_1.14.0               
#> [135] farver_2.1.2                zoo_1.9-0                  
#> [137] jsonlite_2.0.0              BiocParallel_1.38.0        
#> [139] magrittr_2.0.5              ggplotify_0.1.3            
#> [141] dotCall64_1.2               patchwork_1.3.2            
#> [143] Rcpp_1.1.2                  ape_5.8-1                  
#> [145] gdtools_0.5.1               reticulate_1.46.0          
#> [147] stringi_1.8.9               brio_1.1.5                 
#> [149] MASS_7.3-65                 plyr_1.8.9                 
#> [151] pkgbuild_1.4.8              parallel_4.5.1             
#> [153] listenv_1.0.0               ggrepel_0.9.8              
#> [155] deldir_2.0-4                Biostrings_2.78.0          
#> [157] shazam_1.3.2                splines_4.5.1              
#> [159] tensor_1.5.1                hms_1.1.4                  
#> [161] igraph_2.3.3                fastcluster_1.3.0          
#> [163] spatstat.geom_3.8-1         seqinr_4.2-44              
#> [165] markdown_2.0                RcppHNSW_0.7.0             
#> [167] reshape2_1.4.5              stats4_4.5.1               
#> [169] pkgload_1.5.3               evaluate_1.0.5             
#> [171] SeuratObject_5.4.0          tzdb_0.5.0                 
#> [173] foreach_1.5.2               httpuv_1.6.17              
#> [175] RANN_2.6.3                  tidyr_1.3.2                
#> [177] purrr_1.2.2                 polyclip_1.10-7            
#> [179] future_1.75.0               scattermore_1.2            
#> [181] xtable_1.8-8                RSpectra_0.16-2            
#> [183] tidytree_0.4.8              later_1.4.8                
#> [185] viridisLite_0.4.3           tibble_3.3.1               
#> [187] aplot_0.3.1                 memoise_2.0.1              
#> [189] GenomicAlignments_1.46.0    IRanges_2.44.0             
#> [191] cluster_2.1.8.3             globals_0.19.1
```
