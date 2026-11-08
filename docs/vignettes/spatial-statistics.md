Spatial statistics for genes, cells, and clones
================

A map can suggest that a gene clusters, or that a clone stays in one
place. A spatial statistic tells you whether the pattern is stronger
than chance. SIRUP has five tests. Each one answers a different
question:

| Question | Function | Statistic |
|----|----|----|
| Does one variable cluster in space? | `spatialMoran()` | Global Moran’s I |
| Where are the clusters? | `spatialLocalMoran()` | Local Moran’s I |
| Do two variables have the same spatial pattern? | `spatialLee()` | Lee’s L |
| Are spots of type A near spots of type B? | `spatialRipleyCrossK()` | Ripley’s cross K |
| Do the spots of one clone sit next to each other? | `spatialCloneJoinCount()` | Join count |

`spatialCloneSpread()` also measures how far each clone spreads across a
section.

The tests use the spdep and spatstat packages. Each one accepts a SIRUP
object and reads the spot coordinates and image settings from it. The
examples use MG thymus section 2.

## Set up

``` r
library(SIRUP)
library(dplyr)
library(ggplot2)
```

``` r
extdata <- system.file("extdata", package = "SIRUP")
read_receptors <- function(file, locus) {
  airr::read_rearrangement(file.path(extdata, file)) %>%
    extractSpatialIds(
      id = "sequence_id",
      spatial_id = "barcode_10x",
      regex = ",([ACGT]{16})$",
      suffix = "-1"
    ) %>%
    filterAirr(
      locus_val = locus,
      min_v_length = 250,
      min_j_length = 30,
      min_junction_length = 10,
      productive_only = TRUE,
      min_consensus_count = NULL,
      spatial_id = "barcode_10x"
    ) %>%
    filter(!grepl("[^ACGT]", junction))
}

thymus <- readRDS(file.path(extdata, "mg_thymus_slice2_spatial.rds"))

cell_types <- data.table::fread(
  file.path(extdata, "mg_thymus_cell_types.csv.gz"),
  data.table = FALSE
) %>%
  filter(sample == "mg_thymus_slice2")
rownames(cell_types) <- cell_types$barcode
thymus <- Seurat::AddMetaData(thymus, select(cell_types, -sample, -barcode))

bcr_sirup <- createSIRUP(
  spatial = thymus,
  receptors = read_receptors("mg_thymus_slice2_bcr.tsv.gz", "IGH"),
  spatial_id = "barcode_10x",
  image_scale = "lowres",
  swap_coords = FALSE,
  polygons_path = file.path(extdata, "mg_thymus_slice2_gc_polygons.csv")
)
#> filterAirr kept 2195 of 3801 rows (dropped 1606).
bcr_sirup <- updateReceptors(
  bcr_sirup,
  scoper::hierarchicalClones,
  threshold = 0.1,
  summarize_clones = FALSE,
  nproc = 1
)
#> Running defineClonesScoper in bulk mode and only keep heavy chains

tcr_sirup <- createSIRUP(
  spatial = thymus,
  receptors = read_receptors("mg_thymus_slice2_tcr.tsv.gz", "TRB"),
  spatial_id = "barcode_10x",
  image_scale = "lowres",
  swap_coords = FALSE
)
#> filterAirr kept 3718 of 7915 rows (dropped 4197).
tcr_sirup <- updateReceptors(
  tcr_sirup,
  scoper::identicalClones,
  summarize_clones = FALSE,
  nproc = 1
)
#> Running defineClonesScoper in bulk mode and only keep heavy chains
```

`read_receptors()` reads a consensus table, takes the spot barcode from
`sequence_id`, and keeps the high-quality sequences of one locus with
unambiguous junctions. BCR clones come from
`scoper::hierarchicalClones()`. TCRs do not mutate, so TCR clones come
from `scoper::identicalClones()`.

## Neighbors

Most of these statistics compare each spot with its neighbors. You
choose the neighbors in one of two ways:

- `k`: the `k` nearest spots. Use `k = 6` for standard Visium, where
  each spot has six neighbors on a hexagonal grid. Use `k = 8` for
  Visium HD, where each bin has eight neighbors on a square grid. Larger
  values of `k` measure patterns over larger distances.
- `radius`: all spots within a distance. SIRUP measures distances in
  image pixels at the stored `image_scale`, so the same radius in
  micrometers needs a different number at `"lowres"` and at `"hires"`.

`k` does not depend on the image scale, so it is the safer default.

## Global Moran’s I

Moran’s I measures spatial autocorrelation: whether neighboring spots
have similar values.

$$I = \frac{N}{W} \frac{\sum_i \sum_j w_{ij}(x_i - \bar{x})(x_j - \bar{x})}{\sum_i (x_i - \bar{x})^2}$$

Here $x_i$ is the value at spot $i$, $w_{ij}$ is 1 if spot $j$ is a
neighbor of spot $i$ and 0 if not, $N$ is the number of spots, and $W$
is the sum of all $w_{ij}$. $I$ near 0 means no spatial pattern.
Positive $I$ means that similar values cluster. Negative $I$ means that
they alternate, as on a chessboard.

`variable_type` tells SIRUP what `variable` is:

- `"feature"`: a gene, from the assay and layer you give.
- `"metadata"`: a numeric spot metadata column. For a text column, give
  `value` to test the spots that have that value.
- `"receptor"`: 1 for spots with any receptor and 0 for other spots.
  `variable` is only a label here.
- `"clone"`: 1 for spots with the clone that `variable` names.

``` r
aicda <- spatialMoran(
  bcr_sirup,
  variable = "AICDA",
  variable_type = "feature",
  assay = "Spatial",
  layer = "data",
  k = 6
)
aicda$test
#> 
#>  Moran I test under randomisation
#> 
#> data:  dat$value  
#> weights: listw    
#> 
#> Moran I statistic standard deviate = 6.0618, p-value = 6.731e-10
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>      5.944121e-02     -3.123048e-04      9.716874e-05
```

The result is a list. `test` is the `spdep::moran.test()` result: the
statistic, its expected value under no spatial pattern, the variance,
and the p value. For *AICDA*, $I$ is small (0.06) but the p value is
below $10^{-9}$. Visium detects *AICDA* in few spots, and most spots
have a count of zero, so the clustering is real but weak. Read $I$ for
the size of the pattern and the p value for the evidence.

The same function tests cell-type proportions and receptor presence:

``` r
moran_table <- bind_rows(
  tibble(variable = "AICDA expression", I = aicda$test$estimate[[1]], p = aicda$test$p.value),
  {
    res <- spatialMoran(bcr_sirup, variable = "GC B", variable_type = "metadata", k = 6)
    tibble(variable = "GC B cell proportion", I = res$test$estimate[[1]], p = res$test$p.value)
  },
  {
    res <- spatialMoran(bcr_sirup, variable = "BCR", variable_type = "receptor", k = 6)
    tibble(variable = "BCR present", I = res$test$estimate[[1]], p = res$test$p.value)
  },
  {
    res <- spatialMoran(tcr_sirup, variable = "TCR", variable_type = "receptor", k = 6)
    tibble(variable = "TCR present", I = res$test$estimate[[1]], p = res$test$p.value)
  }
)
moran_table
#> # A tibble: 4 × 3
#>   variable                  I         p
#>   <chr>                 <dbl>     <dbl>
#> 1 AICDA expression     0.0594 6.73e- 10
#> 2 GC B cell proportion 0.623  0        
#> 3 BCR present          0.417  0        
#> 4 TCR present          0.332  1.31e-238
```

The GC B cell proportion clusters more strongly ($I$ = 0.62): the
cell2location estimate is smooth and has no zeros from sampling. A p
value of 0 means a value below the smallest number that R can store.

Moran’s I works for binary variables, such as receptor presence, but the
join-count test below suits them better.

## Local Moran’s I

Global Moran’s I gives one number for the whole section. Local Moran’s I
(also called LISA) gives one value for each spot. A spot with a high
positive value has a high (or low) value itself and neighbors with
similar values. Its z-score shows how far it is from the value expected
by chance.

``` r
gc_b_local <- spatialLocalMoran(
  bcr_sirup,
  variable = "GC B",
  variable_type = "metadata",
  k = 6
)
head(gc_b_local[, c("cell", "value", "local_i", "z_i", "p_value")])
#>                 cell  value    local_i       z_i     p_value
#> 1 AAACAAGTATCTCCCA-1 0.0556 0.42119564 1.8775350 0.060444809
#> 2 AAACACCAATAACTGC-1 0.0492 0.02554433 0.1780355 0.858695126
#> 3 AAACAGCTTTCAGAAG-1 0.0128 0.67548360 2.1434191 0.032079462
#> 4 AAACAGGGTCTATATT-1 0.0149 0.44557822 1.5436247 0.122679267
#> 5 AAACAGTGTTCCTGGG-1 0.0251 0.10316582 0.6444111 0.519308914
#> 6 AAACATGGTGAGAGGA-1 0.0915 2.05891309 3.0435870 0.002337758
```

`plotSpatialLocalMoran()` draws the z-scores on the tissue. Red spots
are in clusters of similar values. Blue spots differ from their
neighbors:

``` r
plotSpatialLocalMoran(
  bcr_sirup,
  gc_b_local,
  point_size = 1,
  title = "GC B cell proportion, local Moran's I (z-score)"
) +
  labs(color = "z")
```

<img src="spatial-statistics_files/figure-gfm/local-moran-plot-1.jpeg" alt="" width="100%" />

The clusters of high GC B cell proportion sit where we drew the eGC
polygons. Each spot has its own p value, so correct for multiple testing
(for example with `p.adjust(method = "BH")`) before you call a spot
significant.

## Lee’s L: two variables

Lee’s L combines Pearson correlation with spatial smoothing. It is high
when two variables are high in the same neighborhoods, even if they do
not peak in the same spots.

In MG, patients make antibodies against the acetylcholine receptor
(AChR). *CHRNA1* codes for the AChR α subunit. Thymic myoid cells
express it. Lee’s L tests whether *CHRNA1* and the GC marker *AICDA*
share a spatial pattern in this section. *MS4A1* (CD20) gives a positive
control, because GC B cells express both *MS4A1* and *AICDA*:

``` r
lee_control <- spatialLee(
  bcr_sirup,
  x = "AICDA", y = "MS4A1",
  assay = "Spatial", layer = "data", k = 6
)
lee_chrna1 <- spatialLee(
  bcr_sirup,
  x = "AICDA", y = "CHRNA1",
  assay = "Spatial", layer = "data", k = 6
)

tibble(
  pair = c("AICDA and MS4A1", "AICDA and CHRNA1"),
  L = c(lee_control$L, lee_chrna1$L),
  p = c(lee_control$test$p.value, lee_chrna1$test$p.value)
)
#> # A tibble: 2 × 3
#>   pair                    L        p
#>   <chr>               <dbl>    <dbl>
#> 1 AICDA and MS4A1   0.0975  2.17e-49
#> 2 AICDA and CHRNA1 -0.00171 5.66e- 1
```

The control pair shares a pattern. *AICDA* and *CHRNA1* do not: in this
section, $L$ is near zero and the test is not significant.

The analytical test (`method = "analytical"`) builds a dense matrix with
one entry for each pair of spots. That is fine for Visium. For Visium
HD, use `method = "mc"`, which uses permutations and sparse weights.

## Ripley’s cross K: are A spots near B spots?

Ripley’s cross K counts, for each distance $r$, how many type B points
are within $r$ of an average type A point:

$$K_{AB}(r) = \frac{|W|}{N_A N_B} \sum_{i \in A} \sum_{j \in B} \mathbf{1}(d_{ij} \le r)$$

$|W|$ is the area of the tissue window and $d_{ij}$ is the distance
between points $i$ and $j$. SIRUP compares the observed curve with
curves from random relabeling of the points. If the observed curve is
above the envelope of the random curves, A and B are closer than chance.

T follicular helper (Tfh) cells give B cells the signals that keep a
germinal center going. Here we ask whether spots with GC B cells sit
near spots with Tfh cells. First, label each spot from its cell2location
proportions. `updateSpatial()` runs a function on the Seurat object and
stores the result:

``` r
bcr_sirup <- updateSpatial(bcr_sirup, function(spatial) {
  gc_b <- spatial$`GC B` > 0.10
  tfh <- spatial$Tfh > 0.10
  spatial$gc_tfh_mark <- case_when(
    gc_b & !tfh ~ "GC B",
    tfh & !gc_b ~ "Tfh",
    TRUE ~ NA_character_
  )
  spatial
})

table(getSpatial(bcr_sirup)$gc_tfh_mark, useNA = "ifany")
#> 
#> GC B  Tfh <NA> 
#>   76  211 2916
```

Spots with neither cell type, or with both, get `NA`, and the test
leaves them out. Then run the test. `nsim` sets the number of random
relabelings. Use 99 or more for a final analysis:

``` r
gc_tfh_k <- spatialRipleyCrossK(
  bcr_sirup,
  mark_col = "gc_tfh_mark",
  i = "GC B",
  j = "Tfh",
  nsim = 39,
  max_r = 40
)

plotRipleyCrossK(gc_tfh_k, title = "Cross K: GC B spots to Tfh spots")
```

<img src="spatial-statistics_files/figure-gfm/ripley-1.jpeg" alt="" width="100%" />

The distance $r$ is in image pixels at the stored `image_scale`. To
convert it to micrometers, use the Visium spot size. A spot is 55 µm
wide, and the image object stores its diameter in full-resolution
pixels:

``` r
scale_factors <- getSpatial(bcr_sirup)@images[[1]]@scale.factors
um_per_pixel <- 55 / (scale_factors$spot * scale_factors$lowres)
um_per_pixel
#> [1] 11.7135
max(gc_tfh_k$plot_data$r) * um_per_pixel
#> [1] 468.5401
```

So `max_r = 40` lowres pixels is about 470 µm.

In this section the observed curve is below the envelope at every
distance. Random relabeling keeps the marked spots where they are and
shuffles the labels. A curve below the envelope therefore means that GC
B spots sit next to other GC B spots, and Tfh spots next to other Tfh
spots, more often than next to each other. The two cell types form
separate groups at the 10% threshold. The rule that drops spots with
both cell types adds to the separation. Change the thresholds or the
marks to ask a different question.

## Join counts: does a clone stay together?

The join-count test is the natural test for a clone. Take every spot
with a BCR. Mark the spots that have the clone. Count the neighbor pairs
in which both spots have the clone (the “joins”). If the clone stays in
one place, there are more joins than random labeling gives.

`spatialCloneJoinCount()` runs `spdep::joincount.test()` for each clone
in at least `min_spatial_ids` spots. When the object has stored
polygons, the function runs one test per polygon by default. Set
`use_stored_polygons = FALSE` for one test over the whole section:

``` r
big_clones <- getReceptors(bcr_sirup) %>%
  filter(!is.na(barcode_10x)) %>%
  distinct(clone_id, barcode_10x) %>%
  count(clone_id, sort = TRUE) %>%
  slice_head(n = 5) %>%
  pull(clone_id)

join_count <- spatialCloneJoinCount(
  bcr_sirup,
  clone_val = big_clones,
  use_stored_polygons = FALSE,
  k = 6,
  return_data = TRUE
)

join_count$summary
#>   clone_id         image polygon n_spatial_ids n_clone_spatial_ids
#> 1      598 thymus_slice2     all           797                  11
#> 2      288 thymus_slice2     all           797                  27
#> 3       34 thymus_slice2     all           797                  66
#> 4       64 thymus_slice2     all           797                  10
#> 5      196 thymus_slice2     all           797                  27
#>   observed_join_count expected_join_count   variance statistic      p_value
#> 1                12.5           0.4145729  0.3644266 20.019673 1.855757e-89
#> 2                30.0           2.6457286  2.2476408 18.245749 1.118248e-74
#> 3                67.5          16.1683417 12.6054484 14.457931 1.117067e-47
#> 4                 5.5           0.3391960  0.2987982  9.441226 1.842208e-21
#> 5                15.5           2.6457286  2.2476408  8.574011 4.997114e-18
#>                                   method radius k alternative
#> 1 Join count test under nonfree sampling     NA 6     greater
#> 2 Join count test under nonfree sampling     NA 6     greater
#> 3 Join count test under nonfree sampling     NA 6     greater
#> 4 Join count test under nonfree sampling     NA 6     greater
#> 5 Join count test under nonfree sampling     NA 6     greater
```

The summary has one row per clone. `statistic` is the standardized join
count (a z-score) and `p_value` is the one-sided p value.
`plotSpatialCloneJoinCount()` shows the clone and the joins:

``` r
plotSpatialCloneJoinCount(
  join_count = join_count,
  object = bcr_sirup,
  clone_val = join_count$summary$clone_id[[1]],
  show_neighbors = TRUE
)
```

<img src="spatial-statistics_files/figure-gfm/join-count-plot-1.jpeg" alt="" width="100%" />

The plot needs `return_data = TRUE` in the test.

## Clone spread

`spatialCloneSpread()` describes how far each clone reaches across the
section, from the pairwise distances between its sequences:

``` r
spread <- spatialCloneSpread(bcr_sirup, min_spatial_ids = 3)
spread %>%
  arrange(desc(n_spatial_ids)) %>%
  head()
#> # A tibble: 6 × 5
#>   sample  clone_id n_spatial_ids mean_pairwise_distance max_pairwise_distance
#>   <chr>   <chr>            <int>                  <dbl>                 <dbl>
#> 1 spatial 34                  86                   66.2                 397. 
#> 2 spatial 288                 42                  123.                  304. 
#> 3 spatial 196                 38                  156.                  395. 
#> 4 spatial 739                 16                   21.5                  75.1
#> 5 spatial 810                 14                   80.0                 263. 
#> 6 spatial 907                 14                   54.2                 120.
```

The distances are in image pixels. By default, sequences in the same
spot count as separate points at distance 0. Set
`collapse_within_spot = TRUE` to count each spot once.

## Keep the results

Store results that you want to keep in the object, so they travel with
the data:

``` r
bcr_sirup <- bcr_sirup %>%
  setAnalysis("moran_aicda", aicda) %>%
  setAnalysis("local_moran_gc_b", gc_b_local) %>%
  setAnalysis("lee_aicda_chrna1", lee_chrna1) %>%
  setAnalysis("clone_join_count", join_count) %>%
  setAnalysis("clone_spread", spread)
listAnalyses(bcr_sirup)
#> [1] "moran_aicda"      "local_moran_gc_b" "lee_aicda_chrna1" "clone_join_count"
#> [5] "clone_spread"
```

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
#> [1] sf_1.1-2         scoper_1.5.0     ggplot2_4.0.3    dplyr_1.2.1     
#> [5] SIRUP_0.0.0.9001 testthat_3.3.2  
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
#>  [31] proxy_0.4-29                ggridges_0.5.7             
#>  [33] pbapply_1.7-4               Rsamtools_2.26.0           
#>  [35] systemfonts_1.3.2           yulab.utils_0.2.5          
#>  [37] dbscan_1.2.6                R.utils_2.13.0             
#>  [39] dichromat_2.0-1             kedd_1.0.4                 
#>  [41] parallelly_1.48.0           sessioninfo_1.2.4          
#>  [43] maps_3.4.3                  pals_1.10                  
#>  [45] rstudioapi_0.19.0           generics_0.1.4             
#>  [47] gridGraphics_0.5-1          vroom_1.7.1                
#>  [49] ica_1.0-3                   spatstat.random_3.5-0      
#>  [51] spdep_1.4-2                 Matrix_1.7-3               
#>  [53] S4Vectors_0.48.1            abind_1.4-8                
#>  [55] R.methodsS3_1.8.2           lifecycle_1.0.5            
#>  [57] yaml_2.3.12                 SummarizedExperiment_1.40.0
#>  [59] SparseArray_1.8.0           Rtsne_0.17                 
#>  [61] grid_4.5.1                  alakazam_1.4.3             
#>  [63] promises_1.5.0              crayon_1.5.3               
#>  [65] pwalign_1.6.0               miniUI_0.1.2               
#>  [67] lattice_0.22-7              cowplot_1.2.0              
#>  [69] cigarillo_1.0.0             mapproj_1.2.12             
#>  [71] pillar_1.11.1               knitr_1.51                 
#>  [73] airrblast_0.0.1             GenomicRanges_1.62.1       
#>  [75] boot_1.3-31                 future.apply_1.20.2        
#>  [77] codetools_0.2-20            fastmatch_1.1-8            
#>  [79] wk_0.9.5                    glue_1.8.1                 
#>  [81] ggiraph_0.9.6               ggfun_0.2.1                
#>  [83] spatstat.univar_3.2-0       fontLiberation_0.1.0       
#>  [85] data.table_1.18.4           vctrs_0.7.3                
#>  [87] png_0.1-9                   treeio_1.34.0              
#>  [89] spam_2.11-4                 gtable_0.3.6               
#>  [91] cachem_1.1.0                xfun_0.60                  
#>  [93] S4Arrays_1.10.1             mime_0.13                  
#>  [95] Seqinfo_1.0.0               coda_0.19-4.1              
#>  [97] survival_3.8-9              iterators_1.0.14           
#>  [99] units_0.8-5                 rgl_1.3.36                 
#> [101] ellipsis_0.3.3              fitdistrplus_1.2-6         
#> [103] ROCR_1.0-12                 nlme_3.1-168               
#> [105] ggtree_4.0.5                usethis_3.2.1              
#> [107] bit64_4.8.2                 fontquiver_0.2.1           
#> [109] progress_1.2.3              RcppAnnoy_0.0.23           
#> [111] rprojroot_2.1.1             irlba_2.3.7                
#> [113] airr_2.0.0                  KernSmooth_2.23-26         
#> [115] otel_0.2.0                  colorspace_2.1-3           
#> [117] spData_2.3.5                BiocGenerics_0.56.0        
#> [119] DBI_1.3.0                   ade4_1.7-24                
#> [121] phangorn_2.12.1             tidyselect_1.2.1           
#> [123] bit_4.6.0                   compiler_4.5.1             
#> [125] desc_1.4.3                  fontBitstreamVera_0.1.1    
#> [127] DelayedArray_0.36.1         plotly_4.12.1              
#> [129] scales_1.4.0                classInt_0.4-11            
#> [131] lmtest_0.9-40               quadprog_1.5-8             
#> [133] rappdirs_0.3.4              stringr_1.6.0              
#> [135] digest_0.6.39               goftest_1.2-3              
#> [137] spatstat.utils_3.2-4        dowser_2.5.1.999           
#> [139] rmarkdown_2.31              XVector_0.50.0             
#> [141] base64enc_0.1-6             htmltools_0.5.9            
#> [143] pkgconfig_2.0.3             MatrixGenerics_1.22.0      
#> [145] fastmap_1.2.0               rlang_1.3.0                
#> [147] htmlwidgets_1.6.4           shiny_1.14.0               
#> [149] farver_2.1.2                zoo_1.9-0                  
#> [151] jsonlite_2.0.0              BiocParallel_1.38.0        
#> [153] R.oo_1.27.1                 magrittr_2.0.5             
#> [155] ggplotify_0.1.3             s2_1.1.11                  
#> [157] dotCall64_1.2               patchwork_1.3.2            
#> [159] Rcpp_1.1.2                  ape_5.8-1                  
#> [161] gdtools_0.5.1               reticulate_1.46.0          
#> [163] stringi_1.8.9               brio_1.1.5                 
#> [165] MASS_7.3-65                 plyr_1.8.9                 
#> [167] pkgbuild_1.4.8              parallel_4.5.1             
#> [169] listenv_1.0.0               ggrepel_0.9.8              
#> [171] deldir_2.0-4                Biostrings_2.78.0          
#> [173] shazam_1.3.2                splines_4.5.1              
#> [175] tensor_1.5.1                hms_1.1.4                  
#> [177] igraph_2.3.3                fastcluster_1.3.0          
#> [179] spatstat.geom_3.8-1         seqinr_4.2-44              
#> [181] markdown_2.0                RcppHNSW_0.7.0             
#> [183] reshape2_1.4.5              stats4_4.5.1               
#> [185] pkgload_1.5.3               evaluate_1.0.5             
#> [187] SeuratObject_5.4.0          tzdb_0.5.0                 
#> [189] foreach_1.5.2               httpuv_1.6.17              
#> [191] RANN_2.6.3                  tidyr_1.3.2                
#> [193] purrr_1.2.2                 polyclip_1.10-7            
#> [195] future_1.75.0               scattermore_1.2            
#> [197] lwgeom_0.2-17               xtable_1.8-8               
#> [199] e1071_1.7-17                RSpectra_0.16-2            
#> [201] tidytree_0.4.8              later_1.4.8                
#> [203] class_7.3-23                viridisLite_0.4.3          
#> [205] tibble_3.3.1                aplot_0.3.1                
#> [207] memoise_2.0.1               GenomicAlignments_1.46.0   
#> [209] IRanges_2.44.0              cluster_2.1.8.3            
#> [211] globals_0.19.1
```
