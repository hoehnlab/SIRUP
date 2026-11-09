Get started with SIRUP
================

<div markdown="1"
style="background: #F7E5BC; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<b>Time</b>: ~15 minutes. <b>Data</b>: example MG thymus data in SIRUP.

</div>

Here, we will build a SIRUP object from example data and run a few
simple analyses.

<div markdown="1"
style="background: #9e584a; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #8b3620; color: #ffffff; font-weight: 700;">

<b>Note</b>
</p>

This page starts after read processing. If you have raw or demultiplexed
reads, see [Process demultiplexed reads](processing.md) first.

</div>

## Install SIRUP

Install SIRUP and its dependencies from GitHub with BiocManager:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("hoehnlab/SIRUP")
```

This vignette uses a few of SIRUP’s suggested packages as well. Run the
code below to install any that are missing:

``` r
packages <- c("spatstat.explore", "spatstat.geom", "spatstat.random", "pals")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing)
```

## Look at the example data

SIRUP includes an example spAIRR-seq data set from the thymus of a
patient with myasthenia gravis (MG). The thymus in MG can contain
ectopic germinal centers (eGCs) where B cells undergo affinity
maturation. The data were generated with 10x Visium for spatial
transcriptomics and barcodes, then RAGE-seq + long-read Oxford Nanopore
sequencing to amplify BCRs/TCRs and sequence them.

This vignette uses these files:

| File | Contents |
|----|----|
| `mg_thymus_slice2_spatial.rds` | A spatial Seurat object for MG thymus |
| `mg_thymus_slice2_bcr.tsv.gz` | AIRR rearrangement table of IGH sequences with spatial IDs |
| `mg_thymus_slice2_gc_polygons.csv` | Annotated polygons of the eGCs |
| `mg_thymus_cell_types.csv.gz` | Cell-type proportions for each spot |

The files are in the `extdata` folder of the package. This code finds
the folder:

``` r
extdata <- system.file("extdata", package = "SIRUP")
```

The SIRUP paper describes how the authors made these files. To reduce
file sizes, each Seurat object holds a few hundred genes and a reduced
image.

``` r
library(SIRUP)
library(dplyr)
library(ggplot2)
```

## Load the spatial GEX data and the receptors

Read the Seurat object and the AIRR table:

``` r
thymus <- readRDS(file.path(extdata, "mg_thymus_slice2_spatial.rds"))
thymus
#> An object of class Seurat 
#> 506 features across 3203 samples within 2 assays 
#> Active assay: Spatial (449 features, 0 variable features)
#>  2 layers present: counts, data
#>  1 other assay present: alra
#>  1 spatial field of view present: thymus_slice2

bcrs <- airr::read_rearrangement(file.path(extdata, "mg_thymus_slice2_bcr.tsv.gz"))
bcrs %>%
  select(sequence_id, v_call, j_call, c_call, sequence_alignment) %>%
  head()
#> # A tibble: 6 × 5
#>   sequence_id         v_call                    j_call c_call sequence_alignment
#>   <chr>               <chr>                     <chr>  <chr>  <chr>             
#> 1 9,TGATTCCCGGTTACCT  IGHV2-5*02                IGHJ4… IGHG1  ...ATCACCTTGAAGGA…
#> 2 11,TGATTCCCGGTTACCT IGHV3-30*03,IGHV3-30*06   IGHJ4… IGHG1  .................…
#> 3 2,TGATTCCCGGTTACCT  IGHV5-51*01               IGHJ5… IGHG2  GAGGTGCCGCTGGTGCA…
#> 4 6,TGATTCCCGGTTACCT  IGHV5-51*01               IGHJ5… IGHG2  GAGGTGCCGCTGGTGCA…
#> 5 5,CCTACATTCACAGACG  IGHV3-23*01,IGHV3-23*04,… IGHJ4… IGHG1  .................…
#> 6 3,CCTACATTCACAGACG  IGHV1-2*02,IGHV1-2*04     IGHJ4… IGHG1  CAGGTGCATCTGGTGCA…
```

Each row is one BCR sequence. The processing pipeline wrote the spatial
ID at the end of `sequence_id`, after a comma.

## Match each sequence to a spot

`extractSpatialIds()` is a convenience function to extract the spatial
ID and put it into its own column. The `regex` argument finds the
spatial ID. The `suffix` argument adds `-1`, so that each spatial ID
matches a spot name in the Seurat object.

``` r
bcrs <- extractSpatialIds(
  bcrs,
  id = "sequence_id",
  spatial_id = "barcode_10x",
  regex = ",([ACGT]{16})$",
  suffix = "-1"
)

head(colnames(thymus))
#> [1] "AAACAAGTATCTCCCA-1" "AAACACCAATAACTGC-1" "AAACAGCTTTCAGAAG-1"
#> [4] "AAACAGGGTCTATATT-1" "AAACAGTGTTCCTGGG-1" "AAACATGGTGAGAGGA-1"
mean(bcrs$barcode_10x %in% colnames(thymus))
#> [1] 0.9372528
```

A few sequences have a spatial ID that is not in the Seurat object. The
spot was likely off the tissue or failed the expression quality control.
SIRUP will maintain these receptor sequences in the integrated SIRUP
object, but it will not be able to place them in spatial plots.

## Filter the sequences

Keep only the high-quality sequences. The thresholds used are discussed
in [Process demultiplexed reads](processing.md).

``` r
bcrs <- filterAirr(
  bcrs,
  locus_val = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = NULL,
  spatial_id = "barcode_10x"
)
#> filterAirr kept 2189 of 3793 rows (dropped 1604).
```

## Create a SIRUP object

`createSIRUP()` joins the Seurat object and the receptor table. Tell it
which receptor column holds the spatial ID and which holds the locus.
You can also specify the image scale and coordinate settings when
constructing the object, which are discussed further below.

``` r
sirup <- createSIRUP(
  spatial = thymus,
  receptors = bcrs,
  spatial_id = "barcode_10x",
  locus = "locus",
  image_scale = "lowres",
  swap_coords = FALSE,
  metadata = list(tissue = "MG thymus", section = "slice 2")
)

sirup
#> <SIRUP object>
#> Spatial data:    Seurat, 3,203 spatial IDs, 449 genes, 1 image
#> Receptor data:   tibble, 2,189 receptor sequences x 56 columns
#> Spatial overlap: 817/888 receptor IDs in spatial data (92.0%)
#> Spatial ID field: barcode_10x
#> Analyses:        0

summary(sirup)
#> SIRUP summary
#>   Spatial IDs:      3,203
#>   Receptor records: 2,189
#>   Receptor IDs:     888
#>   Matched IDs:      817 (92.0%)
#>   Matched records:  2,047
#>   Image scale:      lowres
#>   Swap coords:      FALSE
#>   Stored analyses:  0
#> 
#> Per-sample summary
#> # A tibble: 1 × 7
#>   sample  n_spatial_ids n_receptor_records n_receptor_ids n_matched_ids
#>   <chr>           <int>              <int>          <int>         <int>
#> 1 spatial          3203               2189            888           817
#> # ℹ 2 more variables: n_matched_records <int>, id_match_rate <dbl>
```

If a `spatial_id`, `locus`, and/or `clone` column are specified when
constructing the SIRUP object, later functions can find these columns
without including them as arguments.

Two settings control how SIRUP constructs plots over the H&E image:

- `image_scale` is `"lowres"` or `"hires"`. Use the image that tge
  Seurat object stores. `Load10X_Spatial()` stores the low-resolution
  image by default, though it can use the high-resolution image if
  specified in the `image.name` argument.
- `swap_coords` swaps the x and y coordinates. Set it to `TRUE` if the
  spots have the right shape but sit rotated or mirrored over the image.
  This is easiest to determine by plotting the spots as shown below.

``` r
coords <- SeuratObject::GetTissueCoordinates(thymus)
ggplot(coords, aes(x, y)) + geom_point() 
```

![](get-started_files/figure-gfm/check-swapcoords-1.jpeg)

If both arguments are included when constructing the SIRUP object, they
will be stored in the object, and plotting functions will automatically
use them. If they are not specified when constructing the SIRUP object,
both settings should be passed to any plotting functions.

<div markdown="1"
style="background: #9e584a; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #8b3620; color: #ffffff; font-weight: 700;">

<b>Note</b>
</p>

[How the SIRUP object works](sirup-object.md) discusses the SIRUP object
in more detail.

</div>

## Building clones in the SIRUP object

A clone is a group of cells that descend from a common ancestor. Many
SIRUP functions work with clone-level data, but SIRUP does not itself
build clones. The `updateReceptors()` function is a method that can runs
any function on the receptor table and updates the receptor table in the
object with the result. This call runs `hierarchicalClones()` from the
`scoper` package.

``` r
sirup <- updateReceptors(
  sirup,
  scoper::hierarchicalClones,
  threshold = 0.1,
  summarize_clones = FALSE,
  nproc = 1
)
#> Warning in prepare_db(db = db, junction = junction, v_call = v_call, j_call =
#> j_call, : Removed 84 sequences with non ATCG characters.
#> Running defineClonesScoper in bulk mode and only keep heavy chains
sirup
#> <SIRUP object>
#> Spatial data:    Seurat, 3,203 spatial IDs, 449 genes, 1 image
#> Receptor data:   tibble, 2,105 receptor sequences x 57 columns
#> Spatial overlap: 797/868 receptor IDs in spatial data (91.8%)
#> Spatial ID field: barcode_10x
#> Clone ID field:  clone_id
#> Analyses:        0
```

The receptor table now has a `clone_id` column, and SIRUP records it as
the clone field in the object automatically.

## Plot gene expression and receptors on the tissue

`spatialExpressionPlot()` plots gene expression over the H&E image. By
default, it hides spots with zero expression to help the tissue stay
visible, though this can be turned off by setting `zero_tol = -Inf`.

The example object has two assays: `Spatial` (counts and log-normalized
data) and `alra` (ALRA-imputed values), which can help improve
visualization. *AICDA* is a marker for germinal center B cells, and it
is plotted below.

``` r
spatialExpressionPlot(
  sirup,
  features = "AICDA",
  assay = "alra",
  layer = "data",
  title = "AICDA",
  zero_tol = 0,
  point_alpha = 0.5
)
```

![](get-started_files/figure-gfm/expression-1.jpeg)

`plotReceptors()` plots all spots that have a receptor. The `group`
argument is the receptor column that determines the color. All BCRs here
are IGH, so they all are colored the same.

``` r
plotReceptors(
  sirup,
  group = "locus",
  colors = "#0072B2",
  receptor_size = 1.2,
  title = "BCR-positive spots"
)
```

![](get-started_files/figure-gfm/receptors-1.jpeg)

To restrict a receptors plot to a set with only a particular value of
the `group` column, use the `group_val` argument. This plot includes
only the three largest IGH clones.

``` r
top_clones <- getReceptors(sirup) %>%
  count(clone_id, sort = TRUE) %>%
  slice_head(n = 3) %>%
  pull(clone_id)
top_clones
#> [1] "34"  "288" "196"

plotReceptors(
  sirup,
  group = "clone_id",
  group_val = top_clones,
  overlap = "average",
  receptor_size = 1.6,
  title = "Three largest clones"
)
```

![](get-started_files/figure-gfm/clones-1.jpeg)

Because this data was generated with 10X Visium spatial sequencing, each
spot likely contains multiple cells. Resultingly, sequences from
distinct groups may share a spot. By default, `plotReceptors()` gives
the spot a combined label and averages the colors of the two groups. The
`overlap` argument controls this behavior.

`plotReceptorsExpression()` combines the two plots. Spots with gene
expression are filled points. Spots with a receptor have a hollow ring
on top:

``` r
plotReceptorsExpression(
  sirup,
  group = "locus",
  gene = "AICDA",
  assay = "alra",
  layer = "data",
  colors = "#0072B2",
  receptor_size = 1.3,
  feature_size = 2,
  title = "AICDA expression and BCR-positive spots"
)
```

![](get-started_files/figure-gfm/receptors-expression-1.jpeg)

## *In silico* micro-dissection of ROIs

Many spatial questions concern a particular region, such as a germinal
center, a tumor margin, or a lesion. `selectSpatialPolygons()` opens an
app where users can draw polygons over regions of interest using the H&E
image and gene expression (e.g., *AICDA*) or spatial object metadata.

``` r
gc_vertices <- selectSpatialPolygons(
  sirup,
  feature = "AICDA",
  assay = "alra",
  layer = "data",
  plot_image = TRUE,
  save_path = "mg_thymus_slice2_gc_polygons.csv"
)
```

![](figures/polygon-app.jpg)

Click around the edge of a region. Each click snaps to the nearest spot.
**Close polygon** joins the last vertex to the first:

![](figures/polygon-click-1.png)![](figures/polygon-click-2.png)![](figures/polygon-closed.png)

After you click **Done**, the app returns the vertices of every polygon.
It writes them to the CSV file at `save_path` (which can be changed
using the CSV path in the app).

`addPolygons()` attaches polygons to the object. Give it the data frame
from the app (`polygons = gc_vertices`) or the CSV file
(`polygons_path`). SIRUP then determines which polygons, if any, contain
each spot and each receptor. The example data include the annotated eGC
polygons from the SIRUP paper.

``` r
sirup <- addPolygons(
  sirup,
  polygons_path = file.path(extdata, "mg_thymus_slice2_gc_polygons.csv")
)

getPolygons(sirup)
#> # A tibble: 9 × 2
#>   polygon_id vertices        
#>   <chr>      <named list>    
#> 1 1          <tibble [7 × 5]>
#> 2 2          <tibble [4 × 5]>
#> 3 3          <tibble [4 × 5]>
#> 4 4          <tibble [4 × 5]>
#> 5 5          <tibble [6 × 5]>
#> 6 6          <tibble [7 × 5]>
#> 7 7          <tibble [7 × 5]>
#> 8 8          <tibble [4 × 5]>
#> 9 9          <tibble [4 × 5]>

getPolygons(sirup)$vertices[[1]]
#> # A tibble: 7 × 5
#>   spatial_id             x     y image         .sirup_polygon_label
#>   <chr>              <dbl> <dbl> <chr>         <chr>               
#> 1 ATCGACCCAATACAGA-1  300.  327. thymus_slice2 1                   
#> 2 AGTTCCTATTTATGTT-1  307.  327. thymus_slice2 1                   
#> 3 CAGTCTGTATACTGGG-1  310.  320. thymus_slice2 1                   
#> 4 GCGCCGTTCCACGATA-1  307.  314. thymus_slice2 1                   
#> 5 TTATGAATGAAAGGGA-1  300.  314. thymus_slice2 1                   
#> 6 ACATCGATCGTTTACC-1  296.  320. thymus_slice2 1                   
#> 7 ATCGACCCAATACAGA-1  300.  327. thymus_slice2 1
```

Plot the polygons over the H&E image with `plotSpatialPolygons()`.

``` r
plotSpatialPolygons(
  sirup,
  polygon_alpha = 0.3,
  label_outside = TRUE,
  label_offset = 0.02,
  label_color = "white",
  label_leader = FALSE
)
```

![](get-started_files/figure-gfm/plot-spatial-polygons-1.jpeg)

The receptor table now has a `polygon_id` column. Receptors outside
every polygon get `NA`:

``` r
getReceptors(sirup) %>%
  count(polygon_id)
#> # A tibble: 9 × 2
#>   polygon_id     n
#>   <chr>      <int>
#> 1 1             28
#> 2 2             23
#> 3 3             14
#> 4 5             47
#> 5 6              4
#> 6 7             30
#> 7 8              6
#> 8 9              4
#> 9 <NA>        1949
```

## Determine which clones dominate the polygons

`polygonReceptorSummary()` counts the receptors, spots, and clones in
each polygon. It stores the result in the object as an analysis.

``` r
sirup <- polygonReceptorSummary(sirup, polygon_id = "polygon_id")
getAnalysis(sirup, "polygon_receptor_summary")
#> # A tibble: 9 × 5
#>   polygon_id n_receptors n_spatial_ids n_clones top_clone_fraction
#>   <chr>            <int>         <int>    <int>              <dbl>
#> 1 1                   28             6       16              0.143
#> 2 2                   23             3       16              0.174
#> 3 3                   14             3        8              0.286
#> 4 4                    0             0        0             NA    
#> 5 5                   47            12       23              0.128
#> 6 6                    4             4        4              0.25 
#> 7 7                   30             8       12              0.167
#> 8 8                    6             3        4              0.333
#> 9 9                    4             3        4              0.25
```

Clonal dominance measures how much of a region the largest clones fill.
`spatialCloneDominance()` computes the fraction of the receptors in each
polygon that belongs to each clone. `plotSpatialPolygonPies()` draws one
pie chart (or stacked bar chart) for each polygon next to the tissue.

``` r
dominance <- spatialCloneDominance(object = sirup)

dominance$dominance %>%
  group_by(polygon) %>%
  slice_max(proportion, n = 1, with_ties = FALSE) %>%
  select(polygon, clone_id, n_receptors, proportion)
#> # A tibble: 8 × 4
#> # Groups:   polygon [8]
#>   polygon clone_id n_receptors proportion
#>   <chr>   <chr>          <dbl>      <dbl>
#> 1 1       928                4      0.143
#> 2 2       381                4      0.174
#> 3 3       196                4      0.286
#> 4 5       810                6      0.128
#> 5 6       136                1      0.25 
#> 6 7       739                5      0.167
#> 7 8       242                2      0.333
#> 8 9       222                1      0.25

pies <- plotSpatialPolygonPies(
  dominance = dominance,
  object = sirup,
  fill = "clone_id",
  top_n = 5
)
pies$page_1
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text()`).
```

![](get-started_files/figure-gfm/dominance-1.jpeg)

The function returns one plot for each page, with six pies on a page.
The `top_n` argument sets how many clones get their own color in each
pie. The grey `Other` slice holds the remaining clones.

[Draw tissue regions and measure clonal dominance and diversity](regions.md) covers ROIs in more detail.

## Save an analysis in a SIRUP object

`setAnalysis()` stores analysis within the SIRUP object. Stored analyses
can be accessed with `getAnalysis()`.

``` r
sirup <- setAnalysis(sirup, "gc_dominance", dominance)
listAnalyses(sirup)
#> [1] "polygon_receptor_summary" "gc_dominance"

getAnalysis(sirup, "gc_dominance") %>% str()
#> List of 6
#>  $ dominance         :'data.frame':  87 obs. of  11 variables:
#>   ..$ clone_id             : chr [1:87] "928" "720" "859" "17" ...
#>   ..$ v_gene               : chr [1:87] "IGHV4-61*01" "IGHV3-7*04" "IGHV4-34*01" "IGHV5-51*01" ...
#>   ..$ v_gene_family        : chr [1:87] "IGHV4-61" "IGHV3-7" "IGHV4-34" "IGHV5-51" ...
#>   ..$ v_gene_superfamily   : chr [1:87] "IGHV4" "IGHV3" "IGHV4" "IGHV5" ...
#>   ..$ n_receptors          : num [1:87] 4 3 3 2 2 2 2 2 1 1 ...
#>   ..$ proportion           : num [1:87] 0.1429 0.1071 0.1071 0.0714 0.0714 ...
#>   ..$ region_receptor_total: num [1:87] 28 28 28 28 28 28 28 28 28 28 ...
#>   ..$ n_region_spatial_ids : int [1:87] 7 7 7 7 7 7 7 7 7 7 ...
#>   ..$ polygon              : chr [1:87] "1" "1" "1" "1" ...
#>   ..$ polygon_id           : chr [1:87] "1" "1" "1" "1" ...
#>   ..$ polygon_label        : chr [1:87] "1" "1" "1" "1" ...
#>  $ region_spatial_ids:'data.frame':  55 obs. of  3 variables:
#>   ..$ spatial_id: chr [1:55] "ACATCGATCGTTTACC-1" "AGTTCCTATTTATGTT-1" "ATCGACCCAATACAGA-1" "CAGTCTGTATACTGGG-1" ...
#>   ..$ polygon   : chr [1:55] "1" "1" "1" "1" ...
#>   ..$ polygon_id: chr [1:55] "1" "1" "1" "1" ...
#>  $ region_coordinates:'data.frame':  55 obs. of  9 variables:
#>   ..$ cell             : chr [1:55] "ACATCGATCGTTTACC-1" "AGTTCCTATTTATGTT-1" "ATCGACCCAATACAGA-1" "CAGTCTGTATACTGGG-1" ...
#>   ..$ x                : num [1:55] 296 307 300 310 303 ...
#>   ..$ y                : num [1:55] 320 327 327 320 320 ...
#>   ..$ .sirup_y         : num [1:55] -320 -327 -327 -320 -320 ...
#>   ..$ image            : chr [1:55] "thymus_slice2" "thymus_slice2" "thymus_slice2" "thymus_slice2" ...
#>   ..$ .sirup_spatial_id: chr [1:55] "ACATCGATCGTTTACC-1" "AGTTCCTATTTATGTT-1" "ATCGACCCAATACAGA-1" "CAGTCTGTATACTGGG-1" ...
#>   ..$ polygon          : chr [1:55] "1" "1" "1" "1" ...
#>   ..$ polygon_id       : chr [1:55] "1" "1" "1" "1" ...
#>   ..$ polygon_label    : chr [1:55] "1" "1" "1" "1" ...
#>  $ region_receptors  :'data.frame':  156 obs. of  78 variables:
#>   ..$ .sirup_spatial_id : chr [1:156] "ATCGACCCAATACAGA-1" "ATCGACCCAATACAGA-1" "ATCGACCCAATACAGA-1" "ATCGACCCAATACAGA-1" ...
#>   ..$ sequence_id       : chr [1:156] "2,ATCGACCCAATACAGA" "8,ATCGACCCAATACAGA" "1,ATCGACCCAATACAGA" "10,ATCGACCCAATACAGA" ...
#>   ..$ sequence          : chr [1:156] "GAAGTGCAGCTGGTGCAGTCCGGAGCAGAGGTGAAAAAGCCCGGGGAGTCTCTGAGGATCTCCTGTAAGGGTTCTGGATACAGCTTTACCAGCTACTGGATNAGCNGGGTG"| __truncated__ "GAGGTGCAGCTGGTGCAGTCTGGAGCAGAGGTGAAAAAGCCCGGGGAGTCTCTGAAGATCTCCTGTAAGGGTTCTGGATACAGCTTTACCAGCTACTGGATCGGCTGGGTG"| __truncated__ "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGAAGCCTGGGGGGTCCCTTACACTCTCCTGTGCAGCCGCTGGTTTCACTTTCAATAAGGCCTGGATGAACTGGGTC"| __truncated__ "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTTAGACTCTCCTGTGCAGCCTCTGGATTCACTTTCAATTACGCCTGGATGAGCTGGGTC"| __truncated__ ...
#>   ..$ rev_comp          : logi [1:156] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ productive        : logi [1:156] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>   ..$ v_call            : chr [1:156] "IGHV5-10-1*03" "IGHV5-51*01" "IGHV3-15*07" "IGHV3-15*01" ...
#>   ..$ d_call            : chr [1:156] "IGHD4-17*01" "IGHD3-3*01,IGHD3-3*02" "IGHD3-22*01" "IGHD1-26*01" ...
#>   ..$ j_call            : chr [1:156] "IGHJ6*02" "IGHJ4*02" "IGHJ3*02" "IGHJ4*02" ...
#>   ..$ sequence_alignment: chr [1:156] "GAAGTGCAGCTGGTGCAGTCCGGAGCA...GAGGTGAAAAAGCCCGGGGAGTCTCTGAGGATCTCCTGTAAGGGTTCTGGATACAGCTTT............ACCAGCTAC"| __truncated__ "GAGGTGCAGCTGGTGCAGTCTGGAGCA...GAGGTGAAAAAGCCCGGGGAGTCTCTGAAGATCTCCTGTAAGGGTTCTGGATACAGCTTT............ACCAGCTAC"| __truncated__ "GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTGAAGCCTGGGGGGTCCCTTACACTCTCCTGTGCAGCCGCTGGTTTCACTTTC............AATAAGGCC"| __truncated__ "GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTACAGCCTGGGGGGTCCCTTAGACTCTCCTGTGCAGCCTCTGGATTCACTTTC............AATTACGCC"| __truncated__ ...
#>   ..$ germline_alignment: chr [1:156] "GAAGTGCAGCTGGTGCAGTCCGGAGCA...GAGGTGAAAAAGCCCGGGGAGTCTCTGAGGATCTCCTGTAAGGGTTCTGGATACAGCTTT............ACCAGCTAC"| __truncated__ "GAGGTGCAGCTGGTGCAGTCTGGAGCA...GAGGTGAAAAAGCCCGGGGAGTCTCTGAAGATCTCCTGTAAGGGTTCTGGATACAGCTTT............ACCAGCTAC"| __truncated__ "GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTAAAGCCTGGGGGGTCCCTTAGACTCTCCTGTGCAGCCTCTGGTTTCACTTTC............AGTAACGCC"| __truncated__ "GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTAAAGCCTGGGGGGTCCCTTAGACTCTCCTGTGCAGCCTCTGGATTCACTTTC............AGTAACGCC"| __truncated__ ...
#>   ..$ junction          : chr [1:156] "TGTGCGCGCAGCTACGGTGACTACGAGGAGGACTACTACTACTACGGTATGGACGTCTGG" "TGTGCGAGGTTTTGGAGTGGTTACCACGGAACACCTTACTACTTTGACTACTGG" "TGTACCACTCCGCGGCGGACTTATGATAGTGGTGCTTATTATGAACCGGATGCTTTTGATATCTGG" "TGTGCCACGATGGGAGGCAGTGGGAGCTACAAGGTCTACTGG" ...
#>   ..$ junction_aa       : chr [1:156] "CARSYGDYEEDYYYYGMDVW" "CARFWSGYHGTPYYFDYW" "CTTPRRTYDSGAYYEPDAFDIW" "CATMGGSGSYKVYW" ...
#>   ..$ v_cigar           : chr [1:156] "101=1X3=1X185=" "293=" "35=1X19=1X16=1X15=1X3=1X38=1X16=1X5=1X5=1X38=1X44=1X54=" "36=1X51=1X1=1X68=1X70=1X1=1X39=1X21=1X4=" ...
#>   ..$ d_cigar           : chr [1:156] "296S3N13=" "294S10N14=" "312S7N9=1X3=1X6=" "309S5N12=" ...
#>   ..$ j_cigar           : chr [1:156] "316S3N60=" "322S45=" "339S1N49=" "326S10N38=" ...
#>   ..$ stop_codon        : logi [1:156] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ vj_in_frame       : logi [1:156] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>   ..$ locus             : chr [1:156] "IGH" "IGH" "IGH" "IGH" ...
#>   ..$ c_call            : chr [1:156] "IGHD" "IGHM" "IGHA1" "IGHM" ...
#>   ..$ junction_length   : int [1:156] 60 54 66 42 42 42 48 36 51 36 ...
#>   ..$ np1_length        : int [1:156] 5 1 13 10 4 4 10 4 22 4 ...
#>   ..$ np2_length        : int [1:156] 7 14 7 5 2 2 9 0 4 0 ...
#>   ..$ v_sequence_start  : int [1:156] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ v_sequence_end    : int [1:156] 291 293 299 299 296 296 288 296 296 293 ...
#>   ..$ v_germline_start  : int [1:156] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ v_germline_end    : int [1:156] 315 317 317 317 320 320 315 320 317 320 ...
#>   ..$ d_sequence_start  : int [1:156] 297 295 313 310 301 301 299 301 319 298 ...
#>   ..$ d_sequence_end    : int [1:156] 309 308 332 321 313 313 305 312 324 309 ...
#>   ..$ d_germline_start  : int [1:156] 4 11 8 6 7 7 4 7 15 7 ...
#>   ..$ d_germline_end    : int [1:156] 16 24 27 17 19 19 10 18 20 18 ...
#>   ..$ j_sequence_start  : int [1:156] 317 323 340 327 316 316 315 313 329 310 ...
#>   ..$ j_sequence_end    : int [1:156] 376 367 388 364 354 357 361 352 370 349 ...
#>   ..$ j_germline_start  : int [1:156] 4 1 2 11 6 6 2 9 7 9 ...
#>   ..$ j_germline_end    : int [1:156] 63 45 50 48 45 48 48 48 48 48 ...
#>   ..$ v_score           : num [1:156] 449 458 433 442 453 453 430 453 400 445 ...
#>   ..$ v_identity        : num [1:156] 0.993 1 0.963 0.973 0.99 ...
#>   ..$ v_support         : num [1:156] 4.46e-128 6.66e-131 2.25e-123 3.24e-126 1.64e-129 ...
#>   ..$ d_score           : num [1:156] 25.7 27.6 27.6 23.8 25.7 25.7 14.1 23.8 12.2 23.8 ...
#>   ..$ d_identity        : num [1:156] 1 1 0.9 1 1 1 1 1 1 1 ...
#>   ..$ d_support         : num [1:156] 0.005 0.001 0.001 0.017 0.004 0.004 13 0.017 52 0.016 ...
#>   ..$ j_score           : num [1:156] 116 87.2 94.9 73.7 62.2 68 91.1 77.6 69.9 77.6 ...
#>   ..$ j_identity        : num [1:156] 1 1 1 1 0.975 ...
#>   ..$ j_support         : num [1:156] 4.96e-30 2.33e-21 1.19e-23 2.60e-17 7.49e-14 ...
#>   ..$ fwr1              : chr [1:156] "GAAGTGCAGCTGGTGCAGTCCGGAGCA...GAGGTGAAAAAGCCCGGGGAGTCTCTGAGGATCTCCTGTAAGGGTTCT" "GAGGTGCAGCTGGTGCAGTCTGGAGCA...GAGGTGAAAAAGCCCGGGGAGTCTCTGAAGATCTCCTGTAAGGGTTCT" "GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTGAAGCCTGGGGGGTCCCTTACACTCTCCTGTGCAGCCGCT" "GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTACAGCCTGGGGGGTCCCTTAGACTCTCCTGTGCAGCCTCT" ...
#>   ..$ fwr2              : chr [1:156] "ATNAGCNGGGTGCGCCAGATGCCCGGGAAAGGCCTGGAGTGGATGGGGAGG" "ATCGGCTGGGTGCGCCAGATGCCCGGGAAAGGCCTGGAGTGGATGGGGATC" "ATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGACTGGAGTGGGTCGGCCAT" "ATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTGGCCGT" ...
#>   ..$ fwr3              : chr [1:156] "AACTACAGCCCGTCCTTCCAA...GGCCACGTCACCATCTCAGCTGACAAGTCCATCAGCACTGCCTACCTGCAGTGGAGCAGCCTGAAGGCCTCGGACACCGCCATGTATTACTGT" "AGATACAGCCCGTCCTTCCAA...GGCCAGGTCACCATCTCAGCCGACAAGTCCATCAGCACCGCCTACCTGCAGTGGAGCAGCCTGAAGGCCTCGGACACCGCCATGTATTACTGT" "GACTACGCTGCACCCGTGACA...GGCAGATTCACCATCTCAAGAGATGATTCAAAAAACACGCTGTTTCTGCAAATGAACAGCCTGAAAACCGAGGACACAGCCGTGTATTACTGT" "GACTACGCTGCACCCGTGAAA...GGCAGATTCACCATCTCAAGAGATGATTCGAGAAACACGCTGTATCTGCAAATGAACAGCCTGAAAACCGACGACACAGCCGTGTATTACTGT" ...
#>   ..$ fwr4              : chr [1:156] "TGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG" "TGGGGCCAGGGAACCCTGGTCACCGTCTCCT" "TGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG" "TGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG" ...
#>   ..$ cdr1              : chr [1:156] "GGATACAGCTTT............ACCAGCTACTGG" "GGATACAGCTTT............ACCAGCTACTGG" "GGTTTCACTTTC............AATAAGGCCTGG" "GGATTCACTTTC............AATTACGCCTGG" ...
#>   ..$ cdr2              : chr [1:156] "ATTGATCCTAGT......GACTCTTATACC" "ATCTATCCTGGT......GACTCTGATACC" "ATTAGAAGCAGAACTGATGGTGGGACAACA" "ATTAAAAGCGAAACTGATGGTGGGACAACA" ...
#>   ..$ cdr3              : chr [1:156] "GCGCGCAGCTACGGTGACTACGAGGAGGACTACTACTACTACGGTATGGACGTC" "GCGAGGTTTTGGAGTGGTTACCACGGAACACCTTACTACTTTGACTAC" "ACCACTCCGCGGCGGACTTATGATAGTGGTGCTTATTATGAACCGGATGCTTTTGATATC" "GCCACGATGGGAGGCAGTGGGAGCTACAAGGTCTAC" ...
#>   ..$ consensus_count   : int [1:156] 1 1 1 1 1 2 1 1 1 1 ...
#>   ..$ cc_bc_umi         : chr [1:156] "6" "10" "3" "8" ...
#>   ..$ n_productive_reads: chr [1:156] "2" "4" "1" "4" ...
#>   ..$ c_call_freq       : chr [1:156] "1" "1" "1" "1" ...
#>   ..$ barcode_10x       : chr [1:156] "ATCGACCCAATACAGA-1" "ATCGACCCAATACAGA-1" "ATCGACCCAATACAGA-1" "ATCGACCCAATACAGA-1" ...
#>   ..$ v_sequence_length : int [1:156] 291 293 299 299 296 296 288 296 296 293 ...
#>   ..$ j_sequence_length : int [1:156] 60 45 49 38 39 42 47 40 42 40 ...
#>   ..$ clone_id          : chr [1:156] "1" "17" "313" "322" ...
#>   ..$ polygon_id        : chr [1:156] "1" "1" "1" "1" ...
#>   ..$ polygon_1         : logi [1:156] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>   ..$ polygon_2         : logi [1:156] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ polygon_3         : logi [1:156] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ polygon_4         : logi [1:156] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ polygon_5         : logi [1:156] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ polygon_6         : logi [1:156] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ polygon_7         : logi [1:156] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ polygon_8         : logi [1:156] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ polygon_9         : logi [1:156] FALSE FALSE FALSE FALSE FALSE FALSE ...
#>   ..$ cell              : chr [1:156] "ATCGACCCAATACAGA-1" "ATCGACCCAATACAGA-1" "ATCGACCCAATACAGA-1" "ATCGACCCAATACAGA-1" ...
#>   ..$ x                 : num [1:156] 300 300 300 300 300 ...
#>   ..$ y                 : num [1:156] 327 327 327 327 327 ...
#>   ..$ image             : chr [1:156] "thymus_slice2" "thymus_slice2" "thymus_slice2" "thymus_slice2" ...
#>   ..$ .sirup_y          : num [1:156] -327 -327 -327 -327 -327 ...
#>   ..$ polygon           : chr [1:156] "1" "1" "1" "1" ...
#>   ..$ polygon_label     : chr [1:156] "1" "1" "1" "1" ...
#>   ..$ v_gene            : chr [1:156] "IGHV5-10-1*03" "IGHV5-51*01" "IGHV3-15*07" "IGHV3-15*01" ...
#>   ..$ v_gene_family     : chr [1:156] "IGHV5-10-1" "IGHV5-51" "IGHV3-15" "IGHV3-15" ...
#>   ..$ v_gene_superfamily: chr [1:156] "IGHV5" "IGHV5" "IGHV3" "IGHV3" ...
#>  $ polygons          :'data.frame':  47 obs. of  4 variables:
#>   ..$ x       : num [1:47] 300 307 310 307 300 ...
#>   ..$ y       : num [1:47] 327 327 320 314 314 ...
#>   ..$ polygon : chr [1:47] "1" "1" "1" "1" ...
#>   ..$ .sirup_y: num [1:47] -327 -327 -320 -314 -314 ...
#>  $ metadata          :List of 1
#>   ..$ receptor_fields:List of 5
#>   .. ..$ spatial_id: chr "barcode_10x"
#>   .. ..$ clone     : chr "clone_id"
#>   .. ..$ locus     : chr "locus"
#>   .. ..$ count     : NULL
#>   .. ..$ v_call    : chr "v_call"
#>  - attr(*, "class")= chr [1:2] "sirup_spatial_clone_dominance" "list"
```

## Save and read a SIRUP object

`saveSIRUP()` writes the object to one RDS file. `readSIRUP()` reads it
back and validates it.

``` r
path <- file.path(tempdir(), "mg_thymus_slice2.sirup.rds")
saveSIRUP(sirup, path)
sirup <- readSIRUP(path)
sirup
#> <SIRUP object>
#> Spatial data:    Seurat, 3,203 spatial IDs, 449 genes, 1 image
#> Receptor data:   tibble, 2,105 receptor sequences x 67 columns
#> Spatial overlap: 797/868 receptor IDs in spatial data (91.8%)
#> Spatial ID field: barcode_10x
#> Clone ID field:  clone_id
#> Analyses:        2
#> Polygons:        9
```

## Next steps

- To learn more about the structure of the SIRUP object see
  [How the SIRUP object works](sirup-object.md).
- To compare ROIs in the tissue, see
  [Draw tissue regions and measure clonal dominance and diversity](regions.md).
- To build lineage trees and visualize them on the tissue, read
  [Call clones and build lineage trees on the tissue](clones-and-trees.md).
- To test spatial patterns, read
  [Test spatial patterns of genes, cells, and clones](spatial-statistics.md).
- To process your own reads, read
  [Process demultiplexed reads](processing.md).
