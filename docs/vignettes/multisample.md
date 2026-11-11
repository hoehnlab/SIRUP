Analyze several sections with bulk AIRR-seq and blood
================

<div markdown="1"
style="background: #F7E5BC; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #d6b676; color: #ffffff; font-weight: 700;">

<b>Note</b>
</p>

**Time:** ~25 minutes.<br> **Data:** two MG thymus sections, bulk thymus
BCR-seq, and bulk blood BCR-seq, downloaded with
`downloadSIRUPExample()`.<br> **Needs:** SIRUP, scoper, and dowser. The
tree step also needs the IMGT reference sequences.

</div>

Most studies have more than one tissue section. You can also have bulk
AIRR-seq of the same tissue, or blood from other time points. In this
page, you put all of these samples in one SIRUP object. You call clones
on all samples together. You then ask where each clone was seen, compare
clonal dominance in each section, and build a lineage tree that has
blood tips.

A clone is a group of B or T cells that descend from one ancestor cell
and share its receptor rearrangement. If you call clones on all samples
together, one clone ID names the same lineage in every sample.

One SIRUP object can hold all of this:

- `spatial` can be a named list of Seurat objects, one for each section.
- The receptor table holds the rows from every section and every bulk
  sample. A `sample` column tells SIRUP which sample each row belongs
  to.
- Rows whose `sample` does not name a section stay in the table as
  nonspatial rows. SIRUP does not put them on a tissue. Clone calling,
  trees, and your own summaries can use them.

This page follows Figure 5 of the SIRUP paper. The MG thymus data come
from a patient who had a thymectomy. Besides the two Visium sections,
there is bulk BCR-seq of the thymus. There is also bulk BCR-seq of
blood, taken at the thymectomy and 12 months later. The question is
which clones in ectopic germinal centers (eGCs) also appear in the
blood. An eGC is a site in the thymus where B cells proliferate and
mutate their receptors.

## Get the example data

This page uses example files that are not part of the SIRUP package.
This code downloads them. Change `dest` to a folder you want to keep.
The default folder is temporary. The function skips files that are
already there.

``` r
extdata <- SIRUP::downloadSIRUPExample(
  dest = file.path(tempdir(), "SIRUP_example"),
  files = c(
  "mg_thymus_slice1_spatial.rds",
  "mg_thymus_slice2_spatial.rds",
  "mg_thymus_slice1_bcr.tsv.gz",
  "mg_thymus_slice2_bcr.tsv.gz",
  "mg_thymus_slice1_gc_polygons.csv",
  "mg_thymus_slice2_gc_polygons.csv",
  "mg_thymus_bulk_thymus_1.tsv.gz",
  "mg_thymus_bulk_thymus_2.tsv.gz",
  "mg_thymus_bulk_blood_month0_1.tsv.gz",
  "mg_thymus_bulk_blood_month0_2.tsv.gz",
  "mg_thymus_bulk_blood_month12_1.tsv.gz",
  "mg_thymus_bulk_blood_month12_2.tsv.gz"
  )
)
#> ℹ Downloading 'mg_thymus_slice1_spatial.rds'
#> ℹ 'mg_thymus_slice2_spatial.rds' is already in '/var/folders/5l/x6nq60_10t53vqby5_sjxj9r0000gr/T//RtmpqlRSLG/SIRUP_example'.
#> ℹ Downloading 'mg_thymus_slice1_bcr.tsv.gz'
#> ℹ 'mg_thymus_slice2_bcr.tsv.gz' is already in '/var/folders/5l/x6nq60_10t53vqby5_sjxj9r0000gr/T//RtmpqlRSLG/SIRUP_example'.
#> ℹ Downloading 'mg_thymus_slice1_gc_polygons.csv'
#> ℹ 'mg_thymus_slice2_gc_polygons.csv' is already in '/var/folders/5l/x6nq60_10t53vqby5_sjxj9r0000gr/T//RtmpqlRSLG/SIRUP_example'.
#> ℹ Downloading 'mg_thymus_bulk_thymus_1.tsv.gz'
#> ℹ Downloading 'mg_thymus_bulk_thymus_2.tsv.gz'
#> ℹ Downloading 'mg_thymus_bulk_blood_month0_1.tsv.gz'
#> ℹ Downloading 'mg_thymus_bulk_blood_month0_2.tsv.gz'
#> ℹ Downloading 'mg_thymus_bulk_blood_month12_1.tsv.gz'
#> ℹ Downloading 'mg_thymus_bulk_blood_month12_2.tsv.gz'
```

## Load the samples

``` r
library(SIRUP)
library(dplyr)
library(ggplot2)
```

``` r
spatial <- list(
  mg_thymus_slice1 = readRDS(file.path(extdata, "mg_thymus_slice1_spatial.rds")),
  mg_thymus_slice2 = readRDS(file.path(extdata, "mg_thymus_slice2_spatial.rds"))
)

read_bcrs <- function(file, sample, compartment) {
  airr::read_rearrangement(file.path(extdata, file)) %>%
    mutate(sample = sample, compartment = compartment)
}

sections <- bind_rows(
  read_bcrs("mg_thymus_slice1_bcr.tsv.gz", "mg_thymus_slice1", "thymus_slice1"),
  read_bcrs("mg_thymus_slice2_bcr.tsv.gz", "mg_thymus_slice2", "thymus_slice2")
) %>%
  extractSpatialIds(
    id = "sequence_id",
    spatial_id = "barcode_10x",
    regex = ",([ACGT]{16})$",
    suffix = "-1"
  )

bulk <- bind_rows(
  read_bcrs("mg_thymus_bulk_thymus_1.tsv.gz", "bulk_thymus_1", "bulk_thymus"),
  read_bcrs("mg_thymus_bulk_thymus_2.tsv.gz", "bulk_thymus_2", "bulk_thymus"),
  read_bcrs("mg_thymus_bulk_blood_month0_1.tsv.gz", "bulk_blood_month0_1", "blood_month0"),
  read_bcrs("mg_thymus_bulk_blood_month0_2.tsv.gz", "bulk_blood_month0_2", "blood_month0"),
  read_bcrs("mg_thymus_bulk_blood_month12_1.tsv.gz", "bulk_blood_month12_1", "blood_month12"),
  read_bcrs("mg_thymus_bulk_blood_month12_2.tsv.gz", "bulk_blood_month12_2", "blood_month12")
)

bcrs <- bind_rows(sections, bulk) %>%
  mutate(sequence_id = paste(sample, gsub(",", "_", sequence_id), sep = "_")) %>%
  filterAirr(
    locus_val = "IGH",
    min_v_length = 250,
    min_j_length = 30,
    min_junction_length = 10,
    productive_only = TRUE,
    min_consensus_count = NULL
  ) %>%
  filter(!grepl("[^ACGT]", junction))
#> filterAirr kept 5733 of 8849 rows (dropped 3116).

count(bcrs, compartment, sample)
#> # A tibble: 8 × 3
#>   compartment   sample                   n
#>   <chr>         <chr>                <int>
#> 1 blood_month0  bulk_blood_month0_1    182
#> 2 blood_month0  bulk_blood_month0_2    182
#> 3 blood_month12 bulk_blood_month12_1   241
#> 4 blood_month12 bulk_blood_month12_2   108
#> 5 bulk_thymus   bulk_thymus_1          176
#> 6 bulk_thymus   bulk_thymus_2          576
#> 7 thymus_slice1 mg_thymus_slice1      1967
#> 8 thymus_slice2 mg_thymus_slice2      2105
```

`read_bcrs()` reads one AIRR file and adds two columns: `sample` and
`compartment`. After the barcodes come out of the section IDs, the code
makes each `sequence_id` unique and safe for trees. It puts the sample
name in front, because IDs repeat between files. It replaces the comma
with `_`, because dowser renames tree tips that contain a comma.

The section rows get a `sample` that matches the names of the `spatial`
list. The bulk and blood rows have their own sample names and no
barcode.

<div markdown="1"
style="background: #F7E5BC; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #d6b676; color: #ffffff; font-weight: 700;">

<b>Note</b>
</p>

To keep the package small, the bulk files hold a subset of the bulk
data. They hold the blood sequences of clones seen in a section. They
also hold up to 10 bulk thymus sequences of each clone that reaches the
blood.

</div>

## Create a multi-sample object

Give `createSIRUP()` the named list and the name of the `sample` column.

``` r
sirup <- createSIRUP(
  spatial = spatial,
  receptors = bcrs,
  spatial_id = "barcode_10x",
  sample = "sample",
  image_scale = "lowres",
  swap_coords = FALSE
)

sirup
#> <SIRUP object>
#> Spatial data:    Seurat, 6,599 spatial IDs, 449 genes, 2 spatial samples, 2 images
#> Receptor data:   tibble, 5,537 receptor sequences x 58 columns
#> Spatial overlap: 1,567/1,640 receptor IDs in spatial data (95.5%)
#> Samples:         2
#> Spatial ID field: barcode_10x
#> Sample ID field: sample
#> Analyses:        0
spatialSamples(sirup)
#> [1] "mg_thymus_slice1" "mg_thymus_slice2"
```

Visium barcodes repeat between sections. For example,
`AAACAAGTATCTCCCA-1` exists in both. SIRUP matches each receptor to a
spot by the pair of sample and barcode, so the two sections never mix. A
spot is one capture area of the Visium slide.

If the sections need different image settings, give named values:
`image_scale = c(mg_thymus_slice1 = "hires", mg_thymus_slice2 = "lowres")`.

## Get one sample

The accessors take a `sample` argument.

``` r
getSpatial(sirup, "mg_thymus_slice1")
#> An object of class Seurat 
#> 506 features across 3396 samples within 2 assays 
#> Active assay: Spatial (449 features, 0 variable features)
#>  2 layers present: counts, data
#>  1 other assay present: alra
#>  1 spatial field of view present: thymus_slice1
nrow(getReceptors(sirup, "mg_thymus_slice1"))
#> [1] 1967
nrow(getReceptors(sirup))
#> [1] 5537
```

`filterSample()` returns a new SIRUP object with one section. It drops
the nonspatial rows.

``` r
slice1 <- filterSample(sirup, "mg_thymus_slice1")
slice1
#> <SIRUP object>
#> Spatial data:    Seurat, 3,396 spatial IDs, 449 genes, 1 image
#> Receptor data:   tibble, 1,967 receptor sequences x 58 columns
#> Spatial overlap: 770/772 receptor IDs in spatial data (99.7%)
#> Spatial ID field: barcode_10x
#> Analyses:        0
```

`subsetSIRUP()` gives finer control. It can keep a set of spots or a set
of receptor rows. With `keep_nonspatial_receptors = TRUE`, it also keeps
the bulk rows.

`mergeSIRUP()` goes the other way. It combines single-section objects
into a multi-sample object, and it adds a sample column with the names
that you give.

``` r
slice2 <- filterSample(sirup, "mg_thymus_slice2")
merged <- mergeSIRUP(
  slice1, slice2,
  samples = c("mg_thymus_slice1", "mg_thymus_slice2"),
  sample_field = "sample"
)
merged
#> <SIRUP object>
#> Spatial data:    Seurat, 6,599 spatial IDs, 449 genes, 2 spatial samples, 2 images
#> Receptor data:   tibble, 4,072 receptor sequences x 58 columns
#> Spatial overlap: 1,567/1,640 receptor IDs in spatial data (95.5%)
#> Samples:         2
#> Spatial ID field: barcode_10x
#> Sample ID field: sample
#> Analyses:        0
```

## Add polygons for each section

A polygon is a closed outline that selects the spots of one tissue
region. Give `addPolygons()` a named vector of polygon files, one for
each section.

``` r
sirup <- addPolygons(sirup, polygons_path = c(
  mg_thymus_slice1 = file.path(extdata, "mg_thymus_slice1_gc_polygons.csv"),
  mg_thymus_slice2 = file.path(extdata, "mg_thymus_slice2_gc_polygons.csv")
))

getPolygons(sirup)
#> # A tibble: 19 × 3
#>    polygon_id           sample           vertices        
#>    <chr>                <chr>            <named list>    
#>  1 mg_thymus_slice1::1  mg_thymus_slice1 <tibble [6 × 6]>
#>  2 mg_thymus_slice1::10 mg_thymus_slice1 <tibble [6 × 6]>
#>  3 mg_thymus_slice1::2  mg_thymus_slice1 <tibble [4 × 6]>
#>  4 mg_thymus_slice1::3  mg_thymus_slice1 <tibble [5 × 6]>
#>  5 mg_thymus_slice1::4  mg_thymus_slice1 <tibble [5 × 6]>
#>  6 mg_thymus_slice1::5  mg_thymus_slice1 <tibble [6 × 6]>
#>  7 mg_thymus_slice1::6  mg_thymus_slice1 <tibble [8 × 6]>
#>  8 mg_thymus_slice1::7  mg_thymus_slice1 <tibble [5 × 6]>
#>  9 mg_thymus_slice1::8  mg_thymus_slice1 <tibble [4 × 6]>
#> 10 mg_thymus_slice1::9  mg_thymus_slice1 <tibble [5 × 6]>
#> 11 mg_thymus_slice2::1  mg_thymus_slice2 <tibble [7 × 6]>
#> 12 mg_thymus_slice2::2  mg_thymus_slice2 <tibble [4 × 6]>
#> 13 mg_thymus_slice2::3  mg_thymus_slice2 <tibble [4 × 6]>
#> 14 mg_thymus_slice2::4  mg_thymus_slice2 <tibble [4 × 6]>
#> 15 mg_thymus_slice2::5  mg_thymus_slice2 <tibble [6 × 6]>
#> 16 mg_thymus_slice2::6  mg_thymus_slice2 <tibble [7 × 6]>
#> 17 mg_thymus_slice2::7  mg_thymus_slice2 <tibble [7 × 6]>
#> 18 mg_thymus_slice2::8  mg_thymus_slice2 <tibble [4 × 6]>
#> 19 mg_thymus_slice2::9  mg_thymus_slice2 <tibble [4 × 6]>
```

Both files number their eGCs from 1. SIRUP adds the sample to each
polygon ID. `mg_thymus_slice1::1` and `mg_thymus_slice2::1` are
different regions.

## Plot several sections

A spatial plot of a multi-sample object returns a list with one plot for
each section. Set `combine_sample_plots = TRUE` to get one combined
figure. Give `sample_id` to plot one section.

``` r
plotSpatialPolygons(sirup, zoom = TRUE, combine_sample_plots = TRUE)
```

![](multisample_files/figure-gfm/plots-1.jpeg)

``` r
spatialExpressionPlot(
  sirup,
  features = "AICDA",
  assay = "alra",
  layer = "data",
  sample_id = "mg_thymus_slice1",
  title = "AICDA, section 1"
)
```

![](multisample_files/figure-gfm/plot-one-1.jpeg)

## Call clones across all samples

Call clones on the whole receptor table. One clone ID then covers the
sections, the bulk thymus, and the blood.

``` r
sirup <- updateReceptors(
  sirup,
  scoper::hierarchicalClones,
  threshold = 0.1,
  summarize_clones = FALSE,
  nproc = 1
)
dplyr::n_distinct(getReceptors(sirup)$clone_id)
#> [1] 1860
```

## Find where each clone was seen

Each row has a `compartment`. For each clone, the code below counts the
sequences in each compartment. It also records the eGCs that hold the
clone.

``` r
clones <- getReceptors(sirup) %>%
  group_by(clone_id) %>%
  summarise(
    in_section = any(compartment %in% c("thymus_slice1", "thymus_slice2")),
    n_egc = sum(!is.na(polygon_id)),
    egcs = paste(sort(unique(na.omit(polygon_id))), collapse = ", "),
    n_bulk_thymus = sum(compartment == "bulk_thymus"),
    n_blood_month0 = sum(compartment == "blood_month0"),
    n_blood_month12 = sum(compartment == "blood_month12")
  ) %>%
  mutate(
    in_egc = n_egc > 0,
    in_blood = n_blood_month0 + n_blood_month12 > 0
  )

count(clones, in_egc, in_blood)
#> # A tibble: 4 × 3
#>   in_egc in_blood     n
#>   <lgl>  <lgl>    <int>
#> 1 FALSE  FALSE     1575
#> 2 FALSE  TRUE       153
#> 3 TRUE   FALSE      121
#> 4 TRUE   TRUE        11
```

The table shows that most clones stay in the thymus. A few eGC clones
also appear in the blood. The SIRUP paper calls these clones exported.
Fisher’s exact test compares the eGC clones with the other clones of the
sections.

``` r
section_clones <- filter(clones, in_section)
fisher.test(table(section_clones$in_egc, section_clones$in_blood))
#> 
#>  Fisher's Exact Test for Count Data
#> 
#> data:  table(section_clones$in_egc, section_clones$in_blood)
#> p-value = 0.3623
#> alternative hypothesis: true odds ratio is not equal to 1
#> 95 percent confidence interval:
#>  0.625789 2.559146
#> sample estimates:
#> odds ratio 
#>   1.326868
```

For each eGC, count the clones that it holds and the clones that reach
the blood.

``` r
getReceptors(sirup) %>%
  filter(!is.na(polygon_id)) %>%
  distinct(polygon_id, clone_id) %>%
  left_join(select(clones, clone_id, in_blood), by = "clone_id") %>%
  group_by(polygon_id) %>%
  summarise(n_clones = n(), n_exported = sum(in_blood)) %>%
  arrange(desc(n_exported))
#> # A tibble: 17 × 3
#>    polygon_id           n_clones n_exported
#>    <chr>                   <int>      <int>
#>  1 mg_thymus_slice2::5        23          5
#>  2 mg_thymus_slice2::1        16          3
#>  3 mg_thymus_slice1::10       15          1
#>  4 mg_thymus_slice1::2         2          1
#>  5 mg_thymus_slice2::2        16          1
#>  6 mg_thymus_slice2::3         7          1
#>  7 mg_thymus_slice2::6         4          1
#>  8 mg_thymus_slice1::1         6          0
#>  9 mg_thymus_slice1::3         2          0
#> 10 mg_thymus_slice1::4         1          0
#> 11 mg_thymus_slice1::5        10          0
#> 12 mg_thymus_slice1::6        18          0
#> 13 mg_thymus_slice1::7         2          0
#> 14 mg_thymus_slice1::9         5          0
#> 15 mg_thymus_slice2::7        12          0
#> 16 mg_thymus_slice2::8         4          0
#> 17 mg_thymus_slice2::9         4          0
```

The bulk files hold a subset of the bulk data, and the clones come from
that subset. For the full analysis, see Figure 5 of the SIRUP paper.

## Measure clonal dominance in each section

Clonal dominance is the share of the receptors in a region that belong
to its largest clones. The regional functions work on all sections at
once. The dominance table has a `sample` column and the full
`polygon_id`. For a longer treatment, see
[Draw tissue regions and measure clonal dominance and diversity](regions.md).

``` r
dominance <- spatialCloneDominance(object = sirup, locus_val = "IGH")

dominance$dominance %>%
  group_by(polygon_id) %>%
  slice_max(proportion, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sample, polygon_id, clone_id, n_receptors, proportion) %>%
  head()
#> # A tibble: 6 × 5
#>   sample           polygon_id           clone_id n_receptors proportion
#>   <chr>            <chr>                <chr>          <dbl>      <dbl>
#> 1 mg_thymus_slice1 mg_thymus_slice1::1  1270               2      0.222
#> 2 mg_thymus_slice1 mg_thymus_slice1::10 1148               4      0.16 
#> 3 mg_thymus_slice1 mg_thymus_slice1::2  1661               1      0.5  
#> 4 mg_thymus_slice1 mg_thymus_slice1::3  1267               2      0.5  
#> 5 mg_thymus_slice1 mg_thymus_slice1::4  1399               1      1    
#> 6 mg_thymus_slice1 mg_thymus_slice1::5  1385               3      0.231
```

Most plotting functions return a named list with one plot for each
section.

``` r
region_plots <- plotSpatialRegionReceptors(
  dominance = dominance,
  object = sirup,
  top_n = 8,
  point_size = 1.5,
  zoom = TRUE
)
names(region_plots)
#> [1] "mg_thymus_slice1" "mg_thymus_slice2"
region_plots$mg_thymus_slice1
```

![](multisample_files/figure-gfm/region-receptors-1.jpeg)

`plotSpatialPolygonPies()` also works on all sections at once. Give
`sample_id` to draw the pies of one section next to its map.

``` r
pies <- plotSpatialPolygonPies(
  dominance = dominance,
  object = sirup,
  fill = "clone_id",
  top_n = 5,
  sample_id = "mg_thymus_slice1"
)
pies$page_1
#> Warning: Removed 10 rows containing missing values or values outside the scale range
#> (`geom_text()`).
```

![](multisample_files/figure-gfm/pies-1.jpeg)

## Build a tree with blood tips

A lineage tree shows how the sequences of one clone relate to each
other. Each tip of the tree is one sequence. A tree can mix spatial and
bulk sequences. The code below picks an exported eGC clone with blood
sequences from both time points.

``` r
exported <- clones %>%
  filter(in_egc, n_blood_month0 > 0, n_blood_month12 > 0) %>%
  mutate(n_total = n_egc + n_bulk_thymus + n_blood_month0 + n_blood_month12) %>%
  filter(n_total <= 100) %>%
  arrange(desc(n_blood_month0 + n_blood_month12))
exported
#> # A tibble: 4 × 10
#>   clone_id in_section n_egc egcs    n_bulk_thymus n_blood_month0 n_blood_month12
#>   <chr>    <lgl>      <int> <chr>           <int>          <int>           <int>
#> 1 1045     TRUE           2 mg_thy…            10             15               4
#> 2 1673     TRUE           3 mg_thy…             0              2               1
#> 3 868      TRUE           1 mg_thy…             0              1               2
#> 4 85       TRUE           1 mg_thy…            10              1               1
#> # ℹ 3 more variables: in_egc <lgl>, in_blood <lgl>, n_total <int>

tree_clone <- exported$clone_id[[1]]
tree_section <- sub("::.*$", "", exported$egcs[[1]])
tree_section
#> [1] "mg_thymus_slice2"
```

Make the germline of the clone with `dowser::createGermlines()`. The
germline is the sequence of the unmutated ancestor, and the tree starts
there. Then build the tree. As in
[Call clones and build lineage trees on the tissue](clones-and-trees.md), the germline step
needs the IMGT reference sequences. Set the environment variable
`SIRUP_IMGT_DIR` to the folder that holds them. The tree chunks run only
when that folder exists.

Keep `compartment` in `traits`. The tree plot uses it to set the tip
shapes.

``` r
tree_sirup <- updateReceptors(sirup, dplyr::filter, clone_id == tree_clone)
tree_sirup <- updateReceptors(
  tree_sirup,
  dowser::createGermlines,
  references = dowser::readIMGT(imgt_dir),
  id = "sequence_id"
)
tree_sirup <- formatSpatialClones(
  tree_sirup,
  id = "sequence_id",
  traits = c("barcode_10x", "sample", "compartment", "polygon_id")
)
set.seed(1)
tree_sirup <- getSpatialTrees(tree_sirup)
trees <- getAnalysis(tree_sirup, "trees")
```

`plotSpatialCloneTreePairs()` puts the spatial tips on their sections.
Bulk and blood tips appear in the tree only. The argument `shape_by`
sets one shape for each compartment. The argument `color_by` colors the
tips by tree cluster. `clusterTreesTips()` makes the clusters.

``` r
trees <- clusterTreesTips(trees, k = "auto", id = "sequence_id", k_range = 1:6)

pair <- plotSpatialCloneTreePairs(
  trees,
  object = sirup,
  id = "sequence_id",
  color_by = "tree_cluster",
  shape_by = "compartment",
  shape_values = c(
    thymus_slice1 = 21, thymus_slice2 = 21, bulk_thymus = 22,
    blood_month0 = 24, blood_month12 = 23
  ),
  zoom = TRUE,
  polygons = TRUE,
  polygon_alpha = 0,
  polygon_color = "black",
  show_spatial_legend = TRUE
)

names(pair)
```

The list has one element for each clone. If the tips of a clone fall in
more than one section, the list has one element for each clone and
section instead. These elements are named `<clone>::<sample>`. Each one
pairs the full tree with the map of one section. Each element also holds
a `legend` for the tip colors. The function builds it because the call
sets `show_spatial_legend = TRUE`. The code below shows the section that
holds the eGC tips.

``` r
key <- paste0(tree_clone, "::", tree_section)
if (!key %in% names(pair)) key <- tree_clone
section_pair <- pair[[key]]

cowplot::plot_grid(
  cowplot::plot_grid(section_pair$tree, section_pair$spatial_clone, rel_widths = c(1.3, 1)),
  section_pair$legend,
  ncol = 1,
  rel_heights = c(1, section_pair$legend_rel_height)
)
```

Shapes show the compartment of each tip. Circles are tips from a
section. Squares are bulk thymus. Triangles are blood at the thymectomy.
Diamonds are blood 12 months later. The legend shows the colors, which
are the tree clusters from `clusterTreesTips()`. On the map, a dot has
the color of its tip.

## Save the object

`saveSIRUP()` writes the whole object to one file. `readSIRUP()` reads
it back and checks it.

``` r
path <- file.path(tempdir(), "mg_thymus.sirup.rds")
saveSIRUP(sirup, path)
readSIRUP(path)
#> <SIRUP object>
#> Spatial data:    Seurat, 6,599 spatial IDs, 449 genes, 2 spatial samples, 2 images
#> Receptor data:   tibble, 5,537 receptor sequences x 79 columns
#> Spatial overlap: 1,567/1,640 receptor IDs in spatial data (95.5%)
#> Samples:         2
#> Spatial ID field: barcode_10x
#> Sample ID field: sample
#> Clone ID field:  clone_id
#> Analyses:        0
#> Polygons:        19
```

## Next steps

- Call clones and build trees for one section:
  [Call clones and build lineage trees on the tissue](clones-and-trees.md).
- Draw tissue regions and measure clonal dominance and diversity:
  [Draw tissue regions and measure clonal dominance and diversity](regions.md).
- Learn how SIRUP stores samples, polygons, and analyses:
  [How the SIRUP object works](sirup-object.md).
