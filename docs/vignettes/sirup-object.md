How the SIRUP object works
================

> **Time:** about 10 minutes. **Data:** the small synthetic object in
> SIRUP. **Needs:** SIRUP.

A spatial receptor experiment gives you three kinds of data: gene
expression for each spot, receptor sequences, and the regions that you
draw on the tissue. These data must agree with each other. A receptor
must belong to a spot that exists. A clone column must exist in the
receptor table. The SIRUP object keeps all of the data in one place and
checks these rules each time you change it.

This page explains the parts of the object, the rules that SIRUP checks,
and the functions that change the object. It has little code. For a
hands-on start, read [Get started with SIRUP](get-started.md).

The examples use a small synthetic object. It has 25 spots, 10 genes, 5
clones, and 2 polygons. The spots and gene counts are simulated. The
receptor rows are real BCR sequences from a human thymus section,
assigned to synthetic spots.

``` r
library(SIRUP)

x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x
#> <SIRUP object>
#> Spatial data:    Seurat, 25 spatial IDs, 10 genes, 0 images
#> Receptor data:   tibble, 36 receptor sequences x 64 columns
#> Spatial overlap: 19/19 receptor IDs in spatial data (100.0%)
#> Spatial ID field: spatial_id
#> Clone ID field:  clone_id
#> Analyses:        0
#> Polygons:        2
```

## The five parts

A SIRUP object is a list with five parts, in this order:

| Part | Contents | How to get it |
|----|----|----|
| `spatial` | The Seurat object, with expression and the H&E image | `getSpatial(x)` |
| `receptors` | The AIRR receptor table, as a tibble | `getReceptors(x)` |
| `metadata` | The field names and the image settings | `x$metadata` |
| `polygons` | The tissue regions, with one row for each polygon | `getPolygons(x)` |
| `analyses` | Results that you choose to keep | `getAnalysis(x, name)` |

The accessors return ordinary R objects. You can pass them to any
function from another package:

``` r
getSpatial(x)
#> An object of class Seurat 
#> 10 features across 25 samples within 1 assay 
#> Active assay: Spatial (10 features, 0 variable features)
#>  1 layer present: counts
head(getReceptors(x), 3)
#> # A tibble: 3 × 64
#>   sequence_id spatial_id clone_id locus sequence      rev_comp productive v_call
#>   <chr>       <chr>      <chr>    <chr> <chr>         <lgl>    <lgl>      <chr> 
#> 1 seq01       spot21     1        IGH   CAGGTGCATCTG… FALSE    TRUE       IGHV1…
#> 2 seq02       spot19     1        IGH   CAGGTGCATCTG… FALSE    TRUE       IGHV1…
#> 3 seq03       spot16     1        IGH   CAGGTGCATCTG… FALSE    TRUE       IGHV1…
#> # ℹ 56 more variables: d_call <chr>, j_call <chr>, sequence_alignment <chr>,
#> #   germline_alignment <chr>, junction <chr>, junction_aa <chr>, v_cigar <chr>,
#> #   d_cigar <chr>, j_cigar <chr>, stop_codon <lgl>, vj_in_frame <lgl>,
#> #   c_call <chr>, junction_length <int>, np1_length <int>, np2_length <int>,
#> #   v_sequence_start <int>, v_sequence_end <int>, v_germline_start <int>,
#> #   v_germline_end <int>, d_sequence_start <int>, d_sequence_end <int>,
#> #   d_germline_start <int>, d_germline_end <int>, j_sequence_start <int>, …
getPolygons(x)
#> # A tibble: 2 × 2
#>   polygon_id  vertices        
#>   <chr>       <list>          
#> 1 upper_right <tibble [5 × 2]>
#> 2 lower_left  <tibble [5 × 2]>
```

`summary()` gives counts for each part:

``` r
summary(x)
#> SIRUP summary
#>   Spatial IDs:      25
#>   Receptor records: 36
#>   Receptor IDs:     19
#>   Matched IDs:      19 (100.0%)
#>   Matched records:  36
#>   Image scale:      hires
#>   Swap coords:      FALSE
#>   Stored analyses:  0
#>   Polygons:         2
#> 
#> Per-sample summary
#> # A tibble: 1 × 7
#>   sample  n_spatial_ids n_receptor_records n_receptor_ids n_matched_ids
#>   <chr>           <int>              <int>          <int>         <int>
#> 1 spatial            25                 36             19            19
#> # ℹ 2 more variables: n_matched_records <int>, id_match_rate <dbl>
```

## Field names

Your receptor table can use any column names. SIRUP needs to know which
column plays which role. It stores this map in
`x$metadata$receptor_fields`. There are six fields:

| Field | The receptor column that holds… |
|----|----|
| `spatial_id` | the name of the spot that holds the receptor |
| `image` | the name of the Seurat image, if the object has several images |
| `sample` | the sample name, in a multi-sample object |
| `clone` | the clone ID |
| `locus` | the receptor locus, such as `IGH` or `TRB` |
| `count` | the number of reads or UMIs of the sequence |

``` r
x$metadata$receptor_fields
#> $spatial_id
#> [1] "spatial_id"
#> 
#> $image
#> NULL
#> 
#> $sample
#> NULL
#> 
#> $clone
#> [1] "clone_id"
#> 
#> $locus
#> [1] "locus"
#> 
#> $count
#> NULL
```

You set the fields once, in `createSIRUP()`. After that, SIRUP functions
find each column without arguments. To change a field later, give a new
map to `updateReceptors()` with the `receptor_fields` argument. A value
of `NULL` clears an optional field.

## Rules that SIRUP checks

`validateSIRUP()` checks that the parts agree. It gives an error if one
of these rules fails:

1.  The object has class `sirup`.
2.  At least one receptor spatial ID matches a spot name in the Seurat
    object.
3.  Each column that `receptor_fields` names exists in the receptor
    table.
4.  If the `sample` field is set, the column exists in both the receptor
    table and the Seurat metadata.

``` r
validateSIRUP(x, verbose = TRUE)
#> <SIRUP object>
#> Spatial data:    Seurat, 25 spatial IDs, 10 genes, 0 images
#> Receptor data:   tibble, 36 receptor sequences x 64 columns
#> Spatial overlap: 19/19 receptor IDs in spatial data (100.0%)
#> Spatial ID field: spatial_id
#> Clone ID field:  clone_id
#> Analyses:        0
#> Polygons:        2
```

`createSIRUP()`, `updateReceptors()`, `updateSpatial()`, and the other
functions that change the object run this check before they return. A
change that breaks a rule fails at once and leaves your object as it
was. This example gives every receptor a spot name that does not exist:

``` r
tryCatch(
  updateReceptors(x, function(receptors) {
    receptors$spatial_id <- "no-such-spot"
    receptors
  }),
  error = function(e) conditionMessage(e)
)
#> [1] "No receptor IDs match the spatial Seurat object."
```

> **Note.** SIRUP keeps receptor rows that match no spot. Bulk and blood
> sequences are examples. The rules need only one match, so these rows
> stay in the table.

## Image settings

Two settings tell SIRUP how to put spots on the H&E image. Both live in
`x$metadata`:

- `image_scale` is the Seurat image scale, `"lowres"` or `"hires"`.
- `swap_coords` says whether to swap the x and y coordinates. Some older
  Seurat objects store the coordinates in the opposite order to the
  image.

Plots and spatial statistics read these settings when you leave their
own `image_scale` and `swap_coords` arguments at `NULL`. You set them
once.

## Change the data

SIRUP does not wrap every AIRR tool. Three functions run any function on
one part of the object and put the result back:

- `updateReceptors(x, fun, ...)` calls `fun(getReceptors(x), ...)` and
  stores the table that `fun` returns.
- `updateSpatial(x, fun, ...)` calls `fun(getSpatial(x), ...)` and
  stores the Seurat object that `fun` returns.
- `updateAnalysis(x, name, fun, ...)` does the same for a stored
  analysis.

Each function validates the object again before it returns. This call
removes the sequences of one clone:

``` r
x2 <- updateReceptors(x, dplyr::filter, clone_id != "3")
c(before = nrow(getReceptors(x)), after = nrow(getReceptors(x2)))
#> before  after 
#>     36     29
```

This call adds a column to the Seurat metadata:

``` r
x3 <- updateSpatial(x, function(spatial) {
  spatial$gene_rich <- spatial$nFeature_Spatial >= 8
  spatial
})
table(getSpatial(x3)$gene_rich)
#> 
#> FALSE  TRUE 
#>     2    23
```

The same pattern runs germline reconstruction
(`dowser::createGermlines()`) or mutation counting
(`shazam::observedMutations()`) inside the object.
[Call clones and build lineage trees on the tissue](clones-and-trees.md) shows these steps.

## Regions in the object

A region is a polygon on the tissue image. `addPolygons()` stores the
polygons and then marks the polygon of each spot and each receptor. It
adds a `polygon_id` column and one logical `polygon_<ID>` column for
each polygon to both the Seurat metadata and the receptor table:

``` r
grep("^polygon", names(getReceptors(x)), value = TRUE)
#> [1] "polygon_id"          "polygon_upper_right" "polygon_lower_left"
grep("^polygon", names(getSpatial(x)[[]]), value = TRUE)
#> [1] "polygon_id"          "polygon_upper_right" "polygon_lower_left"
```

A spot on the edge of a polygon counts as a member when
`include_boundary` is `TRUE`, which is the default.
[Draw tissue regions and measure clonal dominance and diversity](regions.md) shows how to draw regions and
measure what they hold.

## Several samples

`createSIRUP()` accepts one Seurat object or a named list of Seurat
objects. With a list, the `sample` column of the receptor table names
the list entry of each row. Rows with no matching entry stay in the
table as non-spatial records. Bulk and blood sequences use this feature.
`mergeSIRUP()` joins several objects into one. `subsetSIRUP()` keeps
chosen samples, spots, or receptor rows. `getSpatial()` and
`getReceptors()` return the data of one sample when you give `sample`.
[Analyze several sections with bulk AIRR-seq and blood](multisample.md) has a full example.

## Keep results with the data

`setAnalysis()` stores any R object under a name. Use it for clones,
trees, spatial statistics, and polygon summaries. Stored results stay
with the data that made them, so you do not track separate files:

``` r
x <- setAnalysis(x, "n_clones", 5)
listAnalyses(x)
#> [1] "n_clones"
getAnalysis(x, "n_clones")
#> [1] 5
x <- removeAnalysis(x, "n_clones")
listAnalyses(x)
#> character(0)
```

`saveSIRUP()` validates the object and writes it to one RDS file. The
file also records the versions of SIRUP and SeuratObject that wrote it.
`readSIRUP()` reads the file and does four things. It converts an object
saved under the old package name (SIRAP). It upgrades an older object to
the current layout. It warns you if the package versions differ from the
file. It validates the object.

``` r
path <- saveSIRUP(x, tempfile(fileext = ".sirup.rds"))
x_again <- readSIRUP(path)
x_again
#> <SIRUP object>
#> Spatial data:    Seurat, 25 spatial IDs, 10 genes, 0 images
#> Receptor data:   tibble, 36 receptor sequences x 64 columns
#> Spatial overlap: 19/19 receptor IDs in spatial data (100.0%)
#> Spatial ID field: spatial_id
#> Clone ID field:  clone_id
#> Analyses:        0
#> Polygons:        2
```

## Next steps

- To build an object from your own data, follow
  [Get started with SIRUP](get-started.md).
- To work with regions, read [Draw tissue regions and measure clonal dominance and diversity](regions.md).
- To combine sections, bulk data, and blood, read
  [Analyze several sections with bulk AIRR-seq and blood](multisample.md).
