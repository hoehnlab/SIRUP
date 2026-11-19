How the SIRUP object works
================

<div markdown="1"
style="background: #F7E5BC; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #d6b676; color: #ffffff; font-weight: 700;">

<b>Note</b>
</p>

**Time:** ~10 minutes.<br> **Data:** a small synthetic SIRUP object
included in SIRUP.

</div>

Here, we will explain the parts of a SIRUP object, the rules that SIRUP
checks, and the functions that change the object.

A spatial receptor experiment gives us three kinds of data: gene
expression for each spot, receptor sequences, and regions drawn on the
tissue. These data must agree with each other. For example, a receptor
must belong to a spot that exists, and a clone column must exist in the
receptor table. The SIRUP object keeps all of the data in one place and
checks these rules each time the object changes.

<div markdown="1"
style="background: #F7E5BC; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #d6b676; color: #ffffff; font-weight: 700;">

<b>Note</b>
</p>

This page has little code. To build an object from real data, see
[Get started with SIRUP](get-started.md) first.

</div>

## Load the example object

The examples use a small synthetic object with 25 spots, 10 genes, 5
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

## Look at the five parts of the object

A SIRUP object is a list with five parts, in this order:

| Part | Contents | How to get it |
|----|----|----|
| `spatial` | The Seurat object, with expression and the H&E image | `getSpatial(x)` |
| `receptors` | The AIRR receptor table, as a tibble, or a named list of tables | `getReceptors(x)` |
| `metadata` | The field names and the image settings | `x$metadata` |
| `polygons` | The tissue regions, with one row for each polygon | `getPolygons(x)` |
| `analyses` | Results that are saved in the object | `getAnalysis(x, name)` |

The accessors return ordinary R objects, so they can be passed to any
function from another package.

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

## Set the field names

A receptor table can use any column names, so SIRUP needs to know which
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

The fields are set once, in `createSIRUP()`. After that, SIRUP functions
find each column without extra arguments. To change a field later, give
a new map to `updateReceptors()` with the `receptor_fields` argument. A
value of `NULL` clears an optional field.

## Validate the object

`validateSIRUP()` checks that the parts agree. It gives an error if one
of these rules fails:

1.  The object has class `sirup`.
2.  At least one receptor spatial ID matches a spot name in the Seurat
    object.
3.  Each column that `receptor_fields` names exists in the receptor
    table. If the object has several receptor tables, the spatial ID
    column must be in every table. Each other named column, such as
    `clone_id`, must be in at least one table.
4.  If the sample field is set, the column exists in every receptor
    table and in the Seurat metadata.

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
change that breaks a rule fails right away and leaves the original
object as it was. This example gives every receptor a spot name that
does not exist:

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

<div markdown="1"
style="background: #F7E5BC; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #d6b676; color: #ffffff; font-weight: 700;">

<b>Note</b>
</p>

SIRUP keeps receptor rows that match no spot, such as bulk and blood
sequences. The rules need only one match, so these rows stay in the
table.

</div>

## Set the image settings

Two settings tell SIRUP how to put spots on the H&E image. Both are
stored in `x$metadata`:

- `image_scale` is the Seurat image scale, `"lowres"` or `"hires"`.
- `swap_coords` swaps the x and y coordinates. Some older Seurat objects
  store the coordinates in the opposite order to the image.

Plots and spatial statistics read these settings when their own
`image_scale` and `swap_coords` arguments are `NULL`, so they only need
to be set once.

## Change the data in the object

SIRUP does not wrap every AIRR tool. Instead, three functions run any
function on one part of the object and put the result back:

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

## Add regions to the object

A region is a polygon on the tissue image. `addPolygons()` stores the
polygons and then marks the polygon of each spot and each receptor. It
adds a `polygon_id` column and one logical `polygon_<ID>` column for
each polygon to both the Seurat metadata and the receptor table.

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

## Hold several receptor tables

One object can hold more than one receptor table, e.g., one table of
BCRs and one of TCRs. Give `createSIRUP()` a named list of tables. Each
table is a *slot*. All slots share the spatial data, the polygons, and
the field names. The spatial ID column must be in every table, but a
column such as `clone_id` can be in one table only, e.g., the BCR table.
Clone IDs only need to be unique within one slot.

``` r
bcr <- getReceptors(x)
tcr <- bcr[1:8, ]
tcr$locus <- "TRB"

both <- createSIRUP(
  spatial = getSpatial(x),
  receptors = list(bcr = bcr, tcr = tcr),
  spatial_id = "spatial_id",
  clone = "clone_id",
  image_scale = "hires",
  swap_coords = FALSE,
  polygons = getPolygons(x)
)
receptorSlots(both)
#> [1] "bcr" "tcr"
both
#> <SIRUP object>
#> Spatial data:    Seurat, 25 spatial IDs, 10 genes, 0 images
#> Receptor data:   2 slots
#>   bcr: 36 receptor sequences x 64 columns; 19/19 receptor IDs in spatial data (100.0%)
#>   tcr: 8 receptor sequences x 64 columns; 8/8 receptor IDs in spatial data (100.0%)
#> Spatial ID field: spatial_id
#> Clone ID field:  clone_id
#> Analyses:        0
#> Polygons:        2
```

Functions that use receptors work on one slot at a time, chosen with the
`slot` argument. If the object has several slots, `slot` is required. If
the object has one table, leave it out. That table has the slot name
`"only_slot"`.

``` r
head(getReceptors(both, slot = "tcr"), 3)
#> # A tibble: 3 × 64
#>   sequence_id spatial_id clone_id locus sequence      rev_comp productive v_call
#>   <chr>       <chr>      <chr>    <chr> <chr>         <lgl>    <lgl>      <chr> 
#> 1 seq01       spot21     1        TRB   CAGGTGCATCTG… FALSE    TRUE       IGHV1…
#> 2 seq02       spot19     1        TRB   CAGGTGCATCTG… FALSE    TRUE       IGHV1…
#> 3 seq03       spot16     1        TRB   CAGGTGCATCTG… FALSE    TRUE       IGHV1…
#> # ℹ 56 more variables: d_call <chr>, j_call <chr>, sequence_alignment <chr>,
#> #   germline_alignment <chr>, junction <chr>, junction_aa <chr>, v_cigar <chr>,
#> #   d_cigar <chr>, j_cigar <chr>, stop_codon <lgl>, vj_in_frame <lgl>,
#> #   c_call <chr>, junction_length <int>, np1_length <int>, np2_length <int>,
#> #   v_sequence_start <int>, v_sequence_end <int>, v_germline_start <int>,
#> #   v_germline_end <int>, d_sequence_start <int>, d_sequence_end <int>,
#> #   d_germline_start <int>, d_germline_end <int>, j_sequence_start <int>, …
both <- polygonReceptorSummary(both, slot = "tcr")
listAnalyses(both)
#> [1] "polygon_receptor_summary_tcr"
tryCatch(getReceptors(both), error = function(e) conditionMessage(e))
#> [1] "`slot` is required because this SIRUP object has several receptor slots.\nℹ Slots: bcr, tcr."
```

The polygons apply to every slot. A receptor is in a polygon if it is in
the same sample and its spatial ID lies inside the polygon. The default
analysis names end with the slot name, so results for different slots do
not overwrite each other.

## Hold several samples

`createSIRUP()` accepts one Seurat object or a named list of Seurat
objects, e.g., one for each tissue section. With a list, the column of
the receptor table that `sample_id` names matches each row to a list
entry. Rows with no matching entry stay in the table as non-spatial
records, which is how bulk and blood sequences are stored. Spot barcodes
can repeat between sections, so SIRUP matches each receptor to a spot by
the pair of sample and barcode.

``` r
receptors <- dplyr::bind_rows(
  dplyr::mutate(getReceptors(x), sample = "a"),
  dplyr::mutate(getReceptors(x), sample = "b")
)

two_samples <- createSIRUP(
  spatial = list(a = getSpatial(x), b = getSpatial(x)),
  receptors = receptors,
  spatial_id = "spatial_id",
  clone = "clone_id",
  sample_id = "sample",
  image_scale = "hires",
  swap_coords = FALSE
)
two_samples
#> <SIRUP object>
#> Spatial data:    Seurat, 50 spatial IDs, 10 genes, 2 spatial samples, 0 images
#> Receptor data:   tibble, 72 receptor sequences x 65 columns
#> Spatial overlap: 38/38 receptor IDs in spatial data (100.0%)
#> Samples:         2
#> Spatial ID field: spatial_id
#> Sample ID field: sample
#> Clone ID field:  clone_id
#> Analyses:        0
spatialSamples(two_samples)
#> [1] "a" "b"
```

If the sections need different image settings, give named values, e.g.,
`image_scale = c(a = "hires", b = "lowres")`.

`mergeSIRUP()` goes the other way. It joins several single-sample
objects into one, with the sample names that are given:

``` r
merged <- mergeSIRUP(x, x, samples = c("a", "b"), sample_field = "sample")
spatialSamples(merged)
#> [1] "a" "b"
```

`subsetSIRUP()` keeps chosen samples, spots, or receptor rows.
`getSpatial()` and `getReceptors()` return the data of one sample when
`sample_id` is given. [Analyze multiple spatial sections bulk BCR-seq](multisample.md) has a
full example.

## Save results in the object

`setAnalysis()` stores any R object under a name, e.g., clones, trees,
spatial statistics, and polygon summaries. Stored results stay with the
data that made them, so separate files do not need to be tracked.

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

## Save and read a SIRUP object

`saveSIRUP()` validates the object and writes it to one RDS file. The
file also records the versions of SIRUP and SeuratObject that wrote it.
`readSIRUP()` reads the file and does four things. It converts an object
saved under the old package name (SIRAP), upgrades an older object to
the current layout, warns if the package versions differ from the file,
and validates the object.

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

- To build an object from your own data, see
  [Get started with SIRUP](get-started.md).
- To work with regions, see [Draw tissue regions and measure clonal dominance and diversity](regions.md).
- To combine sections, bulk data, and blood, see
  [Analyze multiple spatial sections bulk BCR-seq](multisample.md).
