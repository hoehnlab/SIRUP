How the SIRUP object works
================

<div markdown="1"
style="background: #F7E5BC; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #d6b676; color: #ffffff; font-weight: 700;">

<b>Note</b>
</p>

**Time:** ~10 minutes.<br> **Data:** A small synthetic SIRUP object
included in the package.

</div>

Here, we will explain the parts of a SIRUP object, the rules that govern
it, and the functions that change the object.

A spatial AIRR-seq (spAIRR) experiment gives us two types of data: gene
expression and adaptive immune receptor sequences measured in spatial
barcodes that map to locations on the tissue. These data need to agree
with each other and be able to interact. The **SIRUP object** keeps the
data in one place and checks these rules.

Because this data is spatially resolved, we are often interested in
studying particular regions of interest (ROIs) on the tissue. The SIRUP
object can also hold such ROIs, which are defined by polygons drawn on
the tissue image.

<div markdown="1"
style="background: #F7E5BC; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #d6b676; color: #ffffff; font-weight: 700;">

<b>Note</b>
</p>

This page has very little code. To build an object from real data, see
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

## Five parts of the SIRUP object

A SIRUP object is an S3 object, and in practice it is a list with five
parts.

| Name | Contents | How to access |
|----|----|----|
| `spatial` | The spatial Seurat object, including gene expression and the H&E image (or named list of spatial Seurat objects) | `getSpatial(x)` |
| `receptors` | The AIRR receptor table, as a tibble (or a named list of tables) | `getReceptors(x)` |
| `metadata` | The important field names (clone_id, sample_id, locus, etc.) and the image settings | `x$metadata` |
| `polygons` | The ROIs that we create, with one row for each polygon | `getPolygons(x)` |
| `analyses` | Analysis results that you want to keep with the object | `getAnalysis(x, name)` |

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

`summary()` describes each part of the object.

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
column has what information. It stores this in
`x$metadata$receptor_fields`. There are six key fields:

| Field | Contents |
|----|----|
| `spatial_id` | identifier of the spatial spot that the receptor is in |
| `image` | the name of the Seurat image if the object has multiple images |
| `sample` | the sample name if the object has multiple samples (you would’ve passed a named list to `spatial` during `createSIRUP()`) |
| `clone` | the clone ID |
| `locus` | the receptor locus, such as `IGH` or `TRB` |
| `count` | the number of reads or UMIs supporting the sequence |

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
can find each column without extra arguments. To change a field later,
give a new map to `updateReceptors()` with the `receptor_fields`
argument. A value of `NULL` clears a field.

## Validate the object

`validateSIRUP()` checks that the parts agree. It gives an error if one
of the below rules fails.

1.  The object has class `sirup`.
2.  At least one receptor spatial ID matches a spot name in the Seurat
    object.
3.  Each column that the metadata `receptor_fields` names exists in the
    receptor table. If the object has several receptor tables, the
    spatial ID column must be in every table. Each other named column,
    such as `clone_id`, must be in at least one table.
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
object as it was.

<div markdown="1"
style="background: #F7E5BC; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #d6b676; color: #ffffff; font-weight: 700;">

<b>Note</b>
</p>

SIRUP keeps receptor rows that match no spot, because they could have
come from another experiment. The validation check only requires that at
least one receptor row matches a spot.

</div>

## Set the image settings

Two settings tell SIRUP how to plot spots on the H&E image. Both are
stored in `x$metadata`:

- `image_scale` is the Seurat image scale, `"lowres"` or `"hires"`.
  Choose based on which image you used to build the spatial Seurat
  object when you called `Seurat::Load10X_Spatial()`
- `swap_coords` swaps the x and y coordinates. Some older Seurat objects
  store the coordinates in the opposite order to the image. You can
  check this by grabbing the tissue coordinates with
  `Seurat::GetTissueCoordinates(getSpatial(x))` and comparing them to
  the image (see [Get started with SIRUP](get-started.md) for an
  example). If the tissue is rotated 90 degrees, set
  `swap_coords = TRUE`.

Plotting and spatial statistics functions read these settings when their
own `image_scale` and `swap_coords` arguments are not provided or are
`NULL`, so they only need to be set when you create the SIRUP object.

## Change the data in the object

Because SIRUP cannot possibly wrap every potential tool, it provides a
few functions that let you run other functions on the data inside the
object and update it. These three functions are:

- `updateReceptors(x, fun, ...)` calls `fun(getReceptors(x), ...)` and
  stores the table that `fun` returns.
- `updateSpatial(x, fun, ...)` calls `fun(getSpatial(x), ...)` and
  stores the Seurat object that `fun` returns.
- `updateAnalysis(x, name, fun, ...)` does the same for a stored
  analysis.

Each function validates the object again before it returns. This code
filters out the sequences in clone 3.

``` r
x2 <- updateReceptors(x, dplyr::filter, clone_id != "3")
c(before = nrow(getReceptors(x)), after = nrow(getReceptors(x2)))
#> before  after 
#>     36     29
```

This code adds a column to the Seurat metadata.

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

Other examples are shown throughout the vignettes. In general, if you
want to update any data in a SIRUP object, you can use one of these
three functions to run your own code on the data and store the result.

## Add regions to the object

See [Annotate tissue ROIs and analyze local clonal structure](regions.md) to learn how to annotate
ROIs and attach them to the SIRUP object.

## Save analyses in the object

`setAnalysis()` stores any R object under a name, e.g., clones, trees,
spatial statistics, and polygon summaries. Stored results then stay with
the SIRUP object for convenience.

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

`saveSIRUP()` validates the object and writes it to a .RDS file. The
file also records the version of SIRUP that wrote it. `readSIRUP()`
reads the file and converts an object saved under the old package name
(SIRAP), upgrades an older object to the current layout, warns if the
package versions differ from the file, and validates the object.

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

## Storing multiple sets of receptors or spatial data in a SIRUP object

### Multiple receptor tables but one spatial sample

One object can hold more than one receptor table, e.g., one table of
BCRs and one of TCRs. Give `createSIRUP()` a named list of AIRR tables.
Each table is a *slot*. All slots share the spatial data, the polygons,
and the field names. The spatial ID column must be in every table, but a
column such as `clone_id` can be in one table only, e.g., the BCR table.
Clone IDs only need to be unique within one slot. See the example below.

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
the object has only one table, you can leave it out.

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
```

### Multiple spatial samples but one receptor table

`createSIRUP()` accepts one spatial sample (i.e., a spatial Seurat
object) or a named list of spatial Seurat objects, e.g., one for each
tissue section. If you pass a list, the `sample_id` column of the
receptor table matches each receptor row to its corresponding spatial
data, so be sure the names in the list match the values in the receptor
table. Rows with no matching entry stay in the table as non-spatial
records, which is how non-spatial receptors are stored. Spot barcodes
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

### Multiple receptor tables and multiple spatial samples

`createSIRUP()` can accept a named list of receptor tables and a named
list of spatial Seurat objects. Each receptor table just needs a sample
column that matches the names of the spatial Seurat objects. Such
objects follow all of the same rules and conventions as discussed in the
previous sections. Functions that use receptors or spatial data work on
one slot or sample at a time, chosen with the `slot` or `sample`
argument.

### Merging and subsetting SIRUP objects

`mergeSIRUP()` joins several single sample objects into one, with the
sample names that are given.

``` r
merged <- mergeSIRUP(x, x, samples = c("a", "b"), sample_field = "sample")
spatialSamples(merged)
#> [1] "a" "b"
```

`subsetSIRUP()` keeps chosen samples, spots, or receptor rows.

## Next steps

- To build a SIRUP object and run a few simple analyses, see
  [Get started with SIRUP](get-started.md).
- To compare ROIs in the tissue, see
  [Annotate tissue ROIs and analyze local clonal structure](regions.md).
- To add more sections or bulk AIRR-seq, see
  [Analyze multiple spatial sections and bulk BCR-seq](multisample.md).
