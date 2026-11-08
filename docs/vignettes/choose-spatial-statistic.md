Choose a spatial statistic
================

> **Time:** about 8 minutes. **Data:** a small synthetic SIRUP object
> bundled with SIRUP. **Needs:** SIRUP with spdep.

This page helps you pick a spatial statistic. Start from the question
you ask, not from the function name. Each section below names one
question. It then gives the function, the input, the way to read the
result, and one common mistake. A *spot* is one capture location of the
spatial assay. A *clone* is a group of receptors that descend from one
ancestor cell.

All tests need a SIRUP object, or a Seurat object with spatial images.
The tests read the spot coordinates from the object. Most tests compare
each spot with its neighbors. You set the neighbors with `k` (the number
of nearest spots) or `radius`.
[Test spatial patterns of genes, cells, and clones](spatial-statistics.md) runs every test on
real data.

## Match the question to the function

| Question | Function | Input |
|----|----|----|
| Is this gene or cell type clustered? | `spatialMoran()` | One gene or numeric metadata column |
| Where are the clusters? | `spatialLocalMoran()` | The same input |
| Do two genes or cell types co-vary? | `spatialLee()` | Two continuous variables |
| Do spots of type A sit near spots of type B? | `spatialRipleyCrossK()` | One metadata column with two labels |
| Does one clone stay together? | `spatialCloneJoinCount()` | Receptor table with clone IDs |
| How far does a clone spread? | `spatialCloneSpread()` | Receptor table with clone IDs |
| How does the repertoire differ between regions? | Polygon summaries | Stored polygons |

## Ask whether a variable clusters

Use `spatialMoran()` to ask if neighboring spots have similar values.
The function computes Moran’s I. This statistic is near 0 when the
values show no spatial pattern. It is positive when similar values sit
together. It is negative when values alternate.

The input is one variable. `variable_type` says what it is: a gene
(`"feature"`), a numeric metadata column (`"metadata"`), receptor
presence (`"receptor"`), or one clone (`"clone"`). For a text metadata
column, give `value` to test one level.

Read $I$ as the size of the pattern. Read the p value as the evidence
for it. A large section gives a tiny p value for a weak pattern.

> **Note.** The common mistake is to read a small p value as a strong
> pattern. Check $I$ first. A sparse gene can have $I$ near 0.06 and a p
> value below $10^{-9}$.

`spatialLocalMoran()` answers the next question: where are the clusters?
It gives one value for each spot, with a z-score and a p value. Each
spot has its own p value, so adjust for multiple testing.
`plotSpatialLocalMoran()` draws the result on the tissue.

This chunk tests one gene on the synthetic object:

``` r
library(SIRUP)
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
moran <- spatialMoran(x, variable = "gene1", k = 4)
moran$test$estimate[[1]]
#> [1] 0.2176226
moran$test$p.value
#> [1] 0.01873459
```

## Ask whether two variables co-vary

Use `spatialLee()` to ask if two continuous variables share a spatial
pattern. Lee’s L combines Pearson correlation with spatial smoothing. L
is high when both variables are high in the same neighborhoods. The
variables do not need to peak in the same spots.

The input is two continuous variables. Each can be a gene or a numeric
metadata column. Give `x_type`, `y_type`, `assay`, and `layer` as
needed. The result holds `L` and a test.

Run a positive control. Pick two genes that you expect to share a
pattern. If the control fails, the data or the neighbor setting may be
wrong.

> **Note.** The common mistake is to run the analytical test on a large
> Visium HD section. It builds a dense matrix with one entry for each
> pair of spots. Use `method = "mc"` for large data.

## Ask whether two spot types sit near each other

Use `spatialRipleyCrossK()` to ask if spots of type A lie near spots of
type B. The function counts the type B spots within each distance of an
average type A spot. It compares the curve with curves from random
relabeling.

The input is one metadata column that holds the labels. Give the two
levels as `i` and `j`. Set `nsim` above 0 to get the envelope. Spots
with `NA` in the column are left out.

Read the curve against the envelope. A curve above the envelope means A
and B sit closer than chance. A curve below the envelope means each type
groups with itself. Use `plotRipleyCrossK()` to draw it.

> **Note.** The common mistake is to set `nsim = 0` and read the curve
> alone. With no simulations, the function returns the observed K only.
> You have no reference for chance.

The distance is in image pixels at the stored `image_scale`. Convert it
before you quote a distance in micrometers.

## Ask whether a clone stays together

Use `spatialCloneJoinCount()` to ask if the spots of one clone sit next
to each other. The function marks the spots that hold the clone. It
counts neighbor pairs in which both spots hold the clone. It then calls
`spdep::joincount.test()`.

The input is a receptor table with clone IDs and spatial IDs. A SIRUP
object holds both. By default, the function tests every clone that
occurs in at least `min_spatial_ids` spots. Give `clone_val` to test
chosen clones. If the object has stored polygons, the function runs one
test for each polygon. Set `use_stored_polygons = FALSE` for one test
over the section.

The summary table has one row for each clone. A high `statistic` (a
z-score) and a small `p_value` mean that the clone stays together. The
default alternative is `"greater"`, so the test looks for more joins
than chance.

> **Note.** The common mistake is to compare the z-scores of clones with
> different sizes. A clone in few spots has little power. Check the
> number of spots next to the statistic.

Moran’s I also accepts one clone with `variable_type = "clone"`. It
treats the clone as a 0/1 value for each spot. The join-count test is
the first choice for clones.

## Ask how far a clone spreads

Use `spatialCloneSpread()` to measure the reach of a clone. It is a
description, not a test. The function reports the mean and the maximum
distance between the points of each clone, for each sample.

``` r
spatialCloneSpread(x)
#> # A tibble: 5 × 5
#>   sample  clone_id n_spatial_ids mean_pairwise_distance max_pairwise_distance
#>   <chr>   <chr>            <int>                  <dbl>                 <dbl>
#> 1 spatial 1                    5                   2.09                  3.16
#> 2 spatial 2                    7                   2.81                  5   
#> 3 spatial 3                    7                   2.47                  5   
#> 4 spatial 4                   10                   2.44                  5   
#> 5 spatial 5                    7                   2.44                  4.47
```

The input is a receptor table with clone IDs and spatial IDs. A clone
needs at least `min_spatial_ids` points (default 2). By default, each
sequence is a point, so sequences in the same spot give distances of 0.
Set `collapse_within_spot = TRUE` to count each spot once. Distances are
in image pixels.

> **Note.** The common mistake is to read a large maximum distance as
> proof that a clone has a wide spread. Two distant spots give a large
> maximum with a clone of two points. Read the mean distance and
> `n_spatial_ids` with it. Use the join-count test to ask if the clone
> clusters.

## Compare regions with polygon summaries

The questions above treat the section as one field. To compare regions
that you drew, use the polygon functions. `polygonReceptorSummary()`
counts receptors, spots, and clones in each polygon.
`spatialCloneDominance()` gives the fraction of receptors in each clone
of each polygon. `polygonDiversity()` gives Hill diversity curves,
rarefied to one size. These functions describe each region. They do not
test a spatial pattern.

> **Note.** The common mistake is to compare the diversity of polygons
> with different receptor counts. `polygonDiversity()` leaves out
> polygons with fewer than `min_n` receptors. Set `min_n` before you
> compare regions.

For the full workflow, see [Draw tissue regions and measure clonal dominance and diversity](regions.md).

## Next steps

- Run every test on real data:
  [Test spatial patterns of genes, cells, and clones](spatial-statistics.md).
- Draw regions and measure their repertoires:
  [Draw tissue regions and measure clonal dominance and diversity](regions.md).
- Call clones before you test them:
  [Call clones and build lineage trees on the tissue](clones-and-trees.md).
