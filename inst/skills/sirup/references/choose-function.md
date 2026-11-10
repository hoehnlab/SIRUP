# Choose a function

Start from the question that the user asks. Then pick the function. Read the usage in `references/functions-*.md` before you call it.

## Spatial patterns

| Question | Function | Input |
|---|---|---|
| Is this gene or cell type clustered? | `spatialMoran()` | One gene or one numeric metadata column |
| Where are the clusters? | `spatialLocalMoran()`, then `plotSpatialLocalMoran()` | The same input |
| Do two genes or cell types co-vary? | `spatialLee()` | Two continuous variables |
| Do spots of type A sit near spots of type B? | `spatialRipleyCrossK()`, then `plotRipleyCrossK()` | One metadata column with two labels |
| Does one clone stay together? | `spatialCloneJoinCount()`, then `plotSpatialCloneJoinCount()` | Receptor table with clone IDs |
| How far does a clone spread? | `spatialCloneSpread()` | Receptor table with clone IDs |

Read each result as follows:

- **Moran's I** is near 0 when the values show no spatial pattern. It is positive when similar values sit together. It is negative when values alternate. Read I as the size of the pattern and the p value as the evidence for it. A large section gives a small p value for a weak pattern, so look at I first. Variable types: a gene (`"feature"`), a numeric metadata column (`"metadata"`), receptor presence (`"receptor"`), or one clone (`"clone"`). For a text metadata column, give `value` to test one level.
- **Local Moran's I** gives one value for each spot, with a z-score and a p value. Each spot has its own p value, so adjust for multiple testing.
- **Lee's L** combines Pearson correlation with spatial smoothing. L is high when both variables are high in the same neighborhoods. The default test is analytical and builds a dense matrix with one entry for each pair of spots. For large data, such as Visium HD, set `method = "mc"` to use permutations. As a positive control, test two genes that you expect to share a pattern.
- **Ripley's cross K** counts the type B spots within each distance of an average type A spot. It compares the curve with curves from random relabeling. Set `nsim` above 0 to get the envelope. A curve above the envelope means that A and B sit closer than chance. A curve below it means that each type groups with itself. Distances are in image pixels at the stored `image_scale`. Convert before you quote micrometers.
- **Join counts** mark the spots that hold the clone and count neighbor pairs in which both spots hold it. The summary table has one row for each clone. A high `statistic` (a z-score) and a small `p_value` mean that the clone stays together. By default the function tests every clone that occurs in at least `min_spatial_ids` spots. Give `clone_val` to test chosen clones. With stored polygons, it runs one test for each polygon. Set `use_stored_polygons = FALSE` for one test over the section. Compare the number of spots of each clone next to the statistic, because a clone in few spots has little power.
- **Clone spread** is a description, not a test. It reports the mean and the maximum distance between the points of each clone, for each sample. By default each sequence is a point, so sequences in one spot give distances of 0. Set `collapse_within_spot = TRUE` to count each spot once. Read the mean distance and `n_spatial_ids` together with the maximum.

The neighbors of a spot come from `k` (the number of nearest spots) or `radius`. Use `k = 6` for the 10x Visium hexagonal grid.

## Regions

| Question | Function |
|---|---|
| How many receptors, spots, and clones does each region hold? | `polygonReceptorSummary()` |
| Which clones dominate each region? | `spatialCloneDominance()`, then `plotSpatialPolygonPies()` or `plotSpatialRegionReceptors()` |
| How diverse is the repertoire of each region? | `polygonDiversity()`, then `summarizePolygonDiversity()` |
| How abundant is each clone in a region? | `polygonAbundance()` |
| How far apart are the regions? | `polygonPairwiseDistance()` |
| Which receptors sit in a zone around a region? | `expandPolygons()`, then `addPolygons()` |
| What do the receptors of each region look like, with my own summary? | `withPolygons()` |

These functions describe each region. They do not test a spatial pattern. `polygonDiversity()` leaves out polygons with fewer than `min_n` receptors (default 30).

## Clones and trees

| Question | Function |
|---|---|
| Which sequences belong to one clone? | `updateReceptors()` with `scoper::hierarchicalClones` (BCRs) or `scoper::identicalClones` (TCRs) |
| How mutated is each sequence? | `updateReceptors()` with `shazam::observedMutations` |
| What is the lineage tree of a clone? | `formatSpatialClones()`, then `getSpatialTrees()` |
| Which tips of a tree are close to each other? | `clusterTreesTips()`, `tipPatristicMatrix()` |
| Where does each tip sit on the tissue? | `plotSpatialCloneTreePairs()` |
| How dense is a clone in 3D? | `plotSpatialCloneDensity3D()` (needs rgl and HTML output) |

## Plots

| Show | Function |
|---|---|
| Gene expression, or a numeric spot value, on the tissue | `spatialExpressionPlot()` |
| Spots that hold a receptor | `plotReceptors()` |
| Expression and receptors together | `plotReceptorsExpression()` |
| The regions | `plotSpatialPolygons()` |
| Several genes, or blended colors | `spatialExpressionMultiPlot()`, `spatialExpressionBlendPlot()` |
| Read-length quality control | `plotLengthQc()` |

`spatialExpressionPlot()` hides spots with zero expression. Set `zero_tol = -Inf` to show them.
