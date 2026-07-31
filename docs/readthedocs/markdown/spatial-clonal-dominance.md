Spatial Clonal Dominance
================

This vignette shows how to draw spatial regions, calculate clone
dominance inside those regions, and make region-level plots for clone
IDs and V-gene annotations.

Clone dominance is the fraction of BCR rows inside a polygon assigned to
each clone. `spatialCloneDominance()` also returns the BCR rows and
spatial coordinates inside each polygon. Its result and the original
SIRAP or Seurat object are the two required inputs to the regional
plotting functions.

By default, this page shows runnable code but does not evaluate it. To
render plots from the full local thymus dataset, render on the cluster
with `SIRAP_RENDER_PLOTS=true`.

## Setup

``` r
library(SIRAP)
library(dplyr)
library(ggplot2)
library(data.table)
library(airr)
library(Seurat)
library(SeuratObject)
library(scoper)

sampleId <- "M5"
sampleLower <- tolower(sampleId)
nproc <- 8

thymusDir <- normalizePath("/dartfs/rc/lab/H/HoehnK/Hunter/Thymus_VisiumHD", mustWork = TRUE)
outDir <- file.path(thymusDir, "SIRAP_walkthrough", sampleId)
finalAirrDir <- file.path(outDir, "03_final_spatial_airr")
plotDir <- file.path(outDir, "plots")
dominanceOutDir <- file.path(plotDir, "clonal_dominance")
dir.create(dominanceOutDir, recursive = TRUE, showWarnings = FALSE)

seuratPath <- file.path(thymusDir, "data", sprintf("Thymus_VisiumHD_%s_8um_hires.rds", sampleId))
finalBCRPath <- file.path(finalAirrDir, sprintf("%s_bcrs_heavy_final_spatial_airr.tsv", sampleLower))
bcrClonesPath <- file.path(outDir, sprintf("%s_bcr_hierarchical_clones.rds", sampleLower))

thymus <- readRDS(seuratPath)
bcrTopBarcode <- airr::read_rearrangement(finalBCRPath)
```

If clones have not already been called, call them with `scoper`.

``` r
if (file.exists(bcrClonesPath)) {
  bcrClones <- readRDS(bcrClonesPath)
} else {
  bcrClones <- scoper::hierarchicalClones(
    filter(bcrTopBarcode, productive),
    threshold = 0.1,
    cell_id = "sequence_id",
    summarize_clones = FALSE,
    nproc = nproc
  )
  saveRDS(bcrClones, bcrClonesPath)
}
```

## Select Polygon Regions

`selectSpatialPolygons()` opens a Shiny app for selecting polygon
vertices over the spatial tissue coordinates. Each click snaps to the
nearest full-resolution spatial barcode, even if `downsample_n` is used
to make the browser plot lighter.

For an HPC interactive node, define the host and port explicitly. Start
the R session on the interactive node, run the code below, then tunnel
the port from your laptop. For example:

``` bash
ssh -N -L 8788:t01.hpcc.dartmouth.edu:8788 f007p0j@discovery.dartmouth.edu
```

Then open <http://127.0.0.1:8788> in your laptop browser. Use
`host = "0.0.0.0"` when Shiny is running on the compute node and the
browser connection arrives through the login node tunnel.

``` r
gcMantlePolygonPath <- file.path(
  dominanceOutDir,
  sprintf("%s_big_center_gc_plus_mantle_polygon_vertices.csv", sampleLower)
)
gcPolygonPath <- file.path(
  dominanceOutDir,
  sprintf("%s_big_gc_polygon_vertices.csv", sampleLower)
)
allGcsPolygonsPath <- file.path(
  dominanceOutDir,
  sprintf("%s_all_gcs_polygons_vertices.csv", sampleLower)
)

dominanceShinyHost <- "0.0.0.0"
dominanceShinyPort <- 8788L
dominanceFeature <- "MS4A1"
dominanceDownsampleN <- 20000
dominancePlotImage <- TRUE
```

Draw or load the GC plus mantle polygon.

``` r
if (file.exists(gcMantlePolygonPath)) {
  gcMantlePolygon <- fread(gcMantlePolygonPath) %>% as.data.frame()
} else {
  message(
    "No saved GC plus mantle polygon found. Starting Shiny on ",
    dominanceShinyHost, ":", dominanceShinyPort, "."
  )

  gcMantlePolygon <- selectSpatialPolygons(
    object = thymus,
    feature = dominanceFeature,
    image_scale = "hires",
    save_path = gcMantlePolygonPath,
    downsample_n = dominanceDownsampleN,
    host = dominanceShinyHost,
    port = dominanceShinyPort,
    plot_image = dominancePlotImage,
    launch_browser = FALSE
  )
}
```

Draw or load a single GC polygon.

``` r
if (file.exists(gcPolygonPath)) {
  gcPolygon <- fread(gcPolygonPath) %>% as.data.frame()
} else {
  message(
    "No saved GC polygon found. Starting Shiny on ",
    dominanceShinyHost, ":", dominanceShinyPort, "."
  )

  gcPolygon <- selectSpatialPolygons(
    object = thymus,
    feature = dominanceFeature,
    image_scale = "hires",
    save_path = gcPolygonPath,
    downsample_n = dominanceDownsampleN,
    host = dominanceShinyHost,
    port = dominanceShinyPort,
    plot_image = dominancePlotImage,
    launch_browser = FALSE
  )
}
```

Draw or load multiple GC polygons. If you draw more than one polygon in
the Shiny app, the returned table has a `polygon` column and downstream
dominance summaries are calculated independently for each polygon.

``` r
if (file.exists(allGcsPolygonsPath)) {
  allGcsPolygons <- fread(allGcsPolygonsPath) %>% as.data.frame()
} else {
  message(
    "No saved all-GC polygons found. Starting Shiny on ",
    dominanceShinyHost, ":", dominanceShinyPort, "."
  )

  allGcsPolygons <- selectSpatialPolygons(
    object = thymus,
    feature = dominanceFeature,
    image_scale = "hires",
    save_path = allGcsPolygonsPath,
    downsample_n = dominanceDownsampleN,
    host = dominanceShinyHost,
    port = dominanceShinyPort,
    plot_image = dominancePlotImage,
    launch_browser = FALSE
  )
}
```

Plot the selected polygons to check that the regions are correct before
using them for dominance summaries.

``` r
allGcsPolygonOverlay <- plotSpatialPolygons(
  object = thymus,
  polygon = allGcsPolygons,
  image_scale = "hires",
  polygon_alpha = 0.2
)

ggsave(
  filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_polygons_overlay.pdf", sampleLower)),
  plot = allGcsPolygonOverlay,
  width = 7,
  height = 6
)
```

## Calculate Clonal Dominance

For this example, the denominator is the whole heavy-chain BCR pool
inside each polygon. `exclude_singletons = TRUE` removes clone IDs
represented by only one BCR row before calculating dominance.

``` r
dominanceBcrs <- bcrClones %>%
  mutate(
    spaced_clone_id = paste0("heavy_", clone_id),
    chain_for_pairing = "heavy"
  )

clonalDominance <- spatialCloneDominance(
  bcrs = dominanceBcrs,
  object = thymus,
  polygon = allGcsPolygons,
  barcode_col = "barcode_8um",
  clone_col = "spaced_clone_id",
  locus_col = "chain_for_pairing",
  locus = "heavy"
)

clonalDominanceNoSingletons <- spatialCloneDominance(
  bcrs = dominanceBcrs,
  object = thymus,
  polygon = allGcsPolygons,
  barcode_col = "barcode_8um",
  clone_col = "spaced_clone_id",
  locus_col = "chain_for_pairing",
  locus = "heavy",
  exclude_singletons = TRUE
)

head(clonalDominance$dominance)
```

The returned object contains:

- `dominance`: clone counts and proportions for each polygon.
- `region_barcodes`: spatial barcodes inside each polygon.
- `region_coordinates`: spatial coordinate rows inside each polygon.
- `region_bcrs`: BCR rows inside each polygon, with plotting coordinates
  and V-gene summaries.
- `polygons`: polygon vertices prepared for plotting.

`dominance` includes allele-level `v_gene`, allele-stripped
`v_gene_family`, and broad `v_gene_superfamily` values for each clone.

## Plot Regional Receptors

`plotSpatialRegionReceptors()` requires the dominance result and its
original SIRAP or Seurat object. It plots the regional receptor rows
already stored by `spatialCloneDominance()` over the tissue image,
without repeating barcode or polygon arguments. In this example those
rows are heavy-chain BCRs. `top_n` keeps the most common displayed
levels and groups lower-frequency levels as `Other`.

``` r
bcrClonePlot <- plotSpatialRegionReceptors(
  clonalDominance,
  object = thymus,
  color_by = "clone_id",
  top_n = 20,
  point_size = 0.5,
  point_alpha = 0.5,
  polygon_alpha = 0.08
)
ggsave(
  filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_bcrs_by_clone_id.pdf", sampleLower)),
  plot = bcrClonePlot,
  width = 8,
  height = 7
)

bcrVGenePlot <- plotSpatialRegionReceptors(
  clonalDominance,
  object = thymus,
  color_by = "v_gene",
  top_n = 20,
  point_size = 0.8,
  point_alpha = 0.75,
  polygon_alpha = 0.08
)
ggsave(
  filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_bcrs_by_v_gene.pdf", sampleLower)),
  plot = bcrVGenePlot,
  width = 8,
  height = 7
)

bcrVFamilyPlot <- plotSpatialRegionReceptors(
  clonalDominance,
  object = thymus,
  color_by = "v_gene_family",
  top_n = 20,
  point_size = 0.8,
  point_alpha = 0.75,
  polygon_alpha = 0.08
)
ggsave(
  filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_bcrs_by_v_gene_family.pdf", sampleLower)),
  plot = bcrVFamilyPlot,
  width = 8,
  height = 7
)

bcrVSuperfamilyPlot <- plotSpatialRegionReceptors(
  clonalDominance,
  object = thymus,
  color_by = "v_gene_superfamily",
  point_size = 0.8,
  point_alpha = 0.75,
  polygon_alpha = 0.08
)
ggsave(
  filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_bcrs_by_v_gene_superfamily.pdf", sampleLower)),
  plot = bcrVSuperfamilyPlot,
  width = 8,
  height = 7
)
```

## Plot Pie Summaries

`plotSpatialPolygonPies()` likewise requires the dominance result and
its original SIRAP or Seurat object. It makes pie summaries directly
from the regional rows and polygons stored in the dominance result. With
`mode = "side_by_side"`, the tissue image and numbered polygons appear
on the left and the pie charts appear on the right. The top levels are
identified within each polygon, then unioned across polygons so a clone
that is top-ranked in any polygon remains named everywhere it appears.
Other levels are grouped as `Other`.

``` r
clonePies <- plotSpatialPolygonPies(
  clonalDominance,
  object = thymus,
  fill_by = "clone_id",
  mode = "side_by_side",
  top_n = 10,
  label_top_n = 5,
  polygon_alpha = 0.2
)
ggsave(
  filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_clone_id_pies.pdf", sampleLower)),
  plot = clonePies,
  width = 12,
  height = 7
)

vFamilyPies <- plotSpatialPolygonPies(
  clonalDominance,
  object = thymus,
  fill_by = "v_gene_family",
  mode = "side_by_side",
  top_n = 10,
  label_top_n = 10,
  polygon_alpha = 0.2
)
ggsave(
  filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_v_gene_family_pies.pdf", sampleLower)),
  plot = vFamilyPies,
  width = 12,
  height = 7
)

vSuperfamilyPies <- plotSpatialPolygonPies(
  clonalDominance,
  object = thymus,
  fill_by = "v_gene_superfamily",
  mode = "side_by_side",
  top_n = 10,
  label_top_n = 10,
  polygon_alpha = 0.2
)
ggsave(
  filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_v_gene_superfamily_pies.pdf", sampleLower)),
  plot = vSuperfamilyPies,
  width = 12,
  height = 7
)
```

The spatially clipped pie mode fills each polygon with the pie slices.
This requires the suggested `sf` package. Labels are usually best left
off for this mode; use `spatial_legend = TRUE` for variables with a
manageable number of levels, such as V-gene superfamily.

``` r
if (requireNamespace("sf", quietly = TRUE)) {
  cloneSpatialPies <- plotSpatialPolygonPies(
    clonalDominance,
    object = thymus,
    fill_by = "clone_id",
    mode = "spatial_clipped",
    top_n = 10,
    label_top_n = 0,
    pie_alpha = 0.9
  )
  ggsave(
    filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_clone_id_spatial_pies.pdf", sampleLower)),
    plot = cloneSpatialPies,
    width = 8,
    height = 7
  )

  vFamilySpatialPies <- plotSpatialPolygonPies(
    clonalDominance,
    object = thymus,
    fill_by = "v_gene_family",
    mode = "spatial_clipped",
    top_n = 10,
    label_top_n = 0,
    pie_alpha = 0.9
  )
  ggsave(
    filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_v_gene_family_spatial_pies.pdf", sampleLower)),
    plot = vFamilySpatialPies,
    width = 8,
    height = 7
  )

  vSuperfamilySpatialPies <- plotSpatialPolygonPies(
    clonalDominance,
    object = thymus,
    fill_by = "v_gene_superfamily",
    mode = "spatial_clipped",
    top_n = 10,
    label_top_n = 0,
    pie_alpha = 0.9,
    spatial_legend = TRUE
  )
  ggsave(
    filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_v_gene_superfamily_spatial_pies.pdf", sampleLower)),
    plot = vSuperfamilySpatialPies,
    width = 8,
    height = 7
  )
} else {
  message("Skipping spatially clipped pie charts because sf is not installed.")
}
```
