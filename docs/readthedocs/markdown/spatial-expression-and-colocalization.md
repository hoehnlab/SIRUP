This vignette shows how to plot gene expression and receptor locations
over a spatial tissue image, then test spatial autocorrelation and
receptor/gene colocalization. The examples use CD3E with TCRs and MS4A1
with BCRs.

For now, this vignette uses the full local thymus dataset from the
walkthrough. We should replace this with a small packaged toy dataset
later.

This read-the-docs version can be rendered in two modes. By default,
code is shown but not evaluated. To generate plots from the full
dataset, render on the cluster with `SIRAP_RENDER_PLOTS=true`.

## Setup

``` r
library(SIRAP)
library(dplyr)
library(ggplot2)
library(cowplot)
library(airr)
library(Seurat)
library(SeuratObject)

sampleId <- "M5"
sampleLower <- tolower(sampleId)

thymusDir <- normalizePath("/dartfs/rc/lab/H/HoehnK/Hunter/Thymus_VisiumHD", mustWork = TRUE)
outDir <- file.path(thymusDir, "SIRAP_walkthrough", sampleId)
finalAirrDir <- file.path(outDir, "03_final_spatial_airr")
plotDir <- file.path(outDir, "plots")
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

seuratPath <- file.path(thymusDir, "data", sprintf("Thymus_VisiumHD_%s_8um_hires.rds", sampleId))
finalBCRPath <- file.path(finalAirrDir, sprintf("%s_bcrs_heavy_final_spatial_airr.tsv", sampleLower))
finalTCRPath <- file.path(finalAirrDir, sprintf("%s_tcrs_beta_final_spatial_airr.tsv", sampleLower))

thymus <- readRDS(seuratPath)
bcrTopBarcode <- airr::read_rearrangement(finalBCRPath)
tcrTopBarcode <- airr::read_rearrangement(finalTCRPath)
```

Add a receptor-type metadata column to the Seurat object. This lets us
use the same object for plotting and spatial statistics.

``` r
bcr_spots <- unique(bcrTopBarcode$barcode_8um)
tcr_spots <- unique(tcrTopBarcode$barcode_8um)

thymus[["sirap_receptor_type"]] <- dplyr::case_when(
  colnames(thymus) %in% bcr_spots & colnames(thymus) %in% tcr_spots ~ "BCR + TCR",
  colnames(thymus) %in% bcr_spots ~ "BCR",
  colnames(thymus) %in% tcr_spots ~ "TCR",
  TRUE ~ "0"
)
```

## Gene Expression Over Tissue

`spatialExpressionPlot()` removes zero-expression points so the tissue
image remains visible.

``` r
cd3ePlot <- spatialExpressionPlot(
  thymus,
  features = "CD3E",
  assay = "Spatial.008um",
  layer = "counts",
  feature_colors = c("grey80", "#d95f02"),
  pt.size.factor = 0.3
)

ms4a1Plot <- spatialExpressionPlot(
  thymus,
  features = "MS4A1",
  assay = "Spatial.008um",
  layer = "counts",
  feature_colors = c("grey80", "#2c7fb8"),
  pt.size.factor = 0.3
)

cowplot::plot_grid(cd3ePlot, ms4a1Plot, ncol = 2)
```

![](/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP/docs/readthedocs/markdown/spatial-expression-and-colocalization_files/figure-gfm/expression-plots-1.png)<!-- -->

## Receptors Over Gene Expression

`plotReceptorExpression()` plots positive gene-expression spots as
filled points and overlays receptor-positive locations as hollow
circles.

``` r
tcrCd3ePlot <- plotReceptorExpression(
  object = thymus,
  receptors = tcrTopBarcode,
  gene = "CD3E",
  barcode_col = "barcode_8um",
  receptor_color = "#d95f02",
  receptor_size = 0.55,
  receptor_stroke = 0.25,
  feature_colors = c("grey80", "#7f2704"),
  feature_size = 0.45,
  title = "CD3E expression with TCR-positive spots",
  receptor_label = "TCR"
)

bcrMs4a1Plot <- plotReceptorExpression(
  object = thymus,
  receptors = bcrTopBarcode,
  gene = "MS4A1",
  barcode_col = "barcode_8um",
  receptor_color = "#2c7fb8",
  receptor_size = 0.55,
  receptor_stroke = 0.25,
  feature_colors = c("grey80", "#08306b"),
  feature_size = 0.45,
  title = "MS4A1 expression with BCR-positive spots",
  receptor_label = "BCR"
)

cowplot::plot_grid(tcrCd3ePlot, bcrMs4a1Plot, ncol = 2)
```

![](/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP/docs/readthedocs/markdown/spatial-expression-and-colocalization_files/figure-gfm/receptor-expression-1.png)<!-- -->

## Moran’s I

Use Moran’s I to test whether a variable is spatially autocorrelated. On
Visium HD data, k-nearest neighbors are often more convenient than
radius-based neighborhoods because the number of spots is large.

``` r
spatialStatsK <- 24

tcrMoran <- spatialMoran(
  thymus,
  variable = "TCR presence",
  variable_type = "receptor",
  receptors = tcrTopBarcode,
  barcode_col = "barcode_8um",
  k = spatialStatsK
)
tcrMoran$test
#> 
#>  Moran I test under randomisation
#> 
#> data:  dat$value  
#> weights: listw    
#> 
#> Moran I statistic standard deviate = 119.27, p-value < 2.2e-16
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>      5.269612e-02     -2.357479e-06      1.952190e-07

cd3eExpressionMoran <- spatialMoran(
  thymus,
  variable = "CD3E",
  variable_type = "feature",
  layer = "counts",
  k = spatialStatsK
)
cd3eExpressionMoran$test
#> 
#>  Moran I test under randomisation
#> 
#> data:  dat$value  
#> weights: listw    
#> 
#> Moran I statistic standard deviate = 18.028, p-value < 2.2e-16
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>      7.962221e-03     -2.357479e-06      1.951858e-07

bcrMoran <- spatialMoran(
  thymus,
  variable = "BCR presence",
  variable_type = "receptor",
  receptors = bcrTopBarcode,
  barcode_col = "barcode_8um",
  k = spatialStatsK
)
bcrMoran$test
#> 
#>  Moran I test under randomisation
#> 
#> data:  dat$value  
#> weights: listw    
#> 
#> Moran I statistic standard deviate = 395.89, p-value < 2.2e-16
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>      1.749105e-01     -2.357479e-06      1.952091e-07

ms4a1ExpressionMoran <- spatialMoran(
  thymus,
  variable = "MS4A1",
  variable_type = "feature",
  layer = "counts",
  k = spatialStatsK
)
ms4a1ExpressionMoran$test
#> 
#>  Moran I test under randomisation
#> 
#> data:  dat$value  
#> weights: listw    
#> 
#> Moran I statistic standard deviate = 398.31, p-value < 2.2e-16
#> alternative hypothesis: greater
#> sample estimates:
#> Moran I statistic       Expectation          Variance 
#>      1.759822e-01     -2.357479e-06      1.952141e-07
```

## Local Moran’s I

Local Moran’s I gives a spot-level statistic that can be plotted back
over the tissue. By default `plotSpatialLocalMoran()` colors the local
Moran z-score.

``` r
cd3eLocalMoran <- spatialLocalMoran(
  thymus,
  variable = "CD3E",
  variable_type = "feature",
  layer = "counts",
  k = spatialStatsK
)

cd3eLocalMoranPlot <- plotSpatialLocalMoran(
  thymus,
  cd3eLocalMoran,
  point_size = 0.55,
  title = paste(sampleId, "CD3E local Moran's I z-score")
)

ms4a1LocalMoran <- spatialLocalMoran(
  thymus,
  variable = "MS4A1",
  variable_type = "feature",
  layer = "counts",
  k = spatialStatsK
)

ms4a1LocalMoranPlot <- plotSpatialLocalMoran(
  thymus,
  ms4a1LocalMoran,
  point_size = 0.55,
  title = paste(sampleId, "MS4A1 local Moran's I z-score")
)

cowplot::plot_grid(cd3eLocalMoranPlot, ms4a1LocalMoranPlot, ncol = 2)
```

![](/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP/docs/readthedocs/markdown/spatial-expression-and-colocalization_files/figure-gfm/local-moran-1.png)<!-- -->

## Ripley’s Cross K

Ripley’s cross K compares the observed spatial relationship between two
marked point types to random relabeling simulations. Here we test
TCR-positive spots near CD3E-positive spots, and BCR-positive spots near
MS4A1-positive spots.

``` r
cd3e_counts <- Seurat::FetchData(thymus, vars = "CD3E", layer = "counts")
ms4a1_counts <- Seurat::FetchData(thymus, vars = "MS4A1", layer = "counts")

thymus[["sirap_cd3e_tcr_mark"]] <- dplyr::case_when(
  colnames(thymus) %in% tcr_spots ~ "TCR",
  cd3e_counts[colnames(thymus), "CD3E"] > 0 ~ "CD3E+",
  TRUE ~ NA_character_
)

thymus[["sirap_ms4a1_bcr_mark"]] <- dplyr::case_when(
  colnames(thymus) %in% bcr_spots ~ "BCR",
  ms4a1_counts[colnames(thymus), "MS4A1"] > 0 ~ "MS4A1+",
  TRUE ~ NA_character_
)
```

``` r
ripleyNsim <- 10 # Increase for real analysis, but this is slow to run on large datasets 

cd3eTcrK <- spatialRipleyCrossK(
  thymus,
  mark_col = "sirap_cd3e_tcr_mark",
  i = "TCR",
  j = "CD3E+",
  nsim = ripleyNsim,
  max_r = 1000
)
#> Generating 10 simulations by evaluating expression  ...
#> 1, 2, 3, 4, 5, 6, 7, 8, 9, 
#> 10.
#> 
#> Done.

cd3eTcrKPlot <- plotRipleyCrossK(
  cd3eTcrK,
  title = paste(sampleId, "Ripley's cross K: TCR to CD3E+"),
  log_y = FALSE
)

ms4a1BcrK <- spatialRipleyCrossK(
  thymus,
  mark_col = "sirap_ms4a1_bcr_mark",
  i = "BCR",
  j = "MS4A1+",
  nsim = ripleyNsim,
  max_r = 1000
)
#> Generating 10 simulations by evaluating expression  ...
#> 1, 2, 3, 4, 5, 6, 7, 8, 9, 
#> 10.
#> 
#> Done.

ms4a1BcrKPlot <- plotRipleyCrossK(
  ms4a1BcrK,
  title = paste(sampleId, "Ripley's cross K: BCR to MS4A1+"),
  log_y = FALSE
)

cowplot::plot_grid(cd3eTcrKPlot, ms4a1BcrKPlot, ncol = 2)
```

![](/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP/docs/readthedocs/markdown/spatial-expression-and-colocalization_files/figure-gfm/ripley-cross-k-1.png)<!-- -->

## Other Spatial Plotting Patterns

Blended expression plots are helpful when two genes form spatial
gradients or partially overlapping niches.

``` r
spatialExpressionBlendPlot(
  thymus,
  features = c("CXCR4", "CXCR5"),
  assay = "Spatial.008um",
  layer = "counts",
  pt.size.factor = 0.2,
  only_blended = TRUE,
  top_left = "blue",
  bottom_left = "grey50",
  bottom_right = "yellow",
  top_right = "green"
)
```

![](/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP/docs/readthedocs/markdown/spatial-expression-and-colocalization_files/figure-gfm/blend-1.png)<!-- -->
