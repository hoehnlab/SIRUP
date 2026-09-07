Processing Visium HD BCR and TCR Data
================

This vignette walks through the Visium HD BCR/TCR processing workflow.
It starts from RAD output, keeps spatial barcodes attached through the
pRESTO/IgBLAST consensus workflow, produces AIRR-style tables with
spatial coordinates, and plots BCR/TCR locations over the tissue.

For now, this vignette uses the full local thymus dataset from the
walkthrough. We should replace this with a small packaged toy dataset
later so the vignette can be run anywhere without large external files.

This processing guide intentionally does not evaluate code during
rendering, even in the read-the-docs version, because the processing
workflow writes files and can launch external pRESTO/IgBLAST commands.

## Setup

``` r
library(SIRAP)
library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(airr)
library(Seurat)
library(SeuratObject)

sampleId <- "M5"
sampleLower <- tolower(sampleId)
nproc <- 8

thymusDir <- normalizePath("../../Thymus_VisiumHD", mustWork = TRUE)
demuxBySample <- c(
  M1 = "022526_Visium_demuxed_bc23",
  M5 = "022526_Visium_demuxed_bc24"
)
demuxDir <- file.path(
  "/dartfs/rc/lab/H/HoehnK/Chinmay/022126_Visium_HD",
  demuxBySample[[sampleId]]
)

rawBCRIgblastPath <- file.path(demuxDir, "igblast_out", "igblast_bcr_1M_example.tsv")
rawTCRIgblastPath <- file.path(demuxDir, "igblast_out", "igblast_tcr_1M_example.tsv")

outDir <- file.path(thymusDir, "SIRAP_walkthrough", sampleId)
prestoInputDir <- file.path(outDir, "01_length_filtered_spatial_fastas")
prestoOutputDir <- file.path(outDir, "02_presto_igblast")
finalAirrDir <- file.path(outDir, "03_final_spatial_airr")
plotDir <- file.path(outDir, "plots")
for (dir in c(prestoInputDir, prestoOutputDir, finalAirrDir, plotDir)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

finalBCRPath <- file.path(finalAirrDir, sprintf("%s_bcrs_heavy_final_spatial_airr.tsv", sampleLower))
finalTCRPath <- file.path(finalAirrDir, sprintf("%s_tcrs_beta_final_spatial_airr.tsv", sampleLower))
```

## Create or Load the Seurat Object

If you need to create the Visium HD Seurat object directly from Space
Ranger output, use `Load10X_Spatial()`. Visium HD loading requires
`arrow`, which can be difficult to install, but SIRAP’s downstream
plotting and analyses do not otherwise require `arrow`. We also find it
easier to create a .csv version of the barcode mappings rather than
relying on the `read_parquet()` function in `arrow` to read the parquet
files directly.

``` r
# Build seurat object
m5_path <- "./analysis/M5_260129/outs/"
thy_m5 <- Load10X_Spatial(
  data.dir = m5_path,
  bin.size = c(8),
  image.name = "tissue_hires_image.png"
)

saveRDS(thy_m5, file.path(thymusDir, "data", "Thymus_VisiumHD_M5_8um_hires.rds"))

# Resave tissue locations as .csvs
m5_barcode_mapping <- read_parquet(file.path(m5_path, "barcode_mappings.parquet"))
write.csv(m5_barcode_mapping, file.path(thymusDir, "M5_barcode_mappings.csv"), 
          row.names = FALSE)

m5_pos_2um <- read_parquet(file.path(
  m5_path, "binned_outputs", "square_002um", "spatial", "tissue_positions.parquet"
  ))
write.csv(m5_pos_2um, file.path(thymusDir, sprintf("%s_tissue_positions_2um.csv", sampleId)),
          row.names = FALSE)

m5_pos_8um <- read_parquet(file.path(
  m5_path, "binned_outputs", "square_008um", "spatial", "tissue_positions.parquet"
  ))
write.csv(m5_pos_8um, file.path(thymusDir, sprintf("%s_tissue_positions_8um.csv", sampleId)),
          row.names = FALSE)

m5_pos_16um <- read_parquet(file.path(
  m5_path, "binned_outputs", "square_016um", "spatial", "tissue_positions.parquet"
  ))
write.csv(m5_pos_16um, file.path(thymusDir, sprintf("%s_tissue_positions_16um.csv", sampleId)),
          row.names = FALSE)
```

For the rest of the vignette, we load the saved object used in the
walkthrough.

``` r
seuratPath <- file.path(thymusDir, "data", sprintf("Thymus_VisiumHD_%s_8um_hires.rds", sampleId))
thymus <- readRDS(seuratPath)
```

## Load Spatial Metadata

The receptor reads carry the 2 um barcode in the sequence identifier.
After pRESTO and IgBLAST, we join through the Visium HD barcode map to
recover the 8 um and 16 um barcodes used for plotting and analysis.

``` r
barcodeMapPath <- file.path(thymusDir, "M5_barcode_mappings.csv")
positionsPath8um <- file.path(thymusDir, sprintf("%s_tissue_positions_8um.csv", sampleId))

barcodeMap <- readBarcodeMapping(barcodeMapPath)

positions8um <- readBarcodeMapping(positionsPath8um) %>%
  mutate(
    barcode_8um = barcode,
    x = pxl_col_in_fullres,
    y = max(pxl_row_in_fullres) - pxl_row_in_fullres
  )

head(barcodeMap)
head(positions8um)
```

## Start from Raw IgBLAST Output

Extract the 2 um spatial barcode before pRESTO and do a first
length/locus filter. `filterAirr()` keeps preprocessing in memory; we
explicitly write AIRR tables and FASTA files only when we need pRESTO
inputs.

``` r
bcrRaw <- data.table::fread(rawBCRIgblastPath)
tcrRaw <- data.table::fread(rawTCRIgblastPath)

bcrRaw <- extractSpatialBarcodes(
  bcrRaw,
  id = "sequence_id",
  barcode = "barcode_2um",
  regex = "[0-9]+_[0-9]+$",
  prefix = "s_002um_",
  suffix = "-1"
)

tcrRaw <- extractSpatialBarcodes(
  tcrRaw,
  id = "sequence_id",
  barcode = "barcode_2um",
  regex = "[0-9]+_[0-9]+$",
  prefix = "s_002um_",
  suffix = "-1"
)

bcrHeavy <- filterAirr(
  bcrRaw,
  locus_val = "IGH",
  min_v_length = 100,
  min_j_length = 10,
  min_junction_length = 1,
  barcode = "barcode_2um",
  add_spatial_sequence_id = TRUE
)

tcrBeta <- filterAirr(
  tcrRaw,
  locus_val = "TRB",
  min_v_length = 100,
  min_j_length = 10,
  min_junction_length = 1,
  barcode = "barcode_2um",
  add_spatial_sequence_id = TRUE
)

spatialReceptorSummary(bcrHeavy, barcode = "barcode_2um")
spatialReceptorSummary(tcrBeta, barcode = "barcode_2um")
```

``` r
bcrHeavySpatialTsv <- file.path(prestoInputDir, sprintf("%s_igblast_bcr_pass_length_spatial_barcode.tsv", sampleLower))
bcrHeavySpatialFasta <- file.path(prestoInputDir, sprintf("%s_igblast_bcr_pass_length_spatial_barcode.fasta", sampleLower))
tcrBetaSpatialTsv <- file.path(prestoInputDir, sprintf("%s_igblast_tcr_pass_length_spatial_barcode.tsv", sampleLower))
tcrBetaSpatialFasta <- file.path(prestoInputDir, sprintf("%s_igblast_tcr_pass_length_spatial_barcode.fasta", sampleLower))

airr::write_rearrangement(bcrHeavy, bcrHeavySpatialTsv)
airr::write_rearrangement(tcrBeta, tcrBetaSpatialTsv)

dowser::dfToFasta(bcrHeavy, file = bcrHeavySpatialFasta, id = "spatial_sequence_id", seq = "sequence_alignment")
dowser::dfToFasta(tcrBeta, file = tcrBetaSpatialFasta, id = "spatial_sequence_id", seq = "sequence_alignment")
```

## Run pRESTO and IgBLAST from R

`processSpatialSample()` wraps the pRESTO/IgBLAST command-line workflow
and prints each command and its output to the R console.

``` r
identityThreshold <- 0.975
consensusFreqThreshold <- 0.51
germlineGlob <- "~/software/germlines/imgt/human/vdj/imgt_human_*.fasta"
igblastBase <- "~/software/igblast"
muscleExec <- "~/software/muscle3.8.31_i86linux64"

bcrPrestoPaths <- processSpatialSample(
  input_fasta = bcrHeavySpatialFasta,
  output_dir = prestoOutputDir,
  receptor = "ig",
  identity = identityThreshold,
  consensus_freq = consensusFreqThreshold,
  germline_glob = germlineGlob,
  igblast_base = igblastBase,
  muscle_exec = muscleExec,
  nproc = nproc
)

tcrPrestoPaths <- processSpatialSample(
  input_fasta = tcrBetaSpatialFasta,
  output_dir = prestoOutputDir,
  receptor = "tr",
  identity = identityThreshold,
  consensus_freq = consensusFreqThreshold,
  germline_glob = germlineGlob,
  igblast_base = igblastBase,
  muscle_exec = muscleExec,
  nproc = nproc
)
```

## Final AIRR Tables with Spatial Coordinates

After pRESTO and IgBLAST, read the consensus `db-pass` files, extract
the barcode from the pRESTO sequence identifier, remap to 8 um/16 um
barcodes, add tissue positions, and collapse obvious spatial bleedover.

``` r
bcrConsensus <- airr::read_rearrangement(bcrPrestoPaths$db_pass)
tcrConsensus <- airr::read_rearrangement(tcrPrestoPaths$db_pass)

bcrConsensus <- extractSpatialBarcodes(
  bcrConsensus,
  id = "sequence_id",
  barcode = "barcode_2um",
  regex = ",([^,]+)"
)

tcrConsensus <- extractSpatialBarcodes(
  tcrConsensus,
  id = "sequence_id",
  barcode = "barcode_2um",
  regex = ",([^,]+)"
)

bcrConsensusFiltered <- filterAirr(
  bcrConsensus,
  locus_val = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = 2,
  max_opening_gap = 2,
  barcode = "barcode_2um"
)

tcrConsensusFiltered <- filterAirr(
  tcrConsensus,
  locus_val = "TRB",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = 2,
  max_opening_gap = 2,
  barcode = "barcode_2um"
)

bcrConsensusFiltered <- bcrConsensusFiltered %>%
  mutate(barcode_2um_swapped = sub("^(.*_)(\\d+)_(\\d+)(-\\d+)$", "\\1\\3_\\2\\4", barcode_2um)) %>%
  left_join(
    barcodeMap %>% select(square_002um, square_008um, square_016um),
    by = c("barcode_2um_swapped" = "square_002um")
  ) %>%
  rename(barcode_8um = square_008um, barcode_16um = square_016um)

tcrConsensusFiltered <- tcrConsensusFiltered %>%
  mutate(barcode_2um_swapped = sub("^(.*_)(\\d+)_(\\d+)(-\\d+)$", "\\1\\3_\\2\\4", barcode_2um)) %>%
  left_join(
    barcodeMap %>% select(square_002um, square_008um, square_016um),
    by = c("barcode_2um_swapped" = "square_002um")
  ) %>%
  rename(barcode_8um = square_008um, barcode_16um = square_016um)

bcrConsensusFilteredPos <- addTissuePositions(
  bcrConsensusFiltered,
  positions8um,
  barcode = "barcode_8um",
  position_barcode = "barcode_8um",
  x = "x",
  y = "y"
) %>% filter(!is.na(x), !is.na(y))

tcrConsensusFilteredPos <- addTissuePositions(
  tcrConsensusFiltered,
  positions8um,
  barcode = "barcode_8um",
  position_barcode = "barcode_8um",
  x = "x",
  y = "y"
) %>% filter(!is.na(x), !is.na(y))

bcrCollapsed <- collapseIdenticalNeighbors(
  data = bcrConsensusFilteredPos,
  tissue_positions = positions8um,
  barcode = "barcode_8um",
  x = "x",
  y = "y",
  sequence = "sequence_alignment",
  consensus_count = "consensus_count",
  num_neighbors_non_boundary = 8L
)

tcrCollapsed <- collapseIdenticalNeighbors(
  data = tcrConsensusFilteredPos,
  tissue_positions = positions8um,
  barcode = "barcode_8um",
  x = "x",
  y = "y",
  sequence = "sequence_alignment",
  consensus_count = "consensus_count",
  num_neighbors_non_boundary = 8L
)

bcrTopBarcode <- collapseTopConsensusPerBarcode(bcrCollapsed, barcode = "barcode_2um", locus = "locus", ties = "all")
tcrTopBarcode <- collapseTopConsensusPerBarcode(tcrCollapsed, barcode = "barcode_2um", locus = "locus", ties = "all")

airr::write_rearrangement(bcrTopBarcode, finalBCRPath)
airr::write_rearrangement(tcrTopBarcode, finalTCRPath)
```

## Plot BCRs and TCRs Over Tissue

`plotReceptors()` and `spatialDiscretePlot()` draw points over the
tissue image rather than over a blank coordinate plane.

``` r
allBCRPlot <- plotReceptors(
  object = thymus,
  receptors = dplyr::mutate(bcrTopBarcode, receptor_type = "BCR"),
  group = "receptor_type",
  group_val = "BCR",
  barcode = "barcode_8um",
  colors = "white",
  receptor_size = 0.5,
  title = paste(sampleId, "all BCR-positive spots")
)

allTCRPlot <- plotReceptors(
  object = thymus,
  receptors = dplyr::mutate(tcrTopBarcode, receptor_type = "TCR"),
  group = "receptor_type",
  group_val = "TCR",
  barcode = "barcode_8um",
  colors = "#d95f02",
  receptor_size = 0.5,
  title = paste(sampleId, "all TCR-positive spots")
)

bcr_spots <- unique(bcrTopBarcode$barcode_8um)
tcr_spots <- unique(tcrTopBarcode$barcode_8um)
thymus[["sirap_receptor_type"]] <- dplyr::case_when(
  colnames(thymus) %in% bcr_spots & colnames(thymus) %in% tcr_spots ~ "BCR + TCR",
  colnames(thymus) %in% bcr_spots ~ "BCR",
  colnames(thymus) %in% tcr_spots ~ "TCR",
  TRUE ~ "0"
)

bcrTCROverlay <- spatialDiscretePlot(
  thymus,
  meta_col = "sirap_receptor_type",
  zero_value = "0",
  image_scale = "hires",
  swap_coords = FALSE,
  point_size = 0.5,
  cols = c("BCR" = "#2c7fb8", "TCR" = "#d95f02", "BCR + TCR" = "white")
) +
  ggtitle(paste(sampleId, "BCR and TCR locations"))

cowplot::plot_grid(allBCRPlot, allTCRPlot, bcrTCROverlay, ncol = 3)
```
