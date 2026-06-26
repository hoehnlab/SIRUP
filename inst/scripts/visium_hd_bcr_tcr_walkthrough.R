###############################################################################
# SIRAP Visium HD BCR/TCR walkthrough
#
# This is a scratchpad-style script for processing and visualizing the thymus
# Visium HD BCR/TCR data. It follows the same logic as the original
# Thymus_VisiumHD processing script, but uses SIRAP functions so we can see
# what the package is doing and where the workflow should feel nicer.
#
# I recommend running this one section at a time.
#
#   setwd("/Volumes/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP")
#   devtools::load_all()
#   source("inst/scripts/visium_hd_bcr_tcr_walkthrough.R")
###############################################################################

# ---- Setup ------------------------------------------------------------------

if (!requireNamespace("SIRAP", quietly = TRUE)) {
  if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".")
  }
}

suppressPackageStartupMessages({
  library(SIRAP)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(data.table)
  library(airr)
})

# ---- Things to change --------------------------------------------------------

# Settings
sampleId <- "M5" # M1 uses demuxed_bc23; M5 uses demuxed_bc24.
sampleLower <- tolower(sampleId)

# Directory and paths for data
thymusDir <- normalizePath("../../Thymus_VisiumHD", mustWork = TRUE)
demuxBySample <- c(M1 = "022526_Visium_demuxed_bc23", M5 = "022526_Visium_demuxed_bc24") # Names of demuxed samples
demuxDir <- file.path("/dartfs/rc/lab/H/HoehnK/Chinmay/022126_Visium_HD", demuxBySample[[sampleId]])

# Using 1M example files from RAD output for speed
rawBCRIgblastPath <- file.path(demuxDir, "igblast_out", "igblast_bcr_1M_example.tsv")
rawTCRIgblastPath <- file.path(demuxDir, "igblast_out", "igblast_tcr_1M_example.tsv")
visiumOutsDir <- file.path(thymusDir, "analysis", sprintf("%s_260129", sampleId), "outs")

# Either specify the parquet file or, as we've done here, a csv that is saved from the parquet so we don't have to use arrow.
barcodeMapPath <- file.path(thymusDir, "M5_barcode_mappings.csv")
positionsPath2um <- file.path(thymusDir, sprintf("%s_tissue_positions_2um.csv", sampleId))
positionsPath8um <- file.path(thymusDir, sprintf("%s_tissue_positions_8um.csv", sampleId))
positionsPath16um <- file.path(thymusDir, sprintf("%s_tissue_positions_16um.csv", sampleId))
seuratPath <- file.path(thymusDir, "data", sprintf("Thymus_VisiumHD_%s_8um_hires.rds", sampleId))

# Set paths for output from this walkthrough.
outDir <- file.path(thymusDir, "SIRAP_walkthrough", sampleId)
prestoInputDir <- file.path(outDir, "01_length_filtered_spatial_fastas")
prestoOutputDir <- file.path(outDir, "02_presto_igblast")
finalAirrDir <- file.path(outDir, "03_final_spatial_airr")
finalBCRPath <- file.path(finalAirrDir, sprintf("%s_bcrs_heavy_final_spatial_airr.tsv", sampleLower))
finalTCRPath <- file.path(finalAirrDir, sprintf("%s_tcrs_beta_final_spatial_airr.tsv", sampleLower))
plotDir <- file.path(outDir, "plots")
for (dir in c(prestoInputDir, prestoOutputDir, finalAirrDir, plotDir)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# Set number of processes for parallel steps
nproc <- 8

# pRESTO/IgBLAST can take a while, so leave this FALSE until the input FASTAs
# look right. When FALSE, the script will look for db-pass TSVs from a prior run.
runPrestoIgblast <- TRUE
runCloneCalling <- TRUE
buildBCRTrees <- TRUE
savePlots <- TRUE
loadSeuratObject <- TRUE
remove_bleed_example <- TRUE
runToyRepair <- TRUE
toyRepairCondaEnv <- Sys.getenv("REPAIR_CONDA_ENV", "repair")
toyRepairCondaBin <- Sys.getenv("CONDA_EXE", "conda")
runClonalDominance <- FALSE


bcrCloneIdsToPlot <- NULL
tcrCloneIdsToPlot <- NULL



# ---- Load the spatial metadata ----------------------------------------------

all(file.exists(c(rawBCRIgblastPath, rawTCRIgblastPath, barcodeMapPath, positionsPath2um, positionsPath8um, positionsPath16um)))

barcodeMap <- readBarcodeMapping(barcodeMapPath)

positions2um <- readBarcodeMapping(positionsPath2um) %>%
  mutate(
    barcode_2um = barcode,
    x = pxl_col_in_fullres,
    y = max(pxl_row_in_fullres) - pxl_row_in_fullres
  )
positions8um <- readBarcodeMapping(positionsPath8um) %>%
  mutate(
    barcode_8um = barcode,
    x = pxl_col_in_fullres,
    y = max(pxl_row_in_fullres) - pxl_row_in_fullres
  )

head(barcodeMap)
head(positions2um)
head(positions8um)

# ---- Start from raw IgBLAST output ------------------------------------------

# Extract the 2 um barcode before pRESTO. We add the 8 um and 16 um
# barcodes after pRESTO by joining to the Visium HD barcode mapping file.

bcrRaw <- data.table::fread(rawBCRIgblastPath)
tcrRaw <- data.table::fread(rawTCRIgblastPath)

dim(bcrRaw)
dim(tcrRaw)

bcrRaw <- extractSpatialBarcode(
  bcrRaw,
  sequence_id_col = "sequence_id",
  barcode_col = "barcode_2um",
  regex = "[0-9]+_[0-9]+$",
  prefix = "s_002um_",
  suffix = "-1"
)

tcrRaw <- extractSpatialBarcode(
  tcrRaw,
  sequence_id_col = "sequence_id",
  barcode_col = "barcode_2um",
  regex = "[0-9]+_[0-9]+$",
  prefix = "s_002um_",
  suffix = "-1"
)

bcrHeavy <- processSpatialAirr(
  bcrRaw,
  loci = "IGH",
  min_v_length = 100,
  min_j_length = 10,
  min_junction_length = 1,
  barcode_col = "barcode_2um",
  add_full_id = TRUE
)
cat("Kept BCR heavy chains: ", nrow(bcrHeavy), "\n")

tcrBeta <- processSpatialAirr(
  tcrRaw,
  loci = "TRB",
  min_v_length = 100,
  min_j_length = 10,
  min_junction_length = 1,
  barcode_col = "barcode_2um",
  add_full_id = TRUE
)
cat("Kept TCR beta chains: ", nrow(tcrBeta), "\n")


# Save the processed BCRs and TCRs as airr.tsvs and FASTAs for pRESTO/IgBLAST
bcrHeavySpatialTsv <- file.path(prestoInputDir, sprintf("%s_igblast_bcr_pass_length_spatial_barcode.tsv", sampleLower))
bcrHeavySpatialFasta <- file.path(prestoInputDir, sprintf("%s_igblast_bcr_pass_length_spatial_barcode.fasta", sampleLower))
tcrBetaSpatialTsv <- file.path(prestoInputDir, sprintf("%s_igblast_tcr_pass_length_spatial_barcode.tsv", sampleLower))
tcrBetaSpatialFasta <- file.path(prestoInputDir, sprintf("%s_igblast_tcr_pass_length_spatial_barcode.fasta", sampleLower))

airr::write_rearrangement(bcrHeavy, bcrHeavySpatialTsv)
airr::write_rearrangement(tcrBeta, tcrBetaSpatialTsv)

dowser::dfToFasta(bcrHeavy, file = bcrHeavySpatialFasta, id = "full_id", seq = "sequence_alignment")
dowser::dfToFasta(tcrBeta, file = tcrBetaSpatialFasta, id = "full_id", seq = "sequence_alignment")

# Summarize pre and post filtering
spatialReceptorSummary(bcrRaw, barcode_col = "barcode_2um")
spatialReceptorSummary(bcrHeavy, barcode_col = "barcode_2um")
spatialReceptorSummary(tcrRaw, barcode_col = "barcode_2um")
spatialReceptorSummary(bcrHeavy, barcode_col = "barcode_2um")

# Plot the distribution of v, j, and junction lengths before and after filtering if requested
if (savePlots) {
  prePrestoQcPlot <- cowplot::plot_grid(
    plotLengthQc(bcrRaw, title = paste(sampleId, "raw IGH reads before filtering")),
    plotLengthQc(bcrHeavy, title = paste(sampleId, "raw IGH reads after filtering")),
    plotLengthQc(tcrRaw, title = paste(sampleId, "raw TRB reads before filtering")),
    plotLengthQc(tcrBeta, title = paste(sampleId, "raw TRB reads after filtering")),
    ncol = 2
  )
  ggsave(file.path(plotDir, sprintf("%s_pre_presto_length_qc.png", sampleLower)), prePrestoQcPlot, width = 10, height = 8, dpi = 300)
}

# Things to inspect:
head(bcrLengthSpatial)
head(tcrLengthSpatial)

# ---- Run pRESTO and IgBLAST from R ------------------------------------------

# This is the command-line part of the original workflow, wrapped in R:
# ClusterSets.py -> ParseHeaders.py -> AlignSets.py -> BuildConsensus.py ->
# AssignGenes.py -> MakeDb.py. The command output prints to the R console.

# Presto and IgBLAST settings and locations
identityThreshold <- 0.975 # How similar sequences have to be to be clustered together by pRESTO within same barcode
consensusFreqThreshold <- 0.51 # What proportion of nucleotides in a sequence cluster need to match to call a consensus at that position
germlineGlob <- "~/software/germlines/imgt/human/vdj/imgt_human_*.fasta"
germlineDir <- path.expand("~/software/germlines/imgt/human/vdj")
igblastBase <- "~/software/igblast"
igblastExec <- "~/software/igblast/bin/igblastn"
muscleExec <- "~/software/muscle3.8.31_i86linux64"

if (runPrestoIgblast) {
  # This function both runs the external presto/igblast commands and returns the paths to the key output files
  bcrPrestoPaths <- processSpatialSample(
    input_fasta = bcrHeavySpatialFasta,
    output_dir = prestoOutputDir,
    receptor = "ig",
    identity = identityThreshold,
    consensus_freq = consensusFreqThreshold,
    germline_glob = germlineGlob,
    igblast_base = igblastBase,
    igblast_exec = igblastExec,
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
    igblast_exec = igblastExec,
    muscle_exec = muscleExec,
    nproc = nproc
  )
} else {
  message("runPrestoIgblast is FALSE, so I am not running the external commands in this pass.")
}

if (!exists("bcrPrestoPaths")) {
  bcrConsensusFasta <- file.path(
    prestoOutputDir,
    sprintf(
      "%s_igblast_bcr_pass_length_spatial_barcode_ident_%s_barcode_cluster_freq_%s_cluster-pass_reheader_align-pass_consensus-pass.fasta",
      sampleLower,
      identityThreshold,
      consensusFreqThreshold
    )
  )
  bcrPrestoPaths <- list(
    consensus = bcrConsensusFasta,
    fmt7 = sub("\\.fasta$", "_igblast.fmt7", bcrConsensusFasta),
    db_pass = sub("\\.fasta$", "_igblast_db-pass.tsv", bcrConsensusFasta)
  )
}
if (!exists("tcrPrestoPaths")) {
  tcrConsensusFasta <- file.path(
    prestoOutputDir,
    sprintf(
      "%s_igblast_tcr_pass_length_spatial_barcode_ident_%s_barcode_cluster_freq_%s_cluster-pass_reheader_align-pass_consensus-pass.fasta",
      sampleLower,
      identity,
      consensusFreq
    )
  )
  tcrPrestoPaths <- list(
    consensus = tcrConsensusFasta,
    fmt7 = sub("\\.fasta$", "_igblast.fmt7", tcrConsensusFasta),
    db_pass = sub("\\.fasta$", "_igblast_db-pass.tsv", tcrConsensusFasta)
  )
}

bcrPrestoPaths$db_pass
tcrPrestoPaths$db_pass

# ---- Continue processing the presto/igblast consensus output -----------------------------

# Now we are back in R. The original script read the db-pass TSVs, kept strong
# productive consensus sequences, recovered the barcode from the pRESTO header,
# swapped the 2 um barcode row/column order, remapped to 8 um/16 um, and kept
# the top consensus sequence per 2 um barcode and locus.

bcrPrestoFullDbPassPath <- "/dartfs/rc/lab/H/HoehnK/Chinmay/022126_Visium_HD/022526_Visium_demuxed_bc24/processed/igblast_bcr_pass_length_ident_0.975_barcode_cluster_freq_0.51_cluster-pass_reheader_align-pass_consensus-pass_igblast_db-pass.tsv"
tcrPrestoFullDbPassPath <- "/dartfs/rc/lab/H/HoehnK/Chinmay/022126_Visium_HD/022526_Visium_demuxed_bc24/processed/igblast_tcr_pass_length_ident_0.975_barcode_cluster_freq_0.51_cluster-pass_reheader_align-pass_consensus-pass_igblast_db-pass.tsv"

all(file.exists(c(bcrPrestoFullDbPassPath, tcrPrestoFullDbPassPath)))

# Load in the consensus output
# Read in the full data instead of this version so that we have more sequences
if (requireNamespace("airr", quietly = TRUE)) {
  bcrConsensus <- airr::read_rearrangement(bcrPrestoFullDbPassPath)
  tcrConsensus <- airr::read_rearrangement(tcrPrestoFullDbPassPath)
} else {
  bcrConsensus <- data.table::fread(bcrPrestoFullDbPassPath)
  tcrConsensus <- data.table::fread(tcrPrestoFullDbPassPath)
}

# Get the spatial barcode back out of the sequence_id
bcrConsensus <- extractSpatialBarcode(
  bcrConsensus,
  sequence_id_col = "sequence_id",
  barcode_col = "barcode_2um",
  regex = ",([^,]+)"
)

tcrConsensus <- extractSpatialBarcode(
  tcrConsensus,
  sequence_id_col = "sequence_id",
  barcode_col = "barcode_2um",
  regex = ",([^,]+)"
)

# Filter the consensus output further based on sequence v_length, j_length, junction length, productivity, and consensus cluster count
# Also don't accept any sequences with an opening gap of more than 2
bcrConsensusFiltered <- processSpatialAirr(
  bcrConsensus,
  loci = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = 2,
  max_opening_gap = 2,
  barcode_col = "barcode_2um"
)
cat("Kept BCR consensus sequences after filtering: ", nrow(bcrConsensusFiltered), "\n")

tcrConsensusFiltered <- processSpatialAirr(
  tcrConsensus,
  loci = "TRB",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = 2,
  max_opening_gap = 2,
  barcode_col = "barcode_2um"
)
cat("Kept TCR consensus sequences after filtering: ", nrow(tcrConsensusFiltered), "\n")

# The barcodes are somehow in the wrong order out of RAD, so swap them
# Also use the barcode map to get the 8 um and 16 um barcodes for each consensus sequence
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

# Add tissue positions to the consensus sequences
bcrConsensusFilteredPos <- addTissuePositions(
  bcrConsensusFiltered,
  positions8um,
  barcode_col = "barcode_8um",
  position_barcode_col = "barcode_8um",
  x_col = "x",
  y_col = "y"
) %>% filter(!is.na(x), !is.na(y))

tcrConsensusFilteredPos <- addTissuePositions(
  tcrConsensusFiltered,
  positions8um,
  barcode_col = "barcode_8um",
  position_barcode_col = "barcode_8um",
  x_col = "x",
  y_col = "y"
) %>% filter(!is.na(x), !is.na(y))

# We want to collapse identical sequences in contiguous clusters together so as to remove bleedover from sequencing
# Do at 8um level, because at 2um level could be within same cell

if (remove_bleed_example) {
  # As an example, grab the 10 most widely expressed heavy chains and plot them in space so we can see if we have little clusters
  bcrConsensusBigSeqs <- bcrConsensusFilteredPos %>%
    group_by(sequence_alignment) %>%
    summarize(n_bc = n(), total_consensus_count = sum(consensus_count)) %>%
    arrange(desc(n_bc), desc(total_consensus_count)) %>%
    slice_head(n = 10) %>%
    pull(sequence_alignment)

  # For each of these sequences, highlight the locations of all the consensus sequences that are identical to it
  plots <- list()
  for (seq in bcrConsensusBigSeqs) {
    seqMatches <- bcrConsensusFilteredPos %>%
      filter(sequence_alignment == seq)
    seqName <- paste0("seq_", which(bcrConsensusBigSeqs == seq))

    if (!is.null(thymus)) {
      seqPlot <- plotReceptorExpression(
        object = thymus,
        receptors = seqMatches,
        barcode_col = "barcode_8um",
        receptor_color = "#2c7fb8",
        receptor_size = 1.2,
        title = paste("Identical sequence:", seq)
      )
      plots[[seqName]] <- seqPlot
      ggsave(file.path(plotDir, sprintf("%s_identical_sequence_%s_bcr_spatial.png", sampleLower, seqName)), seqPlot, width = 6, height = 5, dpi = 300)
      cat("Plotted identical sequence ", seqName, "\n")
    }
  }
  # Try filtering spatially
  # Goal is to collapse identical sequences in contiguous clusters onto the most present location
  # Get the adjacency matrix 

  bcrConsensusFilteredPosSpatialCollapse <- collapseIdenticalNeighbors(
    data = bcrConsensusFilteredPos,
    tissue_positions = positions8um,
    barcode = "barcode_8um",
    x = "x",
    y = "y",
    sequence = "sequence_alignment",
    consensus_count = "consensus_count",
    num_neighbors_non_boundary = 8L
  )

  # Check bcrConsensusBigSeqs again but from this set
  plotsCollapsed <- list()
  for (seq in bcrConsensusBigSeqs) {
    seqMatches <- bcrConsensusFilteredPosSpatialCollapse %>%
      filter(sequence_alignment == seq)
    seqName <- paste0("seq_", which(bcrConsensusBigSeqs == seq))
    
    if (!is.null(thymus)) {
      seqPlot <- plotReceptorExpression(
        object = thymus,
        receptors = seqMatches,
        barcode_col = "barcode_8um",
        receptor_color = "#2c7fb8",
        receptor_size = 1.2,
        title = paste("Identical sequence:", seq)
      )
      plots[[seqName]] <- seqPlot
      ggsave(file.path(plotDir, sprintf("%s_identical_sequence_%s_bcr_spatial_collapsed.png", sampleLower, seqName)), seqPlot, width = 6, height = 5, dpi = 300)
      cat("Plotted identical sequence ", seqName, "\n")
    }
  }
}

bcrConsensusFilteredPosSpatialCollapse <- collapseIdenticalNeighbors(
  data = bcrConsensusFilteredPos,
  tissue_positions = positions8um,
  barcode = "barcode_8um",
  x = "x",
  y = "y",
  sequence = "sequence_alignment",
  consensus_count = "consensus_count",
  num_neighbors_non_boundary = 8L
)

tcrConsensusFilteredPosSpatialCollapse <- collapseIdenticalNeighbors(
  data = tcrConsensusFilteredPos,
  tissue_positions = positions8um,
  barcode = "barcode_8um",
  x = "x",
  y = "y",
  sequence = "sequence_alignment",
  consensus_count = "consensus_count",
  num_neighbors_non_boundary = 8L
)


# Get one consensus sequence per barcode (most present in barcode) at 2um level
bcrTopBarcode <- collapseTopConsensusPerBarcode(bcrConsensusFilteredPosSpatialCollapse, barcode_col = "barcode_2um", locus_col = "locus", ties = "all")
tcrTopBarcode <- collapseTopConsensusPerBarcode(tcrConsensusFilteredPosSpatialCollapse, barcode_col = "barcode_2um", locus_col = "locus", ties = "all")


if (requireNamespace("airr", quietly = TRUE)) {
  airr::write_rearrangement(bcrTopBarcode, finalBCRPath)
  airr::write_rearrangement(tcrTopBarcode, finalTCRPath)
} else {
  data.table::fwrite(bcrTopBarcode, finalBCRPath, sep = "\t")
  data.table::fwrite(tcrTopBarcode, finalTCRPath, sep = "\t")
}

list(
  bcr = spatialReceptorSummary(bcrTopBarcode, barcode_col = "barcode_2um"),
  tcr = spatialReceptorSummary(tcrTopBarcode, barcode_col = "barcode_2um")
)

postPrestoQcPlot <- cowplot::plot_grid(
  plotLengthQc(bcrTopBarcode, title = paste(sampleId, "final BCRs")),
  plotLengthQc(tcrTopBarcode, title = paste(sampleId, "final TCRs")),
  ncol = 1
)
postPrestoQcPlot
if (savePlots) {
  ggsave(file.path(plotDir, sprintf("%s_post_presto_final_length_qc.png", sampleLower)), postPrestoQcPlot, width = 10, height = 8, dpi = 300)
}

# Things to inspect:
head(bcrTopBarcode)
head(tcrTopBarcode)

# ---- Put every BCR and TCR on the tissue slice -------------------------------
# Load thymus data if requested
if (requireNamespace("airr", quietly = TRUE)){
  bcrTopBarcode <- airr::read_rearrangement(finalBCRPath)
  tcrTopBarcode <- airr::read_rearrangement(finalTCRPath)
} else {
  bcrTopBarcode <- data.table::fread(finalBCRPath) %>% as.data.frame()
  tcrTopBarcode <- data.table::fread(finalTCRPath) %>% as.data.frame()
}


thymus <- NULL
if (loadSeuratObject && file.exists(seuratPath)) {
  library(Seurat)
  library(SeuratObject)
  thymus <- readRDS(seuratPath)
}

if (!is.null(thymus)) {
  allBCRPlot <- plotReceptorExpression(
    object = thymus,
    receptors = bcrTopBarcode,
    barcode_col = "barcode_8um",
    receptor_color = "white",
    receptor_size = 0.5,
    title = paste(sampleId, "all BCR-positive spots")
  )
  if (savePlots) {
    ggsave(file.path(plotDir, sprintf("%s_all_bcr_locations.png", sampleLower)), allBCRPlot, width = 6, height = 5, dpi = 300)
  }

  allTCRPlot <- plotReceptorExpression(
    object = thymus,
    receptors = tcrTopBarcode,
    barcode_col = "barcode_8um",
    receptor_color = "#d95f02",
    receptor_size = 0.5,
    title = paste(sampleId, "all TCR-positive spots")
  )
  if (savePlots) {
    ggsave(file.path(plotDir, sprintf("%s_all_tcr_locations.png", sampleLower)), allTCRPlot, width = 6, height = 5, dpi = 300)
  }

  bcr_spots <- unique(bcrTopBarcode$barcode_8um)
  tcr_spots <- unique(tcrTopBarcode$barcode_8um)
  receptor_type <- dplyr::case_when(
    colnames(thymus) %in% bcr_spots & colnames(thymus) %in% tcr_spots ~ "BCR + TCR",
    colnames(thymus) %in% bcr_spots ~ "BCR",
    colnames(thymus) %in% tcr_spots ~ "TCR",
    TRUE ~ "0"
  )
  thymus[["sirap_receptor_type"]] <- receptor_type
  bcrTCROverlay <- spatialDiscretePlot(
    thymus,
    meta_col = "sirap_receptor_type",
    zero_value = "0",
    image.scale = "hires",
    point_size = 0.5,
    cols = c("BCR" = "#2c7fb8", "TCR" = "#d95f02", "BCR + TCR" = "white")
  ) +
    ggtitle(paste(sampleId, "BCR and TCR locations"))
  if (savePlots) {
    ggsave(file.path(plotDir, sprintf("%s_bcr_tcr_locations.png", sampleLower)), bcrTCROverlay, width = 6, height = 5, dpi = 300)
  }
} else {
  message("Skipping receptor spatial plots because the Seurat object was not loaded.")
}


# ---- Add gene expression underneath receptor locations -----------------------

geneOverlayPlots <- list()
genesToShow <- c("CD79A", "CD3E")

if (!is.null(thymus)) {
  for (gene in genesToShow) {
    if (!gene %in% rownames(thymus)) {
      message("Skipping ", gene, ": not found in the Seurat object.")
      next
    }

    bcrGenePlot <- plotReceptorExpression(
      object = thymus,
      receptors = bcrTopBarcode,
      gene = gene,
      barcode_col = "barcode_8um",
      receptor_color = "#2c7fb8",
      receptor_size = 0.5,
      receptor_stroke = 0.25,
      receptor_alpha = 0.8,
      feature_size = 0.5,
      feature_alpha = 0.8,
      title = paste(gene, "expression with BCR spots")
    ) 

    tcrGenePlot <- plotReceptorExpression(
      object = thymus,
      receptors = tcrTopBarcode,
      gene = gene,
      barcode_col = "barcode_8um",
      receptor_color = "#d95f02",
      receptor_size = 0.5,
      receptor_stroke = 0.25,
      receptor_alpha = 0.8,
      feature_size = 0.5,
      feature_alpha = 0.8,
      title = paste(gene, "expression with TCR spots")
    )

    geneOverlayPlots[[gene]] <- cowplot::plot_grid(bcrGenePlot, tcrGenePlot, ncol = 2)
    if (savePlots) {
      ggsave(
        file.path(plotDir, sprintf("%s_%s_receptor_gene_overlay.png", sampleLower, gene)),
        geneOverlayPlots[[gene]],
        width = 11,
        height = 5,
        dpi = 300
      )
    }
  }
} else {
  message("Skipping Seurat gene overlays. Check loadSeuratObject and seuratPath.")
}

if (!is.null(thymus)) {
  blendedGenePlot <- spatialExpressionBlendPlot(
    thymus,
    features = c("CXCR4", "CXCR5"),
    assay = "Spatial.008um",
    layer = "counts",
    image.scale = "hires",
    point_size = 0.2,
    only_blended = TRUE,
    top_left = "blue",
    bottom_left = "grey50",
    bottom_right = "yellow",
    top_right = "green"
  ) +
    ggtitle(paste(sampleId, "CXCR4 + CXCR5"))
  print(blendedGenePlot)
  if (savePlots) {
    ggsave(file.path(plotDir, sprintf("%s_CXCR4_CXCR5_blended_no_zeros.png", sampleLower)), blendedGenePlot, width = 6, height = 5, dpi = 300)
  }
}

# ---- Spatial autocorrelation and colocalization ------------------------------

spatialStatsRadius <- 200
spatialStatsK <- 24 # Two loops around spot
ripleyNsim <- 10 

if (!is.null(thymus)) {
  bcrMoran <- spatialMoran(
    thymus,
    variable = "BCR presence",
    variable_type = "receptor",
    receptors = bcrTopBarcode,
    barcode_col = "barcode_8um",
    # radius = spatialStatsRadius # Unadvisable to use the radius for visium HD data since nSpots is so high
    k = spatialStatsK # Instead use k-nearest neighbors
  )
  bcrMoran$test

  if ("CD79A" %in% rownames(thymus)) {
    cd79aExpressionMoran <- spatialMoran(
      thymus,
      variable = "CD79A",
      variable_type = "feature",
      layer = "counts",
      k = spatialStatsK
    )
    cd79aExpressionMoran$test

    cd79aPositiveMoran <- spatialMoran(
      thymus,
      variable = "CD79A",
      variable_type = "feature",
      positive = TRUE,
      layer = "counts",
      k = spatialStatsK
    )
    cd79aPositiveMoran$test

    aicdaLocalMoran <- spatialLocalMoran( # slower than global Moran's I but gives spot-level values we can plot
      thymus,
      variable = "AICDA",
      variable_type = "feature",
      layer = "counts",
      k = spatialStatsK
    )
    aicdaLocalMoranPlot <- plotSpatialLocalMoran(
      thymus,
      aicdaLocalMoran,
      point_size = 0.55,
      title = paste(sampleId, "AICDA local Moran's I z-score")
    )
    
    if (savePlots) {
      ggsave(file.path(plotDir, sprintf("%s_AICDA_local_moran.png", sampleLower)), aicdaLocalMoranPlot, width = 6, height = 5, dpi = 300)
    }
  }

  if ("seurat_clusters" %in% colnames(thymus@meta.data)) {
    clusterToTest <- as.character(thymus@meta.data$seurat_clusters[!is.na(thymus@meta.data$seurat_clusters)][1])
    clusterMoran <- spatialMoran(
      thymus,
      variable = "seurat_clusters",
      variable_type = "metadata",
      value = clusterToTest,
      radius = spatialStatsRadius
    )
    clusterMoran$test
  }

  bcrTcrColocalizationPlot <- spatialDiscretePlot(
    thymus,
    meta_col = "sirap_receptor_type",
    zero_value = "0",
    cols = c("BCR" = "#2c7fb8", "TCR" = "#d95f02", "BCR + TCR" = "white"),
    point_size = 0.5
  ) +
    ggtitle(paste(sampleId, "BCR, TCR, and dual-positive spots"))
  if (savePlots) {
    ggsave(file.path(plotDir, sprintf("%s_bcr_tcr_k_marks.png", sampleLower)), bcrTcrColocalizationPlot, width = 6, height = 5, dpi = 300)
  }

  bcrTcrK <- spatialRipleyCrossK(
    thymus,
    mark_col = "sirap_receptor_type",
    i = "BCR",
    j = "TCR",
    nsim = ripleyNsim,
    max_r = 1000
  )
  bcrTcrKPlot <- plotRipleyCrossK(bcrTcrK, title = paste(sampleId, "Ripley's cross K: BCR to TCR"))
  if (savePlots) {
      ggsave(file.path(plotDir, sprintf("%s_bcr_tcr_ripley_cross_k.png", sampleLower)), bcrTcrKPlot, width = 6, height = 4, dpi = 300)
  }



W_thy2 <- owin(xrange = c(900, 11100), yrange = c(700, 11000))
spatial_data_thy2 <- spatial_data %>% filter(patient == "m5") %>% mutate(cell_type = factor(cell_type))
thy2K <- spatialRipleyCrossK(
  thy_merged_harmonized_spatial_only,
  mark_col = "cell_type",
  i = "GC B",
  j = "CHRNA1+ mTEC",
  nsim = 99,
  max_r = 2500
)
thy2KPlot <- plotRipleyCrossK(thy2K, title = paste(sampleId, "Ripley's cross K: GC B to CHRNA1+ mTEC"))
ggsave(file.path(plotDir, sprintf("%s_thy2_gc_b_chrna1_mtec_ripley_cross_k.png", sampleLower)), thy2KPlot, width = 6, height = 4, dpi = 300)

kcross_data_thy2_ppp <- ppp(x = spatial_data_thy2$y, y = -spatial_data_thy2$x, marks = spatial_data_thy2$cell_type, window = W_thy2)
kcross_thy2 <- Kcross(kcross_data_thy2_ppp, i = "GC B", j = "CHRNA1+ mTEC")

# Test this against a simulated version
shuffle_thy2 <- expression(rlabel(kcross_data_thy2_ppp))
kcross_thy2_sim <- envelope(kcross_data_thy2_ppp, Kcross, nsim = 99, simulate = shuffle_thy2, i = "GC B", j = "CHRNA1+ mTEC")

# Check against other TECs as a comparison
kcross_thy2_mtecs <- Kcross(kcross_data_thy2_ppp, i = "GC B", j = "mTEC")
shuffle_thy2_mtecs <- expression(rlabel(kcross_data_thy2_ppp))
kcross_thy2_mtecs_sim <- envelope(kcross_data_thy2_ppp, Kcross, nsim = 99, simulate = shuffle_thy2_mtecs, i = "GC B", j = "mTEC")

kcross_thy2_ctecs <- Kcross(kcross_data_thy2_ppp, i = "GC B", j = "cTEC")
shuffle_thy2_ctecs <- expression(rlabel(kcross_data_thy2_ppp))
kcross_thy2_ctecs_sim <- envelope(kcross_data_thy2_ppp, Kcross, nsim = 99, simulate = shuffle_thy2_ctecs, i = "GC B", j = "cTEC")

pdf("GEX_plots/kcross_GCB_CHRNA1+mTECs_thy2.pdf")
plot(kcross_thy2)
plot(kcross_thy2_sim)
plot(kcross_thy2_mtecs_sim)
plot(kcross_thy2_ctecs_sim)
dev.off()


  if ("BCR + TCR" %in% thymus$sirap_receptor_type) {
    bcrDualPositiveK <- spatialRipleyCrossK(
      thymus,
      mark_col = "sirap_receptor_type",
      i = "BCR",
      j = "BCR + TCR",
      nsim = ripleyNsim
    )
    bcrDualPositiveKPlot <- plotRipleyCrossK(bcrDualPositiveK, title = paste(sampleId, "Ripley's cross K: BCR to BCR + TCR"))
    if (savePlots) {
      ggsave(file.path(plotDir, sprintf("%s_bcr_dual_positive_ripley_cross_k.png", sampleLower)), bcrDualPositiveKPlot, width = 6, height = 4, dpi = 300)
    }
  }
}

# ---- Call clones -------------------------------------------------------------

bcrClonesPath <- file.path(outDir, sprintf("%s_bcr_hierarchical_clones.rds", sampleLower))
tcrClonesPath <- file.path(outDir, sprintf("%s_tcr_hierarchical_clones.rds", sampleLower))
bcrTreesPath <- file.path(outDir, sprintf("%s_bcr_pml_trees.rds", sampleLower))

buildClones <- FALSE # Set to TRUE to build clones, which can take a while. When FALSE, will look for RDS files from a prior run.
if (buildClones) {
  if (!requireNamespace("scoper", quietly = TRUE)) {
    stop("Install scoper to call clones.")
  }

  bcrClones <- scoper::hierarchicalClones(
    filter(bcrTopBarcode, productive),
    threshold = 0.1,
    cell_id = "sequence_id",
    summarize_clones = FALSE,
    nproc = nproc
  )

  tcrClones <- scoper::identicalClones(
    filter(tcrTopBarcode, productive),
    cell_id = "sequence_id",
    summarize_clones = FALSE,
    nproc = nproc
  )

  saveRDS(bcrClones, bcrClonesPath)
  saveRDS(tcrClones, tcrClonesPath)
} else if (file.exists(bcrClonesPath) && file.exists(tcrClonesPath)) {
  bcrClones <- readRDS(bcrClonesPath)
  tcrClones <- readRDS(tcrClonesPath)
} else {
  stop("Clone RDS files not found. Set runCloneCalling <- TRUE and rerun this section.")
}

bcrCloneCounts <- bcrClones %>%
  filter(!is.na(clone_id)) %>%
  count(clone_id, sort = TRUE)
head(bcrCloneCounts, 20)

tcrCloneCounts <- tcrClones %>%
  filter(!is.na(clone_id)) %>%
  count(clone_id, sort = TRUE)
head(tcrCloneCounts, 20)

# ---- Pick a couple of clones and plot them ----------------------------------
# Grab a few clones with 30-50 sequences
bcrCloneIds <- bcrCloneCounts %>%
  filter(n >= 30, n <= 50) %>%
  slice_head(n = 2) %>%
  pull(clone_id)
tcrCloneIds <- tcrCloneCounts %>%
  filter(n >= 30, n <= 50) %>%
  slice_head(n = 2) %>%
  pull(clone_id)


if (!is.null(thymus)) {
  bcrClonePlots <- lapply(
    bcrCloneIds,
    function(cloneId) {
      plotReceptorExpression(
        object = thymus,
        receptors = filter(bcrClones, .data$clone_id == .env$cloneId),
        barcode_col = "barcode_8um",
        receptor_color = "#2c7fb8",
        receptor_size = 1.7,
        title = paste("BCR clone", cloneId)
      )
    }
  )
  selectedBCRClonePlot <- cowplot::plot_grid(plotlist = bcrClonePlots, ncol = 2)
  selectedBCRClonePlot
  if (savePlots) {
    ggsave(file.path(plotDir, sprintf("%s_selected_bcr_clones_spatial.png", sampleLower)), selectedBCRClonePlot, width = 10, height = 5, dpi = 300)
  }
} else {
  message("Skipping BCR clone plots because the Seurat object was not loaded.")
}

if (!is.null(thymus)) {
  tcrClonePlots <- lapply(
    tcrCloneIds,
    function(cloneId) {
      plotReceptorExpression(
        object = thymus,
        receptors = filter(tcrClones, .data$clone_id == .env$cloneId),
        barcode_col = "barcode_8um",
        receptor_color = "#d95f02",
        receptor_size = 1.7,
        title = paste("TCR clone", cloneId)
      )
    }
  )
  selectedTCRClonePlot <- cowplot::plot_grid(plotlist = tcrClonePlots, ncol = 2)
  selectedTCRClonePlot
  if (savePlots) {
    ggsave(file.path(plotDir, sprintf("%s_selected_tcr_clones_spatial.png", sampleLower)), selectedTCRClonePlot, width = 10, height = 5, dpi = 300)
  }
} else {
  message("Skipping TCR clone plots because the Seurat object was not loaded.")
}

# ---- Build BCR trees and match tips back to space ----------------------------
buildBCRTrees <- FALSE
if (buildBCRTrees) {
  if (!requireNamespace("dowser", quietly = TRUE)) {
    stop("Install dowser to build BCR lineage trees.")
  }
  if (!requireNamespace("alakazam", quietly = TRUE)) {
    stop("Install alakazam to read IMGT germlines for dowser tree building.")
  }

  refs <- dowser::readIMGT(germlineDir)

  bcrTreeInput <- bcrClones %>%
    filter(clone_id %in% bcrCloneIds, !is.na(clone_id))

  bcrTreeInput <- dowser::createGermlines(bcrTreeInput, refs, id = "sequence_id")

  bcrFormattedClones <- dowser::formatClones(
    bcrTreeInput,
    id = "sequence_id",
    germ = "germline_alignment",
    traits = c("barcode_8um"),
    columns = c("x", "y", "barcode_8um"),
    filterstop = FALSE,
    collapse = TRUE,
    add_count = "consensus_count" %in% colnames(bcrTreeInput),
    num_fields = if ("consensus_count" %in% colnames(bcrTreeInput)) "consensus_count" else NULL
  )
  bcrTrees <- dowser::getTrees(bcrFormattedClones, nproc = nproc)

  saveRDS(bcrTrees, bcrTreesPath)
} else if (file.exists(bcrTreesPath)) {
  bcrTrees <- readRDS(bcrTreesPath)
} else {
  stop("BCR tree RDS file not found. Set buildBCRTrees <- TRUE and rerun this section.")
}

# ---- Cluster BCR tree tips into genetically similar subtrees -----------------

# This cuts the k - 1 longest branches in each tree, which gives k connected
# subtrees. We can then color the tree and spatial plot by subtree instead of
# giving every tip a unique color.
treeTipClusterK <- 6

bcrTreesClustered <- clusterTreesTips(
  bcrTrees,
  k = treeTipClusterK,
  id_col = "sequence_id",
  cluster_col = "tree_cluster"
)
bcrTreesClustered$tree_cluster_k

if (is.null(thymus)) {
  message("Skipping tree/spatial pairs because the Seurat object was not loaded.")
} else {
  pairedTreeSpatialPlots <- plotSpatialCloneTreePairs(
    bcrTrees,
    object = thymus,
    barcode_col = "barcode_8um",
    clone_ids = bcrCloneIds,
    id_col = "sequence_id",
    treescale = 0.01,
    zoom = FALSE
  )

  pairedTreeSpatialClusterPlots <- plotSpatialCloneTreePairs(
    bcrTreesClustered,
    object = thymus,
    barcode_col = "barcode_8um",
    clone_ids = bcrCloneIds,
    id_col = "sequence_id",
    color_by = "tree_cluster",
    treescale = 0.01,
    zoom = FALSE
  )
}

if (savePlots && exists("pairedTreeSpatialPlots")) {
  spatialCloneTreePairsToPDF(
    pairedTreeSpatialPlots,
    file.path(plotDir, sprintf("%s_bcr_tree_spatial_pairs.pdf", sampleLower)),
    width = 12,
    height = 6
  )
}

if (savePlots && exists("pairedTreeSpatialClusterPlots")) {
  spatialCloneTreePairsToPDF(
    pairedTreeSpatialClusterPlots,
    file.path(plotDir, sprintf("%s_bcr_tree_spatial_pairs_tree_clusters.pdf", sampleLower)),
    width = 12,
    height = 6
  )
}


# ---- Toy SPaCeD heavy/light pairing example ---------------------------------

# This small example is intentionally explicit. Heavy clone 1 is near light
# clone 1, and heavy clone 2 is near light clone 2, even though the paired
# clones do not occupy exactly the same barcode. For Visium HD, use 8um
# barcodes as these native spatial bins rather than adding another
# coordinate-binning step in SIRAP.
toyPositions <- data.frame(
  barcode = paste0("toy_", 1:6),
  x = c(0, 2, 4, 100, 102, 104),
  y = c(0, 1, 0, 100, 101, 100)
)

toyHeavyClones <- data.frame(
  barcode = c("toy_1", "toy_2", "toy_3", "toy_4", "toy_5"),
  clone_id = c("1", "1", "1", "2", "2")
)

toyLightClones <- data.frame(
  barcode = c("toy_2", "toy_3", "toy_5", "toy_6"),
  clone_id = c("1", "1", "2", "2")
)

toyHeavyPoints <- buildSpacedCloneTable(
  toyHeavyClones,
  positions = toyPositions,
  clone_prefix = "heavy_"
)

toyLightPoints <- buildSpacedCloneTable(
  toyLightClones,
  positions = toyPositions,
  clone_prefix = "light_"
)

toySpatialAffinity <- spacedSpatialAffinity(
  heavy_points = toyHeavyPoints,
  light_points = toyLightPoints,
  min_x_obs = 1
)
toySpatialAffinity

toyRepairInput <- buildRepairInputMatrix(toyHeavyPoints, toyLightPoints)
toyRepairInput

toyRepairDir <- file.path(outDir, "toy_spaced_repair")
toyRepairInputPath <- file.path(toyRepairDir, "toy_spaced_repair_input.tsv")
toyRepairOutputDir <- file.path(toyRepairDir, "repair_out")
dir.create(toyRepairDir, recursive = TRUE, showWarnings = FALSE)
dir.create(toyRepairOutputDir, recursive = TRUE, showWarnings = FALSE)
writeRepairInputMatrix(toyRepairInput, toyRepairInputPath)

toyRepairCommand <- repairAnalyzeCommand(
  input_path = toyRepairInputPath,
  output_dir = toyRepairOutputDir,
  a_pattern = "^heavy",
  b_pattern = "^light",
  mxo = 1
)
toyRepairCommand

toyRepairCondaCommand <- paste(
  shQuote(toyRepairCondaBin),
  "run",
  "-n",
  shQuote(toyRepairCondaEnv),
  toyRepairCommand
)
toyRepairCondaCommand

if (runToyRepair) {
  checkSpacedRepair(conda_env = toyRepairCondaEnv, conda_exe = toyRepairCondaBin)
}

if (runToyRepair) {
  toyRepairStatus <- system(toyRepairCondaCommand)
  if (toyRepairStatus != 0) {
    stop("Toy REPAIR command failed with exit status ", toyRepairStatus, ": ", toyRepairCondaCommand)
  }
}

toyRepairModelFiles <- list.files(
  toyRepairOutputDir,
  pattern = "_model_matrix\\.tsv$",
  full.names = TRUE
)
if (length(toyRepairModelFiles) != 1L) {
  stop(
    "Expected exactly one REPAIR model matrix in ", toyRepairOutputDir,
    ". Found: ", paste(basename(toyRepairModelFiles), collapse = ", "),
    ". Make sure `", toyRepairCondaEnv, "` is the modified REPAIR conda environment used by SPaCeD."
  )
}
toySpacedPairs <- pairHeavyLightSpaced(
  spatial_affinity = toySpatialAffinity,
  mapping_matrix = toyRepairModelFiles[[1]],
  weights = c(0, 0.5, 1),
  top_n_light_chains = 5
)
toySpacedPairs


# ---- Spatial clonal dominance inside selected polygons ----------------------

# This section draws or loads polygon vertices, calculates BCR clone dominance
# inside each polygon, and writes the same diagnostic plots used in the M5
# clonal dominance inspection script. It is off by default because selecting
# polygons can launch a Shiny app.
if (runClonalDominance) {
  if (is.null(thymus)) {
    stop("Load the Seurat object before running clonal dominance.")
  }
  if (!exists("bcrClones")) {
    stop("Run or load BCR clone calling before running clonal dominance.")
  }

  dominanceOutDir <- file.path(plotDir, "clonal_dominance")
  dir.create(dominanceOutDir, recursive = TRUE, showWarnings = FALSE)

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

  # Edit these directly when running interactively on the HPC. Use host
  # 0.0.0.0 when tunneling from the login node to an interactive compute node.
  dominanceShinyHost <- "0.0.0.0"
  dominanceShinyPort <- 8788L
  dominanceFeature <- "MS4A1"
  dominanceDownsampleN <- 20000
  dominancePlotImage <- TRUE

  if (file.exists(gcMantlePolygonPath)) {
    gcMantlePolygon <- fread(gcMantlePolygonPath) %>% as.data.frame()
  } else {
    message(
      "No saved GC plus mantle polygon found. Starting Shiny on ",
      dominanceShinyHost, ":", dominanceShinyPort, "."
    )
    message("From your laptop, tunnel this port and open http://127.0.0.1:", dominanceShinyPort)

    gcMantlePolygon <- selectSpatialPolygons(
      object = thymus,
      feature = dominanceFeature,
      image.scale = "hires",
      save_path = gcMantlePolygonPath,
      downsample_n = dominanceDownsampleN,
      host = dominanceShinyHost,
      port = dominanceShinyPort,
      plot_image = dominancePlotImage,
      launch.browser = FALSE
    )
  }

  if (file.exists(gcPolygonPath)) {
    gcPolygon <- fread(gcPolygonPath) %>% as.data.frame()
  } else {
    message(
      "No saved GC polygon found. Starting Shiny on ",
      dominanceShinyHost, ":", dominanceShinyPort, "."
    )
    message("From your laptop, tunnel this port and open http://127.0.0.1:", dominanceShinyPort)

    gcPolygon <- selectSpatialPolygons(
      object = thymus,
      feature = dominanceFeature,
      image.scale = "hires",
      save_path = gcPolygonPath,
      downsample_n = dominanceDownsampleN,
      host = dominanceShinyHost,
      port = dominanceShinyPort,
      plot_image = dominancePlotImage,
      launch.browser = FALSE
    )
  }

  if (file.exists(allGcsPolygonsPath)) {
    allGcsPolygons <- fread(allGcsPolygonsPath) %>% as.data.frame()
  } else {
    message(
      "No saved all-GC polygons found. Starting Shiny on ",
      dominanceShinyHost, ":", dominanceShinyPort, "."
    )
    message("From your laptop, tunnel this port and open http://127.0.0.1:", dominanceShinyPort)

    allGcsPolygons <- selectSpatialPolygons(
      object = thymus,
      feature = dominanceFeature,
      image.scale = "hires",
      save_path = allGcsPolygonsPath,
      downsample_n = dominanceDownsampleN,
      host = dominanceShinyHost,
      port = dominanceShinyPort,
      plot_image = dominancePlotImage,
      launch.browser = FALSE
    )
  }

  allGcsPolygonOverlay <- plotSpatialPolygons(
    object = thymus,
    polygon = allGcsPolygons,
    image.scale = "hires",
    polygon_alpha = 0.2
  )
  ggsave(
    filename = file.path(dominanceOutDir, sprintf("%s_all_gcs_polygons_overlay.pdf", sampleLower)),
    plot = allGcsPolygonOverlay,
    width = 7,
    height = 6
  )

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

  fwrite(
    clonalDominance$dominance,
    file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_clonal_dominance.tsv", sampleLower)),
    sep = "\t"
  )
  fwrite(
    clonalDominanceNoSingletons$dominance,
    file.path(dominanceOutDir, sprintf("%s_all_gcs_heavy_clonal_dominance_no_singletons.tsv", sampleLower)),
    sep = "\t"
  )

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
}




# ---- Done -------------------------------------------------------------------
