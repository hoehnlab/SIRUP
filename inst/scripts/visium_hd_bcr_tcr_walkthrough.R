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

# Notes while testing:
# - TODO: Let's assume people have the airr package
# - TODO:
# - TODO:

# ---- Things to change --------------------------------------------------------

# Settings
sampleId <- "M5" # M1 uses demuxed_bc23; M5 uses demuxed_bc24.
sampleLower <- tolower(sampleId)

# Directory and paths for data
thymusDir <- normalizePath("../../Thymus_VisiumHD", mustWork = TRUE)
demuxBySample <- c(M1 = "022526_Visium_demuxed_bc23", M5 = "022526_Visium_demuxed_bc24") # Names of demuxed samples
demuxDir <- file.path("/dartfs/rc/lab/H/HoehnK/Chinmay/022126_Visium_HD", demuxBySample[[sampleId]])

# Using 1M example files from RAD output for speed
rawBcrIgblastPath <- file.path(demuxDir, "igblast_out", "igblast_bcr_1M_example.tsv")
rawTcrIgblastPath <- file.path(demuxDir, "igblast_out", "igblast_tcr_1M_example.tsv")
visiumOutsDir <- file.path(thymusDir, "analysis", sprintf("%s_260129", sampleId), "outs")

# Either specify the parquet file or, as we've done here, a csv that is saved from the parquet so we don't have to use arrow.
barcodeMapPath <- file.path(thymusDir, "M5_barcode_mappings.csv")
positionsPath2um <- file.path(thymusDir, sprintf("%s_tissue_positions_2um.csv", sampleId))
positionsPath8um <- file.path(thymusDir, sprintf("%s_tissue_positions_8um.csv", sampleId))
positionsPath16um <- file.path(thymusDir, sprintf("%s_tissue_positions_16um.csv", sampleId))
seuratPath <- file.path(thymusDir, "data", sprintf("Thymus_VisiumHD_%s_8um_hires_sr_processed.rds", sampleId))

# Set paths for output from this walkthrough.
outDir <- file.path(thymusDir, "SIRAP_walkthrough", sampleId)
prestoInputDir <- file.path(outDir, "01_length_filtered_spatial_fastas")
prestoOutputDir <- file.path(outDir, "02_presto_igblast")
finalAirrDir <- file.path(outDir, "03_final_spatial_airr")
plotDir <- file.path(outDir, "plots")
for (dir in c(prestoInputDir, prestoOutputDir, finalAirrDir, plotDir)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# Set number of processes for parallel steps
nproc <- 4

# pRESTO/IgBLAST can take a while, so leave this FALSE until the input FASTAs
# look right. When FALSE, the script will look for db-pass TSVs from a prior run.
runPrestoIgblast <- TRUE
runCloneCalling <- TRUE
buildBcrTrees <- TRUE
savePlots <- TRUE
loadSeuratObject <- TRUE
remove_bleed_example <- TRUE


bcrCloneIdsToPlot <- NULL
tcrCloneIdsToPlot <- NULL



# ---- Load the spatial metadata ----------------------------------------------

all(file.exists(c(rawBcrIgblastPath, rawTcrIgblastPath, barcodeMapPath, positionsPath2um, positionsPath8um, positionsPath16um)))

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

bcrRaw <- data.table::fread(rawBcrIgblastPath)
tcrRaw <- data.table::fread(rawTcrIgblastPath)

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

# Notes:
# - TODO: Should processSpatialAirr() print a drop-count summary?
# Yes
# - TODO: Should full_id include barcode16um too?
#  No, I think we only include other barcodes after presto and igblast because only visiumHD has multiple of our datasets
# TODO: Change name from runPrestoIgblastSample to processSpatialSample

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
  # Try dry_run = TRUE first to check the commands and paths without actually running
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

# Notes:
# - TODO: Should runPrestoIgblastSample() check all external executables first?
# No, but it should print helpful error messages if they are not found when it tries to run them.

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

# TODO: @Codex update this section and clean it up
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

    seqPlot <- plotReceptorOverlay(
      tissue_positions = positions8um,
      receptors = seqMatches,
      receptor_color = "#2c7fb8",
      receptor_size = 1.2,
      title = paste("Identical sequence:", seq)
    )
    plots[[seqName]] <- seqPlot
    ggsave(file.path(plotDir, sprintf("%s_identical_sequence_%s_bcr_spatial.png", sampleLower, seqName)), seqPlot, width = 6, height = 5, dpi = 300)
    cat("Plotted identical sequence ", seqName, "\n")
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
    
    seqPlot <- plotReceptorOverlay(
      tissue_positions = positions8um,
      receptors = seqMatches,
      receptor_color = "#2c7fb8",
      receptor_size = 1.2,
      title = paste("Identical sequence:", seq)
    )
    plots[[seqName]] <- seqPlot
    ggsave(file.path(plotDir, sprintf("%s_identical_sequence_%s_bcr_spatial_collapsed.png", sampleLower, seqName)), seqPlot, width = 6, height = 5, dpi = 300)
    cat("Plotted identical sequence ", seqName, "\n")
  }
}

# TODO: @Codex move collapseIdenticalNeighbors to thhe processing.R file
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


finalBcrPath <- file.path(finalAirrDir, sprintf("%s_bcrs_heavy_final_spatial_airr.tsv", sampleLower))
finalTcrPath <- file.path(finalAirrDir, sprintf("%s_tcrs_beta_final_spatial_airr.tsv", sampleLower))
if (requireNamespace("airr", quietly = TRUE)) {
  airr::write_rearrangement(bcrTopBarcode, finalBcrPath)
  airr::write_rearrangement(tcrTopBarcode, finalTcrPath)
} else {
  data.table::fwrite(bcrTopBarcode, finalBcrPath, sep = "\t")
  data.table::fwrite(tcrTopBarcode, finalTcrPath, sep = "\t")
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

# Notes:
# - TODO: Should the post-pRESTO block become a single high-level function?
# Answer: No, I think those are all different enough steps and not hard to do so we should let people do theme
# - TODO: Should collapseTopConsensusPerBarcode() have built-in diagnostics for ties?
# Answer: Sure

# ---- Put every BCR and TCR on the tissue slice -------------------------------
# TODO: @Codex update plotReceptorOverlay to use Seurat's spatial plotting so that we can actually see the tissue underneath
allBcrPlot <- plotReceptorOverlay(
  tissue_positions = positions8um,
  receptors = bcrTopBarcode,
  receptor_color = "#2c7fb8",
  receptor_size = 1.2,
  title = paste(sampleId, "all BCR-positive spots")
)
if (savePlots) {
  ggsave(file.path(plotDir, sprintf("%s_all_bcr_locations.png", sampleLower)), allBcrPlot, width = 6, height = 5, dpi = 300)
}


allTcrPlot <- plotReceptorOverlay(
  tissue_positions = positions8um,
  receptors = tcrTopBarcode,
  receptor_color = "#d95f02",
  receptor_size = 1.2,
  title = paste(sampleId, "all TCR-positive spots")
)
if (savePlots) {
  ggsave(file.path(plotDir, sprintf("%s_all_tcr_locations.png", sampleLower)), allTcrPlot, width = 6, height = 5, dpi = 300)
}

bcrTcrOverlay <- ggplot(positions8um, aes(x = x, y = y)) +
  geom_point(color = "grey88", size = 0.55, alpha = 0.75) +
  geom_point(data = bcrTopBarcode, aes(x = x, y = y), inherit.aes = FALSE, shape = 21, fill = NA, color = "#2c7fb8", stroke = 0.45, size = 1.4) +
  geom_point(data = tcrTopBarcode, aes(x = x, y = y), inherit.aes = FALSE, shape = 21, fill = NA, color = "#d95f02", stroke = 0.45, size = 1.4) +
  coord_fixed() +
  scale_y_reverse() +
  theme_minimal() +
  labs(title = paste(sampleId, "BCR and TCR locations"), x = "x", y = "y")
if (savePlots) {
  ggsave(file.path(plotDir, sprintf("%s_bcr_tcr_locations.png", sampleLower)), bcrTcrOverlay, width = 6, height = 5, dpi = 300)
}

# Notes:
# - TODO: Should plotReceptorOverlay() support BCR and TCR layers at the same time?
# Answer: Yes, but I think it's easier to do that with ggplot directly since then we can control the colors and legend better. We can add a helper that makes it easy to plot two sets of receptors with ggplot though.

# ---- Add gene expression underneath receptor locations -----------------------

thymus <- NULL
geneOverlayPlots <- list()
genesToShow <- c("CD79A", "CD3E")

if (loadSeuratObject && file.exists(seuratPath) && requireNamespace("SeuratObject", quietly = TRUE)) {
  thymus <- readRDS(seuratPath)

  for (gene in genesToShow) {
    if (!gene %in% rownames(thymus)) {
      message("Skipping ", gene, ": not found in the Seurat object.")
      next
    }

    expr <- SeuratObject::FetchData(thymus, vars = gene)
    expr$barcode_8um <- rownames(expr)

    tissueGene <- positions %>%
      left_join(expr, by = "barcode_8um") %>%
      mutate(expression = .data[[gene]])

    bcrGenePlot <- plotReceptorOverlay(
      tissue_positions = tissueGene,
      receptors = bcr,
      gene_col = "expression",
      receptor_color = "#2c7fb8",
      receptor_size = 1.25,
      title = paste(gene, "expression with BCR spots")
    )

    tcrGenePlot <- plotReceptorOverlay(
      tissue_positions = tissueGene,
      receptors = tcr,
      gene_col = "expression",
      receptor_color = "#d95f02",
      receptor_size = 1.25,
      title = paste(gene, "expression with TCR spots")
    )

    geneOverlayPlots[[gene]] <- cowplot::plot_grid(bcrGenePlot, tcrGenePlot, ncol = 2)
    print(geneOverlayPlots[[gene]])
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
  message("Skipping gene overlays. Check loadSeuratObject, seuratPath, and SeuratObject.")
}



# Notes:
# - TODO: Should SIRAP include a helper that fetches Seurat expression and joins tissue positions?
# No, we should just already be working with seurat


# @Codex I haven't looked at the below yet, so make updates as you see fit, but that's why I haven't left any comments below here
# ---- Call clones -------------------------------------------------------------

bcrClonesPath <- file.path(outDir, sprintf("%s_bcr_hierarchical_clones.rds", sampleLower))
tcrClonesPath <- file.path(outDir, sprintf("%s_tcr_hierarchical_clones.rds", sampleLower))
bcrTreesPath <- file.path(outDir, sprintf("%s_bcr_pml_trees.rds", sampleLower))

if (runCloneCalling) {
  if (!requireNamespace("scoper", quietly = TRUE)) {
    stop("Install scoper to call clones.")
  }

  bcrClones <- scoper::hierarchicalClones(
    filter(bcr, productive %in% c(TRUE, "TRUE", "T")),
    threshold = 0.1,
    cell_id = "sequence_id",
    summarize_clones = FALSE,
    nproc = nproc
  )

  tcrClones <- scoper::hierarchicalClones(
    filter(tcr, productive %in% c(TRUE, "TRUE", "T")),
    threshold = 0.1,
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
  message("Clone RDS files not found. Set runCloneCalling <- TRUE and rerun this section.")
}

if (exists("bcrClones")) {
  bcrCloneCounts <- bcrClones %>%
    filter(!is.na(clone_id)) %>%
    count(clone_id, sort = TRUE)
  head(bcrCloneCounts, 20)
}

if (exists("tcrClones")) {
  tcrCloneCounts <- tcrClones %>%
    filter(!is.na(clone_id)) %>%
    count(clone_id, sort = TRUE)
  head(tcrCloneCounts, 20)
}

# ---- Pick a couple of clones and plot them ----------------------------------

chooseBigClones <- function(clones, cloneIds = NULL, n = 2) {
  if (!is.null(cloneIds)) {
    return(as.character(cloneIds))
  }
  clones %>%
    filter(!is.na(clone_id)) %>%
    count(clone_id, sort = TRUE) %>%
    slice_head(n = n) %>%
    pull(clone_id) %>%
    as.character()
}

if (exists("bcrClones")) {
  bcrCloneIds <- chooseBigClones(bcrClones, bcrCloneIdsToPlot, n = 2)
  bcrCloneIds
}

if (exists("tcrClones")) {
  tcrCloneIds <- chooseBigClones(tcrClones, tcrCloneIdsToPlot, n = 2)
  tcrCloneIds
}

# After inspecting clone counts, set these manually if you want:
# bcrCloneIdsToPlot <- c("123", "456")
# tcrCloneIdsToPlot <- c("789", "101")

if (exists("bcrClones") && exists("bcrCloneIds")) {
  bcrClonePlots <- lapply(
    bcrCloneIds,
    function(cloneId) {
      plotReceptorOverlay(
        tissue_positions = positions,
        receptors = filter(bcrClones, .data$clone_id == .env$cloneId),
        receptor_color = "#2c7fb8",
        receptor_size = 1.7,
        title = paste("BCR clone", cloneId)
      )
    }
  )
  selectedBcrClonePlot <- cowplot::plot_grid(plotlist = bcrClonePlots, ncol = 2)
  selectedBcrClonePlot
  if (savePlots) {
    ggsave(file.path(plotDir, sprintf("%s_selected_bcr_clones_spatial.png", sampleLower)), selectedBcrClonePlot, width = 10, height = 5, dpi = 300)
  }
}

if (exists("tcrClones") && exists("tcrCloneIds")) {
  tcrClonePlots <- lapply(
    tcrCloneIds,
    function(cloneId) {
      plotReceptorOverlay(
        tissue_positions = positions,
        receptors = filter(tcrClones, .data$clone_id == .env$cloneId),
        receptor_color = "#d95f02",
        receptor_size = 1.7,
        title = paste("TCR clone", cloneId)
      )
    }
  )
  selectedTcrClonePlot <- cowplot::plot_grid(plotlist = tcrClonePlots, ncol = 2)
  selectedTcrClonePlot
  if (savePlots) {
    ggsave(file.path(plotDir, sprintf("%s_selected_tcr_clones_spatial.png", sampleLower)), selectedTcrClonePlot, width = 10, height = 5, dpi = 300)
  }
}

# ---- Build BCR trees and match tips back to space ----------------------------

if (buildBcrTrees) {
  if (!exists("bcrClones") || !exists("bcrCloneIds")) {
    stop("Run clone calling/loading and clone selection before tree building.")
  }
  if (!requireNamespace("dowser", quietly = TRUE)) {
    stop("Install dowser to build BCR lineage trees.")
  }
  if (!requireNamespace("alakazam", quietly = TRUE)) {
    stop("Install alakazam to read IMGT germlines for dowser tree building.")
  }

  refs <- alakazam::readIMGT(germlineDir)

  bcrTreeInput <- bcrClones %>%
    filter(clone_id %in% bcrCloneIds, !is.na(clone_id)) %>%
    mutate(
      sequence_id_tree = make.unique(as.character(sequence_id)),
      c_call_simple = case_when(
        stringr::str_detect(c_call, "IGHM") ~ "IgM",
        stringr::str_detect(c_call, "IGHD") ~ "IgD",
        stringr::str_detect(c_call, "IGHG1") ~ "IgG1",
        stringr::str_detect(c_call, "IGHG2") ~ "IgG2",
        stringr::str_detect(c_call, "IGHG3") ~ "IgG3",
        stringr::str_detect(c_call, "IGHG4") ~ "IgG4",
        stringr::str_detect(c_call, "IGHA1") ~ "IgA1",
        stringr::str_detect(c_call, "IGHA2") ~ "IgA2",
        stringr::str_detect(c_call, "IGHE") ~ "IgE",
        TRUE ~ "unknown"
      )
    ) %>%
    assignTipColors(clone_ids = bcrCloneIds, id_col = "sequence_id_tree")

  bcrTreeInput <- dowser::createGermlines(bcrTreeInput, refs, id = "sequence_id_tree")

  bcrFormattedClones <- dowser::formatClones(
    bcrTreeInput,
    id = "sequence_id_tree",
    traits = c("tip_color_id", "tip_color", "c_call_simple", "barcode_8um"),
    columns = c("x", "y", "barcode_8um"),
    filterstop = FALSE,
    collapse = TRUE,
    add_count = "consensus_count" %in% colnames(bcrTreeInput),
    num_fields = if ("consensus_count" %in% colnames(bcrTreeInput)) "consensus_count" else NULL
  )
  bcrTrees <- dowser::getTrees(bcrFormattedClones, build = "pml", nproc = nproc)
  bcrTreeResult <- list(formatted_clones = bcrFormattedClones, trees = bcrTrees)

  saveRDS(bcrTreeResult, bcrTreesPath)
} else if (file.exists(bcrTreesPath)) {
  bcrTreeResult <- readRDS(bcrTreesPath)
}

if (exists("bcrTreeResult") && exists("bcrTreeInput")) {
  treePalette <- unique(bcrTreeInput$tip_color)
  names(treePalette) <- unique(bcrTreeInput$tip_color_id)

  bcrTreePlots <- dowser::plotTrees(
    bcrTreeResult$trees,
    tips = "tip_color_id",
    palette = treePalette,
    tipsize = if ("consensus_count" %in% colnames(bcrTreeInput)) "consensus_count" else NULL
  )

  pairedTreeSpatialPlots <- list()
  for (i in seq_along(bcrCloneIds)) {
    cloneId <- bcrCloneIds[[i]]
    spatialPlot <- plotColoredCloneSpatial(
      bcrTreeInput,
      clone_id = cloneId,
      tissue_positions = positions
    )

    pairedTreeSpatialPlots[[i]] <- cowplot::plot_grid(
      bcrTreePlots[[i]] + ggplot2::ggtitle(paste("BCR clone", cloneId, "tree")),
      spatialPlot + ggplot2::ggtitle(paste("BCR clone", cloneId, "matched tissue colors")),
      ncol = 2,
      rel_widths = c(1.2, 1)
    )
    print(pairedTreeSpatialPlots[[i]])
  }

  if (savePlots) {
    pdf(file.path(plotDir, sprintf("%s_bcr_tree_spatial_pairs.pdf", sampleLower)), width = 12, height = 6)
    for (page in pairedTreeSpatialPlots) {
      print(page)
    }
    dev.off()
  }
}

# Notes:
# - TODO: Should plotTreeSpatialPairs() become a package function?

# ---- Done -------------------------------------------------------------------

message("Walkthrough objects are loaded in the current R session.")
message("Output directory: ", outDir)
