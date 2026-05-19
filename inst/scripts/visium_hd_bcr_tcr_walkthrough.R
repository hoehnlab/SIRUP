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
})

# Notes while testing:
# - TODO:
# - TODO:
# - TODO:

# ---- Things to change --------------------------------------------------------

sampleId <- "M5" # M1 uses demuxed_bc23; M5 uses demuxed_bc24.
sampleLower <- tolower(sampleId)

thymusDir <- normalizePath("../../Thymus_VisiumHD", mustWork = TRUE)
demuxBySample <- c(M1 = "022526_Visium_demuxed_bc23", M5 = "022526_Visium_demuxed_bc24")
demuxDir <- file.path("/dartfs/rc/lab/H/HoehnK/Chinmay/022126_Visium_HD", demuxBySample[[sampleId]])

rawBcrIgblastPath <- file.path(demuxDir, "igblast_out", "igblast_bcr.tsv")
rawTcrIgblastPath <- file.path(demuxDir, "igblast_out", "igblast_tcr.tsv")
barcodeMapPath <- file.path(thymusDir, sprintf("%s_barcode_mappings.csv", sampleId))
positionsPath <- file.path(thymusDir, sprintf("%s_tissue_positions_8um.csv", sampleId))
seuratPath <- file.path(thymusDir, "data", sprintf("Thymus_VisiumHD_%s_8um_hires_sr_processed.rds", sampleId))

outDir <- file.path(thymusDir, "SIRAP_walkthrough", sampleId)
prestoInputDir <- file.path(outDir, "01_length_filtered_spatial_fastas")
prestoOutputDir <- file.path(outDir, "02_presto_igblast")
finalAirrDir <- file.path(outDir, "03_final_spatial_airr")
plotDir <- file.path(outDir, "plots")
for (dir in c(prestoInputDir, prestoOutputDir, finalAirrDir, plotDir)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

nproc <- 4
genesToShow <- c("CD79A", "CD3E")

# pRESTO/IgBLAST can take a while, so leave this FALSE until the input FASTAs
# look right. When FALSE, the script will look for db-pass TSVs from a prior run.
runPrestoIgblast <- FALSE
runCloneCalling <- FALSE
buildBcrTrees <- FALSE
savePlots <- TRUE
loadSeuratObject <- TRUE

identity <- 0.975
consensusFreq <- 0.51
germlineGlob <- "~/software/germlines/imgt/human/vdj/imgt_human_*.fasta"
germlineDir <- path.expand("~/software/germlines/imgt/human/vdj")
igblastBase <- "~/software/igblast"
igblastExec <- "~/software/igblast/bin/igblastn"
muscleExec <- "~/software/muscle3.8.31_i86linux64"

bcrCloneIdsToPlot <- NULL
tcrCloneIdsToPlot <- NULL

# ---- Small helpers for this walkthrough -------------------------------------

checkFiles <- function(paths) {
  fileStatus <- data.frame(path = paths, exists = file.exists(paths))
  print(fileStatus)
  if (!all(fileStatus$exists)) {
    stop("Some input files are missing. Edit the paths above, then rerun this section.")
  }
  invisible(fileStatus)
}

expectedPrestoPaths <- function(input_fasta, output_dir, identity = 0.975, consensus_freq = 0.51) {
  stem <- tools::file_path_sans_ext(basename(input_fasta))
  consensus <- file.path(
    output_dir,
    sprintf(
      "%s_ident_%s_barcode_cluster_freq_%s_cluster-pass_reheader_align-pass_consensus-pass.fasta",
      stem,
      identity,
      consensus_freq
    )
  )
  list(
    consensus = consensus,
    fmt7 = sub("\\.fasta$", "_igblast.fmt7", consensus),
    db_pass = sub("\\.fasta$", "_igblast_db-pass.tsv", consensus)
  )
}

savePlot <- function(plot, filename, width = 8, height = 6) {
  if (savePlots) {
    ggsave(file.path(plotDir, filename), plot, width = width, height = height, dpi = 300)
  }
  invisible(plot)
}

# ---- Load the spatial metadata ----------------------------------------------

checkFiles(c(rawBcrIgblastPath, rawTcrIgblastPath, barcodeMapPath, positionsPath))

barcodeMap <- readBarcodeMapping(barcodeMapPath)

positions <- readBarcodeMapping(positionsPath) |>
  mutate(
    barcode_8um = barcode,
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres
  )

head(barcodeMap)
head(positions)

# ---- Start from raw IgBLAST output ------------------------------------------

# The original script first kept IGH and TRB reads with V >= 100, J >= 10, and
# a real junction. Then it pulled the 2 um barcode from sequence_id, mapped it
# to 8 um and 16 um barcodes, and wrote both AIRR TSVs and FASTAs for pRESTO.

bcrRaw <- readSpatialAirr(rawBcrIgblastPath)
tcrRaw <- readSpatialAirr(rawTcrIgblastPath)

dim(bcrRaw)
dim(tcrRaw)

bcrLengthSpatial <- preprocessSpatialAirr(
  input = bcrRaw,
  barcode_map = barcodeMap,
  loci = "IGH",
  out_tsv = file.path(prestoInputDir, sprintf("%s_igblast_bcr_pass_length_spatial_barcode.tsv", sampleLower)),
  out_fasta = file.path(prestoInputDir, sprintf("%s_igblast_bcr_pass_length_spatial_barcode.fasta", sampleLower)),
  min_v_length = 100,
  min_j_length = 10,
  min_junction_length = 1,
  sequence_col = "sequence_alignment"
)

tcrLengthSpatial <- preprocessSpatialAirr(
  input = tcrRaw,
  barcode_map = barcodeMap,
  loci = "TRB",
  out_tsv = file.path(prestoInputDir, sprintf("%s_igblast_tcr_pass_length_spatial_barcode.tsv", sampleLower)),
  out_fasta = file.path(prestoInputDir, sprintf("%s_igblast_tcr_pass_length_spatial_barcode.fasta", sampleLower)),
  min_v_length = 100,
  min_j_length = 10,
  min_junction_length = 1,
  sequence_col = "sequence_alignment"
)

spatialReceptorSummary(bcrLengthSpatial, barcode_col = "barcode_2um")
spatialReceptorSummary(tcrLengthSpatial, barcode_col = "barcode_2um")

prePrestoQcPlot <- cowplot::plot_grid(
  plotLengthQc(bcrLengthSpatial, title = paste(sampleId, "raw IGH reads after length/spatial filtering")),
  plotLengthQc(tcrLengthSpatial, title = paste(sampleId, "raw TRB reads after length/spatial filtering")),
  ncol = 1
)
prePrestoQcPlot
savePlot(prePrestoQcPlot, sprintf("%s_pre_presto_length_qc.png", sampleLower), width = 10, height = 8)

# Things to inspect:
head(bcrLengthSpatial)
head(tcrLengthSpatial)

# Notes:
# - TODO: Should preprocessSpatialAirr() print a drop-count summary?
# - TODO: Should full_id include barcode16um too?

# ---- Run pRESTO and IgBLAST from R ------------------------------------------

# This is the command-line part of the original workflow, wrapped in R:
# ClusterSets.py -> ParseHeaders.py -> AlignSets.py -> BuildConsensus.py ->
# AssignGenes.py -> MakeDb.py. The command output prints to the R console.

bcrLengthSpatialFasta <- file.path(prestoInputDir, sprintf("%s_igblast_bcr_pass_length_spatial_barcode.fasta", sampleLower))
tcrLengthSpatialFasta <- file.path(prestoInputDir, sprintf("%s_igblast_tcr_pass_length_spatial_barcode.fasta", sampleLower))

if (runPrestoIgblast) {
  bcrPrestoPaths <- runPrestoIgblastSample(
    input_fasta = bcrLengthSpatialFasta,
    output_dir = prestoOutputDir,
    receptor = "ig",
    identity = identity,
    consensus_freq = consensusFreq,
    germline_glob = germlineGlob,
    igblast_base = igblastBase,
    igblast_exec = igblastExec,
    muscle_exec = muscleExec,
    nproc = nproc
  )

  tcrPrestoPaths <- runPrestoIgblastSample(
    input_fasta = tcrLengthSpatialFasta,
    output_dir = prestoOutputDir,
    receptor = "tr",
    identity = identity,
    consensus_freq = consensusFreq,
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
  bcrPrestoPaths <- expectedPrestoPaths(bcrLengthSpatialFasta, prestoOutputDir, identity, consensusFreq)
}
if (!exists("tcrPrestoPaths")) {
  tcrPrestoPaths <- expectedPrestoPaths(tcrLengthSpatialFasta, prestoOutputDir, identity, consensusFreq)
}

bcrPrestoPaths$db_pass
tcrPrestoPaths$db_pass

# Notes:
# - TODO: Should SIRAP expose expectedPrestoPaths() as a real package function?
# - TODO: Should runPrestoIgblastSample() check all external executables first?

# ---- Polish the pRESTO/IgBLAST consensus output -----------------------------

# Now we are back in R. The original script read the db-pass TSVs, kept strong
# productive consensus sequences, recovered the barcode from the pRESTO header,
# swapped the 2 um barcode row/column order, remapped to 8 um/16 um, and kept
# the top consensus sequence per 2 um barcode and locus.

checkFiles(c(bcrPrestoPaths$db_pass, tcrPrestoPaths$db_pass))

bcrConsensus <- readSpatialAirr(bcrPrestoPaths$db_pass)
tcrConsensus <- readSpatialAirr(tcrPrestoPaths$db_pass)

bcrConsensusFiltered <- filterReceptorReads(
  bcrConsensus,
  loci = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 11,
  productive_only = TRUE,
  min_consensus_count = 2,
  max_opening_gap = 2
)

tcrConsensusFiltered <- filterReceptorReads(
  tcrConsensus,
  loci = "TRB",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 11,
  productive_only = TRUE,
  min_consensus_count = 2,
  max_opening_gap = 2
)

bcrSpatialConsensus <- addSpatialBarcodes(bcrConsensusFiltered, barcodeMap, mode = "presto_header", swap = TRUE)
tcrSpatialConsensus <- addSpatialBarcodes(tcrConsensusFiltered, barcodeMap, mode = "presto_header", swap = TRUE)

bcrTopBarcode <- collapseTopConsensusPerBarcode(bcrSpatialConsensus, barcode_col = "barcode_2um", locus_col = "locus", ties = "all")
tcrTopBarcode <- collapseTopConsensusPerBarcode(tcrSpatialConsensus, barcode_col = "barcode_2um", locus_col = "locus", ties = "all")

bcr <- addTissuePositions(
  bcrTopBarcode,
  positions,
  barcode_col = "barcode_8um",
  position_barcode_col = "barcode_8um",
  x_col = "x",
  y_col = "y"
) |>
  filter(!is.na(x), !is.na(y))

tcr <- addTissuePositions(
  tcrTopBarcode,
  positions,
  barcode_col = "barcode_8um",
  position_barcode_col = "barcode_8um",
  x_col = "x",
  y_col = "y"
) |>
  filter(!is.na(x), !is.na(y))

finalBcrPath <- file.path(finalAirrDir, sprintf("%s_bcrs_heavy_final_spatial_airr.tsv", sampleLower))
finalTcrPath <- file.path(finalAirrDir, sprintf("%s_tcrs_beta_final_spatial_airr.tsv", sampleLower))
writeSpatialOutputs(bcr, out_tsv = finalBcrPath)
writeSpatialOutputs(tcr, out_tsv = finalTcrPath)

list(
  bcr = spatialReceptorSummary(bcr, barcode_col = "barcode_8um"),
  tcr = spatialReceptorSummary(tcr, barcode_col = "barcode_8um")
)

postPrestoQcPlot <- cowplot::plot_grid(
  plotLengthQc(bcr, title = paste(sampleId, "final BCRs")),
  plotLengthQc(tcr, title = paste(sampleId, "final TCRs")),
  ncol = 1
)
postPrestoQcPlot
savePlot(postPrestoQcPlot, sprintf("%s_post_presto_final_length_qc.png", sampleLower), width = 10, height = 8)

# Things to inspect:
head(bcr)
head(tcr)

# Notes:
# - TODO: Should the post-pRESTO block become a single high-level function?
# - TODO: Should collapseTopConsensusPerBarcode() have built-in diagnostics for ties?

# ---- Put every BCR and TCR on the tissue slice -------------------------------

allBcrPlot <- plotReceptorOverlay(
  tissue_positions = positions,
  receptors = bcr,
  receptor_color = "#2c7fb8",
  receptor_size = 1.2,
  title = paste(sampleId, "all BCR-positive spots")
)
allBcrPlot

allTcrPlot <- plotReceptorOverlay(
  tissue_positions = positions,
  receptors = tcr,
  receptor_color = "#d95f02",
  receptor_size = 1.2,
  title = paste(sampleId, "all TCR-positive spots")
)
allTcrPlot

bcrTcrOverlay <- ggplot(positions, aes(x = x, y = y)) +
  geom_point(color = "grey88", size = 0.55, alpha = 0.75) +
  geom_point(data = bcr, aes(x = x, y = y), inherit.aes = FALSE, shape = 21, fill = NA, color = "#2c7fb8", stroke = 0.45, size = 1.4) +
  geom_point(data = tcr, aes(x = x, y = y), inherit.aes = FALSE, shape = 21, fill = NA, color = "#d95f02", stroke = 0.45, size = 1.4) +
  coord_fixed() +
  scale_y_reverse() +
  theme_minimal() +
  labs(title = paste(sampleId, "BCR and TCR locations"), x = "x", y = "y")
bcrTcrOverlay

savePlot(allBcrPlot, sprintf("%s_all_bcr_locations.png", sampleLower), width = 6, height = 5)
savePlot(allTcrPlot, sprintf("%s_all_tcr_locations.png", sampleLower), width = 6, height = 5)
savePlot(bcrTcrOverlay, sprintf("%s_bcr_tcr_locations.png", sampleLower), width = 6, height = 5)

# Notes:
# - TODO: Should plotReceptorOverlay() support BCR and TCR layers at the same time?

# ---- Add gene expression underneath receptor locations -----------------------

thymus <- NULL
geneOverlayPlots <- list()

if (loadSeuratObject && file.exists(seuratPath) && requireNamespace("SeuratObject", quietly = TRUE)) {
  thymus <- readRDS(seuratPath)

  for (gene in genesToShow) {
    if (!gene %in% rownames(thymus)) {
      message("Skipping ", gene, ": not found in the Seurat object.")
      next
    }

    expr <- SeuratObject::FetchData(thymus, vars = gene)
    expr$barcode_8um <- rownames(expr)

    tissueGene <- positions |>
      left_join(expr, by = "barcode_8um") |>
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
    savePlot(geneOverlayPlots[[gene]], sprintf("%s_%s_receptor_gene_overlay.png", sampleLower, gene), width = 11, height = 5)
  }
} else {
  message("Skipping gene overlays. Check loadSeuratObject, seuratPath, and SeuratObject.")
}

# Notes:
# - TODO: Should SIRAP include a helper that fetches Seurat expression and joins tissue positions?

# ---- Call clones -------------------------------------------------------------

bcrClonesPath <- file.path(outDir, sprintf("%s_bcr_hierarchical_clones.rds", sampleLower))
tcrClonesPath <- file.path(outDir, sprintf("%s_tcr_hierarchical_clones.rds", sampleLower))
bcrTreesPath <- file.path(outDir, sprintf("%s_bcr_pml_trees.rds", sampleLower))

if (runCloneCalling) {
  bcrClones <- assignSpatialClones(
    bcr,
    threshold = 0.1,
    productive_only = TRUE,
    cell_id = "sequence_id",
    nproc = nproc
  )

  tcrClones <- assignSpatialClones(
    tcr,
    threshold = 0.1,
    productive_only = TRUE,
    cell_id = "sequence_id",
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
  bcrCloneCounts <- bcrClones |>
    filter(!is.na(clone_id)) |>
    count(clone_id, sort = TRUE)
  head(bcrCloneCounts, 20)
}

if (exists("tcrClones")) {
  tcrCloneCounts <- tcrClones |>
    filter(!is.na(clone_id)) |>
    count(clone_id, sort = TRUE)
  head(tcrCloneCounts, 20)
}

# Notes:
# - TODO: Should assignSpatialClones() have better defaults for spatial data?

# ---- Pick a couple of clones and plot them ----------------------------------

chooseBigClones <- function(clones, cloneIds = NULL, n = 2) {
  if (!is.null(cloneIds)) {
    return(as.character(cloneIds))
  }
  clones |>
    filter(!is.na(clone_id)) |>
    count(clone_id, sort = TRUE) |>
    slice_head(n = n) |>
    pull(clone_id) |>
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
  savePlot(selectedBcrClonePlot, sprintf("%s_selected_bcr_clones_spatial.png", sampleLower), width = 10, height = 5)
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
  savePlot(selectedTcrClonePlot, sprintf("%s_selected_tcr_clones_spatial.png", sampleLower), width = 10, height = 5)
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

  bcrTreeInput <- bcrClones |>
    filter(clone_id %in% bcrCloneIds, !is.na(clone_id)) |>
    mutate(
      sequence_id_tree = make.unique(as.character(sequence_id)),
      c_call_simple = simplifyIghConstantCall(c_call)
    ) |>
    assignTipColors(clone_ids = bcrCloneIds, id_col = "sequence_id_tree")

  bcrTreeInput <- dowser::createGermlines(bcrTreeInput, refs, id = "sequence_id_tree")

  bcrTreeResult <- buildDowserTrees(
    bcrTreeInput,
    clone_ids = bcrCloneIds,
    id = "sequence_id_tree",
    traits = c("tip_color_id", "tip_color", "c_call_simple", "barcode_8um"),
    columns = c("x", "y", "barcode_8um"),
    method = "pml",
    nproc = nproc
  )

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
# - TODO: Should buildDowserTrees() handle createGermlines()?
# - TODO: Should plotTreeSpatialPairs() become a package function?

# ---- Done -------------------------------------------------------------------

message("Walkthrough objects are loaded in the current R session.")
message("Output directory: ", outDir)
