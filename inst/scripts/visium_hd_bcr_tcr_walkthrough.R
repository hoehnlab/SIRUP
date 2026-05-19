###############################################################################
# SIRAP interactive Visium HD BCR/TCR walkthrough
#
# This is intentionally an interactive script, not a command-line pipeline.
# Open it in RStudio/Positron/VS Code, run one section at a time, inspect the
# objects that are created, and add notes where a workflow step feels awkward.
#
# Suggested startup from the package root:
#
#   setwd("/Volumes/HoehnK/Hunter/spatial_bcr_tcr_package")
#   devtools::load_all()
#   source("inst/scripts/visium_hd_bcr_tcr_walkthrough.R")
#
# Or run line-by-line / section-by-section, which is the intended use.
###############################################################################

# ---- 0. Package setup --------------------------------------------------------

# If you are developing the package locally, run devtools::load_all() once before
# this script. This block tries to be helpful if you forgot.
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

# Notes / wishlist while testing:
# - TODO:
# - TODO:
# - TODO:

# ---- 1. Settings to edit -----------------------------------------------------

thymusDir <- normalizePath("../Thymus_VisiumHD", mustWork = TRUE)
sampleId <- "M5" # Try "M1" for the other slice.

outDir <- file.path(thymusDir, "SIRAP_walkthrough")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

genesToShow <- c("CD79A", "CD3E")
nproc <- 4

# Interactive toggles. Keep expensive steps off until the earlier objects look
# right. If cached clone/tree RDS files exist, the script will load them.
savePlots <- TRUE
loadSeuratObject <- TRUE
runCloneCalling <- FALSE
buildBcrTrees <- FALSE

germlineDir <- path.expand("~/software/germlines/imgt/human/vdj")

# Leave NULL to choose the two largest clones after clone calling. Once you see
# clone counts, set these manually and rerun the plotting/tree sections.
bcrCloneIdsToPlot <- NULL
tcrCloneIdsToPlot <- NULL

# ---- 2. Paths ----------------------------------------------------------------

sampleLower <- tolower(sampleId)

bcrPath <- file.path(thymusDir, "data", sprintf("%s_bcrs_heavy_final_swapped_002.tsv", sampleLower))
tcrPath <- file.path(thymusDir, "data", sprintf("%s_tcrs_beta_final_swapped_002.tsv", sampleLower))
positionsPath <- file.path(thymusDir, sprintf("%s_tissue_positions_8um.csv", sampleId))
seuratPath <- file.path(thymusDir, "data", sprintf("Thymus_VisiumHD_%s_8um_hires_sr_processed.rds", sampleId))

bcrClonesPath <- file.path(outDir, sprintf("%s_bcr_hierarchical_clones.rds", sampleLower))
tcrClonesPath <- file.path(outDir, sprintf("%s_tcr_hierarchical_clones.rds", sampleLower))
bcrTreesPath <- file.path(outDir, sprintf("%s_bcr_pml_trees.rds", sampleLower))

pathsToCheck <- c(bcrPath, tcrPath, positionsPath)
if (!all(file.exists(pathsToCheck))) {
  print(data.frame(path = pathsToCheck, exists = file.exists(pathsToCheck)))
  stop("One or more required input files are missing. Edit thymusDir/sampleId above.")
}

# ---- 3. Load BCR/TCR AIRR tables and tissue coordinates ----------------------

bcr <- readSpatialAirr(bcrPath)
tcr <- readSpatialAirr(tcrPath)
positions <- readBarcodeMapping(positionsPath)

positions <- positions |>
  mutate(
    barcode_8um = barcode,
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres
  )

bcr <- bcr |>
  left_join(positions |> select(barcode_8um, x, y), by = "barcode_8um") |>
  filter(!is.na(x), !is.na(y))

tcr <- tcr |>
  left_join(positions |> select(barcode_8um, x, y), by = "barcode_8um") |>
  filter(!is.na(x), !is.na(y))

spatialSummary <- list(
  bcr = spatialReceptorSummary(bcr, barcode_col = "barcode_8um"),
  tcr = spatialReceptorSummary(tcr, barcode_col = "barcode_8um")
)
spatialSummary

# Things to inspect interactively:
head(bcr)
head(tcr)
head(positions)

# Notes:
# - TODO: Are barcode/coordinate column names intuitive enough?
# - TODO: Should addTissuePositions() do this join automatically for this case?

# ---- 4. Plot all BCR and TCR locations --------------------------------------

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

if (savePlots) {
  ggsave(file.path(outDir, sprintf("%s_all_bcr_locations.png", sampleLower)), allBcrPlot, width = 6, height = 5, dpi = 300)
  ggsave(file.path(outDir, sprintf("%s_all_tcr_locations.png", sampleLower)), allTcrPlot, width = 6, height = 5, dpi = 300)
  ggsave(file.path(outDir, sprintf("%s_bcr_tcr_locations.png", sampleLower)), bcrTcrOverlay, width = 6, height = 5, dpi = 300)
}

# Notes:
# - TODO: Should plotReceptorOverlay() support two receptor layers directly?
# - TODO: Are hollow circles readable at 8um spot density?

# ---- 5. Gene expression overlays with hollow receptor circles ----------------

thymus <- NULL
geneOverlayPlots <- list()

if (loadSeuratObject && file.exists(seuratPath) && requireNamespace("SeuratObject", quietly = TRUE)) {
  thymus <- readRDS(seuratPath)

  for (gene in genesToShow) {
    if (!gene %in% rownames(thymus)) {
      message("Skipping ", gene, ": not found in Seurat object.")
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

    if (savePlots) {
      ggsave(
        file.path(outDir, sprintf("%s_%s_receptor_gene_overlay.png", sampleLower, gene)),
        geneOverlayPlots[[gene]],
        width = 11,
        height = 5,
        dpi = 300
      )
    }
  }
} else {
  message("Skipping Seurat gene overlays. Set loadSeuratObject <- TRUE and check seuratPath if needed.")
}

# Notes:
# - TODO: Should we add a package function to fetch expression and join positions?
# - TODO: Should gene expression overlays use Seurat image coordinates instead?

# ---- 6. Clone calling or cached clone loading --------------------------------

if (runCloneCalling) {
  message("Calling BCR clones with scoper::hierarchicalClones().")
  bcrClones <- assignSpatialClones(
    bcr,
    threshold = 0.1,
    productive_only = TRUE,
    cell_id = "sequence_id",
    nproc = nproc
  )

  message("Calling TCR clones with scoper::hierarchicalClones().")
  tcrClones <- assignSpatialClones(
    tcr,
    threshold = 0.1,
    productive_only = TRUE,
    cell_id = "sequence_id",
    nproc = nproc
  )

  saveRDS(bcrClones, bcrClonesPath)
  saveRDS(tcrClones, tcrClonesPath)
} else {
  if (file.exists(bcrClonesPath) && file.exists(tcrClonesPath)) {
    bcrClones <- readRDS(bcrClonesPath)
    tcrClones <- readRDS(tcrClonesPath)
  } else {
    message("Clone RDS files not found. Set runCloneCalling <- TRUE and rerun this section.")
  }
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
# - TODO: Should assignSpatialClones() cache automatically?
# - TODO: Does the default threshold/cell_id feel right for spatial long reads?

# ---- 7. Pick clones to visualize --------------------------------------------

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

# After inspecting clone counts, you can set these manually, then rerun from here:
# bcrCloneIdsToPlot <- c("123", "456")
# tcrCloneIdsToPlot <- c("789", "101")

# ---- 8. Plot selected BCR/TCR clones on the slice ----------------------------

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
  cowplot::plot_grid(plotlist = bcrClonePlots, ncol = 2)

  if (savePlots) {
    ggsave(
      file.path(outDir, sprintf("%s_selected_bcr_clones_spatial.png", sampleLower)),
      cowplot::plot_grid(plotlist = bcrClonePlots, ncol = 2),
      width = 10,
      height = 5,
      dpi = 300
    )
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
  cowplot::plot_grid(plotlist = tcrClonePlots, ncol = 2)

  if (savePlots) {
    ggsave(
      file.path(outDir, sprintf("%s_selected_tcr_clones_spatial.png", sampleLower)),
      cowplot::plot_grid(plotlist = tcrClonePlots, ncol = 2),
      width = 10,
      height = 5,
      dpi = 300
    )
  }
}

# Notes:
# - TODO: Should there be a clone-specific wrapper around plotReceptorOverlay()?
# - TODO: Should clone spot sizes scale by consensus_count?

# ---- 9. Optional BCR trees with matched spatial tip colors -------------------

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
  if (!dir.exists(germlineDir)) {
    stop("germlineDir does not exist: ", germlineDir)
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
  message("Loaded cached BCR trees: ", bcrTreesPath)
}

# Run this after building or loading bcrTreeResult.
if (exists("bcrTreeResult") && exists("bcrTreeInput")) {
  treePalette <- unique(bcrTreeInput$tip_color)
  names(treePalette) <- unique(bcrTreeInput$tip_color_id)

  if ("consensus_count" %in% colnames(bcrTreeInput)) {
    bcrTreePlots <- dowser::plotTrees(
      bcrTreeResult$trees,
      tips = "tip_color_id",
      palette = treePalette,
      tipsize = "consensus_count"
    )
  } else {
    bcrTreePlots <- dowser::plotTrees(
      bcrTreeResult$trees,
      tips = "tip_color_id",
      palette = treePalette
    )
  }

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
    pdf(file.path(outDir, sprintf("%s_bcr_tree_spatial_pairs.pdf", sampleLower)), width = 12, height = 6)
    for (page in pairedTreeSpatialPlots) {
      print(page)
    }
    dev.off()
  }
}

# Notes:
# - TODO: Should buildDowserTrees() handle createGermlines() too?
# - TODO: Should plotTreeSpatialPairs() become a package function?
# - TODO: How should we choose colors when a clone has many tips?

# ---- 10. End -----------------------------------------------------------------

message("Interactive walkthrough objects are loaded in the current R session.")
message("Output directory: ", outDir)
