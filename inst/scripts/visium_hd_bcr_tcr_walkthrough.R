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

sample_id <- "M5" # M1 uses demuxed_bc23; M5 uses demuxed_bc24.
sample_lower <- tolower(sample_id)

thymus_dir <- normalizePath("../../Thymus_VisiumHD", mustWork = TRUE)
demux_by_sample <- c(M1 = "022526_Visium_demuxed_bc23", M5 = "022526_Visium_demuxed_bc24")
demux_dir <- file.path("/dartfs/rc/lab/H/HoehnK/Chinmay/022126_Visium_HD", demux_by_sample[[sample_id]])

raw_bcr_igblast_path <- file.path(demux_dir, "igblast_out", "igblast_bcr.tsv")
raw_tcr_igblast_path <- file.path(demux_dir, "igblast_out", "igblast_tcr.tsv")
barcode_map_path <- file.path(thymus_dir, sprintf("%s_barcode_mappings.csv", sample_id))
positions_path <- file.path(thymus_dir, sprintf("%s_tissue_positions_8um.csv", sample_id))
seurat_path <- file.path(thymus_dir, "data", sprintf("Thymus_VisiumHD_%s_8um_hires_sr_processed.rds", sample_id))

out_dir <- file.path(thymus_dir, "SIRAP_walkthrough", sample_id)
presto_input_dir <- file.path(out_dir, "01_length_filtered_spatial_fastas")
presto_output_dir <- file.path(out_dir, "02_presto_igblast")
final_airr_dir <- file.path(out_dir, "03_final_spatial_airr")
plot_dir <- file.path(out_dir, "plots")
for (dir in c(presto_input_dir, presto_output_dir, final_airr_dir, plot_dir)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

nproc <- 4
genes_to_show <- c("CD79A", "CD3E")

# pRESTO/IgBLAST can take a while, so leave this FALSE until the input FASTAs
# look right. When FALSE, the script will look for db-pass TSVs from a prior run.
run_presto_igblast <- FALSE
run_clone_calling <- FALSE
build_bcr_trees <- FALSE
save_plots <- TRUE
load_seurat_object <- TRUE

identity <- 0.975
consensus_freq <- 0.51
germline_glob <- "~/software/germlines/imgt/human/vdj/imgt_human_*.fasta"
germline_dir <- path.expand("~/software/germlines/imgt/human/vdj")
igblast_base <- "~/software/igblast"
igblast_exec <- "~/software/igblast/bin/igblastn"
muscle_exec <- "~/software/muscle3.8.31_i86linux64"

bcr_clone_ids_to_plot <- NULL
tcr_clone_ids_to_plot <- NULL

# ---- Small helpers for this walkthrough -------------------------------------

check_files <- function(paths) {
  file_status <- data.frame(path = paths, exists = file.exists(paths))
  print(file_status)
  if (!all(file_status$exists)) {
    stop("Some input files are missing. Edit the paths above, then rerun this section.")
  }
  invisible(file_status)
}

expected_presto_paths <- function(input_fasta, output_dir, identity = 0.975, consensus_freq = 0.51) {
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

save_plot <- function(plot, filename, width = 8, height = 6) {
  if (save_plots) {
    ggsave(file.path(plot_dir, filename), plot, width = width, height = height, dpi = 300)
  }
  invisible(plot)
}

# ---- Load the spatial metadata ----------------------------------------------

check_files(c(raw_bcr_igblast_path, raw_tcr_igblast_path, barcode_map_path, positions_path))

barcode_map <- read_barcode_mapping(barcode_map_path)

positions <- read_barcode_mapping(positions_path) |>
  mutate(
    barcode_8um = barcode,
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres
  )

head(barcode_map)
head(positions)

# ---- Start from raw IgBLAST output ------------------------------------------

# The original script first kept IGH and TRB reads with V >= 100, J >= 10, and
# a real junction. Then it pulled the 2 um barcode from sequence_id, mapped it
# to 8 um and 16 um barcodes, and wrote both AIRR TSVs and FASTAs for pRESTO.

bcr_raw <- read_spatial_airr(raw_bcr_igblast_path)
tcr_raw <- read_spatial_airr(raw_tcr_igblast_path)

dim(bcr_raw)
dim(tcr_raw)

bcr_length_spatial <- preprocess_spatial_airr(
  input = bcr_raw,
  barcode_map = barcode_map,
  loci = "IGH",
  out_tsv = file.path(presto_input_dir, sprintf("%s_igblast_bcr_pass_length_spatial_barcode.tsv", sample_lower)),
  out_fasta = file.path(presto_input_dir, sprintf("%s_igblast_bcr_pass_length_spatial_barcode.fasta", sample_lower)),
  min_v_length = 100,
  min_j_length = 10,
  min_junction_length = 1,
  sequence_col = "sequence_alignment"
)

tcr_length_spatial <- preprocess_spatial_airr(
  input = tcr_raw,
  barcode_map = barcode_map,
  loci = "TRB",
  out_tsv = file.path(presto_input_dir, sprintf("%s_igblast_tcr_pass_length_spatial_barcode.tsv", sample_lower)),
  out_fasta = file.path(presto_input_dir, sprintf("%s_igblast_tcr_pass_length_spatial_barcode.fasta", sample_lower)),
  min_v_length = 100,
  min_j_length = 10,
  min_junction_length = 1,
  sequence_col = "sequence_alignment"
)

spatial_receptor_summary(bcr_length_spatial, barcode_col = "barcode_2um")
spatial_receptor_summary(tcr_length_spatial, barcode_col = "barcode_2um")

pre_presto_qc_plot <- cowplot::plot_grid(
  plot_length_qc(bcr_length_spatial, title = paste(sample_id, "raw IGH reads after length/spatial filtering")),
  plot_length_qc(tcr_length_spatial, title = paste(sample_id, "raw TRB reads after length/spatial filtering")),
  ncol = 1
)
pre_presto_qc_plot
save_plot(pre_presto_qc_plot, sprintf("%s_pre_presto_length_qc.png", sample_lower), width = 10, height = 8)

# Things to inspect:
head(bcr_length_spatial)
head(tcr_length_spatial)

# Notes:
# - TODO: Should preprocess_spatial_airr() print a drop-count summary?
# - TODO: Should full_id include barcode16um too?

# ---- Run pRESTO and IgBLAST from R ------------------------------------------

# This is the command-line part of the original workflow, wrapped in R:
# ClusterSets.py -> ParseHeaders.py -> AlignSets.py -> BuildConsensus.py ->
# AssignGenes.py -> MakeDb.py. The command output prints to the R console.

bcr_length_spatial_fasta <- file.path(presto_input_dir, sprintf("%s_igblast_bcr_pass_length_spatial_barcode.fasta", sample_lower))
tcr_length_spatial_fasta <- file.path(presto_input_dir, sprintf("%s_igblast_tcr_pass_length_spatial_barcode.fasta", sample_lower))

if (run_presto_igblast) {
  bcr_presto_paths <- run_presto_igblast_sample(
    input_fasta = bcr_length_spatial_fasta,
    output_dir = presto_output_dir,
    receptor = "ig",
    identity = identity,
    consensus_freq = consensus_freq,
    germline_glob = germline_glob,
    igblast_base = igblast_base,
    igblast_exec = igblast_exec,
    muscle_exec = muscle_exec,
    nproc = nproc
  )

  tcr_presto_paths <- run_presto_igblast_sample(
    input_fasta = tcr_length_spatial_fasta,
    output_dir = presto_output_dir,
    receptor = "tr",
    identity = identity,
    consensus_freq = consensus_freq,
    germline_glob = germline_glob,
    igblast_base = igblast_base,
    igblast_exec = igblast_exec,
    muscle_exec = muscle_exec,
    nproc = nproc
  )
} else {
  message("run_presto_igblast is FALSE, so I am not running the external commands in this pass.")
}

if (!exists("bcr_presto_paths")) {
  bcr_presto_paths <- expected_presto_paths(bcr_length_spatial_fasta, presto_output_dir, identity, consensus_freq)
}
if (!exists("tcr_presto_paths")) {
  tcr_presto_paths <- expected_presto_paths(tcr_length_spatial_fasta, presto_output_dir, identity, consensus_freq)
}

bcr_presto_paths$db_pass
tcr_presto_paths$db_pass

# Notes:
# - TODO: Should SIRAP expose expected_presto_paths() as a real package function?
# - TODO: Should run_presto_igblast_sample() check all external executables first?

# ---- Polish the pRESTO/IgBLAST consensus output -----------------------------

# Now we are back in R. The original script read the db-pass TSVs, kept strong
# productive consensus sequences, recovered the barcode from the pRESTO header,
# swapped the 2 um barcode row/column order, remapped to 8 um/16 um, and kept
# the top consensus sequence per 2 um barcode and locus.

check_files(c(bcr_presto_paths$db_pass, tcr_presto_paths$db_pass))

bcr_consensus <- read_spatial_airr(bcr_presto_paths$db_pass)
tcr_consensus <- read_spatial_airr(tcr_presto_paths$db_pass)

bcr_consensus_filtered <- filter_receptor_reads(
  bcr_consensus,
  loci = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 11,
  productive_only = TRUE,
  min_consensus_count = 2,
  max_opening_gap = 2
)

tcr_consensus_filtered <- filter_receptor_reads(
  tcr_consensus,
  loci = "TRB",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 11,
  productive_only = TRUE,
  min_consensus_count = 2,
  max_opening_gap = 2
)

bcr_spatial_consensus <- add_spatial_barcodes(bcr_consensus_filtered, barcode_map, mode = "presto_header", swap = TRUE)
tcr_spatial_consensus <- add_spatial_barcodes(tcr_consensus_filtered, barcode_map, mode = "presto_header", swap = TRUE)

bcr_top_barcode <- collapse_top_consensus_per_barcode(bcr_spatial_consensus, barcode_col = "barcode_2um", locus_col = "locus", ties = "all")
tcr_top_barcode <- collapse_top_consensus_per_barcode(tcr_spatial_consensus, barcode_col = "barcode_2um", locus_col = "locus", ties = "all")

bcr <- add_tissue_positions(
  bcr_top_barcode,
  positions,
  barcode_col = "barcode_8um",
  position_barcode_col = "barcode_8um",
  x_col = "x",
  y_col = "y"
) |>
  filter(!is.na(x), !is.na(y))

tcr <- add_tissue_positions(
  tcr_top_barcode,
  positions,
  barcode_col = "barcode_8um",
  position_barcode_col = "barcode_8um",
  x_col = "x",
  y_col = "y"
) |>
  filter(!is.na(x), !is.na(y))

final_bcr_path <- file.path(final_airr_dir, sprintf("%s_bcrs_heavy_final_spatial_airr.tsv", sample_lower))
final_tcr_path <- file.path(final_airr_dir, sprintf("%s_tcrs_beta_final_spatial_airr.tsv", sample_lower))
write_spatial_outputs(bcr, out_tsv = final_bcr_path)
write_spatial_outputs(tcr, out_tsv = final_tcr_path)

list(
  bcr = spatial_receptor_summary(bcr, barcode_col = "barcode_8um"),
  tcr = spatial_receptor_summary(tcr, barcode_col = "barcode_8um")
)

post_presto_qc_plot <- cowplot::plot_grid(
  plot_length_qc(bcr, title = paste(sample_id, "final BCRs")),
  plot_length_qc(tcr, title = paste(sample_id, "final TCRs")),
  ncol = 1
)
post_presto_qc_plot
save_plot(post_presto_qc_plot, sprintf("%s_post_presto_final_length_qc.png", sample_lower), width = 10, height = 8)

# Things to inspect:
head(bcr)
head(tcr)

# Notes:
# - TODO: Should the post-pRESTO block become a single high-level function?
# - TODO: Should collapse_top_consensus_per_barcode() have built-in diagnostics for ties?

# ---- Put every BCR and TCR on the tissue slice -------------------------------

all_bcr_plot <- plot_receptor_overlay(
  tissue_positions = positions,
  receptors = bcr,
  receptor_color = "#2c7fb8",
  receptor_size = 1.2,
  title = paste(sample_id, "all BCR-positive spots")
)
all_bcr_plot

all_tcr_plot <- plot_receptor_overlay(
  tissue_positions = positions,
  receptors = tcr,
  receptor_color = "#d95f02",
  receptor_size = 1.2,
  title = paste(sample_id, "all TCR-positive spots")
)
all_tcr_plot

bcr_tcr_overlay <- ggplot(positions, aes(x = x, y = y)) +
  geom_point(color = "grey88", size = 0.55, alpha = 0.75) +
  geom_point(data = bcr, aes(x = x, y = y), inherit.aes = FALSE, shape = 21, fill = NA, color = "#2c7fb8", stroke = 0.45, size = 1.4) +
  geom_point(data = tcr, aes(x = x, y = y), inherit.aes = FALSE, shape = 21, fill = NA, color = "#d95f02", stroke = 0.45, size = 1.4) +
  coord_fixed() +
  scale_y_reverse() +
  theme_minimal() +
  labs(title = paste(sample_id, "BCR and TCR locations"), x = "x", y = "y")
bcr_tcr_overlay

save_plot(all_bcr_plot, sprintf("%s_all_bcr_locations.png", sample_lower), width = 6, height = 5)
save_plot(all_tcr_plot, sprintf("%s_all_tcr_locations.png", sample_lower), width = 6, height = 5)
save_plot(bcr_tcr_overlay, sprintf("%s_bcr_tcr_locations.png", sample_lower), width = 6, height = 5)

# Notes:
# - TODO: Should plot_receptor_overlay() support BCR and TCR layers at the same time?

# ---- Add gene expression underneath receptor locations -----------------------

thymus <- NULL
gene_overlay_plots <- list()

if (load_seurat_object && file.exists(seurat_path) && requireNamespace("SeuratObject", quietly = TRUE)) {
  thymus <- readRDS(seurat_path)

  for (gene in genes_to_show) {
    if (!gene %in% rownames(thymus)) {
      message("Skipping ", gene, ": not found in the Seurat object.")
      next
    }

    expr <- SeuratObject::FetchData(thymus, vars = gene)
    expr$barcode_8um <- rownames(expr)

    tissue_gene <- positions |>
      left_join(expr, by = "barcode_8um") |>
      mutate(expression = .data[[gene]])

    bcr_gene_plot <- plot_receptor_overlay(
      tissue_positions = tissue_gene,
      receptors = bcr,
      gene_col = "expression",
      receptor_color = "#2c7fb8",
      receptor_size = 1.25,
      title = paste(gene, "expression with BCR spots")
    )

    tcr_gene_plot <- plot_receptor_overlay(
      tissue_positions = tissue_gene,
      receptors = tcr,
      gene_col = "expression",
      receptor_color = "#d95f02",
      receptor_size = 1.25,
      title = paste(gene, "expression with TCR spots")
    )

    gene_overlay_plots[[gene]] <- cowplot::plot_grid(bcr_gene_plot, tcr_gene_plot, ncol = 2)
    print(gene_overlay_plots[[gene]])
    save_plot(gene_overlay_plots[[gene]], sprintf("%s_%s_receptor_gene_overlay.png", sample_lower, gene), width = 11, height = 5)
  }
} else {
  message("Skipping gene overlays. Check load_seurat_object, seurat_path, and SeuratObject.")
}

# Notes:
# - TODO: Should SIRAP include a helper that fetches Seurat expression and joins tissue positions?

# ---- Call clones -------------------------------------------------------------

bcr_clones_path <- file.path(out_dir, sprintf("%s_bcr_hierarchical_clones.rds", sample_lower))
tcr_clones_path <- file.path(out_dir, sprintf("%s_tcr_hierarchical_clones.rds", sample_lower))
bcr_trees_path <- file.path(out_dir, sprintf("%s_bcr_pml_trees.rds", sample_lower))

if (run_clone_calling) {
  bcr_clones <- assign_spatial_clones(
    bcr,
    threshold = 0.1,
    productive_only = TRUE,
    cell_id = "sequence_id",
    nproc = nproc
  )

  tcr_clones <- assign_spatial_clones(
    tcr,
    threshold = 0.1,
    productive_only = TRUE,
    cell_id = "sequence_id",
    nproc = nproc
  )

  saveRDS(bcr_clones, bcr_clones_path)
  saveRDS(tcr_clones, tcr_clones_path)
} else if (file.exists(bcr_clones_path) && file.exists(tcr_clones_path)) {
  bcr_clones <- readRDS(bcr_clones_path)
  tcr_clones <- readRDS(tcr_clones_path)
} else {
  message("Clone RDS files not found. Set run_clone_calling <- TRUE and rerun this section.")
}

if (exists("bcr_clones")) {
  bcr_clone_counts <- bcr_clones |>
    filter(!is.na(clone_id)) |>
    count(clone_id, sort = TRUE)
  head(bcr_clone_counts, 20)
}

if (exists("tcr_clones")) {
  tcr_clone_counts <- tcr_clones |>
    filter(!is.na(clone_id)) |>
    count(clone_id, sort = TRUE)
  head(tcr_clone_counts, 20)
}

# Notes:
# - TODO: Should assign_spatial_clones() have better defaults for spatial data?

# ---- Pick a couple of clones and plot them ----------------------------------

choose_big_clones <- function(clones, clone_ids = NULL, n = 2) {
  if (!is.null(clone_ids)) {
    return(as.character(clone_ids))
  }
  clones |>
    filter(!is.na(clone_id)) |>
    count(clone_id, sort = TRUE) |>
    slice_head(n = n) |>
    pull(clone_id) |>
    as.character()
}

if (exists("bcr_clones")) {
  bcr_clone_ids <- choose_big_clones(bcr_clones, bcr_clone_ids_to_plot, n = 2)
  bcr_clone_ids
}

if (exists("tcr_clones")) {
  tcr_clone_ids <- choose_big_clones(tcr_clones, tcr_clone_ids_to_plot, n = 2)
  tcr_clone_ids
}

# After inspecting clone counts, set these manually if you want:
# bcr_clone_ids_to_plot <- c("123", "456")
# tcr_clone_ids_to_plot <- c("789", "101")

if (exists("bcr_clones") && exists("bcr_clone_ids")) {
  bcr_clone_plots <- lapply(
    bcr_clone_ids,
    function(clone_id) {
      plot_receptor_overlay(
        tissue_positions = positions,
        receptors = filter(bcr_clones, .data$clone_id == .env$clone_id),
        receptor_color = "#2c7fb8",
        receptor_size = 1.7,
        title = paste("BCR clone", clone_id)
      )
    }
  )
  selected_bcr_clone_plot <- cowplot::plot_grid(plotlist = bcr_clone_plots, ncol = 2)
  selected_bcr_clone_plot
  save_plot(selected_bcr_clone_plot, sprintf("%s_selected_bcr_clones_spatial.png", sample_lower), width = 10, height = 5)
}

if (exists("tcr_clones") && exists("tcr_clone_ids")) {
  tcr_clone_plots <- lapply(
    tcr_clone_ids,
    function(clone_id) {
      plot_receptor_overlay(
        tissue_positions = positions,
        receptors = filter(tcr_clones, .data$clone_id == .env$clone_id),
        receptor_color = "#d95f02",
        receptor_size = 1.7,
        title = paste("TCR clone", clone_id)
      )
    }
  )
  selected_tcr_clone_plot <- cowplot::plot_grid(plotlist = tcr_clone_plots, ncol = 2)
  selected_tcr_clone_plot
  save_plot(selected_tcr_clone_plot, sprintf("%s_selected_tcr_clones_spatial.png", sample_lower), width = 10, height = 5)
}

# ---- Build BCR trees and match tips back to space ----------------------------

if (build_bcr_trees) {
  if (!exists("bcr_clones") || !exists("bcr_clone_ids")) {
    stop("Run clone calling/loading and clone selection before tree building.")
  }
  if (!requireNamespace("dowser", quietly = TRUE)) {
    stop("Install dowser to build BCR lineage trees.")
  }
  if (!requireNamespace("alakazam", quietly = TRUE)) {
    stop("Install alakazam to read IMGT germlines for dowser tree building.")
  }

  refs <- alakazam::readIMGT(germline_dir)

  bcr_tree_input <- bcr_clones |>
    filter(clone_id %in% bcr_clone_ids, !is.na(clone_id)) |>
    mutate(
      sequence_id_tree = make.unique(as.character(sequence_id)),
      c_call_simple = simplify_igh_constant_call(c_call)
    ) |>
    assign_tip_colors(clone_ids = bcr_clone_ids, id_col = "sequence_id_tree")

  bcr_tree_input <- dowser::createGermlines(bcr_tree_input, refs, id = "sequence_id_tree")

  bcr_tree_result <- build_dowser_trees(
    bcr_tree_input,
    clone_ids = bcr_clone_ids,
    id = "sequence_id_tree",
    traits = c("tip_color_id", "tip_color", "c_call_simple", "barcode_8um"),
    columns = c("x", "y", "barcode_8um"),
    method = "pml",
    nproc = nproc
  )

  saveRDS(bcr_tree_result, bcr_trees_path)
} else if (file.exists(bcr_trees_path)) {
  bcr_tree_result <- readRDS(bcr_trees_path)
}

if (exists("bcr_tree_result") && exists("bcr_tree_input")) {
  tree_palette <- unique(bcr_tree_input$tip_color)
  names(tree_palette) <- unique(bcr_tree_input$tip_color_id)

  bcr_tree_plots <- dowser::plotTrees(
    bcr_tree_result$trees,
    tips = "tip_color_id",
    palette = tree_palette,
    tipsize = if ("consensus_count" %in% colnames(bcr_tree_input)) "consensus_count" else NULL
  )

  paired_tree_spatial_plots <- list()
  for (i in seq_along(bcr_clone_ids)) {
    clone_id <- bcr_clone_ids[[i]]
    spatial_plot <- plot_colored_clone_spatial(
      bcr_tree_input,
      clone_id = clone_id,
      tissue_positions = positions
    )

    paired_tree_spatial_plots[[i]] <- cowplot::plot_grid(
      bcr_tree_plots[[i]] + ggplot2::ggtitle(paste("BCR clone", clone_id, "tree")),
      spatial_plot + ggplot2::ggtitle(paste("BCR clone", clone_id, "matched tissue colors")),
      ncol = 2,
      rel_widths = c(1.2, 1)
    )
    print(paired_tree_spatial_plots[[i]])
  }

  if (save_plots) {
    pdf(file.path(plot_dir, sprintf("%s_bcr_tree_spatial_pairs.pdf", sample_lower)), width = 12, height = 6)
    for (page in paired_tree_spatial_plots) {
      print(page)
    }
    dev.off()
  }
}

# Notes:
# - TODO: Should build_dowser_trees() handle createGermlines()?
# - TODO: Should plotTreeSpatialPairs() become a package function?

# ---- Done -------------------------------------------------------------------

message("Walkthrough objects are loaded in the current R session.")
message("Output directory: ", out_dir)
