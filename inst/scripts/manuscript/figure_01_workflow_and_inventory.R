# Figure 1: SIRAP workflow and data overview

library(SIRAP)
library(dplyr)
library(ggplot2)
library(airr)
library(Seurat)
library(SeuratObject)

sample_id <- Sys.getenv("SIRAP_SAMPLE_ID", "M5")
sample_lower <- tolower(sample_id)

thymus_dir <- Sys.getenv(
  "SIRAP_THYMUS_DIR",
  "/dartfs/rc/lab/H/HoehnK/Hunter/Thymus_VisiumHD"
)
walkthrough_dir <- file.path(thymus_dir, "SIRAP_walkthrough", sample_id)
final_airr_dir <- file.path(walkthrough_dir, "03_final_spatial_airr")
manuscript_out_dir <- file.path(thymus_dir, "SIRAP_manuscript", "figure_01", sample_id)
dir.create(manuscript_out_dir, recursive = TRUE, showWarnings = FALSE)

seurat_path <- file.path(thymus_dir, "data", sprintf("Thymus_VisiumHD_%s_8um_hires.rds", sample_id))
final_bcr_path <- file.path(final_airr_dir, sprintf("%s_bcrs_heavy_final_spatial_airr.tsv", sample_lower))
final_tcr_path <- file.path(final_airr_dir, sprintf("%s_tcrs_beta_final_spatial_airr.tsv", sample_lower))

stopifnot(file.exists(seurat_path), file.exists(final_bcr_path), file.exists(final_tcr_path))

thymus <- readRDS(seurat_path)
bcr_top_barcode <- airr::read_rearrangement(final_bcr_path)
tcr_top_barcode <- airr::read_rearrangement(final_tcr_path)

receptor_summary <- tibble::tibble(
  sample_id = sample_id,
  receptor = c("BCR heavy", "TCR beta"),
  n_rows = c(nrow(bcr_top_barcode), nrow(tcr_top_barcode)),
  n_spatial_barcodes = c(
    dplyr::n_distinct(bcr_top_barcode$barcode_8um),
    dplyr::n_distinct(tcr_top_barcode$barcode_8um)
  ),
  n_productive = c(
    sum(bcr_top_barcode$productive %in% TRUE, na.rm = TRUE),
    sum(tcr_top_barcode$productive %in% TRUE, na.rm = TRUE)
  )
)

utils::write.table(
  receptor_summary,
  file.path(manuscript_out_dir, "receptor_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

bcr_spots <- unique(bcr_top_barcode$barcode_8um)
tcr_spots <- unique(tcr_top_barcode$barcode_8um)

thymus[["sirap_receptor_type"]] <- dplyr::case_when(
  colnames(thymus) %in% bcr_spots & colnames(thymus) %in% tcr_spots ~ "BCR + TCR",
  colnames(thymus) %in% bcr_spots ~ "BCR",
  colnames(thymus) %in% tcr_spots ~ "TCR",
  TRUE ~ "0"
)

if ("MS4A1" %in% rownames(thymus)) {
  bcr_ms4a1_plot <- plotReceptorExpression(
    object = thymus,
    receptors = bcr_top_barcode,
    gene = "MS4A1",
    barcode_col = "barcode_8um",
    receptor_color = "#2c7fb8",
    receptor_size = 0.55,
    receptor_stroke = 0.25,
    feature_colors = c("grey80", "#08306b"),
    feature_size = 0.45,
    title = paste(sample_id, "MS4A1 expression with BCR-positive spots"),
    receptor_label = "BCR"
  )
  ggsave(file.path(manuscript_out_dir, "bcr_ms4a1_overlay.pdf"), bcr_ms4a1_plot, width = 7, height = 6)
}

if ("CD3E" %in% rownames(thymus)) {
  tcr_cd3e_plot <- plotReceptorExpression(
    object = thymus,
    receptors = tcr_top_barcode,
    gene = "CD3E",
    barcode_col = "barcode_8um",
    receptor_color = "#d95f02",
    receptor_size = 0.55,
    receptor_stroke = 0.25,
    feature_colors = c("grey80", "#7f2704"),
    feature_size = 0.45,
    title = paste(sample_id, "CD3E expression with TCR-positive spots"),
    receptor_label = "TCR"
  )
  ggsave(file.path(manuscript_out_dir, "tcr_cd3e_overlay.pdf"), tcr_cd3e_plot, width = 7, height = 6)
}
