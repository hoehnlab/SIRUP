# Figure 2: ectopic thymic GC clonality

library(SIRAP)
library(dplyr)
library(ggplot2)
library(data.table)
library(airr)
library(Seurat)
library(SeuratObject)
library(scoper)

sample_id <- Sys.getenv("SIRAP_SAMPLE_ID", "M5")
sample_lower <- tolower(sample_id)
nproc <- as.integer(Sys.getenv("SIRAP_NPROC", "8"))

thymus_dir <- Sys.getenv(
  "SIRAP_THYMUS_DIR",
  "/dartfs/rc/lab/H/HoehnK/Hunter/Thymus_VisiumHD"
)
walkthrough_dir <- file.path(thymus_dir, "SIRAP_walkthrough", sample_id)
final_airr_dir <- file.path(walkthrough_dir, "03_final_spatial_airr")
manuscript_out_dir <- file.path(thymus_dir, "SIRAP_manuscript", "figure_02", sample_id)
dir.create(manuscript_out_dir, recursive = TRUE, showWarnings = FALSE)

seurat_path <- file.path(thymus_dir, "data", sprintf("Thymus_VisiumHD_%s_8um_hires.rds", sample_id))
final_bcr_path <- file.path(final_airr_dir, sprintf("%s_bcrs_heavy_final_spatial_airr.tsv", sample_lower))
bcr_clones_path <- file.path(walkthrough_dir, sprintf("%s_bcr_hierarchical_clones.rds", sample_lower))
all_gcs_polygons_path <- file.path(
  walkthrough_dir,
  "plots",
  "clonal_dominance",
  sprintf("%s_all_gcs_polygons_vertices.csv", sample_lower)
)

stopifnot(file.exists(seurat_path), file.exists(final_bcr_path))

thymus <- readRDS(seurat_path)
bcr_top_barcode <- airr::read_rearrangement(final_bcr_path)

if (file.exists(bcr_clones_path)) {
  bcr_clones <- readRDS(bcr_clones_path)
} else {
  bcr_clones <- scoper::hierarchicalClones(
    filter(bcr_top_barcode, productive),
    threshold = 0.1,
    cell_id = "sequence_id",
    summarize_clones = FALSE,
    nproc = nproc
  )
  saveRDS(bcr_clones, bcr_clones_path)
}

clone_counts <- bcr_clones |>
  filter(!is.na(clone_id)) |>
  count(clone_id, sort = TRUE)

utils::write.table(
  clone_counts,
  file.path(manuscript_out_dir, "bcr_clone_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!file.exists(all_gcs_polygons_path)) {
  stop(
    "Missing GC polygon file: ", all_gcs_polygons_path,
    "\nCreate polygons with selectSpatialPolygons() before running Figure 2."
  )
}

all_gcs_polygons <- data.table::fread(all_gcs_polygons_path) |> as.data.frame()

polygon_overlay <- plotSpatialPolygons(
  object = thymus,
  polygon = all_gcs_polygons,
  image.scale = "hires",
  polygon_alpha = 0.2
)
ggsave(file.path(manuscript_out_dir, "all_gcs_polygons_overlay.pdf"), polygon_overlay, width = 7, height = 6)

bcr_for_region <- bcr_clones |>
  mutate(
    spaced_clone_id = paste0("heavy_", clone_id),
    chain_for_pairing = "heavy"
  )

dominance <- spatialCloneDominance(
  bcrs = bcr_for_region,
  object = thymus,
  polygon = all_gcs_polygons,
  barcode_col = "barcode_8um",
  clone_col = "spaced_clone_id",
  locus_col = "chain_for_pairing",
  locus = "heavy"
)

utils::write.table(
  dominance$dominance,
  file.path(manuscript_out_dir, "gc_clone_dominance_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

join_count <- spatialCloneJoinCount(
  object = thymus,
  receptors = bcr_for_region,
  clone_col = "spaced_clone_id",
  barcode_col = "barcode_8um",
  polygon = all_gcs_polygons,
  min_barcodes = 3,
  k = 6,
  style = "B"
)

utils::write.table(
  join_count$summary,
  file.path(manuscript_out_dir, "gc_clone_join_count_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
