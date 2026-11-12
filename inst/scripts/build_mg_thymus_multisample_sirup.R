#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Demo: a MULTI-SAMPLE SIRUP object from the two MG thymus Visium slices.
#
# Shows the new machinery: `spatial` is a NAMED LIST of Seurat objects keyed by
# the receptor `sample` field, with ONE combined AIRR (receptor) table. Builds the
# object, (re)builds clones across both slices in one shared clone space, pulls
# receptors/spatial data for a chosen sample, makes a couple of plots that either
# select one section (images=) or facet across both, and saves the object.
#
# Run with the package loaded:  devtools::load_all("SIRUP")  (or library(SIRUP)).
# ---------------------------------------------------------------------------

setwd("/Volumes/HoehnK/Hunter/spatial_bcr_tcr_package/")

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(scoper); library(ggplot2); library(Seurat)
})
devtools::load_all("SIRUP")

repo <- "/Volumes/HoehnK/Hunter/spatial_bcr_tcr_package"
out_dir <- file.path(repo, "analysis", "fig_3_bcr_multisample", "demo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Inputs: the two slice Seurats + the MG thymus receptor (AIRR) table ----
# Spatial is a NAMED list; the names are the sample ids the receptors map to.
seurat_paths <- c(
  mg_thymus_slice1 = file.path(repo, "analysis/fig_2_gex/01_processing/objects/thymus_slice1_processed_seurat.rds"),
  mg_thymus_slice2 = file.path(repo, "analysis/fig_2_gex/01_processing/objects/thymus_slice2_processed_seurat.rds")
)
spatial_1 <- readRDS(seurat_paths["mg_thymus_slice1"])
spatial_2 <- readRDS(seurat_paths["mg_thymus_slice2"])

# The pooled MG thymus AIRR already carries a `sample` column
# (mg_thymus_slice1 / mg_thymus_slice2). Derive the 10x barcode so spatial_id
# matches each Seurat's colnames (ACGT...-1). Drop any existing clone_id so we can
# demo clone calling below.
airr <- readRDS(file.path(repo, "analysis/fig_3_bcr/01_clone_building/analysis/clones_thymus_spatial_only_igh.rds")) %>%
  as.data.frame() %>%
  select(-any_of("clone_id")) %>%
  mutate(barcode_10x = paste0(str_extract(sequence_id, "^[ACGT]+"), "-1"))


# --- 2. Build the multi-sample SIRUP object --------------------
# The multi-sample `spatial` argument is a NAMED LIST keyed by the receptor
# `sample` field; we build it inline from the two slice objects. `sample = "sample"`
# routes each receptor to its section's Seurat, and shared 16bp barcodes across
# slices are disambiguated internally by the composite (sample, spatial_id) key.
# No clone field yet -- we call clones below and attach them.
obj <- createSIRUP(
  spatial     = list(mg_thymus_slice1 = spatial_1, mg_thymus_slice2 = spatial_2),
  receptors   = airr,
  spatial_id  = "barcode_10x",
  sample_id      = "sample",
  locus       = "locus",
  count       = "consensus_count",
  image_scale = "hires",   # or per-sample: list(mg_thymus_slice1 = "hires", mg_thymus_slice2 = "hires")
  swap_coords = FALSE
)

obj
# <SIRUP object>
#  Spatial data:    Seurat, 6,599 spatial IDs, 36,601 genes, 2 spatial samples, 2 images
# Receptor data:   tibble, 5,312 receptor sequences x 62 columns
# Spatial overlap: 1,471/1,530 receptor IDs in spatial data (96.1%)
# Spatial ID field: barcode_10x
# Analyses:        0
summary(obj)

# --- 3. Explore the machinery: per-sample accessors (work pre-clones) ----------
message("spatialSamples(obj): ", paste(spatialSamples(obj), collapse = ", "))

slice1_seurat <- getSpatial(obj, "mg_thymus_slice1")   # one section's Seurat
message("slice1 spots: ", ncol(slice1_seurat))

all_receptors    <- getReceptors(obj)                   # one combined AIRR table
slice2_receptors <- getReceptors(obj, "mg_thymus_slice2") # one sample's receptor rows
message("receptors: ", nrow(all_receptors), " total; ", nrow(slice2_receptors), " in slice 2")

# spatialCoordinateTable() flattens all sections, each row tagged by `sample`.
coord_tbl <- spatialCoordinateTable(obj)
print(dplyr::count(coord_tbl, sample))

# --- 4. Build clones across BOTH slices, then attach to the object -------------
# Call IGH clones on the object's combined receptor table (one shared clone space);
# tissue_seq_id keeps every spatial sequence.

# Attach clone_id back and register it as the object's clone field.
obj <- updateReceptors(
  obj,
  hierarchicalClones, 
  threshold = 0.1,
  cell_id = "tissue_seq_id",
  summarize_clones = FALSE, 
  nproc = 4
)


# --- 5. Add the GC polygons (section-tagged) and run dominance -----------------
poly_csv <- c(
  mg_thymus_slice1 = file.path(repo, "analysis/fig_2_gex/02_polygons/analysis/gc_polygons/mg_thymus_slice1/mg_thymus_slice1_gc_polygons_vertices.csv"),
  mg_thymus_slice2 = file.path(repo, "analysis/fig_2_gex/02_polygons/analysis/gc_polygons/mg_thymus_slice2/mg_thymus_slice2_gc_polygons_vertices.csv")
)
# `addPolygons()` now accepts named polygon CSV paths for multi-sample objects.
# The names tag each CSV by sample as it is loaded, so polygons from slice 1 are
# only applied to `mg_thymus_slice1` and polygons from slice 2 are only applied
# to `mg_thymus_slice2`. If polygon IDs are reused across slices, SIRUP stores
# composite IDs internally (for example `mg_thymus_slice1::GC1`) to keep the
# geometry and annotations separate.
obj <- addPolygons(obj, polygons_path = poly_csv)

# Keep the dominance result inside the object (stored as an analysis), and grab a
# handle for inspection.
obj <- setAnalysis(obj, "gc_clonal_dominance",
                   spatialCloneDominance(obj, locus_val = "IGH"))
dominance <- getAnalysis(obj, "gc_clonal_dominance")
print(summary(obj))

# --- 6. Plots: one section (images=) or faceted across both -------------------
p_one  <- plotReceptorsExpression(obj, gene = "AICDA", group = "locus",
                                  assay = "alra", layer = "data",
                                  images = "mg_thymus_slice1")      # single section
p_both <- plotReceptorsExpression(obj, gene = "AICDA", group = "locus",
                                  assay = "alra", layer = "data",
                                  combine_sample_plots = TRUE)      # combined figure for saving
p_both_separate <- plotReceptorsExpression(obj, gene = "AICDA", group = "locus",
                                           assay = "alra", layer = "data",
                                           combine_sample_plots = FALSE) 
p_poly <- plotSpatialPolygons(obj, label = TRUE,
                              combine_sample_plots = TRUE)          # one combined figure (both sections)
p_poly_separate <- plotSpatialPolygons(obj, label = TRUE,
                                       combine_sample_plots = FALSE) # list: one ggplot per section

ggsave(file.path(out_dir, "demo_aicda_slice1.png"), p_one,  width = 6, height = 5, dpi = 200)
ggsave(file.path(out_dir, "demo_aicda_both.png"),  p_both, width = 10, height = 5, dpi = 200)
ggsave(file.path(out_dir, "demo_gc_polygons.png"), p_poly, width = 10, height = 5, dpi = 200)

# --- 7. Save the object -------------------------------------------------------
saveSIRUP(obj, file.path(out_dir, "mg_thymus_multisample_sirup.rds"))
message("Saved object + demo plots to: ", out_dir)
