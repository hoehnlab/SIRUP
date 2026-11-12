# Build the example data shipped in inst/extdata and used by the vignettes.
#
# Source: the two myasthenia gravis (MG) thymus Visium sections analysed in the
# SIRUP manuscript, their UMI-supported BCR (IGH) and TCR (TRB) consensus
# tables, matched bulk BCR-seq from thymus and blood, cell2location cell-type
# proportions, and the hand-drawn ectopic germinal centre polygons.
#
# Receptor files are the pipeline outputs themselves. The script copies them,
# or copies a subset of their lines, without changing any row:
#   * spatial BCRs: the full `*_umi_support_2.tsv` consensus tables of both
#     sections;
#   * spatial TCRs: the full TRB consensus table of section 2;
#   * bulk BCRs: lines of the formatted bulk AIRR tables. The manuscript clone
#     assignments choose the lines: every blood row of a clone seen in a
#     section, and up to 10 bulk thymus rows of each clone that reaches blood;
#   * reads: lines of the RAD IgBLAST read table of section 2, enriched for
#     reads that pass the read-level filter.
#
# Spatial objects are slimmed: a marker panel plus variable genes, Spatial and
# ALRA assays, and the H&E reduced to Space Ranger lowres size (polygon
# vertices are rescaled to match).
#
# Run on Discovery (needs ~40 GB of memory):
#   module load R/4.5.1; Rscript data-raw/make_example_data.R

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(Seurat)
  library(SeuratObject)
})

repo <- "/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package"
devtools::load_all(file.path(repo, "SIRUP"), quiet = TRUE)
out_dir <- file.path(repo, "SIRUP", "inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20261012)

umi_dir <- file.path(repo, "analysis/processing_umi/output")
umi_file <- function(dataset, chain) {
  list.files(file.path(umi_dir, dataset, chain), pattern = "_umi_support_2[.]tsv$", full.names = TRUE)
}
bcr_path <- file.path(repo, "analysis/fig_5_mg_umi/01_build_object/analysis/mg_thymus_multisample.rds")
pooled_clones_path <- file.path(repo, "analysis/fig_5_mg_umi/01_build_object/analysis/pooled_clones_from_raw.rds")
bulk_dir <- "/dartfs/rc/lab/H/HoehnK/Hunter/Thymus_Spatial/BulkBCRs"
bulk_files <- c(
  bulk_thymus_1 = "17_formatted_airr.tsv", bulk_thymus_2 = "45_formatted_airr.tsv",
  bulk_blood_month0_1 = "21_formatted_airr.tsv", bulk_blood_month0_2 = "25_formatted_airr.tsv",
  bulk_blood_month12_1 = "23_formatted_airr.tsv", bulk_blood_month12_2 = "33_formatted_airr.tsv"
)
c2l_dir <- file.path(repo, "analysis/fig_2_gex/02_deconvolution/analysis/cell2location_analysis/results")
polygon_dir <- file.path(repo, "analysis/fig_2_gex/01_build_objects/inputs/gc_polygons")
raw_reads_path <- file.path(repo, "analysis/processing_10_8/rad_igblast_data/KH02_BCR.airr.tsv")

samples <- c("mg_thymus_slice1", "mg_thymus_slice2")
short_names <- c(mg_thymus_slice1 = "thymus_slice1", mg_thymus_slice2 = "thymus_slice2")

# Copy a text file, or selected lines of it, to a gzipped file. Lines are
# written unchanged; `keep` is a logical vector over the data lines.
copy_lines <- function(from, to, keep = NULL) {
  lines <- readLines(from)
  if (!is.null(keep)) lines <- c(lines[[1]], lines[-1][keep])
  con <- gzfile(file.path(out_dir, to), "w")
  writeLines(lines, con)
  close(con)
  message(to, ": ", length(lines) - 1L, " rows")
}

bcr_obj <- readSIRUP(bcr_path)

# ---- 1. Gene panel -------------------------------------------------------
marker_genes <- c(
  # B cells, germinal centre B cells, plasma cells
  "MS4A1", "CD19", "CD79A", "CR2", "FCER2", "AICDA", "BCL6", "RGS13", "LMO2",
  "MKI67", "TOP2A", "IGHM", "IGHD", "IGHG1", "IGHG3", "IGHA1", "JCHAIN",
  "MZB1", "XBP1", "IRF4", "PRDM1", "CXCR4", "CD83",
  # T cells and T follicular helper cells
  "CD3E", "CD4", "CD8A", "FOXP3", "CCR7", "CXCL13", "CXCR5", "PDCD1", "ICOS",
  "GZMK", "NKG7", "CCL5", "PTPRC",
  # Thymic cortex, medulla and epithelium
  "CD1A", "CD1E", "RAG1", "PSMB11", "PRSS16", "AIRE", "KRT5", "KRT14",
  "KRT15", "CCL19", "CCL21", "CHRNA1", "EPCAM",
  # Type-I interferon response
  "IFI44", "IFIT3", "MX1", "IRF9", "ISG15",
  # Myeloid and stromal
  "CD68", "LYZ", "COL1A1"
)

variable_genes <- character()
for (s in samples) {
  sp <- getSpatial(bcr_obj, s)
  DefaultAssay(sp) <- "Spatial"
  sp <- FindVariableFeatures(sp, selection.method = "vst", nfeatures = 300, verbose = FALSE)
  variable_genes <- union(variable_genes, VariableFeatures(sp))
}
all_genes <- rownames(getSpatial(bcr_obj, samples[[1]])[["Spatial"]])
marker_genes <- intersect(marker_genes, all_genes)
keep_genes <- union(marker_genes, intersect(variable_genes, all_genes))
message("genes kept: ", length(keep_genes), " (", length(marker_genes), " markers)")

# ---- 2. Slim Seurat objects with a lowres image --------------------------
for (s in samples) {
  sp <- getSpatial(bcr_obj, s)
  image_name <- Images(sp)[[1]]

  spatial_counts <- LayerData(sp, assay = "Spatial", layer = "counts")[keep_genes, ]
  spatial_data <- LayerData(sp, assay = "Spatial", layer = "data")[keep_genes, ]
  alra_data <- LayerData(sp, assay = "alra", layer = "data")[marker_genes, ]
  alra_data <- as(round(as.matrix(alra_data), 3), "dgCMatrix")

  slim <- CreateSeuratObject(counts = spatial_counts, assay = "Spatial", project = short_names[[s]])
  slim <- SetAssayData(slim, assay = "Spatial", layer = "data", new.data = spatial_data)
  slim[["alra"]] <- CreateAssay5Object(data = alra_data)
  DefaultAssay(slim) <- "Spatial"
  slim$nCount_Spatial <- sp$nCount_Spatial[colnames(slim)]
  slim$nFeature_Spatial <- sp$nFeature_Spatial[colnames(slim)]
  slim$seurat_clusters <- sp$seurat_clusters[colnames(slim)]

  # Reduce the hires raster to Space Ranger lowres size by block averaging.
  img <- sp@images[[image_name]]
  hires <- img@image
  shrink <- img@scale.factors$hires / img@scale.factors$lowres
  row_group <- ceiling(seq_len(dim(hires)[[1]]) / shrink)
  col_group <- ceiling(seq_len(dim(hires)[[2]]) / shrink)
  lowres <- array(0, dim = c(max(row_group), max(col_group), 3))
  for (channel in 1:3) {
    m <- rowsum(hires[, , channel], row_group) / as.vector(table(row_group))
    m <- t(rowsum(t(m), col_group) / as.vector(table(col_group)))
    lowres[, , channel] <- m
  }
  img@image <- round(lowres * 255) / 255
  slim[[image_name]] <- img
  message(s, ": ", ncol(slim), " spots, image ", paste(dim(img@image), collapse = " x "))

  saveRDS(slim, file.path(out_dir, paste0(s, "_spatial.rds")), compress = "xz")
}

# ---- 3. Cell-type proportions -------------------------------------------
c2l <- list()
for (s in samples) {
  props <- fread(file.path(c2l_dir, short_names[[s]], "cell2location_proportions.csv"))
  setnames(props, 1, "barcode")
  props <- props[barcode %in% colnames(getSpatial(bcr_obj, s))]
  numeric_cols <- setdiff(names(props), "barcode")
  props[, (numeric_cols) := lapply(.SD, round, 4), .SDcols = numeric_cols]
  props[, sample := s]
  c2l[[s]] <- props
}
c2l <- rbindlist(c2l, use.names = TRUE)
setcolorder(c2l, c("sample", "barcode"))
fwrite(c2l, file.path(out_dir, "mg_thymus_cell_types.csv.gz"))
message("cell-type table: ", nrow(c2l), " spots x ", ncol(c2l) - 2, " cell types")

# ---- 4. GC polygons rescaled to lowres ----------------------------------
for (s in samples) {
  poly <- fread(file.path(polygon_dir, s, paste0(s, "_gc_polygons_vertices.csv")))
  sf <- getSpatial(bcr_obj, s)@images[[1]]@scale.factors
  poly[, x := round(x * sf$lowres / sf$hires, 3)]
  poly[, y := round(y * sf$lowres / sf$hires, 3)]
  fwrite(poly, file.path(out_dir, paste0(s, "_gc_polygons.csv")))
}

# ---- 5. Spatial BCR and TCR consensus tables (copied unchanged, except the 8 rows below) ----
copy_lines(umi_file("kh01", "bcr"), "mg_thymus_slice1_bcr.tsv.gz")
# Section 2 BCRs: leave out the rows whose junction has a character other than
# A, C, G, T, N, or ?. These are 8 sequences with a gap ("-"). scoper cannot read them
# and stops with an error, so the vignettes could not call clones on this table.
bcr2_file <- umi_file("kh02", "bcr")
bcr2_fields <- strsplit(readLines(bcr2_file), "\t", fixed = TRUE)
bcr2_junction_col <- match("junction", bcr2_fields[[1]])
bcr2_junction <- vapply(bcr2_fields[-1], function(x) x[[bcr2_junction_col]], "")
copy_lines(bcr2_file, "mg_thymus_slice2_bcr.tsv.gz", keep = !grepl("[^ACGTN?]", bcr2_junction))
copy_lines(umi_file("kh02", "tcr"), "mg_thymus_slice2_tcr.tsv.gz")

# ---- 6. Bulk BCR lines of clones shared with the sections ----------------
# The manuscript clone calls only choose which lines to copy.
pooled <- readRDS(pooled_clones_path)
spatial_clones <- unique(pooled$clone_id[pooled$sample %in% samples])
blood <- pooled[grepl("blood", pooled$tissue_sample) & pooled$clone_id %in% spatial_clones, ]
bulk_thymus <- pooled[grepl("bulk_thymus", pooled$tissue_sample) & pooled$clone_id %in% unique(blood$clone_id), ]
bulk_thymus <- bulk_thymus %>% group_by(clone_id) %>% slice_sample(n = 10) %>% ungroup()
chosen <- bind_rows(blood, bulk_thymus)
for (s in names(bulk_files)) {
  path <- file.path(bulk_dir, bulk_files[[s]])
  ids <- fread(path, select = "sequence_id", showProgress = FALSE)$sequence_id
  keep_ids <- chosen$sequence_id[chosen$tissue_sample == s]
  copy_lines(path, paste0("mg_thymus_", s, ".tsv.gz"), keep = ids %in% keep_ids)
}

# ---- 7. Read-level IgBLAST lines, copied unchanged -----------------------
# Only ~0.4% of reads pass the read-level filter, so a plain random sample has
# too few passing reads for the QC plots. Copy up to 750 lines that pass the
# IGH filter used in the vignette and 1,250 that do not. fread() with `nrows`
# on the full 34 GB file segfaults (quoted fields with embedded tabs), so
# read the first 600,000 lines through `head`.
set.seed(20261012)
raw_cols <- c(
  "sequence_id", "locus", "productive", "v_sequence_start", "v_sequence_end",
  "j_sequence_start", "j_sequence_end", "junction_length", "sequence_alignment",
  "cell_barcode"
)
head_cmd <- paste("head -n 600001", shQuote(raw_reads_path))
raw <- fread(cmd = head_cmd, select = raw_cols, nThread = 1, showProgress = FALSE)
passing <- filterAirr(
  mutate(as.data.frame(raw), barcode = cell_barcode),
  locus_val = "IGH", min_v_length = 100, min_j_length = 10,
  min_junction_length = 1, productive_only = FALSE,
  min_consensus_count = NULL, spatial_id = "barcode"
)
is_passing <- raw$sequence_id %in% passing$sequence_id
message("read-level pass rate in the first 600,000 rows: ", round(100 * mean(is_passing), 2), "%")
rows <- sort(c(
  sample(which(is_passing), min(sum(is_passing), 750)),
  sample(which(!is_passing), 1250)
))
idx_file <- tempfile()
writeLines(as.character(rows + 1L), idx_file)
system(paste(
  head_cmd, "| awk 'NR==FNR {keep[$1]; next} FNR == 1 || (FNR in keep)'", shQuote(idx_file), "- |",
  "gzip >", shQuote(file.path(out_dir, "mg_thymus_slice2_reads.tsv.gz"))
))
message("mg_thymus_slice2_reads.tsv.gz: ", length(rows), " rows")

# ---- 8. Check: the slim objects reproduce the source polygon membership --
read_section <- function(file, sample) {
  airr::read_rearrangement(file.path(out_dir, file)) %>%
    extractSpatialIds(id = "sequence_id", spatial_id = "barcode_10x",
                      regex = ",([ACGT]{16})$", suffix = "-1") %>%
    mutate(sample = sample)
}
check <- createSIRUP(
  spatial = list(
    mg_thymus_slice1 = readRDS(file.path(out_dir, "mg_thymus_slice1_spatial.rds")),
    mg_thymus_slice2 = readRDS(file.path(out_dir, "mg_thymus_slice2_spatial.rds"))
  ),
  receptors = bind_rows(
    read_section("mg_thymus_slice1_bcr.tsv.gz", "mg_thymus_slice1"),
    read_section("mg_thymus_slice2_bcr.tsv.gz", "mg_thymus_slice2")
  ),
  spatial_id = "barcode_10x",
  sample_id = "sample",
  image_scale = "lowres",
  swap_coords = FALSE
)
check <- addPolygons(check, polygons_path = c(
  mg_thymus_slice1 = file.path(out_dir, "mg_thymus_slice1_gc_polygons.csv"),
  mg_thymus_slice2 = file.path(out_dir, "mg_thymus_slice2_gc_polygons.csv")
))
for (s in samples) {
  old_ids <- getSpatial(bcr_obj, s)$polygon_id
  new_ids <- getSpatial(check, s)$polygon_id[names(old_ids)]
  message(s, ": polygon membership identical = ", identical(unname(old_ids), unname(new_ids)))
}
print(check)
for (f in list.files(out_dir, pattern = "[.]tsv[.]gz$", full.names = TRUE)) {
  message(basename(f), ": valid AIRR = ", airr::validate_rearrangement(airr::read_rearrangement(f)))
}
print(file.info(list.files(out_dir, full.names = TRUE))[, "size", drop = FALSE])
