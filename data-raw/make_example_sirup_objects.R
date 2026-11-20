# Build the saved SIRUP object in example-data/ that most guides load with
# readSIRUP(). Only the "Get started", "How the SIRUP object works", and
# "multisample" guides show how to build an object. The other guides start from
# this file.
#
# mg_thymus_slice2_brighter.sirup.rds is MG thymus section 2, with the brightened H&E image
# (see brighten_example_images.R). It has a BCR slot and a TCR
# slot, the cell2location cell-type proportions in the Seurat metadata, and BCR
# clones. It has no polygons.
#
# The code follows the "Get started" guide and uses the example files from
# make_example_data.R.
#
# Run it from the SIRUP package folder:
#   Rscript data-raw/make_example_sirup_objects.R

suppressPackageStartupMessages({
  library(dplyr)
  devtools::load_all(".", quiet = TRUE)
})

example_dir <- "example-data"

# --- Section 2 ---------------------------------------------------------------

thymus <- readRDS(file.path(example_dir, "mg_thymus_slice2_spatial_brighter.rds"))

cell_types <- data.table::fread(
  file.path(example_dir, "mg_thymus_cell_types.csv.gz"),
  data.table = FALSE
) %>%
  filter(sample == "mg_thymus_slice2")
rownames(cell_types) <- cell_types$barcode
thymus <- Seurat::AddMetaData(thymus, select(cell_types, -sample, -barcode))

bcrs <- airr::read_rearrangement(file.path(example_dir, "mg_thymus_slice2_bcr.tsv.gz")) %>%
  extractSpatialIds(
    id = "sequence_id",
    spatial_id = "barcode_10x",
    regex = ",([ACGT]{16})$",
    suffix = "-1"
  ) %>%
  mutate(sequence_id = gsub(",", "_", sequence_id)) %>%
  filterAirr(
    locus_val = "IGH",
    min_v_length = 250,
    min_j_length = 30,
    min_junction_length = 10,
    productive_only = TRUE,
    min_consensus_count = NULL,
    spatial_id = "barcode_10x"
  ) %>%
  filter(!grepl("[^ACGT]", junction))

tcrs <- airr::read_rearrangement(file.path(example_dir, "mg_thymus_slice2_tcr.tsv.gz")) %>%
  extractSpatialIds(
    id = "sequence_id",
    spatial_id = "barcode_10x",
    regex = ",([ACGT]{16})$",
    suffix = "-1"
  ) %>%
  mutate(sequence_id = gsub(",", "_", sequence_id)) %>%
  filterAirr(
    locus_val = "TRB",
    min_v_length = 250,
    min_j_length = 30,
    min_junction_length = 10,
    productive_only = TRUE,
    min_consensus_count = NULL,
    spatial_id = "barcode_10x"
  ) %>%
  filter(!grepl("[^ACGT]", junction))

sirup <- createSIRUP(
  spatial = thymus,
  receptors = list(BCR = bcrs, TCR = tcrs),
  spatial_id = "barcode_10x",
  image_scale = "lowres",
  swap_coords = FALSE
)
sirup <- updateReceptors(
  sirup,
  scoper::hierarchicalClones,
  slot = "BCR",
  threshold = 0.1,
  summarize_clones = FALSE,
  nproc = 1
)
saveSIRUP(sirup, file.path(example_dir, "mg_thymus_slice2_brighter.sirup.rds"))
