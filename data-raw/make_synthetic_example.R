# Build the small synthetic SIRUP object used in the function examples.
#
# The object has no real data. It has:
#   * 25 spots on a 5 x 5 grid, with coordinates in the metadata columns `x`
#     and `y` (no tissue image);
#   * 10 genes, gene1 to gene10, with simulated counts. gene1 is high in the
#     upper right and gene2 is high in the lower left. The other genes have
#     random counts with no spatial pattern;
#   * a metadata column `cell_type` with the values "A" and "B", set at random;
#   * 5 clones (clone_id "1" to "5"). Each clone has 1 to 10 sequences, drawn
#     with a discrete uniform distribution. Each sequence goes to a random spot;
#   * real BCR sequences. Each synthetic clone gets the sequences of one real
#     clone from the myasthenia gravis thymus section 2 in inst/extdata. The
#     script keeps every AIRR column of these rows. It changes only the
#     sequence ID, the spatial ID, and the clone ID. The spots stay synthetic;
#   * 2 polygons: "upper_right" and "lower_left".
#
# Run from the package root:
#   Rscript data-raw/make_synthetic_example.R
#
# The germline step needs the IMGT human references. Set the environment
# variable SIRUP_IMGT_DIR to the folder that holds the V, D, and J FASTA files.

suppressPackageStartupMessages({
  library(SeuratObject)
  library(dplyr)
  library(tibble)
})
devtools::load_all(".", quiet = TRUE)
set.seed(20261106)

# ---- 1. Spots on a 5 x 5 grid ---------------------------------------------
# y grows downward, as in an image row index. Small y is the top of the plot.
spots <- sprintf("spot%02d", 1:25)
grid <- expand.grid(x = 1:5, y = 1:5)
grid$spot <- spots

# ---- 2. Gene counts --------------------------------------------------------
# u goes from 0 (lower left) to 1 (upper right).
u <- ((grid$x - 1) + (5 - grid$y)) / 8
gene_means <- matrix(NA_real_, nrow = 10, ncol = 25)
gene_means[1, ] <- 0.5 + 8 * u
gene_means[2, ] <- 0.5 + 8 * (1 - u)
gene_means[3:10, ] <- runif(8, min = 0.5, max = 4)
counts <- matrix(
  rpois(10 * 25, lambda = as.vector(gene_means)),
  nrow = 10,
  dimnames = list(paste0("gene", 1:10), spots)
)
counts <- methods::as(counts, "CsparseMatrix")

spatial <- CreateSeuratObject(counts = counts, assay = "Spatial")
spatial$x <- grid$x
spatial$y <- grid$y
spatial$cell_type <- sample(c("A", "B"), size = 25, replace = TRUE)

# ---- 3. Clones and receptor table ------------------------------------------
n_sequences <- sample(1:10, size = 5, replace = TRUE)
receptors <- tibble(
  clone_id = rep(as.character(1:5), times = n_sequences),
  spatial_id = sample(spots, size = sum(n_sequences), replace = TRUE),
  locus = "IGH"
)
receptors$sequence_id <- sprintf("seq%02d", seq_len(nrow(receptors)))

# ---- 3b. Real sequences for each synthetic clone ----------------------------
# Read the real table and keep the high-quality IGH sequences (the same
# thresholds as the main vignette). Then call clones with scoper.
real <- airr::read_rearrangement("inst/extdata/mg_thymus_slice2_bcr.tsv.gz")
real <- extractSpatialIds(
  real,
  id = "sequence_id",
  spatial_id = "barcode_10x",
  regex = ",([ACGT]{16})$",
  suffix = "-1"
)
real <- filterAirr(
  real,
  locus_val = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = NULL,
  spatial_id = "barcode_10x"
)
real <- filter(real, !grepl("[^ACGT]", junction))
real <- scoper::hierarchicalClones(
  real,
  threshold = 0.1,
  summarize_clones = FALSE,
  nproc = 1
)

# Use the real clones with at least 10 different sequences, smallest first.
# Synthetic clone i takes the real clone number i in this list. The function
# keeps the n most abundant different sequences of that clone.
real_sizes <- real %>%
  group_by(clone_id) %>%
  summarize(n_unique = n_distinct(sequence_alignment), .groups = "drop") %>%
  filter(n_unique >= 10) %>%
  arrange(n_unique, clone_id)
real_clone_ids <- head(real_sizes$clone_id, 5)

real_rows <- list()
for (i in 1:5) {
  rows <- real %>%
    filter(clone_id == real_clone_ids[[i]]) %>%
    arrange(desc(consensus_count), sequence_id) %>%
    distinct(sequence_alignment, .keep_all = TRUE) %>%
    head(n_sequences[[i]])
  real_rows[[i]] <- rows
}
real_rows <- bind_rows(real_rows)

# Put the real rows in the synthetic table. The synthetic sequence ID, spatial
# ID, and clone ID replace the real values.
real_rows$sequence_id <- receptors$sequence_id
real_rows$spatial_id <- receptors$spatial_id
real_rows$clone_id <- receptors$clone_id
real_rows$barcode_10x <- NULL
receptors <- real_rows %>%
  select(sequence_id, spatial_id, clone_id, locus, everything())

# One germline for each clone, from the full IMGT human references. The
# folder that comes with dowser has too few alleles. Set SIRUP_IMGT_DIR to the
# folder with the human V, D, and J FASTA files.
imgt_dir <- Sys.getenv(
  "SIRUP_IMGT_DIR",
  unset = "/Volumes/HoehnK/Hunter/germlines/imgt/human/vdj"
)
references <- dowser::readIMGT(imgt_dir)
receptors <- dowser::createGermlines(
  receptors,
  references = references,
  id = "sequence_id",
  nproc = 1
)
receptors <- tibble::as_tibble(receptors)

# ---- 4. Polygons -----------------------------------------------------------
# Each polygon covers a 2 x 2 block of spots.
polygons <- tibble(
  polygon_id = c("upper_right", "lower_left"),
  vertices = list(
    data.frame(x = c(3.5, 5.5, 5.5, 3.5), y = c(0.5, 0.5, 2.5, 2.5)),
    data.frame(x = c(0.5, 2.5, 2.5, 0.5), y = c(3.5, 3.5, 5.5, 5.5))
  )
)

# ---- 5. SIRUP object -------------------------------------------------------
sirup_synthetic <- createSIRUP(
  spatial = spatial,
  receptors = receptors,
  spatial_id = "spatial_id",
  clone = "clone_id",
  locus = "locus",
  image_scale = "hires",
  swap_coords = FALSE,
  polygons = polygons
)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveSIRUP(sirup_synthetic, "inst/extdata/sirup_synthetic.rds")
print(sirup_synthetic)
print(table(receptors$clone_id))
