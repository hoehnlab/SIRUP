# SIRUP workflows

Each recipe lists what it needs, gives code, and says what to check. The code runs on the example data that comes with SIRUP. For the user's data, change the file names and check the column names.

Recipes:

- A. Build an object and call clones
- B. Regions: clonal dominance and diversity
- C. Lineage trees on the tissue
- D. Spatial statistics
- E. Several samples, bulk data, and blood
- F. Process reads into an AIRR table

All recipes start with this setup:

```r
library(SIRUP)
library(dplyr)
extdata <- system.file("extdata", package = "SIRUP")
```

The example data are a myasthenia gravis (MG) thymus, in two 10x Visium sections:

| File in `extdata` | Holds |
|---|---|
| `mg_thymus_slice2_spatial_brighter.rds` | A Seurat object for section 2 (a few hundred genes, a low-resolution image) |
| `mg_thymus_slice2_bcr.tsv.gz` | IGH consensus sequences of section 2, in an AIRR table |
| `mg_thymus_slice2_tcr.tsv.gz` | TRB consensus sequences of section 2 |
| `mg_thymus_slice2_gc_polygons.csv` | Outlines of the ectopic germinal centers (eGCs) of section 2 |
| `mg_thymus_cell_types.csv.gz` | Cell-type proportions for each spot |
| `mg_thymus_slice1_*` | The same files for section 1 |
| `mg_thymus_bulk_*.tsv.gz` | Bulk BCR sequences of the thymus and of blood |

## A. Build an object and call clones

Needs: airr, Seurat, scoper.

```r
thymus <- readRDS(file.path(extdata, "mg_thymus_slice2_spatial_brighter.rds"))
bcrs <- airr::read_rearrangement(file.path(extdata, "mg_thymus_slice2_bcr.tsv.gz"))

# The pipeline writes the spot barcode at the end of sequence_id, after a comma.
bcrs <- extractSpatialIds(
  bcrs, id = "sequence_id", spatial_id = "barcode_10x",
  regex = ",([ACGT]{16})$", suffix = "-1"
)
mean(bcrs$barcode_10x %in% colnames(thymus))     # the share of sequences that match a spot

bcrs <- filterAirr(
  bcrs, locus_val = "IGH", min_v_length = 250, min_j_length = 30,
  min_junction_length = 10, productive_only = TRUE,
  min_consensus_count = NULL, spatial_id = "barcode_10x"
)

sirup <- createSIRUP(
  spatial = thymus, receptors = bcrs, spatial_id = "barcode_10x",
  locus = "locus", image_scale = "lowres", swap_coords = FALSE,
  metadata = list(tissue = "MG thymus", section = "slice 2")
)
sirup <- updateReceptors(
  sirup, scoper::hierarchicalClones,
  threshold = 0.1, summarize_clones = FALSE, nproc = 1
)
summary(sirup)
```

Check: the match share is high (a few sequences fall outside the tissue). `summary()` shows the receptor count. The receptor table has a `clone_id` column. scoper prints a warning when it removes sequences with an `N` in the junction. Report the count that it removed.

Plots:

```r
# Gene expression over the H&E image. The alra assay gives smoother maps.
spatialExpressionPlot(sirup, features = "AICDA", assay = "alra", layer = "data", zero_tol = 0)

# Spots that hold a receptor, colored by a receptor column.
plotReceptors(sirup, group = "locus", colors = "#0072B2", receptor_size = 1.2)

# The three largest clones.
top_clones <- getReceptors(sirup) %>% count(clone_id, sort = TRUE) %>%
  slice_head(n = 3) %>% pull(clone_id)
plotReceptors(sirup, group = "clone_id", group_val = top_clones, receptor_size = 1.6)

# Expression and receptors in one plot.
plotReceptorsExpression(sirup, group = "locus", gene = "AICDA", assay = "alra",
                        layer = "data", colors = "#0072B2")
```

A Visium spot holds several cells, so two selected clones can share a spot. `plotReceptors()` then averages the colors. The `overlap` argument changes this.

## B. Regions: clonal dominance and diversity

Needs: a polygon CSV file from the user (see `references/object.md`, "Regions").

```r
sirup <- addPolygons(
  sirup, polygons_path = file.path(extdata, "mg_thymus_slice2_gc_polygons.csv")
)
getPolygons(sirup)
plotSpatialPolygons(sirup, polygon_alpha = 0.3)

# Receptors, spots, and clones in each polygon. The result is stored in the object.
sirup <- polygonReceptorSummary(sirup, polygon_id = "polygon_id")
getAnalysis(sirup, "polygon_receptor_summary")

# Clonal dominance: the fraction of each polygon's receptors that each clone holds.
dominance <- spatialCloneDominance(object = sirup, locus_val = "IGH")
dominance$dominance %>% arrange(polygon, desc(proportion)) %>% head()
pies <- plotSpatialPolygonPies(dominance = dominance, object = sirup,
                               fill = "clone_id", top_n = 5)
pies$page_1
sirup <- setAnalysis(sirup, "gc_dominance", dominance)
```

Other region functions: `expandPolygons()` makes a zone around a polygon. `withPolygons()` runs a function on the receptors of each polygon. `polygonDiversity()` and `polygonAbundance()` estimate clonal diversity and abundance. `polygonPairwiseDistance()` measures the distances between polygon centers. Read their usage in `references/functions-spatial-polygons.md`. Both functions resample all polygons to the same number of receptors, the size of the smallest polygon that has at least `min_n` receptors (default 10), and leave out smaller polygons. Set `min_n` before you compare regions.

Check: `getPolygons(sirup)` has one row for each polygon. The receptor table has a `polygon_id` column, with `NA` for receptors outside every polygon.

## C. Lineage trees on the tissue

Needs: dowser, shazam, scoper, and the IMGT reference folder (the user downloads it and sets `SIRUP_IMGT_DIR`). Tree building is slow for large clones. Ask the user first.

Tree tip names cannot contain a comma, so replace the comma in `sequence_id` before you create the object:

```r
bcrs <- mutate(bcrs, sequence_id = gsub(",", "_", sequence_id))
# ... createSIRUP() and the clone call from recipe A, then addPolygons() from recipe B.

imgt_dir <- path.expand(Sys.getenv("SIRUP_IMGT_DIR"))
stopifnot(dir.exists(imgt_dir))
references <- dowser::readIMGT(imgt_dir)
sirup <- updateReceptors(
  sirup, dowser::createGermlines,
  references = references, id = "sequence_id", nproc = 1
)

# The SIRUP paper also removed likely PCR chimeras after germline reconstruction.
sirup <- updateReceptors(sirup, function(receptors) {
  chimera <- shazam::slideWindowDb(
    receptors, sequenceColumn = "sequence_alignment",
    germlineColumn = "germline_alignment_d_mask", mutThresh = 6, windowSize = 10
  )
  receptors[!chimera, ]
})

# Choose expanded clones. This keeps clones with 10 to 60 sequences and one sequence
# in a polygon at least.
clone_sizes <- getReceptors(sirup) %>%
  group_by(clone_id) %>%
  summarise(n_sequences = n(), n_in_polygon = sum(!is.na(polygon_id))) %>%
  filter(n_sequences >= 10, n_sequences <= 60, n_in_polygon > 0) %>%
  arrange(desc(n_in_polygon))
tree_clones <- head(clone_sizes$clone_id, 3)

# Build the trees on a small object that holds only these clones.
tree_sirup <- updateReceptors(sirup, dplyr::filter, clone_id %in% tree_clones)
tree_sirup <- formatSpatialClones(
  tree_sirup, analysis_name = "formatted_clones", id = "sequence_id",
  germ = "germline_alignment_d_mask", traits = c("barcode_10x", "polygon_id")
)
set.seed(1)    # tree building is random; the seed makes the trees repeatable
tree_sirup <- getSpatialTrees(
  tree_sirup, f_clones_analysis = "formatted_clones",
  analysis_name = "trees", build = "pratchet"
)
trees <- getAnalysis(tree_sirup, "trees")

# Show each tree next to its spots on the tissue.
clustered <- clusterTreesTips(trees, k = 4, id = "sequence_id", cluster = "tree_cluster")
pairs <- plotSpatialCloneTreePairs(
  clustered, object = sirup, clone_val = tree_clones, id = "sequence_id",
  color_by = "tree_cluster", zoom = TRUE, polygons = TRUE,
  polygon_alpha = 0, polygon_color = "black", show_spatial_legend = TRUE
)
cowplot::plot_grid(
  cowplot::plot_grid(pairs[[1]]$tree, pairs[[1]]$spatial_clone, rel_widths = c(1.2, 1)),
  pairs[[1]]$legend, ncol = 1, rel_heights = c(1, pairs[[1]]$legend_rel_height)
)
```

`plotSpatialCloneTreePairs()` builds the legend only when `show_spatial_legend = TRUE`. `tipPatristicMatrix()` gives the distances between tips in a tree. `spatialCloneTreePairsToPDF()` writes one PDF page for each clone.

Check: `trees` holds one tree for each chosen clone. The plot shows each tip as a circle on the map when the tip has a spot in this section.

## D. Spatial statistics

Needs: spdep (installed with SIRUP). Ripley's cross K needs the spatstat packages. Choose the statistic with `references/choose-function.md`.

```r
# Is a gene clustered? Moran's I. k is the number of nearest spots that count as neighbors.
aicda <- spatialMoran(sirup, variable = "AICDA", variable_type = "feature",
                      assay = "Spatial", layer = "data", k = 6)
aicda$test            # I is in $test$estimate; the p value is in $test$p.value

# Where are the clusters? Local Moran's I. One row for each spot.
local <- spatialLocalMoran(sirup, variable = "AICDA", variable_type = "feature",
                           assay = "Spatial", layer = "data", k = 6)
plotSpatialLocalMoran(sirup, local, point_size = 1)

# Do two genes co-vary? Lee's L.
lee <- spatialLee(sirup, x = "AICDA", y = "MS4A1", assay = "Spatial", layer = "data", k = 6)
lee$L; lee$test$p.value

# Does one clone stay together? Join counts, for chosen clones.
jc <- spatialCloneJoinCount(sirup, clone_val = top_clones, use_stored_polygons = FALSE,
                            k = 6, return_data = TRUE)
jc$summary

# How far does a clone spread? A description, not a test.
spatialCloneSpread(sirup, min_spatial_ids = 3)
```

`variable_type` is `"feature"` (a gene), `"metadata"` (a numeric Seurat metadata column), `"receptor"` (receptor present or not), or `"clone"` (one clone). For a text metadata column, give `value` to test one level.

For Ripley's cross K, first mark the spots in the Seurat metadata with `updateSpatial()`. Then call `spatialRipleyCrossK(sirup, mark_col = "mark", i = "A", j = "B", nsim = 39, max_r = 40)`. With `nsim = 0` the function returns only the observed curve. Distances are in image pixels at the stored `image_scale`.

Read the result as follows. Moran's I measures the size of the pattern, and the p value measures the evidence. A large section gives a small p value for a weak pattern. Look at I first. Local Moran's I gives one p value for each spot, so adjust for multiple testing, for example with `p.adjust(method = "BH")`.

Store results with `setAnalysis()`.

## E. Several samples, bulk data, and blood

Needs: a named list of Seurat objects, and one receptor table with a `sample` column.

```r
# spatial: a named list, for example list(mg_thymus_slice1 = seurat1, mg_thymus_slice2 = seurat2)
# receptors: one table. Section rows have a barcode column that matches a spot.
# Bulk and blood rows have a sample name that is not in the list, and no spot.
sirup <- createSIRUP(
  spatial = spatial_list, receptors = receptor_table,
  spatial_id = "barcode_10x", sample_id = "sample",
  image_scale = "lowres", swap_coords = FALSE
)
spatialSamples(sirup)
```

Call clones on all samples at once, so one clone ID names the same lineage in every sample. Make each `sequence_id` unique across samples, for example by pasting the sample name in front. Use `getReceptors(sirup)` to count, for each clone, the sequences in each compartment (section, bulk, blood).

Alternatively, build one object for each section and join them with `mergeSIRUP(a, b, samples = c("a", "b"))`. `subsetSIRUP(x, sample_id = "a")` keeps one sample.

## F. Process reads into an AIRR table

Needs: demultiplexed reads, and these tools: IgBLAST, pRESTO, vsearch, MUSCLE. The R packages rrad (demultiplexing) and airrblast (runs IgBLAST) are also needed. Do not install the tools without the user's approval. Run the pipeline on a compute node.

```r
# 1. Demultiplex (rrad package). Write a FASTA file with the barcode and UMI in the read name.
demux <- rrad::rad_demux(layout = "five_prime", fastq = "reads.fq.gz",
                         out_dir = "prepared", output = "sample")
rrad::rad_reformat(demux, output_fastq = "prepared/sample-igblast.fasta",
                   to_fasta = TRUE, header_format = "presto")

# 2. Build consensus sequences. Use receptor = "tr" for TCRs, in a separate call.
result <- processSpatialSample(
  input = "prepared/sample-igblast.fasta", output_dir = "processed/bcr",
  receptor = "ig", organism = "human",
  igblast_base = "path/to/igblast", igblast_refs = "path/to/igblast_refs",
  imgt_dir = "path/to/imgt", muscle_exec = "path/to/muscle3.8.31",
  identity = 1.0, consensus_freq = 0.51, consensus_min_count = 2, nproc = 8
)
result$db_pass
```

Parameters from the SIRUP paper: `identity = 1.0` for BCRs, so that SIRUP does not merge sequences that differ by one mutation. Use `identity = 0.99` for TCRs. `consensus_freq = 0.51` gave the most high-quality sequences for both chains.

After processing, extract the spot barcode from the sequence ID, add quality measures, and filter:

```r
bcrs <- extractSpatialIds(consensus, id = "sequence_id", spatial_id = "barcode_10x",
                          regex = ",([ACGT]{16})$", suffix = "-1")
bcrs <- addSequenceQualityMetrics(bcrs, sequence = "sequence", alignment = "sequence_alignment")
bcrs <- filterAirr(bcrs, locus_val = "IGH", min_v_length = 250, min_j_length = 30,
                   min_junction_length = 10, productive_only = TRUE,
                   min_consensus_count = NULL, max_gap_length = 15,
                   spatial_id = "barcode_10x")
```

Use `max_gap_length = 15` for IGH and 24 for TRB. Then continue with recipe A.

For Visium HD, reads carry a 2 µm barcode. `readSpatialIdMapping()` reads the Space Ranger file that maps 2 µm barcodes to 8 µm and 16 µm bins. Join the mapping to the receptor table, load the 8 µm bins with `Seurat::Load10X_Spatial(bin.size = 8)`, and give `spatial_id = "barcode_8um"` to `createSIRUP()`. The processing page on the documentation site explains each step.
