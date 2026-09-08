# SIRAP

Spatial Immune Receptor Analysis and Processing: how sweet is that?

`SIRAP` is an R package-in-progress for long-read spatial immune receptor
sequencing data. Its processing pipeline starts from prepared, demultiplexed
FASTA reads with barcode annotations and runs a first IgBLAST pass, read
filtering, clustering within each barcode, consensus generation, and a second
IgBLAST pass on the consensus sequences. It collects reusable pieces from the
existing thymus BCR/TCR analysis scripts into package functions for:

- processing prepared FASTA reads through IgBLAST and pRESTO
- AIRR TSV reading/writing
- receptor length and productivity QC
- Visium HD 2 um to 8 um/16 um barcode mapping
- adding tissue coordinates (`barcode`, `x`, `y`)
- generating pRESTO/IgBLAST shell scripts for barcode-cluster consensus calling
- collapsing consensus calls to top sequences per barcode
- validated SIRAP objects pairing a spatial Seurat object with a tibble-backed
  AIRR-style receptor table
- lightweight clone/spatial visualization helpers
- spatial polygon selection and clonal dominance summaries

## Development install

```r
install.packages("remotes")
remotes::install_local("/path/to/SIRAP")
```

Processing requires IgBLAST, its reference databases, the pRESTO/Change-O tools,
vsearch, and muscle. Configure their paths for your installation; see
`?processSpatialSample` and `?runIgblast`.

## Prepare FASTA, then process it in SIRAP

Demultiplexing and read-layout/whitelist setup happen upstream of SIRAP.
The standalone [rrad package](https://github.com/indianewok/rrad) is one way
to prepare those reads. It is listed in `Suggests`, with its GitHub source
pinned to `v1.2.1`. Install it separately when using the preparation example
below; SIRAP does not require rrad to process an existing FASTA.

```r
# Separate upstream preparation with rrad:
demux <- rrad::rad_demux(
  layout = "five_prime", # choose the layout matching your library
  fastq = "reads.fq.gz",
  out_dir = "prepared",
  output = "sample"
)
rrad::rad_reformat(
  demux,
  output_fastq = "prepared/sample-igblast.fasta",
  to_fasta = TRUE,
  header_format = "presto"
)

# SIRAP begins at the prepared FASTA and runs both IgBLAST passes:
result <- SIRAP::processSpatialSample(
  input = "prepared/sample-igblast.fasta",
  output_dir = "processed",
  receptor = "ig",
  igblast_db_dir = "/path/to/igblast/database",
  igblast_base = "/path/to/igblast",
  germline_glob = "/path/to/germlines/imgt/human/vdj/imgt_human_*.fasta"
)
```

The preparation writes headers such as `read1|BARCODE=BC01|UMI=ACGT`, keeping
barcode information inside the identifier retained by IgBLAST. Ordinary FASTA
input runs the first annotation pass automatically. Use
`input_type = "spatial_fasta"` only for a FASTA already annotated, filtered,
and formatted for within-barcode consensus; that entry point skips the first
IgBLAST pass.

See [Preparing FASTA with rrad for SIRAP](docs/readthedocs/markdown/demultiplexing-with-rad.md)
for the input contract, whitelist downloads, and Visium HD preparation, and
[Custom read layouts and whitelists](docs/readthedocs/markdown/custom-read-layouts-with-rad.md)
for custom upstream preparation.

## Minimal preprocessing example

```r
library(SIRAP)

m1_bcr_heavy <- data.table::fread("igblast_bcr.tsv")
m1_bcr_heavy <- extractSpatialIds(
  m1_bcr_heavy,
  id = "sequence_id",
  spatial_id = "barcode_2um",
  regex = "[0-9]+_[0-9]+$",
  prefix = "s_002um_",
  suffix = "-1"
)
m1_bcr_heavy <- filterAirr(
  m1_bcr_heavy,
  locus_val = "IGH",
  spatial_id = "barcode_2um",
  add_spatial_sequence_id = TRUE
)
```

## Run pRESTO/IgBLAST consensus calling

```r
processSpatialSample(
  input = "m1/igblast_bcr_pass_length_spatial_id.fasta",
  input_type = "spatial_fasta", # already annotated and filtered
  output_dir = "m1/processed",
  receptor = "ig",
  germline_glob = "/path/to/germlines/imgt/human/vdj/imgt_human_*.fasta"
)
```

## Post-consensus filtering sketch

```r
processed <- airr::read_rearrangement("processed_db-pass.tsv")
processed <- extractSpatialIds(
  processed,
  id = "sequence_id",
  spatial_id = "barcode_2um",
  regex = ",([^,]+)"
)
processed <- filterAirr(
  processed,
  locus_val = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = 2,
  spatial_id = "barcode_2um"
)
bc_map <- readSpatialIdMapping("barcode_mappings.parquet")
processed <- processed |>
  dplyr::mutate(barcode_2um_swapped = sub("^(.*_)(\\d+)_(\\d+)(-\\d+)$", "\\1\\3_\\2\\4", barcode_2um)) |>
  dplyr::left_join(
    bc_map |> dplyr::select(square_002um, square_008um, square_016um),
    by = c("barcode_2um_swapped" = "square_002um")
  ) |>
  dplyr::rename(barcode_8um = square_008um, barcode_16um = square_016um)
top_per_barcode <- collapseTopConsensusPerSpatialId(processed, ties = "all")
```

## Pair spatial and receptor data

`createSIRAP()` records the relationship between a spatial Seurat object and
an AIRR-style receptor table. The receptor component is stored as a tibble, and
unmatched bulk or other nonspatial receptor rows are retained.

```r
sirap <- createSIRAP(
  spatial = seurat_object,
  receptors = cloned_receptors,
  spatial_id = "barcode_8um",
  clone = "clone_id",
  image_scale = "hires",
  swap_coords = FALSE
)

sirap
getReceptors(sirap)
listAnalyses(sirap)

# Work with the receptor table directly when an external function needs it.
clone_counts <- dplyr::count(
  getReceptors(sirap),
  clone_id,
  sort = TRUE
)

# Store a derived result with the data that produced it.
sirap <- setAnalysis(sirap, "clone_counts", clone_counts)
getAnalysis(sirap, "clone_counts")

plotReceptors(
  object = sirap,
  group = "clone_id"
)

spatialExpressionPlot(
  object = sirap,
  features = "MS4A1"
)

spatialReceptorSummary(sirap)
```

Read-only plotting, spatial-analysis, and receptor-summary functions accept a
SIRAP object directly. Functions continue to accept separate Seurat and
receptor objects, so existing workflows remain supported.

The component helpers work with functions from SIRAP, other packages, or your
own code:

- `updateReceptors()`, `updateSpatial()`, and `updateAnalysis()` put the result
  back into a returned SIRAP object.
- `getReceptors()`, `getSpatial()`, and `getAnalysis()` return stored data
  when a function should inspect it without modifying the SIRAP object.
- `setAnalysis()`, `getAnalysis()`, `listAnalyses()`, and `removeAnalysis()`
  manage arbitrary derived results, including Dowser clones or trees and
  spatial statistics.

This means clone calling or germline creation can be run directly against the
stored receptor table, for example
`updateReceptors(sirap, scoper::hierarchicalClones, ...)`. A conventional
`clone_id` column added by the function is recorded automatically.

## Save and load a SIRAP object

`saveSIRAP()` writes a validated object to a single RDS file and `readSIRAP()`
reads it back. Unlike a bare `saveRDS()`/`readRDS()`, this pair validates on
write, records the SIRAP and Seurat versions used, upgrades older objects to the
current schema, warns on version drift, and re-validates on read.

```r
saveSIRAP(sirap, "m1.sirap.rds")
sirap <- readSIRAP("m1.sirap.rds")
```
