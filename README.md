# SIRAP

Spatial Immune Receptor Analysis and Processing: how sweet is that?

`SIRAP` is an R package-in-progress for long-read spatial immune receptor
sequencing data. It is the R companion to the
[rad](https://github.com/indianewok/rad) command-line demultiplexer, and
collects reusable pieces from the existing thymus BCR/TCR analysis scripts
into package functions for:

- driving the `rad` CLI (`radPrep`, `radDemux`, `radReformat`,
  `radList`, `radModify`, `radScanWl`) from R
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

The CLI wrappers require the `rad` binary on `PATH`, or pointed at
explicitly via `options(SIRAP.radPath = "/path/to/rad")`, the
`SIRAP_RAD_PATH` environment variable, or the `rad =` argument on each
wrapper. See `?radPath` for the resolution order.

## Drive the rad CLI from R

For an end-to-end walkthrough — locating/building `rad`, listing layouts and
whitelists, and running `prep -> scan-wl -> demux -> reformat` on rad's public
test data (including the Visium HD spatial path) — see the
[Demultiplexing Long Reads with RAD from R](docs/readthedocs/markdown/demultiplexing-with-rad.md)
guide. To describe a custom or non-standard read structure of your own — the
layout CSV schema, orientation handling, and validating it before demux — see
[Defining Custom Read Layouts for RAD](docs/readthedocs/markdown/custom-read-layouts-with-rad.md).

```r
library(SIRAP)

# Inspect a built-in layout
radPrep(layout = "five_prime", read_layout = TRUE)

# Layout + position map from a real FASTQ
radPrep(
  layout       = "five_prime",
  position_map = TRUE,
  fastq        = "reads.fq.gz",
  output       = "run/demo",
  threads      = 8
)

# Full demux pipeline
radDemux(
  layout  = "five_prime",
  fastq   = "reads.fq.gz",
  output  = "demo",
  out_dir = "run",
  threads = 8
)
```

## Minimal preprocessing example

```r
library(SIRAP)

m1_bcr_heavy <- data.table::fread("igblast_bcr.tsv")
m1_bcr_heavy <- extractSpatialBarcode(
  m1_bcr_heavy,
  sequence_id_col = "sequence_id",
  barcode_col = "barcode_2um",
  regex = "[0-9]+_[0-9]+$",
  prefix = "s_002um_",
  suffix = "-1"
)
m1_bcr_heavy <- processSpatialAirr(
  m1_bcr_heavy,
  loci = "IGH",
  barcode_col = "barcode_2um",
  add_full_id = TRUE
)
```

## Run pRESTO/IgBLAST consensus calling

```r
processSpatialSample(
  input_fasta = "m1/igblast_bcr_pass_length_spatial_barcode.fasta",
  output_dir = "m1/processed",
  receptor = "ig",
  germline_glob = "/path/to/germlines/imgt/human/vdj/imgt_human_*.fasta"
)
```

## Post-consensus filtering sketch

```r
processed <- airr::read_rearrangement("processed_db-pass.tsv")
processed <- extractSpatialBarcode(
  processed,
  sequence_id_col = "sequence_id",
  barcode_col = "barcode_2um",
  regex = ",([^,]+)"
)
processed <- processSpatialAirr(
  processed,
  loci = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = 2,
  barcode_col = "barcode_2um"
)
bc_map <- readBarcodeMapping("barcode_mappings.parquet")
processed <- processed |>
  dplyr::mutate(barcode_2um_swapped = sub("^(.*_)(\\d+)_(\\d+)(-\\d+)$", "\\1\\3_\\2\\4", barcode_2um)) |>
  dplyr::left_join(
    bc_map |> dplyr::select(square_002um, square_008um, square_016um),
    by = c("barcode_2um_swapped" = "square_002um")
  ) |>
  dplyr::rename(barcode_8um = square_008um, barcode_16um = square_016um)
top_per_barcode <- collapseTopConsensusPerBarcode(processed, ties = "all")
```

## Pair spatial and receptor data

`createSIRAP()` records the relationship between a spatial Seurat object and
an AIRR-style receptor table. The receptor component is stored as a tibble, and
unmatched bulk or other nonspatial receptor rows are retained.

```r
sirap <- createSIRAP(
  spatial = seurat_object,
  receptors = cloned_receptors,
  barcode_col = "barcode_8um",
  clone_col = "clone_id"
)

sirap
receptorTable(sirap)

plotReceptors(
  object = sirap,
  group_col = "clone_id"
)
```

Functions continue to accept separate Seurat and receptor objects, so existing
workflows remain supported.
