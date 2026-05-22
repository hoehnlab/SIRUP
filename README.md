# SIRAP

Spatial Immune Receptor Analysis and Processing: how sweet is that?

`SIRAP` is an R package-in-progress for long-read spatial immune receptor
sequencing data. It collects reusable pieces from the existing thymus BCR/TCR
analysis scripts into package functions for:

- receptor length and productivity QC
- Visium HD 2 um to 8 um/16 um barcode mapping
- adding tissue coordinates (`barcode`, `x`, `y`)
- generating pRESTO/IgBLAST shell scripts for barcode-cluster consensus calling
- collapsing consensus calls to top sequences per barcode
- lightweight clone/spatial visualization helpers

## Development install

```r
install.packages("remotes")
remotes::install_local("/Volumes/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP")
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

## Generate the pRESTO/IgBLAST consensus script

```r
writePrestoIgblastScript(
  input_fastas = c("m1/igblast_bcr_pass_length_spatial_barcode.fasta",
                   "m5/igblast_bcr_pass_length_spatial_barcode.fasta"),
  output_db_paths = c("m1/processed/final_db-pass.tsv",
                      "m5/processed/final_db-pass.tsv"),
  script_path = "inst/scripts/process_bcr_heavy.sh",
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
