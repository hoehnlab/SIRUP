# SIRAP

Spatial Immune Receptor Analysis and Processing: how sweet is that?

`SIRAP` is an R package-in-progress for long-read spatial immune receptor
sequencing data. It collects reusable pieces from the existing thymus BCR/TCR
analysis scripts into package functions for:

- AIRR TSV reading/writing
- receptor length and productivity QC
- Visium HD 2 um to 8 um/16 um barcode mapping
- adding tissue coordinates (`barcode`, `x`, `y`)
- generating pRESTO/IgBLAST shell scripts for barcode-cluster consensus calling
- collapsing consensus calls to top sequences per barcode
- clone calling helpers
- lightweight clone/spatial visualization helpers

## Development install

```r
install.packages("remotes")
remotes::install_local("/Volumes/HoehnK/Hunter/spatial_bcr_tcr_package")
```

## Minimal preprocessing example

```r
library(SIRAP)

bc_map <- readBarcodeMapping("M1_barcode_mappings.csv")

m1_bcr_heavy <- preprocessSpatialAirr(
  input = "igblast_bcr.tsv",
  barcode_map = bc_map,
  loci = "IGH",
  out_tsv = "igblast_bcr_pass_length_spatial_barcode.tsv",
  out_fasta = "igblast_bcr_pass_length_spatial_barcode.fasta"
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
processed <- readSpatialAirr("processed_db-pass.tsv")
processed <- filterReceptorReads(
  processed,
  loci = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = 2
)
processed <- addSpatialBarcodes(processed, bc_map, mode = "presto_header", swap = TRUE)
top_per_barcode <- collapseTopConsensusPerBarcode(processed, ties = "all")
```
