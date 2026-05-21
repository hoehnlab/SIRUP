# SIRAP

Spatial Immune Receptor Analysis and Processing: how sweet is that?

`SIRAP` is an R package-in-progress for long-read spatial immune receptor
sequencing data. It is the R companion to the
[rad](https://github.com/indianewok/rad) command-line demultiplexer, and
collects reusable pieces from the existing thymus BCR/TCR analysis scripts
into package functions for:

- driving the `rad` CLI (`rad_prep`, `rad_demux`, `rad_reformat`,
  `rad_list`, `rad_modify`, `rad_scan_wl`) from R
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
remotes::install_local("/path/to/SIRAP")
```

The CLI wrappers require the `rad` binary on `PATH`, or pointed at
explicitly via `options(SIRAP.rad_path = "/path/to/rad")`, the
`SIRAP_RAD_PATH` environment variable, or the `rad =` argument on each
wrapper. See `?rad_path` for the resolution order.

## Drive the rad CLI from R

```r
library(SIRAP)

# Inspect a built-in layout
rad_prep(layout = "five_prime", read_layout = TRUE)

# Layout + position map from a real FASTQ
rad_prep(
  layout       = "five_prime",
  position_map = TRUE,
  fastq        = "reads.fq.gz",
  output       = "run/demo",
  threads      = 8
)

# Full demux pipeline
rad_demux(
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

bc_map <- read_barcode_mapping("M1_barcode_mappings.csv")

m1_bcr_heavy <- preprocess_spatial_airr(
  input = "igblast_bcr.tsv",
  barcode_map = bc_map,
  loci = "IGH",
  out_tsv = "igblast_bcr_pass_length_spatial_barcode.tsv",
  out_fasta = "igblast_bcr_pass_length_spatial_barcode.fasta"
)
```

## Generate the pRESTO/IgBLAST consensus script

```r
write_presto_igblast_script(
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
processed <- read_spatial_airr("processed_db-pass.tsv")
processed <- filter_receptor_reads(
  processed,
  loci = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = 2
)
processed <- add_spatial_barcodes(processed, bc_map, mode = "presto_header", swap = TRUE)
top_per_barcode <- collapse_top_consensus_per_barcode(processed, ties = "all")
```
