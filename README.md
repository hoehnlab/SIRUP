# SIRAP

Spatial Immune Receptor Analysis and Processing: how sweet is that?

`SIRAP` is an R package-in-progress for long-read spatial immune receptor
sequencing data. It includes the [RAD](https://github.com/indianewok/rad)
demultiplexing engine through its required `rrad` dependency, and
collects reusable pieces from the existing thymus BCR/TCR analysis scripts
into package functions for:

- zero-setup native RAD whitelist discovery, demultiplexing, and reformatting
  with `rad_scan_wl()`, `rad_demux()`, and `rad_reformat()`
- optional compatibility access to the external `rad` CLI (`radPrep`,
  `radDemux`, `radReformat`, `radList`, `radModify`, `radScanWl`)
- AIRR TSV reading/writing
- receptor length and productivity QC
- Visium HD 2 um to 8 um/16 um barcode mapping
- adding tissue coordinates (`barcode`, `x`, `y`)
- a resumable two-pass IgBLAST workflow: read-level V/D/J assignment followed
  by spatial-barcode/sequence-cluster consensus building and consensus-level
  V/D/J assignment
- collapsing consensus calls to top sequences per barcode
- validated SIRAP objects pairing a spatial Seurat object with a tibble-backed
  AIRR-style receptor table
- lightweight clone/spatial visualization helpers
- spatial polygon selection and clonal dominance summaries

## Private installation

Both SIRAP and its required `rrad` dependency are currently private. Before
installing, your GitHub account must have read access to both repositories.
Store the credential in the system credential manager; never put a token in an
R script, repository URL, or committed `.Renviron` file.

```r
install.packages(c("gitcreds", "remotes"))
gitcreds::gitcreds_set() # one-time interactive credential setup

# Installs SIRAP and its pinned rrad v1.2.1 dependency.
remotes::install_github("hjm2217/SIRAP", dependencies = NA)
```

If you already have a SIRAP checkout, replace the last call with
`remotes::install_local("/path/to/SIRAP", dependencies = NA)`. Access failures
for a private repository commonly appear as repository-not-found or 404
errors. The credential used by a direct GitHub installation must be able to
read both repositories.

The embedded engine currently supports macOS and Linux. Linux builds use
OpenMP when R exposes it; current macOS builds run the demultiplexer serially
and record one effective thread even if a larger value was requested.

The native `rad_scan_wl()`, `rad_demux()`, and `rad_reformat()` paths use the
RAD core compiled into the required `rrad` dependency; users do not install or
configure a separate executable. Large built-in 10x whitelists are not fetched
during package installation. The first native run that needs one downloads it
to R's per-user cache, verifies its pinned size and SHA-256 checksum, and reuses
it thereafter. `rad_download_whitelists()` is available only when a machine or
cluster cache must be prepared before going offline.

The older camel-case CLI wrappers remain available as a compatibility path.
They require the `rad` binary on `PATH`, or pointed at explicitly via
`options(SIRAP.rad_path = "/path/to/rad")`, the `SIRAP_RAD_PATH` environment
variable, or the `rad =` argument. See `?radPath` for the resolution order.

## Run RAD natively from R

`rad_demux()` performs a compatibility handshake before every run and returns
a structured result containing the upstream commit, patched-engine and
resource digests, normalized configuration, artifacts, statistics, and SIRAP
provenance. It recalibrates layout caches by default; opt into reuse only with
`reuse_cache = TRUE` when the existing files belong to the same input.
Choose either `kit` or `global_whitelist` as the reference barcode source, not
both. Automatic whitelist discovery and `custom_whitelist` are likewise
alternatives for the true/sample barcode set.

```r
library(SIRAP)

rad_core_info()

run <- rad_demux(
  layout  = "five_prime",
  fastq   = "reads.fq.gz",
  output  = "demo",
  out_dir = "run",
  threads = 8
)

run$artifacts
run$stats
run$provenance$core_handshake
```

Whitelist discovery and post-demultiplexing reformatting use the same
in-process engine and compatibility handshake:

```r
scan <- rad_scan_wl(
  input = "reads.fq.gz",
  output_prefix = "run/cells",
  adapter_seq = "ACGTACGT",
  barcode_length = 16L,
  whitelist = "10x_3v3"
)

formatted <- rad_reformat(
  fastq = run,
  output_fastq = "run/demo-igblast.fasta",
  to_fasta = TRUE,
  header_format = "presto"
)
```

That conversion is streamed and produces canonical identifiers such as
`read-F|BARCODE=AACC|UMI=TGCA`; FASTQ qualities are not copied because IgBLAST
accepts FASTA, not FASTQ.

Ordinary runs need no setup call. To populate the verified cache explicitly
before an offline batch, use `rad_download_whitelists()`; a shared cache can be
selected with `options(rrad.cache_dir = "/shared/path")`.

To deliberately use an externally installed CLI instead, select it explicitly:

```r
run <- rad_demux(
  layout = "five_prime",
  fastq = "reads.fq.gz",
  backend = "cli",
  rad = "/path/to/rad"
)
```

The compatibility adapter accepts RAD CLI 1.0.2 only. Because that legacy
binary does not reliably propagate every parser, gzip-writer, or cache error,
SIRAP validates the complete input and output FASTQ streams, reports its cache
policy as uncontrolled, and rejects CLI debug mode. These extra safety passes
apply only when `backend = "cli"`; the native path remains the recommended one.

## Drive the legacy rad CLI from R

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

## RAD to spatial AIRR: the two IgBLAST passes

```r
library(SIRAP)

processed <- processSpatialSample(
  input_fasta = run, # the rad_demux() result can be passed directly
  output_dir = "m1/processed",
  receptor = "ig",
  igblast_base = "/path/to/igblast",
  germline_glob = "/path/to/germlines/imgt/human/vdj/imgt_human_*.fasta"
)

processed$artifacts$read_igblast_airr # pass 1: V/D/J for each RAD read
processed$artifacts$spatial_fasta     # filtered, pRESTO-ready reads
processed$artifacts$consensus_fasta   # spatial barcode/sequence-cluster consensus
processed$artifacts$final_airr        # pass 2: consensus V/D/J calls
```

The first pass is not redundant: it supplies each read's V/D/J calls and
aligned receptor sequence before spatial grouping. The second pass annotates
the consensus sequence produced by pRESTO. Both use the same `igblast_base`
installation and conventional germline databases. An existing read-level AIRR table
can resume the pipeline with `input_type = "airr"`; an already prepared
spatial FASTA uses `input_type = "spatial_fasta"` and skips pass 1. A bare
FASTA path otherwise defaults to the full two-pass route, so pass 1 is not
silently skipped; gzip-compressed RAD FASTA is normalized automatically.

The read-level pass streams the FASTA into complete, approximately
base-balanced chunks and runs up to four IgBLAST processes concurrently.
`nproc` is the total CPU budget, so `nproc = 12` uses four processes with
three IgBLAST threads each. Small inputs automatically retain the original
single-process path, and chunk results are merged back in input order. Direct
`runIgblast()` callers can set `max_processes = 1` to force monolithic
execution or supply `chunk_bases` when benchmarking alternative chunk sizes;
inside `processSpatialSample()`, use `igblast_max_processes` to tune or disable
the worker pool.

## Post-consensus filtering sketch

```r
processed <- airr::read_rearrangement("processed_db-pass.tsv")
processed$barcode_2um <- processed$barcode
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
