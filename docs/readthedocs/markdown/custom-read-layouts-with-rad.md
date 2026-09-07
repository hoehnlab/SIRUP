Defining Custom Read Layouts for RAD
================

The whole point of [`rad`](https://github.com/indianewok/rad) is that it
is *read-structure agnostic*: instead of hard-coding one chemistry, you
hand it a small CSV that describes how your reads are laid out —
primers, barcodes, UMIs, the cDNA — and it figures out the alignment and
extraction from there. The [demultiplexing
guide](demultiplexing-with-rad.html) drove rad on its **bundled**
layouts (`five_prime`, `three_prime`, `sctagger`, `visium_hd`, …). This
guide is the other half: how to write your *own* layout for a wonky or
custom format and run the same pipeline on it through SIRAP.

SIRAP exposes two deliberately separate RAD interfaces. The recommended
snake-case functions — `rad_demux()`, `rad_scan_wl()`, and
`rad_reformat()` — delegate to the RAD core compiled into the required
`rrad` package. They need no separately downloaded `rad` executable. The
preserved camel-case functions — `radDemux()`, `radScanWl()`,
`radReformat()`, and the CLI-only helpers `radList()`, `radPrep()`,
`radModify()`, and `radRun()` — launch an external RAD CLI for
compatibility with older workflows.

This guide uses the native snake-case interface for data processing.
Sections that demonstrate a CLI-only layout inspection or registry
command are labeled accordingly. The code is not evaluated while the
guide is rendered because the examples need a FASTQ matching the custom
structure and write processing artifacts; copy the chunks into a live R
session against your own reads.

## When You Need a Custom Layout

Reach for a custom layout when none of the bundled keys describe your
reads:

  - a non-10x chemistry, or a 10x variant with an extra/odd element;
  - targeted or amplicon libraries with custom primers flanking the
    insert;
  - multi-round combinatorial barcodes (SPLiT-seq-like) beyond the
    bundled `splitseq` template;
  - anything you would describe to a colleague as “it’s mostly 5’,
    but…”.

If your reads are *close* to a bundled layout, don’t start from a blank
file — copy the closest bundled CSV and edit it (see [Start From a
Bundled Template](#start-from-a-bundled-template)).

## Setup

``` r
library(SIRAP)

# Everything in this guide writes under one scratch directory.
runDir <- file.path(tempdir(), "rad_layout")
dir.create(runDir, recursive = TRUE, showWarnings = FALSE)

nproc <- 4

# Confirm that SIRAP and rrad agree on the embedded RAD core and resources.
rad_core_info()
```

No external binary setup is needed for the native examples. Only the
sections explicitly labeled **legacy CLI** require `radPath()` to find
an executable; see the [demultiplexing
guide](demultiplexing-with-rad.html#point-r-at-the-rad-binary-legacy-cli-only)
for its resolution rules and `updateRad()` for maintaining a source
checkout.

## The Anatomy of a Read Layout

rad models a read as an **ordered list of elements**. Two kinds matter:

  - **static** elements are fixed sequences rad *aligns* — adapters,
    primers, the TSO, a poly-T stretch. They are the anchors.
  - **variable** elements are the slots rad *extracts* — the cell
    barcode, the UMI, and the cDNA `read` itself. rad reads them out of
    the windows between the static anchors.

At `demux` time rad aligns the static anchors, then pulls each variable
element from its position relative to those anchors, and (for barcode
elements) corrects it against a whitelist. So a good layout is really a
good description of “what fixed sequence sits where, and what variable
thing comes between.”

## The Layout CSV Schema

A layout is a CSV with one row per element. The columns rad expects:

| Column            | Meaning                                                                                                                                                                |
| ----------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`              | A unique name for the element (`forw_primer`, `barcode`, `umi`, `read`, …).                                                                                            |
| `seq`             | The literal sequence — required for **static** elements, left blank for variable ones.                                                                                 |
| `expected_length` | Length in bases — required for **variable** elements (barcode 16, UMI 10/12), blank for static (length comes from `seq`).                                              |
| `type`            | `static` or `variable`.                                                                                                                                                |
| `class`           | What the element *is*: `barcode`, `umi`, `read`, and adapter/primer-style anchors. The `barcode`/`umi`/`read` classes drive extraction.                                |
| `whitelist`       | For a `barcode` element, the kit key (or path) whose whitelist corrects it — e.g. `10x_5v1`. Blank for non-barcode rows.                                               |
| `direction`       | *Optional.* `forward_only` / `reverse_only`; omit it to get forward + an auto-generated reverse complement (see [Orientation](#orientation-and-the-direction-column)). |

A minimal 5’-style layout, written as a data frame and saved to disk:

``` r
layout <- data.frame(
  id              = c("forw_primer",            "barcode",  "umi",      "read", "rev_primer"),
  seq             = c("CTACACGACGCTCTTCCGATCT", "",         "",         "",     "GTACTCTGCGTTGATACCACTGCTT"),
  expected_length = c("",                       "16",       "10",       "",     ""),
  type            = c("static",                 "variable", "variable", "variable", "static"),
  class           = c("",                       "barcode",  "umi",      "read", ""),
  whitelist       = c("",                       "10x_5v1",  "",         "",     ""),
  stringsAsFactors = FALSE
)

layoutPath <- file.path(runDir, "my_layout.csv")
write.csv(layout, layoutPath, row.names = FALSE, quote = FALSE)
```

That describes: forward primer → 16 bp barcode (corrected against the
`10x_5v1` whitelist) → 10 bp UMI → the cDNA read → reverse primer. rad’s
own [layouts-and-whitelists
guide](https://github.com/indianewok/rad/blob/main/docs/layouts-and-whitelists.md)
documents the schema and has more worked examples.

## Start From a Bundled Template

Hand-writing a layout from scratch is error-prone. The faster path is to
copy a layout bundled with `rrad` that is close to yours and edit it. The
installed layout directory can be inspected without a RAD executable:

``` r
layoutDir <- system.file("rad", "resources", "read_layout", package = "rrad")
list.files(layoutDir)
```

Copy the closest CSV into your scratch directory and edit the copy —
change a barcode length, swap a primer, or add an element:

``` r
tmplPath <- file.path(layoutDir, "five_prime_read_layout.csv")

myLayout <- file.path(runDir, "five_prime_custom.csv")
file.copy(tmplPath, myLayout, overwrite = TRUE)

# ...edit myLayout in your editor or with data.table::fread()/fwrite()...
```

The built-in keys accepted by `rad_demux()` include `five_prime`,
`three_prime`, `sctagger`, `splitseq`, `visium`, `visium_hd`, and
`nanopore_rapid_bc`. If you already maintain the external CLI, the
legacy `radList(layouts = TRUE)` call prints the CLI’s registry,
including any of its local overrides; native `rrad` runs intentionally
use their packaged resources or the explicit file path you supply.

## Validate the Layout (legacy CLI helper)

Every native `rad_demux()` call validates the custom CSV before
processing and fails before publishing output if it is malformed. There
is not yet a separate native dry-run function that prints the parsed
layout. If you have the legacy CLI installed,
`radPrep(read_layout = TRUE)` provides that optional inspection: it
normalizes the layout, injects its internal start/stop sentinels,
auto-generates reverse-complement elements, and prints the ordered
structure. Nothing is written and no FASTQ is needed.

``` r
radPrep(layout = layoutPath, read_layout = TRUE)
```

Read the output back against your intent: are the elements in the right
order, are the static sequences correct, and did the variable elements
land with the lengths you expected? A malformed CSV (a missing column, a
`variable` row with no `expected_length`) surfaces here, not three steps
later inside demux.

## Check the Layout Against Real Reads

The native demultiplexer calculates the layout and position map as part
of the run. A bounded pilot checks how the layout lands on real data
while keeping the same engine and result contract used for the full run:

``` r
myFastq <- "my_reads.fq.gz"  # reads with the structure your layout describes

pilot <- rad_demux(
  layout = layoutPath,
  fastq = myFastq,
  out_dir = runDir,
  output = "my_pilot",
  max_reads = 50000L,
  threads = nproc
)

pilot$artifacts[c("layout", "position_map")]
```

If the position map shows your barcode or UMI consistently landing a few
bases off, that is the signal to adjust an `expected_length` or a primer
`seq` and re-run — cheaper to fix here than after a full demux. Native
runs recalibrate these artifacts by default; set `reuse_cache = TRUE`
only when they are known to belong to the same input and layout. The
legacy CLI alternative is
`radPrep(position_map = TRUE, fastq = myFastq, output = "...")`.

## Orientation and the `direction` Column

Long reads arrive in both orientations. By default — when you omit the
`direction` column entirely — rad treats every row as forward and
**auto- generates the reverse-complement** counterpart for you, so a
single forward description covers reads sequenced either way. That is
what the example above relies on.

Add a `direction` column only when an element is one-sided. Mark
specific rows `forward_only` or `reverse_only`:

``` r
oneSided <- data.frame(
  id              = c("adapter_rc",            "barcode"),
  seq             = c("AGATCGGAAGAGCGTCGTGTAG", ""),
  expected_length = c("",                       "16"),
  type            = c("static",                 "variable"),
  class           = c("adapter",                "barcode"),
  whitelist       = c("",                       "10x_5v1"),
  direction       = c("reverse_only",           "forward_only"),
  stringsAsFactors = FALSE
)
```

Use this when a primer or adapter genuinely only appears on one strand;
otherwise leave `direction` out and let rad build the reverse
complements.

## Whitelists for Barcode Elements

The `whitelist` column on a `barcode` row tells rad which barcode set to
correct against. You have two options:

1.  **A bundled kit key.** Put a registered key (`10x_3v3`,
    `10x_5v1`, `10x_Vis_V1`..`V5`, `splitseq_bc1/2`,
    `visium_hd_bc1/2`, …) in the `whitelist` column. Note the
    underscores — RAD expects `10x_3v3`, not `10xv3`. Native calls
    resolve packaged resources automatically and download a large
    built-in whitelist on first use when needed. If an external CLI is
    installed, legacy `radList(whitelists = TRUE)` prints its registry.

    ``` r
    # Optional: prepare all downloadable native whitelists for offline use.
    rad_download_whitelists()
    ```

2.  **Your own whitelist.** If your barcodes aren’t a known kit, build
    one from the reads with native `rad_scan_wl()` — it scans for
    barcodes downstream of an adapter/primer and writes the selected
    whitelist and scan statistics. It also supports paired BC1+BC2
    scans. See the [demultiplexing guide](demultiplexing-with-rad.html)
    for a worked call. Then point the layout’s `whitelist` column or
    `rad_demux(custom_whitelist = ...)` at
    `scan$artifacts$whitelist`.

    ``` r
    scan <- rad_scan_wl(
      input = myFastq,
      output_prefix = file.path(runDir, "my_cells"),
      adapter_seq = "CTACACGACGCTCTTCCGATCT",
      barcode_length = 16L,
      threads = nproc
    )
    ```

## Register the Layout as a Reusable Key (legacy CLI only)

Native runs accept an explicit layout CSV path directly and intentionally
do not consult the external CLI’s per-user registry. Keeping that path in
an R variable or configuration file is the most reproducible approach.
If you deliberately use the legacy CLI, `radModify()` can register a
short key for its camel-case wrappers:

``` r
radModify(layout_key = "my_kit", add = layoutPath)
radList(layouts = TRUE)$stdout   # `my_kit` now appears in the listing
```

``` r
radModify(layout_key = "my_kit", remove = TRUE)  # unregister
```

## Demultiplex With Your Layout

Once the layout validates, pass its CSV path to native `rad_demux()`
exactly where the demultiplexing guide used a bundled key. The barcode
element already names its whitelist, so a minimal run needs only the
layout, a FASTQ, an output directory, and a prefix:

``` r
nativeRun <- rad_demux(
  layout  = layoutPath,
  fastq   = myFastq,
  out_dir = runDir,
  output  = "custom",
  threads = nproc
)
```

This writes `custom.fq.gz` (corrected, `CB:Z`/`UB:Z`-tagged reads) plus
the whitelist summaries under `runDir`, and returns their paths, run
statistics, and the verified core identity. The native options from the
[demultiplexing guide](demultiplexing-with-rad.html) apply unchanged,
including `bc_correction_mode`, `joint_bc_mode = "strict"` for paired
barcodes, `custom_whitelist`, and `write_debug = TRUE`.

## Continue Natively, or Use the Legacy One-Call Wrapper

Native results compose directly. For example, this streams the tagged
FASTQ to a pRESTO/IgBLAST-ready FASTA without loading its records into a
data frame:

``` r
igblastInput <- rad_reformat(
  fastq = nativeRun,
  output_fastq = file.path(runDir, "custom-igblast.fasta"),
  to_fasta = TRUE,
  header_format = "presto",
  threads = nproc
)
```

There is no native `rad_run()` convenience function: call `rad_demux()`
and `rad_reformat()` explicitly so both structured results and their
provenance stay available. If an existing workflow intentionally uses
the external CLI, legacy `radRun()` chains `radPrep()` → `radDemux()` →
optional `radReformat()`:

``` r
radRun(
  layout  = layoutPath,
  fastq   = myFastq,
  out_dir = runDir,
  output  = "custom",
  threads = nproc
)
```

## seqspec YAML (Newer rad Builds)

The native engine also accepts a [Pachter-lab
seqspec](https://github.com/pachterlab/seqspec) YAML in place of the
layout CSV. Pass the `.yaml` path directly to
`rad_demux(layout = ...)`:

``` r
yamlRun <- rad_demux(
  layout = "my_assay.seqspec.yaml",
  fastq = myFastq,
  out_dir = runDir,
  output = "my_assay",
  threads = nproc
)
```

This path is newer than the CSV schema; if you already maintain a
seqspec for an assay it saves re-describing the structure. The legacy
CLI can optionally inspect it with
`radPrep(layout = "my_assay.seqspec.yaml", read_layout = TRUE)` when the
installed CLI build supports seqspec.

## Where to Next

  - [Demultiplexing Long Reads with RAD](demultiplexing-with-rad.html) —
    the full `prep → scan-wl → demux → reformat` pipeline on bundled
    layouts, including the Visium HD paired-barcode case.
  - [Visium HD processing](visium-hd-processing.html) — pick up from RAD
    output: spatial barcodes through pRESTO/IgBLAST to spatial AIRR
    tables.
  - rad’s own docs: [layouts and
    whitelists](https://github.com/indianewok/rad/blob/main/docs/layouts-and-whitelists.md),
    [CLI
    reference](https://github.com/indianewok/rad/blob/main/docs/cli-reference.md),
    [installation](https://github.com/indianewok/rad/blob/main/docs/installation.md).
