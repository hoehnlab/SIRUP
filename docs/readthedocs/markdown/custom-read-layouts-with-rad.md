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
custom format and run the same pipeline on it, all from R through
SIRAP’s `rad*()` wrappers.

As with the other rad guide, SIRAP does not reimplement rad — every
wrapper here shells out to the same `rad` binary you would call from a
terminal. The wrappers used below are `radList()`, `radPrep()`,
`radModify()`, and `radDemux()`, with pointers to `radScanWl()` for
building the whitelists your layout references.

This guide does not evaluate code during rendering: the
layout-validation steps shell out to the external `rad` binary, and the
demux steps need a FASTQ that matches your custom structure. Copy the
chunks into a live R session against your own reads.

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

# Confirm R can see the rad binary (resolution order: explicit `rad =`,
# options(SIRAP.rad_path=), $SIRAP_RAD_PATH, then `rad` on your PATH).
radPath()
#> "/Users/you/code/rad/build/rad"
```

If rad is not on your `PATH`, point SIRAP at it once per session — see
the [demultiplexing
guide](demultiplexing-with-rad.html#point-r-at-the-rad-binary) for the
full resolution rules and `updateRad()` for building rad from source.

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
copy a bundled layout that is close to yours and edit it. `radList()`
prints each registered layout key *and the CSV it resolves to* — start
there:

``` r
radList(layouts = TRUE)$stdout
#> [layouts]
#> five_prime          .../resources/read_layout/five_prime_read_layout.csv
#> three_prime         .../resources/read_layout/three_prime_read_layout.csv
#> sctagger            .../resources/read_layout/sctagger_sim_read_layout.csv
#> splitseq            .../resources/read_layout/splitseq_read_layout.csv
#> visium              .../resources/read_layout/visium_three_prime_read_layout.csv
#> visium_hd           .../resources/read_layout/visium_hd_read_layout.csv
#> nanopore_rapid_bc   .../resources/read_layout/nanopore_bulk_rapid_bc_read_layout.csv
```

Grab the path for the closest key, copy it into your scratch dir, and
edit the copy — change a barcode length, swap a primer, add an element:

``` r
# Pull the resolved path for `five_prime` out of the listing and copy it.
listing  <- radList(layouts = TRUE)$stdout
tmplPath <- sub(".*\\s", "", grep("five_prime", strsplit(listing, "\n")[[1]],
                                  value = TRUE)[1])

myLayout <- file.path(runDir, "five_prime_custom.csv")
file.copy(tmplPath, myLayout, overwrite = TRUE)

# ...edit myLayout in your editor or with data.table::fread()/fwrite()...
```

## Validate the Layout

Before touching any reads, have rad parse the CSV and show you how it
interpreted each element. `radPrep(read_layout = TRUE)` does exactly
this — it normalizes the layout, injects its internal start/stop
sentinels, auto-generates reverse- complement elements, and prints the
ordered structure. Nothing is written and no FASTQ is needed, so it is a
fast syntax check:

``` r
radPrep(layout = layoutPath, read_layout = TRUE)
```

Read the output back against your intent: are the elements in the right
order, are the static sequences correct, and did the variable elements
land with the lengths you expected? A malformed CSV (a missing column, a
`variable` row with no `expected_length`) surfaces here, not three steps
later inside demux.

## Check the Layout Against Real Reads

`read_layout` only checks the *description*. To see how the layout
actually lands on your data, add `position_map = TRUE`. rad samples
reads, estimates where each element falls, and reports misalignment
statistics. This mode needs a FASTQ and an output base path; rad writes
`<output>_layout.csv` and `<output>_position_map.csv`, which `demux`
reuses on later runs to the same prefix.

``` r
myFastq <- "my_reads.fq.gz"  # reads with the structure your layout describes

radPrep(
  layout       = layoutPath,
  position_map = TRUE,
  fastq        = myFastq,
  output       = file.path(runDir, "my_prep"),
  threads      = nproc
)
```

If the position map shows your barcode or UMI consistently landing a few
bases off, that is the signal to adjust an `expected_length` or a primer
`seq` and re-run — cheaper to fix here than after a full demux.

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

1.  **A bundled kit key.** `radList(whitelists = TRUE)` prints the
    registered kits (`10x_3v3`, `10x_5v1`, `10x_Vis_V1`..`V5`,
    `splitseq_bc1/2`, `visium_hd_bc1/2`, …). Put the key in the
    `whitelist` column. Note the underscores — rad expects `10x_3v3`,
    not `10xv3`.
    
    ``` r
    radList(whitelists = TRUE)$stdout
    ```

2.  **Your own whitelist.** If your barcodes aren’t a known kit, build
    one from the reads with `radScanWl()` — it scans for the barcodes
    sitting downstream of an adapter/primer and writes a `.csv`/`.txt`
    whitelist, with a split BC1+BC2 mode for paired barcodes. See the
    [demultiplexing
    guide](demultiplexing-with-rad.html#discover-barcodes-de-novo-optional)
    for a worked `radScanWl()` run. Then point the layout’s `whitelist`
    column (or `radDemux(custom_whitelist = ...)`) at the file you
    generated.

## Register the Layout as a Reusable Key (Optional)

Passing the CSV path around works everywhere, but if you reuse a layout
a lot you can register it under a short key with `radModify()`, so
`radList()` and the `layout =` argument accept the name:

``` r
radModify(layout_key = "my_kit", add = layoutPath)
radList(layouts = TRUE)$stdout   # `my_kit` now appears in the listing
```

``` r
radModify(layout_key = "my_kit", remove = TRUE)  # unregister
```

**Caveat (current build):** rad’s registry edits are *process-local* — a
key added in one `rad` invocation does not persist to the next, separate
one. For scripts and reproducible pipelines, prefer passing the layout
**CSV path** directly to `layout =` (as every chunk below does); treat
`radModify()` as an interactive convenience.

## Demultiplex With Your Layout

Once the layout validates, it is a drop-in for the `layout =` argument
of every rad wrapper — pass your CSV path exactly where the demux guide
passed a bundled key. The barcode element already names its whitelist,
so a minimal run needs only the layout, a FASTQ, an output directory,
and a prefix:

``` r
radDemux(
  layout  = layoutPath,
  fastq   = myFastq,
  out_dir = runDir,
  output  = "custom",
  threads = nproc
)
```

This writes `custom.fq.gz` (corrected, `CB:Z`/`UB:Z`-tagged reads) plus
the whitelist summaries under `runDir`. From here everything in the
[demultiplexing guide](demultiplexing-with-rad.html#demultiplex) applies
unchanged — `bc_correction_mode`, `joint_bc_mode = "strict"` for paired
barcodes, `custom_whitelist =`, `write_dbg = TRUE`, and the
`radReformat()` / `radRun()` steps all take your layout the same way. To
chain prep → demux → reformat in one call:

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

Recent rad builds also accept a [Pachter-lab
seqspec](https://github.com/pachterlab/seqspec) YAML in `prep`/`demux`
in place of the layout CSV — hand the `.yaml` path to `layout =` and the
wrappers pass it straight through:

``` r
radPrep(layout = "my_assay.seqspec.yaml", read_layout = TRUE)
```

This path is newer than the CSV schema; if you already maintain a
seqspec for an assay it saves re-describing the structure. Confirm your
`rad` is recent enough (`updateRad()`), and check rad’s docs for the
current state of seqspec support.

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
