###############################################################################
# SIRAP RAD walkthrough for 10x Visium HD 3' data
#
# This script is a section-by-section walkthrough for using SIRAP's R wrappers
# around Chinmay's `rad` command line tool. It is written for the BCR/TCR
# RAGE-seq data nested under:
#
#   /Volumes/HoehnK/data/Hoehn/260120_10X_Visium_HD_3p/RAGEseq_Hoehn_Feb2026
#
# The sibling `fastq/` and `analysis/` directories under the project root are
# the spatial RNA-seq data. We use `analysis/` only as Space Ranger context for
# spatial coordinate annotation; the receptor reads come from the RAGE-seq
# subdirectory.
#
# The goal is to make the RAD part of the workflow readable and repeatable:
#
#   1. point SIRAP at the `rad` binary
#   2. inspect the layouts and whitelist kits RAD knows about
#   3. check the Visium HD read layout on a small subset of reads
#   4. run RAD demultiplexing for each sample
#   5. optionally reformat RAD output FASTQs for downstream tools
#   6. optionally annotate RAD barcode-pair counts with Space Ranger coordinates
#
# Run this file one section at a time the first time through. Once the settings
# are correct, set dryRun <- FALSE and run the sample loop.
#
#   setwd("/Volumes/HoehnK/Hunter/spatial_bcr_tcr_package/SIRAP")
#   pkgload::load_all()
#   source("inst/scripts/visium_hd_3p_rad_walkthrough.R")
###############################################################################

# ---- Setup ------------------------------------------------------------------

sirapLoadError <- NULL
sirapFunctionNames <- c(
  "radPath", "radSrcPath", "radList", "radHelp", "radPrep", "radScanWl",
  "radDemux", "radRun", "radReformat", "buildAnnotatedPairs"
)

sirapFunctionsAvailable <- function() {
  all(vapply(sirapFunctionNames, exists, logical(1), mode = "function"))
}

if (!requireNamespace("SIRAP", quietly = TRUE) && !sirapFunctionsAvailable()) {
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    tryCatch(
      pkgload::load_all("."),
      error = function(e) {
        sirapLoadError <<- conditionMessage(e)
      }
    )
  } else if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
    tryCatch(
      devtools::load_all("."),
      error = function(e) {
        sirapLoadError <<- conditionMessage(e)
      }
    )
  }
}

if (!requireNamespace("SIRAP", quietly = TRUE) && !sirapFunctionsAvailable()) {
  stop(
    "SIRAP is not installed or loaded. Run this script from the SIRAP package ",
    "directory with pkgload::load_all(), or install SIRAP first.",
    if (!is.null(sirapLoadError)) paste0("\n\nLoad error was:\n", sirapLoadError) else "",
    call. = FALSE
  )
}

if (!is.null(sirapLoadError) && sirapFunctionsAvailable()) {
  warning(
    "SIRAP's R functions loaded, but the compiled shared object did not load. ",
    "The RAD wrapper sections in this walkthrough should still be usable. ",
    "Rebuild/reinstall SIRAP before running sections that call compiled code.\n\n",
    "Load error was:\n",
    sirapLoadError,
    call. = FALSE
  )
}

if (requireNamespace("SIRAP", quietly = TRUE) && !sirapFunctionsAvailable()) {
  for (fn in sirapFunctionNames) {
    assign(fn, get(fn, envir = asNamespace("SIRAP"), inherits = FALSE), envir = environment())
  }
}

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Things to change --------------------------------------------------------

# Keep this TRUE until the sample table, RAD binary path, layout, and whitelist
# choices all look right. In dry-run mode, the script prints what it would do
# without starting long RAD jobs.
dryRun <- TRUE

# RAD jobs can be large. Start with 4-8 threads on a laptop or login node; use
# more only on a machine where that is appropriate.
nThreads <- 8L

# This is the root directory Hunter pointed us to. Most subdirectories here are
# spatial RNA-seq outputs. The BCR/TCR receptor reads are under `rageSeqDir`.
dataDir <- "/Volumes/HoehnK/data/Hoehn/260120_10X_Visium_HD_3p"
spaceRangerAnalysisDir <- file.path(dataDir, "analysis")
rageSeqDir <- file.path(dataDir, "RAGEseq_Hoehn_Feb2026")
rageSeqRunDir <- file.path(
  rageSeqDir,
  "Run1_021626",
  "20260216_1440_P2S-00525-A_PBE80350_8f821e29"
)
rageSeqFastqPassDir <- file.path(rageSeqRunDir, "fastq_pass")

# Write RAD outputs into the same project directory, but keep them separated
# from the raw FASTQs and Space Ranger outputs.
radOutRoot <- file.path(dataDir, "rad_walkthrough")

# Layout tells RAD how to interpret bases in the read. For this Visium HD 3'
# project, `visium` is the expected built-in layout key. Confirm it appears in
# the `radList(layouts = TRUE)` output below before running the full jobs.
layoutKey <- "visium"

# Whitelist kits control barcode correction. RAD's canonical kit names are
# shown by `radList(whitelists = TRUE)`. Visium HD uses two spatial ID
# halves, so the two-part whitelist scan uses both of these.
bc1WhitelistKit <- "visium_hd_bc1"
bc2WhitelistKit <- "visium_hd_bc2"

# `radDemux()` currently accepts one `kit` argument. Leave this NULL unless
# Chinmay's RAD instructions for this layout say to pass one of the whitelist
# kits directly to demux. `radScanWl()` below is the section that explicitly
# uses both Visium HD whitelist halves.
demuxKit <- NULL

# `radScanWl()` searches for spatial IDs relative to an adapter/primer
# sequence. This sequence is experiment/layout-specific, so the walkthrough
# skips that section until this value is filled in. Ask Chinmay or check the RAD
# layout CSV if this section is needed.
adapterSeq <- NA_character_

# Use a small number of reads when checking layout and whitelist settings. This
# keeps the first pass fast and catches obvious mistakes before the full run.
qcMaxReads <- 100000L

# Full run setting. Keep NULL for all reads; set to a number such as 1000000L
# when testing the full demux command on a smaller subset.
demuxMaxReads <- NULL

# This sample table is the main place a student should edit when new datasets
# arrive. The receptor run was already demultiplexed by Nanopore barcode. The
# `barcode_alignment_*.tsv` file shows barcode23 and barcode24 are the large
# receptor libraries; these also match the older bc23/bc24 sample convention
# used in our BCR/TCR walkthrough.
#
# `receptor_fastq_dir` points to many FASTQ chunks for each sample. The helper
# `prepareFastqForRad()` below concatenates those chunks into one gzipped FASTQ
# per sample before RAD runs, because RAD's R wrapper accepts one FASTQ path.
#
# `spaceranger_dir` still points to the spatial RNA-seq Space Ranger outputs.
# Those files are not the BCR/TCR reads; they are used later to add Visium HD
# coordinates and in-tissue flags to RAD barcode-pair output.
samples <- data.table::data.table(
  sample_id = c("M1", "M5"),
  nanopore_barcode = c("barcode23", "barcode24"),
  receptor_fastq_dir = file.path(
    rageSeqFastqPassDir,
    c("barcode23", "barcode24")
  ),
  spaceranger_dir = file.path(
    spaceRangerAnalysisDir,
    c("M1_260129/outs/binned_outputs", "M5_260129/outs/binned_outputs")
  )
)

samples[, rad_out_dir := file.path(radOutRoot, sample_id)]
samples[, output_prefix := paste0(tolower(sample_id), "_visium_hd_3p")]
samples[, receptor_fastq := file.path(rad_out_dir, paste0(output_prefix, "_rage_pass.fastq.gz"))]

# ---- Small helper functions --------------------------------------------------

section <- function(title) {
  cat("\n\n# ---- ", title, " ", paste(rep("-", max(1, 72 - nchar(title))), collapse = ""), "\n", sep = "")
}

say <- function(...) {
  message(sprintf(...))
}

checkFiles <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(label, " missing:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}

checkDirs <- function(paths, label) {
  missing <- paths[!dir.exists(paths)]
  if (length(missing) > 0L) {
    stop(label, " missing:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}

fastqChunks <- function(fastq_dir) {
  sort(list.files(fastq_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE))
}

prepareFastqForRad <- function(fastq_dir, out_fastq, overwrite = FALSE) {
  chunks <- fastqChunks(fastq_dir)
  if (length(chunks) == 0L) {
    stop("No .fastq.gz chunks found in: ", fastq_dir, call. = FALSE)
  }

  if (file.exists(out_fastq) && !isTRUE(overwrite)) {
    say("Using existing concatenated FASTQ: %s", out_fastq)
    return(out_fastq)
  }

  dir.create(dirname(out_fastq), recursive = TRUE, showWarnings = FALSE)

  # Gzip permits concatenated members, so joining .fastq.gz chunks byte-for-byte
  # creates a valid gzipped FASTQ without decompressing the full dataset.
  ok <- file.copy(chunks[1], out_fastq, overwrite = TRUE)
  if (!ok) stop("Could not create: ", out_fastq, call. = FALSE)

  if (length(chunks) > 1L) {
    con_out <- file(out_fastq, open = "ab")
    on.exit(close(con_out), add = TRUE)
    for (chunk in chunks[-1]) {
      con_in <- file(chunk, open = "rb")
      repeat {
        buf <- readBin(con_in, what = "raw", n = 1024L * 1024L)
        if (length(buf) == 0L) break
        writeBin(buf, con_out)
      }
      close(con_in)
    }
  }

  out_fastq
}

runOrPrint <- function(label, expr, dry_run = dryRun) {
  section(label)
  if (isTRUE(dry_run)) {
    say("dryRun is TRUE, so this section is not running yet.")
    say("Review the code in this block, then set dryRun <- FALSE when ready.")
    return(invisible(NULL))
  }
  force(expr)
}

printRadResult <- function(result) {
  if (is.null(result)) return(invisible(NULL))
  if (!is.null(result$command)) {
    cat("Command:\n")
    cat(paste(shQuote(result$command), collapse = " "), "\n\n")
  }
  if (!is.null(result$stdout) && nzchar(result$stdout)) {
    cat(result$stdout, "\n")
  }
  invisible(result)
}

# Find a RAD CSV that looks like a valid barcode-pair count table. This is used
# after demux/scan-wl because RAD output filenames can differ by version.
findValidPairsCsv <- function(sample_out_dir) {
  csvs <- list.files(
    sample_out_dir,
    pattern = "\\.csv(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(csvs) == 0L) return(NA_character_)

  for (path in csvs) {
    dt <- tryCatch(data.table::fread(path, nrows = 20L), error = function(e) NULL)
    if (is.null(dt)) next
    needed <- c("bc1", "bc2", "pair", "count")
    if (all(needed %in% names(dt))) return(path)
  }

  NA_character_
}

# ---- Confirm inputs ----------------------------------------------------------

section("Confirm input files")

stopifnot(dir.exists(dataDir))
stopifnot(dir.exists(spaceRangerAnalysisDir))
stopifnot(dir.exists(rageSeqDir))
stopifnot(dir.exists(rageSeqFastqPassDir))

checkDirs(samples$receptor_fastq_dir, "RAGE-seq receptor FASTQ directories")
checkDirs(samples$spaceranger_dir, "Space Ranger binned output directories")

samples[, receptor_fastq_chunks := lengths(lapply(receptor_fastq_dir, fastqChunks))]
if (any(samples$receptor_fastq_chunks == 0L)) {
  stop("At least one receptor FASTQ directory has no .fastq.gz chunks.", call. = FALSE)
}

print(samples)

# ---- Locate RAD and inspect its registered resources -------------------------

section("Locate RAD")

# `radPath()` resolves the command line binary. It checks, in order:
#   1. an explicit path passed to radPath()
#   2. options(SIRAP.rad_path = "/path/to/rad")
#   3. Sys.getenv("SIRAP_RAD_PATH")
#   4. a `rad` executable on PATH
#
# If this fails, install/build RAD and then either put it on PATH or set:
#
#   options(SIRAP.rad_path = "/absolute/path/to/rad")
#
radBin <- tryCatch(
  radPath(),
  error = function(e) {
    if (isTRUE(dryRun)) {
      say("RAD binary was not found, but dryRun is TRUE so the walkthrough can still be reviewed.")
      say("RAD lookup error: %s", conditionMessage(e))
      return(NA_character_)
    }
    stop(e)
  }
)

if (!is.na(radBin)) {
  say("Using RAD binary: %s", radBin)
}

# These listings are worth looking at before every new project. The layout list
# confirms that `layoutKey` exists. The whitelist list confirms the exact kit
# names; RAD is strict about underscores and capitalization.
section("List RAD layouts")
if (is.na(radBin)) {
  say("Skipping radList(layouts = TRUE) until RAD is installed/configured.")
} else {
  layouts <- radList(layouts = TRUE, capture = TRUE, rad = radBin)
  printRadResult(layouts)
}

section("List RAD whitelist kits")
if (is.na(radBin)) {
  say("Skipping radList(whitelists = TRUE) until RAD is installed/configured.")
} else {
  whitelists <- radList(whitelists = TRUE, capture = TRUE, rad = radBin)
  printRadResult(whitelists)
}

# `radHelp()` is a thin wrapper around `rad help <topic>`. It is useful when a
# student wants to see the raw command line documentation without leaving R.
section("Show RAD help for the main commands")
if (is.na(radBin)) {
  say("Skipping radHelp() until RAD is installed/configured.")
} else {
  printRadResult(radHelp("prep", rad = radBin))
  printRadResult(radHelp("demux", rad = radBin))
  printRadResult(radHelp("reformat", rad = radBin))
}

# ---- Optional: update the RAD binary ----------------------------------------

# `updateRad()` pulls the RAD source tree and rebuilds the binary. It is useful
# when Chinmay has pushed fixes upstream. It is intentionally not run by default
# because it changes software state and can take a few minutes.
#
# Uncomment after confirming `radSrcPath()` finds the source clone:
#
#   radSrcPath()
#   updateRad(jobs = nThreads)

# ---- Layout check with radPrep() --------------------------------------------

runOrPrint("Check read layout and position map with radPrep()", {
  for (i in seq_len(nrow(samples))) {
    sample <- samples[i]
    dir.create(sample$rad_out_dir, recursive = TRUE, showWarnings = FALSE)
    receptorFastq <- prepareFastqForRad(sample$receptor_fastq_dir, sample$receptor_fastq)

    # `radPrep()` wraps `rad prep`.
    #
    # read_layout = TRUE asks RAD to print/write how it interprets the layout:
    # which bases are barcode, UMI, spacer, and biological sequence.
    #
    # position_map = TRUE asks RAD to sample reads and estimate where expected
    # sequence elements appear. This is a quick sanity check that we selected
    # the right layout and FASTQ before launching the full demultiplexing run.
    prep <- radPrep(
      layout = layoutKey,
      read_layout = TRUE,
      position_map = TRUE,
      fastq = receptorFastq,
      output = file.path(sample$rad_out_dir, paste0(sample$output_prefix, "_prep_qc")),
      max_reads = qcMaxReads,
      threads = nThreads,
      verbose = TRUE,
      capture = FALSE,
      rad = radBin
    )

    printRadResult(prep)
  }
})

# Expected outputs from this section:
#
#   <sample_out>/<prefix>_prep_qc_layout.csv
#   <sample_out>/<prefix>_prep_qc_position_map.csv
#
# Open the position map CSV if RAD reports unexpected offsets or low match
# rates. A bad position map usually means the wrong read, wrong layout, or wrong
# chemistry assumptions.

# ---- Optional whitelist scan with radScanWl() --------------------------------

runOrPrint("Generate/check Visium HD barcode-pair whitelist with radScanWl()", {
  if (is.na(adapterSeq) || !nzchar(adapterSeq)) {
    say("Skipping radScanWl(): adapterSeq is not set.")
  } else {
    for (i in seq_len(nrow(samples))) {
      sample <- samples[i]
      dir.create(sample$rad_out_dir, recursive = TRUE, showWarnings = FALSE)
      receptorFastq <- prepareFastqForRad(sample$receptor_fastq_dir, sample$receptor_fastq)

      # `radScanWl()` wraps `rad scan-wl`.
      #
      # This command scans reads around an adapter/primer sequence and counts
      # valid barcode pairs. For Visium HD, we pass the BC1 and BC2 whitelist
      # kit names so RAD can validate both halves of the spatial ID.
      scan <- radScanWl(
        adapter_seq = adapterSeq,
        input = receptorFastq,
        output_prefix = file.path(sample$rad_out_dir, paste0(sample$output_prefix, "_scan_wl")),
        max_reads = demuxMaxReads,
        bc1_whitelist = bc1WhitelistKit,
        bc2_whitelist = bc2WhitelistKit,
        threads = nThreads,
        verbose = TRUE,
        capture = FALSE,
        rad = radBin
      )

      printRadResult(scan)
    }
  }
})

# ---- Full demultiplexing with radDemux() -------------------------------------

runOrPrint("Run full RAD demultiplexing with radDemux()", {
  for (i in seq_len(nrow(samples))) {
    sample <- samples[i]
    dir.create(sample$rad_out_dir, recursive = TRUE, showWarnings = FALSE)
    receptorFastq <- prepareFastqForRad(sample$receptor_fastq_dir, sample$receptor_fastq)

    # `radDemux()` wraps `rad demux`, the main RAD processing command.
    #
    # Important arguments:
    #   layout: tells RAD how to parse the read.
    #   fastq: the concatenated RAGE-seq receptor FASTQ for this sample.
    #   kit: optional whitelist kit for barcode correction; left NULL here
    #        because Visium HD has two barcode halves and scan-wl handles both.
    #   output: filename prefix for all outputs from this sample.
    #   out_dir: directory where RAD should write this sample's outputs.
    #   max_reads: NULL means process all reads; set a number while testing.
    #   write_dbg: TRUE writes extra debug/metrics files that are useful for
    #              checking a new layout or chemistry.
    demux <- radDemux(
      layout = layoutKey,
      fastq = receptorFastq,
      kit = demuxKit,
      max_reads = demuxMaxReads,
      output = sample$output_prefix,
      out_dir = sample$rad_out_dir,
      log_file = file.path(sample$rad_out_dir, paste0(sample$output_prefix, ".rad.log")),
      write_dbg = TRUE,
      threads = nThreads,
      verbose = TRUE,
      capture = FALSE,
      rad = radBin
    )

    printRadResult(demux)
  }
})

# After demux finishes, inspect each sample directory:
#
#   list.files(samples$rad_out_dir[1], recursive = TRUE)
#
# The exact filenames can vary by RAD version, but useful outputs usually
# include a demultiplexed/tagged FASTQ, CSV count tables, debug files, and a log.

# ---- Optional one-call pipeline with radRun() --------------------------------

# `radRun()` is a convenience wrapper that chains:
#
#   radPrep() -> radDemux() -> optionally radReformat()
#
# It is nice once settings are stable, but the walkthrough uses separate calls
# above so a new student can see what each step does. This block is the compact
# version for later.
#
# runOrPrint("Run prep + demux as one pipeline with radRun()", {
#   for (i in seq_len(nrow(samples))) {
#     sample <- samples[i]
#     radRun(
#       layout = layoutKey,
#       fastq = prepareFastqForRad(sample$receptor_fastq_dir, sample$receptor_fastq),
#       out_dir = sample$rad_out_dir,
#       output = sample$output_prefix,
#       kit = demuxKit,
#       threads = nThreads,
#       position_map = TRUE,
#       do_reformat = FALSE,
#       max_reads = demuxMaxReads,
#       write_dbg = TRUE,
#       verbose = TRUE,
#       capture = FALSE,
#       rad = radBin
#     )
#   }
# })

# ---- Optional FASTQ reformatting with radReformat() --------------------------

runOrPrint("Reformat RAD FASTQ outputs with radReformat()", {
  for (i in seq_len(nrow(samples))) {
    sample <- samples[i]

    # By default `radRun()` expects demux to write <out_dir>/<output>.fq.gz.
    # We use the same convention here, but skip gracefully if the RAD version
    # wrote its FASTQ under a different name.
    demuxFastq <- file.path(sample$rad_out_dir, paste0(sample$output_prefix, ".fq.gz"))
    if (!file.exists(demuxFastq)) {
      say("Skipping reformat for %s; no FASTQ found at %s", sample$sample_id, demuxFastq)
      next
    }

    # `radReformat()` wraps `rad reformat`.
    #
    # reformat_header = TRUE converts headers to a compact QNAME_CB_UB form,
    # which many downstream tools prefer.
    #
    # coordinate = TRUE asks RAD to rewrite Visium HD coordinate information
    # using its default coordinate mode. `bin_size = 2` keeps 2 um resolution.
    #
    # split_bc = FALSE keeps one output file. Set TRUE only if you really need
    # one FASTQ per barcode; that can create many files.
    reformattedDir <- file.path(sample$rad_out_dir, "reformatted")
    dir.create(reformattedDir, recursive = TRUE, showWarnings = FALSE)

    reformat <- radReformat(
      fastq = demuxFastq,
      out_dir = reformattedDir,
      split_bc = FALSE,
      reformat_header = TRUE,
      reformat_delim = "_",
      coordinate = TRUE,
      bin_size = 2L,
      threads = nThreads,
      verbose = TRUE,
      capture = FALSE,
      rad = radBin
    )

    printRadResult(reformat)
  }
})

# ---- Optional coordinate annotation with buildAnnotatedPairs() ----------------

runOrPrint("Annotate RAD barcode-pair counts with Space Ranger coordinates", {
  # `buildAnnotatedPairs()` is not a RAD CLI wrapper. It is an SIRAP helper that
  # takes RAD's valid barcode-pair counts and joins them to:
  #
  #   1. RAD's Visium HD coordinate reference
  #   2. Space Ranger 2 um tissue positions
  #   3. Space Ranger filtered barcode sets at 2 um, 8 um, and 16 um
  #
  # The output is a practical bridge table for spatial QC and plotting.
  radSrc <- radSrcPath()
  radCoords <- file.path(radSrc, "resources", "wl", "visium_hd_coordinates.csv.gz")
  checkFiles(radCoords, "RAD Visium HD coordinate reference")

  for (i in seq_len(nrow(samples))) {
    sample <- samples[i]
    validPairs <- findValidPairsCsv(sample$rad_out_dir)

    if (is.na(validPairs)) {
      say("No valid barcode-pair CSV found for %s under %s", sample$sample_id, sample$rad_out_dir)
      next
    }

    annotatedPath <- file.path(
      sample$rad_out_dir,
      paste0(sample$output_prefix, "_annotated_pairs.csv.gz")
    )

    annotated <- buildAnnotatedPairs(
      valid_pairs = validPairs,
      rad_coords = radCoords,
      spaceranger_dir = sample$spaceranger_dir,
      out_csv = annotatedPath,
      verbose = TRUE
    )

    say("Wrote annotated pairs for %s to %s", sample$sample_id, annotatedPath)
    print(utils::head(annotated))
  }
})

# ---- Minimal QC checklist ----------------------------------------------------

section("QC checklist")

cat(
  "\nAfter running with dryRun <- FALSE, check these before handing outputs onward:\n",
  "\n1. RAD logs exist and do not end in an error:",
  "\n   file.path(samples$rad_out_dir, paste0(samples$output_prefix, '.rad.log'))",
  "\n\n2. The prep position-map CSVs show the expected barcode/UMI positions.",
  "\n\n3. Demux output directories contain non-empty FASTQ and/or CSV outputs.",
  "\n\n4. If buildAnnotatedPairs() ran, the annotated CSV has non-zero rows and",
  "\n   sensible `in_tissue`, `in_sr_filtered_2um`, and `rad_valid_16um` rates.",
  "\n\n5. Keep dryRun <- TRUE in committed/shared copies of this script so a",
  "\n   student has to intentionally turn on long-running processing.\n",
  sep = ""
)
