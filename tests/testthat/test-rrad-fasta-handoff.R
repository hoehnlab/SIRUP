# These tests verify the file contract and orchestration only. External
# IgBLAST, pRESTO, and Change-O calls are deliberately mocked.

.rrad_handoff_airr <- function(ids, barcodes = NULL, umis = NULL) {
  n <- length(ids)
  out <- data.frame(
    sequence_id = ids,
    sequence = rep(paste0("NNN", strrep("ACGT", 30), "NNN"), n),
    sequence_alignment = rep(strrep("ACGT", 30), n),
    v_sequence_start = rep(1L, n),
    v_sequence_end = rep(100L, n),
    j_sequence_start = rep(101L, n),
    j_sequence_end = rep(120L, n),
    junction_length = rep(12L, n),
    productive = rep("T", n),
    complete_vdj = rep("T", n),
    stringsAsFactors = FALSE
  )
  if (!is.null(barcodes)) out$barcode <- barcodes
  if (!is.null(umis)) out$umi <- umis
  out
}

test_that("current rrad FASTA headers and legacy tags give the same metadata", {
  fixture <- test_path("fixtures", "rrad-handoff", "canonical.fasta")
  lines <- readLines(fixture)
  ids <- sub("^>", "", lines[startsWith(lines, ">")])
  expect_equal(ids, c(
    "read-1|BARCODE=AACC|UMI=TGCA",
    "read-2|BARCODE=TTGA"
  ))
  canonical <- parseRadSequenceIds(data.frame(sequence_id = ids))
  legacy <- parseRadSequenceIds(data.frame(sequence_id = c(
    "read-1|umi:TGCA|barcode:AACC", "read-2|barcode:TTGA"
  )))
  expect_equal(canonical$sequence_id, c("read-1", "read-2"))
  expect_equal(canonical$barcode, c("AACC", "TTGA"))
  expect_equal(canonical$umi, c("TGCA", NA_character_))
  expect_equal(canonical, legacy)
  expect_equal(
    parseRadSequenceIds(data.frame(sequence_id = ids), strip = FALSE)$sequence_id,
    ids
  )
  mixed_case <- parseRadSequenceIds(data.frame(
    sequence_id = "mixed|barcode=AACC|uMi=TGCA"
  ))
  expect_equal(mixed_case$barcode, "AACC")
  expect_equal(mixed_case$umi, "TGCA")
})

test_that("header parsing preserves tabulated metadata and rejects conflicts", {
  airr <- data.frame(
    sequence_id = c("bare", "tagged|BARCODE=TTGA|UMI=CGTA"),
    barcode = c("AACC", "TTGA"),
    umi = c("TGCA", NA_character_),
    stringsAsFactors = FALSE
  )
  parsed <- parseRadSequenceIds(airr)
  expect_equal(parsed$barcode, c("AACC", "TTGA"))
  expect_equal(parsed$umi, c("TGCA", "CGTA"))
  expect_equal(parsed$sequence_id, c("bare", "tagged"))

  airr$barcode[[2L]] <- "DIFFERENT"
  expect_error(parseRadSequenceIds(airr))
  airr$barcode[[2L]] <- "TTGA"
  airr$umi[[2L]] <- "DIFFERENT"
  expect_error(parseRadSequenceIds(airr))
})

test_that("ambiguous or malformed metadata is rejected before grouping", {
  bad_ids <- c(
    "read|BARCODE=AACC|barcode=TTGA",
    "read|BARCODE=AACC|barcode:AACC",
    "read|BARCODE=AACC|UMI=TGCA|umi=TGCA",
    "read|BARCODE=AA CC",
    "read|BARCODE=AA,CC",
    "read|BARCODE=AA=CC"
  )
  for (id in bad_ids) {
    expect_error(parseRadSequenceIds(data.frame(sequence_id = id)), info = id)
  }
})

test_that("rrad FASTA enters the read annotation pass before all six consensus commands", {
  fixture <- test_path("fixtures", "rrad-handoff", "canonical.fasta")
  work <- withr::local_tempdir(pattern = "sirap handoff ")
  germline <- file.path(work, "germline reference.fasta")
  writeLines(c(">reference", strrep("ACGT", 30)), germline)
  observed <- new.env(parent = emptyenv())
  observed$steps <- character()
  observed$commands <- list()
  testthat::local_mocked_bindings(
    runIgblast = function(input, output_tsv, ...) {
      observed$steps <- c(observed$steps, "read_igblast")
      expect_equal(normalizePath(input), normalizePath(fixture))
      lines <- readLines(input)
      ids <- sub("^>", "", lines[startsWith(lines, ">")])
      data.table::fwrite(.rrad_handoff_airr(ids), output_tsv, sep = "\t")
      invisible(output_tsv)
    },
    .sirap_check_tools = function(...) invisible(TRUE),
    .sirap_run_command = function(exec, args, ...) {
      observed$steps <- c(observed$steps, exec)
      observed$commands[[exec]] <- args
      invisible(0L)
    },
    .package = "SIRAP"
  )
  result <- processSpatialSample(
    input = fixture,
    output_dir = work,
    receptor = "tr",
    igblast_db_dir = work,
    germline_glob = germline,
    nproc = 2
  )
  expect_equal(observed$steps, c(
    "read_igblast", "ClusterSets.py", "ParseHeaders.py", "AlignSets.py",
    "BuildConsensus.py", "AssignGenes.py", "MakeDb.py"
  ))
  expect_equal(result$n_raw, 2L)
  expect_equal(result$n_filtered, 2L)
  filtered <- data.table::fread(result$filtered_airr)
  expect_equal(filtered$barcode, c("AACC", "TTGA"))
  expect_equal(filtered$umi[[1L]], "TGCA")
  spatial <- readLines(result$spatial_fasta)
  expect_match(spatial[[1L]], "barcode=AACC", ignore.case = TRUE)
  expect_match(spatial[[3L]], "barcode=TTGA", ignore.case = TRUE)
  expect_equal(spatial[c(2L, 4L)], filtered$sequence_alignment)
  expect_false(any(spatial[c(2L, 4L)] == filtered$sequence))
  second_pass <- observed$commands[["AssignGenes.py"]]
  expect_equal(second_pass[match("-s", second_pass) + 1L], result$consensus)
  expect_equal(second_pass[match("--loci", second_pass) + 1L], "tr")
  expect_false(fixture %in% second_pass)
})

test_that("AIRR input keeps existing metadata and skips the raw read pass", {
  work <- withr::local_tempdir()
  germline <- file.path(work, "reference.fasta")
  writeLines(c(">reference", strrep("ACGT", 30)), germline)
  testthat::local_mocked_bindings(
    runIgblast = function(...) stop("AIRR input must skip read-level IgBLAST"),
    .sirap_check_tools = function(...) invisible(TRUE),
    .sirap_run_command = function(...) invisible(0L),
    .package = "SIRAP"
  )
  result <- processSpatialSample(
    input = .rrad_handoff_airr("plain-read", "AACC", "TGCA"),
    output_dir = work,
    germline_glob = germline
  )
  filtered <- data.table::fread(result$filtered_airr)
  expect_null(result$raw_airr)
  expect_equal(filtered$barcode, "AACC")
  expect_equal(filtered$umi, "TGCA")
})

test_that("missing barcodes in surviving reads fail before consensus", {
  work <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    .sirap_check_tools = function(...) stop("Consensus must not start"),
    .sirap_run_command = function(...) stop("Consensus must not start"),
    .package = "SIRAP"
  )
  for (missing_barcode in c(NA_character_, "")) {
    expect_error(
      processSpatialSample(
        input = .rrad_handoff_airr(
          c("good", "missing"), c("AACC", missing_barcode)
        ),
        output_dir = work,
        germline_glob = file.path(work, "unused.fasta")
      ),
      "[Bb]arcode|[Ss]patial"
    )
  }
})

test_that("reads removed by quality filtering need no spatial barcode", {
  work <- withr::local_tempdir()
  germline <- file.path(work, "reference.fasta")
  writeLines(c(">reference", strrep("ACGT", 30)), germline)
  testthat::local_mocked_bindings(
    .sirap_check_tools = function(...) invisible(TRUE),
    .sirap_run_command = function(...) invisible(0L),
    .package = "SIRAP"
  )
  airr <- .rrad_handoff_airr(c("good", "discard"), c("AACC", NA_character_))
  airr$productive[[2L]] <- "F"
  result <- processSpatialSample(
    input = airr,
    output_dir = work,
    germline_glob = germline
  )
  expect_equal(result$n_raw, 2L)
  expect_equal(result$n_filtered, 1L)
  expect_equal(data.table::fread(result$filtered_airr)$barcode, "AACC")
})
