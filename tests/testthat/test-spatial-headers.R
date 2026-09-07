test_that("parseRadSequenceIds reads canonical pRESTO RAD metadata", {
  input <- data.frame(
    sequence_id = c(
      "read-F|BARCODE=AACC|UMI=TGCA",
      "read-R|UMI=CCCC|BARCODE=TTGA"
    ),
    stringsAsFactors = FALSE
  )

  out <- parseRadSequenceIds(input, require_barcode = TRUE)

  expect_identical(out$sequence_id, c("read-F", "read-R"))
  expect_identical(out$rad_read_id, c("read-F", "read-R"))
  expect_identical(out$barcode, c("AACC", "TTGA"))
  expect_identical(out$umi, c("TGCA", "CCCC"))
})

test_that("parseRadSequenceIds accepts native and legacy headers at ingress", {
  input <- data.frame(
    sequence_id = c(
      "native-F\tCR:Z:raw CB:Z:AACC TS:A:+ UB:Z:TGCA",
      "coord-R BC:Z:s_002um_00001_00002-1 CB:Z:AA-TT RG:Z:2um TS:A:-",
      "legacy|umi:GGGG|barcode:CCCC"
    ),
    stringsAsFactors = FALSE
  )

  out <- parseRadSequenceIds(input, require_barcode = TRUE)

  expect_identical(out$sequence_id, c("native-F", "coord-R", "legacy"))
  expect_identical(
    out$barcode,
    c("AACC", "s_002um_00001_00002-1", "CCCC")
  )
  expect_identical(out$umi, c("TGCA", NA_character_, "GGGG"))
  expect_identical(out$raw_barcode, c("raw", NA_character_, NA_character_))
  expect_identical(out$rad_cb, c("AACC", "AA-TT", NA_character_))
  expect_identical(out$rad_ts, c("+", "-", NA_character_))
})

test_that("RAD header parsing rejects missing barcodes and duplicate fields", {
  expect_error(
    parseRadSequenceIds(
      data.frame(sequence_id = "read-only"),
      require_barcode = TRUE
    ),
    "Missing spatial barcode"
  )
  expect_error(
    parseRadSequenceIds(data.frame(
      sequence_id = "read|BARCODE=AACC|barcode=TTGA"
    )),
    "Duplicate pRESTO field"
  )
})

test_that("spatial FASTA writer emits the canonical pRESTO grammar", {
  work <- withr::local_tempdir(pattern = "sirap-spatial-header-")
  path <- file.path(work, "spatial.fasta")
  airr <- data.frame(
    sequence_id = c("read-F", "read-R"),
    barcode = c("AACC", "TTGA"),
    umi = c("TGCA", NA_character_),
    sequence_alignment = c("ACGT", "TGCA"),
    stringsAsFactors = FALSE
  )

  result <- .sirap_write_spatial_fasta(airr, path)

  expect_identical(result, normalizePath(path))
  expect_identical(readLines(path), c(
    ">read-F|BARCODE=AACC|UMI=TGCA", "ACGT",
    ">read-R|BARCODE=TTGA", "TGCA"
  ))
})

test_that("spatial FASTA headers reject ambiguity and duplicate IDs", {
  expect_error(
    .sirap_format_spatial_headers("read|bad", "AACC"),
    "delimiter or whitespace"
  )
  expect_error(
    .sirap_format_spatial_headers("read", "AA,CC"),
    "delimiter or whitespace"
  )
  expect_error(
    .sirap_format_spatial_headers(c("read", "read"), c("AACC", "TTGA")),
    "must be unique"
  )
})
