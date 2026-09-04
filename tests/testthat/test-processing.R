test_that("gap length helpers work", {
  expect_equal(.sirap_max_gap_length("AA...TT..G"), 3L)
  expect_equal(.sirap_opening_gap_length("...A.."), 3L)
  expect_equal(.sirap_opening_gap_length("A.."), 0L)
})

test_that("filterAirr filters by locus and lengths", {
  airr <- data.frame(
    locus = c("IGH", "IGK", "IGH"),
    v_sequence_start = c(1, 1, 1),
    v_sequence_end = c(120, 120, 50),
    j_sequence_start = c(1, 1, 1),
    j_sequence_end = c(20, 20, 20),
    junction_length = c(12, 12, 12),
    productive = c(TRUE, TRUE, TRUE)
  )
  out <- filterAirr(airr, locus_val = "IGH", productive_only = TRUE)
  expect_equal(nrow(out), 1L)
  expect_equal(out$locus, "IGH")
})

test_that("filterAirr handles character productive and complete_vdj", {
  airr <- data.frame(
    v_sequence_start = rep(1, 4),
    v_sequence_end = rep(120, 4),
    j_sequence_start = rep(1, 4),
    j_sequence_end = rep(20, 4),
    junction_length = rep(12, 4),
    productive = c("T", "F", "T", "T"),
    complete_vdj = c("T", "T", "F", NA),
    stringsAsFactors = FALSE
  )
  out <- filterAirr(airr, productive_only = TRUE, complete_vdj = FALSE)
  expect_equal(nrow(out), 3L)

  out <- filterAirr(airr, productive_only = TRUE, complete_vdj = TRUE)
  expect_equal(nrow(out), 1L)
  expect_equal(out$complete_vdj, "T")

  # complete_vdj = TRUE is a no-op when the column is absent.
  out <- filterAirr(airr[, setdiff(colnames(airr), "complete_vdj")], complete_vdj = TRUE)
  expect_equal(nrow(out), 3L)
})

test_that("filterAirr builds spatial FASTA headers that writeFastqas round-trips", {
  airr <- data.frame(
    sequence_id = c("read1", "read2"),
    sequence = c("ACGTACGT", "TTTTAAAA"),
    v_sequence_start = c(1, 1),
    v_sequence_end = c(120, 120),
    j_sequence_start = c(1, 1),
    j_sequence_end = c(20, 20),
    junction_length = c(12, 12),
    barcode = c("BC01", "BC02"),
    stringsAsFactors = FALSE
  )
  out <- filterAirr(
    airr,
    min_consensus_count = NULL,
    spatial_id = "barcode",
    add_spatial_sequence_id = TRUE,
    spatial_id_key = "barcode"
  )
  expect_equal(out$spatial_sequence_id, c(
    "sequence_id=read1|barcode=BC01",
    "sequence_id=read2|barcode=BC02"
  ))

  fa <- file.path(withr::local_tempdir(), "spatial.fasta")
  writeFastqas(
    data.frame(id = out$spatial_sequence_id, seq = out$sequence, stringsAsFactors = FALSE),
    fn = fa,
    type = "fa",
    header_col = "id"
  )
  lines <- readLines(fa)
  expect_equal(lines[1], ">sequence_id=read1|barcode=BC01")
  expect_equal(lines[3], ">sequence_id=read2|barcode=BC02")
})
