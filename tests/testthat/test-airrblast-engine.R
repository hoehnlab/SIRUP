# These tests run processSpatialSample() on a spatial FASTA. The pRESTO commands and
# the two airrblast functions are replaced, so no IgBLAST is needed. They check what
# SIRUP does with the airrblast result: the header fields, the pass filter, and the
# arguments it passes on.

.run_with_airrblast <- function(assigned, receptor = "ig", ...) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  fasta <- file.path(dir, "s.fasta")
  writeLines(c(">a", "ACGT"), fasta)
  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    assignGenes = function(file, locus, ...) {
      seen$assign_locus <- locus
      assigned
    },
    addGaps = function(db, gapdb, locus, ...) {
      seen$gapdb <- gapdb
      seen$gaps_locus <- locus
      db
    },
    .package = "airrblast",
    .env = parent.frame()
  )
  testthat::local_mocked_bindings(
    .sirup_check_tools = function(...) invisible(TRUE),
    .sirup_run_command = function(...) invisible(0L),
    .package = "SIRUP",
    .env = parent.frame()
  )
  result <- suppressMessages(processSpatialSample(
    fasta, input_type = "spatial_fasta", output_dir = dir, receptor = receptor, ...
  ))
  list(
    table = data.table::fread(result$db_pass, sep = "\t", colClasses = "character",
                              na.strings = "", data.table = FALSE),
    seen = seen
  )
}

.assigned_rows <- function(sequence_id, junction = rep("TGT", length(sequence_id)),
                           productive = rep(TRUE, length(sequence_id))) {
  data.frame(
    sequence_id = sequence_id, v_call = "V", j_call = "J", productive = productive,
    sequence_alignment = "ACGT", junction = junction, c_call = "TRBC1*01",
    stringsAsFactors = FALSE
  )
}

test_that("pRESTO header fields become columns, as MakeDb.py --extended did", {
  # The barcode in the first header contains a comma, as in real consensus headers.
  assigned <- .assigned_rows(c(
    "1,AAA|CONSCOUNT=3|CC_BC_UMI=4|N_PRODUCTIVE_READS=3|C_CALL=TRBC2*01|C_CALL_FREQ=1.0",
    "2,CCC|CONSCOUNT=1"
  ), productive = c(TRUE, FALSE))
  out <- .run_with_airrblast(assigned, germline_glob = "/g/imgt/human/vdj/imgt_human_*.fasta")$table
  expect_equal(out$sequence_id, c("1,AAA", "2,CCC"))
  expect_equal(out$consensus_count, c("3", "1"))
  expect_equal(out$cc_bc_umi, c("4", NA))
  expect_equal(out$n_productive_reads, c("3", NA))
  expect_equal(out$c_call_freq, c("1.0", NA))
  # A header value replaces the IgBLAST value. A row without that key keeps it.
  expect_equal(out$c_call, c("TRBC2*01", "TRBC1*01"))
  # Logical columns are written as T/F, as MakeDb.py wrote them.
  expect_equal(out$productive, c("T", "F"))
})

test_that("plain IDs and malformed headers are left alone", {
  plain <- .run_with_airrblast(.assigned_rows(c("a", "b")),
                               germline_glob = "/g/imgt/human/vdj/x_*.fasta")$table
  expect_equal(plain$sequence_id, c("a", "b"))
  expect_false("consensus_count" %in% names(plain))

  expect_warning(
    bad <- .run_with_airrblast(.assigned_rows("a|NOEQUALS"),
                               germline_glob = "/g/imgt/human/vdj/x_*.fasta")$table,
    "without `=`"
  )
  expect_equal(bad$sequence_id, "a|NOEQUALS")
})

test_that("rows without a junction, V call, J call, alignment or productive flag are dropped", {
  assigned <- .assigned_rows(c("keep", "no_junction", "no_flag"),
                             junction = c("TGT", NA, "TGT"),
                             productive = c(TRUE, TRUE, NA))
  out <- .run_with_airrblast(assigned, germline_glob = "/g/imgt/human/vdj/x_*.fasta")$table
  expect_equal(out$sequence_id, "keep")
})

test_that("the IMGT folder comes from germline_glob unless imgt_dir is given", {
  glob <- "/g/imgt/human/vdj/imgt_human_*.fasta"
  inferred <- .run_with_airrblast(.assigned_rows("a"), germline_glob = glob)$seen
  expect_equal(inferred$gapdb, "/g/imgt")
  given <- .run_with_airrblast(.assigned_rows("a"), germline_glob = glob, imgt_dir = "/other/imgt")$seen
  expect_equal(given$gapdb, "/other/imgt")

  expect_error(processSpatialSample("x.fasta", input_type = "spatial_fasta", germline_glob = "/g/other/*.fasta"),
               "Cannot infer")
  expect_error(processSpatialSample("x.fasta", input_type = "spatial_fasta",
                                    germline_glob = "/g/imgt/mouse/vdj/*.fasta"), "Cannot infer")
  expect_error(processSpatialSample("x.fasta", input_type = "spatial_fasta"), "imgt_dir")
})

test_that("both airrblast functions get the locus name they accept", {
  tr <- .run_with_airrblast(.assigned_rows("a"), receptor = "tr",
                            germline_glob = "/g/imgt/human/vdj/x_*.fasta")$seen
  expect_equal(tr$assign_locus, "TR")
  expect_equal(tr$gaps_locus, "TR")
  ig <- .run_with_airrblast(.assigned_rows("a"), receptor = "ig",
                            germline_glob = "/g/imgt/human/vdj/x_*.fasta")$seen
  expect_equal(ig$assign_locus, "Ig")
  expect_equal(ig$gaps_locus, "Ig")
})
