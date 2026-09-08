equivalence_script <- test_path("..", "..", "inst", "scripts", "check_rrad_equivalence.R")
if (!file.exists(equivalence_script)) {
  equivalence_script <- system.file("scripts", "check_rrad_equivalence.R", package = "SIRAP")
}
equivalence_helpers <- new.env(parent = globalenv())
sys.source(equivalence_script, equivalence_helpers)

test_that("equivalence ignores row ordering but keeps values and multiplicity", {
  compare <- equivalence_helpers$equivalence_compare
  x <- data.frame(sequence_id = c("a", "a", "b"), barcode = c("AA", "AA", "BB"),
                  v_call = c("V1", "V1", "V2"), consensus_count = c(1L, 1L, 2L),
                  sequence_alignment = c(".AC-G", ".AC-G", "TTGC"))
  expect_true(compare(x, x[c(3, 1, 2), rev(names(x))])$pass)
  expect_false(compare(x, x[-1, ])$pass)
  expect_false(compare(x, rbind(x, x[1, ]))$pass)
  changed <- x
  changed$barcode[1] <- "CC"
  expect_false(compare(x, changed)$pass)
  changed <- x
  changed$v_call[1] <- "V3"
  expect_false(compare(x, changed)$pass)
  changed <- x
  changed$sequence_alignment[1] <- "ACG"
  expect_false(compare(x, changed)$pass)
  expect_false(compare(x, x[, -2])$pass)
  expect_false(compare(x[FALSE, ], x[FALSE, ])$pass)
})

test_that("FASTA sampling preserves complete multiline records", {
  path <- tempfile(fileext = ".fasta")
  writeLines(c(">one|BARCODE=AC|UMI=GT", "AC", "GT", ">two|BARCODE=TG", "CA", ">three|BARCODE=CA", "TG"), path)
  sampled <- equivalence_helpers$equivalence_read_fasta(path, 2L)
  expect_identical(sampled$sequence, c("ACGT", "CA"))
  parsed <- equivalence_helpers$equivalence_parse_tags(sampled)
  expect_identical(parsed$sequence_id, c("one", "two"))
  expect_identical(parsed$barcode, c("AC", "TG"))
  expect_identical(parsed$umi, c("GT", NA_character_))
  unlink(path)
})

test_that("equivalence CLI requires flags and bounded integer thread counts", {
  parse <- equivalence_helpers$equivalence_options
  required <- c("--input", "reads.fasta", "--output-dir", "new_directory",
                "--igblast-base", "/tools/igblast", "--germline-glob", "/ref/*.fasta",
                "--muscle", "/tools/muscle")
  expect_identical(parse(c(required, "--threads", "2"))$threads, "2")
  invalid_flag <- required
  invalid_flag[1] <- "input"
  expect_error(parse(invalid_flag), "must start with --")
  for (threads in c("0", "-1", "1.5", "2147483648", "99999999999999999999")) {
    expect_error(parse(c(required, "--threads", threads)), "positive integer")
  }
})
