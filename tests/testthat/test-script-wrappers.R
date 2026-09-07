test_that("paf2csv uses the shared RAD header codec", {
  work <- withr::local_tempdir(pattern = "sirap-paf-header-")
  paf <- file.path(work, "reads.paf")
  output <- file.path(work, "output")
  writeLines(c(
    paste(c(
      "read-1|BARCODE=AACC|UMI=TGCA", 100, 0, 100, "+",
      "target", 100, 0, 100, 100, 100, 60
    ), collapse = "\t"),
    paste(c(
      "read-2|umi:CCCC|barcode:TTGA", 100, 0, 100, "+",
      "target", 100, 0, 100, 100, 100, 60
    ), collapse = "\t")
  ), paf)

  result <- paf2csv(paf, output, include_rest = FALSE)
  processed <- data.table::fread(result$output_processed_csv)

  expect_identical(processed$id, c("read-1", "read-2"))
  expect_identical(processed$barcode, c("AACC", "TTGA"))
  expect_identical(processed$umi, c("TGCA", "CCCC"))
})
