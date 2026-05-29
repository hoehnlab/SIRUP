test_that("extractSpatialBarcode extracts raw sequence suffixes", {
  airr <- data.frame(sequence_id = c("read_12_34", "read_10_20"))
  out <- extractSpatialBarcode(
    airr,
    barcode_col = "barcode_2um",
    regex = "[0-9]+_[0-9]+$",
    prefix = "s_002um_",
    suffix = "-1"
  )
  expect_equal(out$barcode_2um, c("s_002um_12_34-1", "s_002um_10_20-1"))
})

test_that("extractSpatialBarcode uses a capture group when present", {
  reads <- data.frame(sequence_id = "read,s_002um_12_34-1")
  out <- extractSpatialBarcode(reads, barcode_col = "barcode_2um", regex = ",([^,]+)")
  expect_equal(out$barcode_2um, "s_002um_12_34-1")
})
