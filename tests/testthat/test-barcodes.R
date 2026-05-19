test_that("swapVisiumHdBarcode swaps row and column", {
  expect_equal(
    swapVisiumHdBarcode("s_002um_12_34-1"),
    "s_002um_34_12-1"
  )
})

test_that("addSpatialBarcodes maps raw sequence suffixes", {
  airr <- data.frame(sequence_id = c("read_12_34", "read_10_20"))
  mapping <- data.frame(
    square_002um = c("s_002um_12_34-1", "s_002um_10_20-1"),
    square_008um = c("s_008um_3_8-1", "s_008um_2_5-1"),
    square_016um = c("s_016um_2_4-1", "s_016um_1_3-1")
  )
  out <- addSpatialBarcodes(airr, mapping)
  expect_equal(out$barcode_8um, mapping$square_008um)
})
