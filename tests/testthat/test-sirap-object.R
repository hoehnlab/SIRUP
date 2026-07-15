make_test_sirap <- function(barcode_col = "barcode", clone_col = "clone_id") {
  skip_if_not_installed("SeuratObject")
  spatial <- SeuratObject::CreateSeuratObject(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 4,
      dimnames = list("gene", paste0("cell", 1:4))
    )
  )
  spatial$x <- c(0, 1, 2, 3)
  spatial$y <- c(0, 1, 2, 3)
  receptors <- data.frame(
    barcode_value = c("cell1", "cell2", "bulk1"),
    clone_value = c("A", "B", "A"),
    locus = c("IGH", "IGH", "IGH"),
    stringsAsFactors = FALSE
  )
  names(receptors)[names(receptors) == "barcode_value"] <- barcode_col
  names(receptors)[names(receptors) == "clone_value"] <- clone_col
  createSIRAP(
    spatial = spatial,
    receptors = receptors,
    barcode_col = barcode_col,
    clone_col = clone_col
  )
}

test_that("createSIRAP creates a validated tibble-backed container", {
  x <- make_test_sirap()

  expect_s3_class(x, "sirap")
  expect_true(isSIRAP(x))
  expect_s3_class(x$receptors, "tbl_df")
  expect_identical(names(x), c("spatial", "receptors", "schema", "metadata"))
  expect_identical(x$schema$barcode, "barcode")
  expect_identical(x$schema$clone, "clone_id")
  expect_identical(x$schema$locus, "locus")
  expect_identical(receptorTable(x), x$receptors)
  expect_identical(spatialObject(x), x$spatial)
  expect_invisible(validateSIRAP(x))
})

test_that("createSIRAP retains unmatched nonspatial receptors", {
  x <- make_test_sirap()
  info <- summary(x)

  expect_equal(nrow(receptorTable(x)), 3)
  expect_equal(info$n_receptor_barcodes, 3)
  expect_equal(info$n_matched_barcodes, 2)
  expect_equal(info$n_matched_records, 2)
  expect_equal(info$barcode_match_rate, 2 / 3)
})

test_that("default optional schema fields are NULL when absent", {
  skip_if_not_installed("SeuratObject")
  spatial <- SeuratObject::CreateSeuratObject(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )
  x <- createSIRAP(spatial, data.frame(barcode = "cell1"))

  expect_null(x$schema$clone)
  expect_null(x$schema$locus)
  expect_error(
    createSIRAP(spatial, data.frame(barcode = "cell1"), clone_col = "missing"),
    "missing column"
  )
})

test_that("createSIRAP rejects receptor tables without spatial overlap", {
  skip_if_not_installed("SeuratObject")
  spatial <- SeuratObject::CreateSeuratObject(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )

  expect_error(
    createSIRAP(spatial, data.frame(barcode = "bulk1")),
    "No receptor barcodes match"
  )
  expect_error(
    createSIRAP(spatial, data.frame(barcode = "cell1"), metadata = list("unnamed")),
    "named list"
  )
})

test_that("SIRAP schema supplies receptor inputs to spatial analyses", {
  x <- make_test_sirap(barcode_col = "spot", clone_col = "clonal_group")

  variable <- spatialVariableData(
    object = x,
    variable = "receptor_present",
    variable_type = "receptor"
  )
  expect_equal(variable$value[match(paste0("cell", 1:4), variable$cell)], c(1, 1, 0, 0))

  polygon <- data.frame(
    x = c(-1, 4, 4, -1, -1),
    y = c(-1, -1, 4, 4, -1)
  )
  dominance <- spatialCloneDominance(x, polygon = polygon)
  expect_s3_class(dominance, "sirap_spatial_clone_dominance")
  expect_equal(sort(dominance$dominance$clone_id), c("A", "B"))
  expect_equal(sort(dominance$region_bcrs$.sirap_barcode), c("cell1", "cell2"))
})
