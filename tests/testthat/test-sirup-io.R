make_io_test_sirup <- function() {
  skip_if_not_installed("SeuratObject")
  spatial <- sirup_test_seurat(
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
    spatial_id = c("cell1", "cell2", "bulk1"),
    clone_id = c("A", "B", "A"),
    locus = c("IGH", "IGH", "IGH"),
    stringsAsFactors = FALSE
  )
  createSIRUP(
    spatial = spatial,
    receptors = receptors,
    spatial_id = "spatial_id",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
}

test_that("saveSIRUP/readSIRUP round-trips every component", {
  x <- make_io_test_sirup()
  file <- tempfile(fileext = ".sirup.rds")

  expect_invisible(saveSIRUP(x, file))
  expect_true(file.exists(file))

  y <- readSIRUP(file)
  expect_s3_class(y, "sirup")
  expect_identical(names(y), c("spatial", "receptors", "metadata", "polygons", "analyses"))
  expect_identical(getReceptors(y), getReceptors(x))
  expect_identical(y$metadata, x$metadata)
  expect_identical(getPolygons(y), getPolygons(x))
  expect_identical(y$analyses, x$analyses)
  expect_identical(colnames(y$spatial), colnames(x$spatial))
})

test_that("readSIRUP reads objects saved under the old SIRAP name", {
  x <- make_io_test_sirup()
  legacy <- x
  legacy$analyses$dominance <- structure(
    data.frame(clone_id = "A"),
    class = c("sirap_spatial_clone_dominance", "data.frame")
  )
  class(legacy) <- c("sirap", "list")
  attr(legacy, "sirap_pkg_version") <- as.character(utils::packageVersion("SIRUP"))
  file <- tempfile(fileext = ".sirap.rds")
  saveRDS(legacy, file)

  y <- readSIRUP(file)
  expect_s3_class(y, "sirup")
  expect_false(inherits(y, "sirap"))
  expect_s3_class(y$analyses$dominance, "sirup_spatial_clone_dominance")
  expect_identical(getReceptors(y), getReceptors(x))
})

test_that("saveSIRUP validates the object before writing", {
  x <- make_io_test_sirup()
  x$receptors[["spatial_id"]] <- rep(NA_character_, nrow(x$receptors))
  expect_error(saveSIRUP(x, tempfile(fileext = ".sirup.rds")))
})

test_that("readSIRUP re-validates and rejects a corrupted object", {
  x <- make_io_test_sirup()
  # Break the spatial join without going through saveSIRUP's validation.
  x$receptors[["spatial_id"]] <- c("nomatch1", "nomatch2", "nomatch3")
  file <- tempfile(fileext = ".sirup.rds")
  saveRDS(x, file)
  expect_error(readSIRUP(file), "match")
})

test_that("readSIRUP upgrades a legacy object missing polygons/analyses", {
  x <- make_io_test_sirup()
  legacy <- structure(
    list(spatial = x$spatial, receptors = x$receptors, metadata = x$metadata),
    class = c("sirup", "list")
  )
  file <- tempfile(fileext = ".sirup.rds")
  saveRDS(legacy, file)

  y <- readSIRUP(file)
  expect_identical(names(y), c("spatial", "receptors", "metadata", "polygons", "analyses"))
  expect_equal(nrow(getPolygons(y)), 0)
  expect_identical(y$analyses, list())
})

test_that("readSIRUP warns when the saving package version differs", {
  x <- make_io_test_sirup()
  file <- tempfile(fileext = ".sirup.rds")
  saveSIRUP(x, file)

  stamped <- readRDS(file)
  attr(stamped, "sirup_pkg_version") <- "0.0.0.0000"
  saveRDS(stamped, file)

  expect_warning(readSIRUP(file), "SIRUP")
})

test_that("readSIRUP errors on a non-SIRUP file", {
  file <- tempfile(fileext = ".rds")
  saveRDS(list(1, 2, 3), file)
  expect_error(readSIRUP(file), "SIRUP object")
})

test_that("readSIRUP errors when the file is missing", {
  expect_error(readSIRUP(tempfile(fileext = ".rds")), "does not exist")
})
