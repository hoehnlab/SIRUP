make_io_test_sirap <- function() {
  skip_if_not_installed("SeuratObject")
  spatial <- sirap_test_seurat(
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
  createSIRAP(
    spatial = spatial,
    receptors = receptors,
    spatial_id = "spatial_id",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
}

test_that("saveSIRAP/readSIRAP round-trips every component", {
  x <- make_io_test_sirap()
  file <- tempfile(fileext = ".sirap.rds")

  expect_invisible(saveSIRAP(x, file))
  expect_true(file.exists(file))

  y <- readSIRAP(file)
  expect_s3_class(y, "sirap")
  expect_identical(names(y), c("spatial", "receptors", "metadata", "polygons", "analyses"))
  expect_identical(getReceptors(y), getReceptors(x))
  expect_identical(y$metadata, x$metadata)
  expect_identical(getPolygons(y), getPolygons(x))
  expect_identical(y$analyses, x$analyses)
  expect_identical(colnames(y$spatial), colnames(x$spatial))
})

test_that("saveSIRAP validates the object before writing", {
  x <- make_io_test_sirap()
  x$receptors[["spatial_id"]] <- rep(NA_character_, nrow(x$receptors))
  expect_error(saveSIRAP(x, tempfile(fileext = ".sirap.rds")))
})

test_that("readSIRAP re-validates and rejects a corrupted object", {
  x <- make_io_test_sirap()
  # Break the spatial join without going through saveSIRAP's validation.
  x$receptors[["spatial_id"]] <- c("nomatch1", "nomatch2", "nomatch3")
  file <- tempfile(fileext = ".sirap.rds")
  saveRDS(x, file)
  expect_error(readSIRAP(file), "match")
})

test_that("readSIRAP upgrades a legacy object missing polygons/analyses", {
  x <- make_io_test_sirap()
  legacy <- structure(
    list(spatial = x$spatial, receptors = x$receptors, metadata = x$metadata),
    class = c("sirap", "list")
  )
  file <- tempfile(fileext = ".sirap.rds")
  saveRDS(legacy, file)

  y <- readSIRAP(file)
  expect_identical(names(y), c("spatial", "receptors", "metadata", "polygons", "analyses"))
  expect_equal(nrow(getPolygons(y)), 0)
  expect_identical(y$analyses, list())
})

test_that("readSIRAP warns when the saving package version differs", {
  x <- make_io_test_sirap()
  file <- tempfile(fileext = ".sirap.rds")
  saveSIRAP(x, file)

  stamped <- readRDS(file)
  attr(stamped, "sirap_pkg_version") <- "0.0.0.0000"
  saveRDS(stamped, file)

  expect_warning(readSIRAP(file), "SIRAP")
})

test_that("readSIRAP errors on a non-SIRAP file", {
  file <- tempfile(fileext = ".rds")
  saveRDS(list(1, 2, 3), file)
  expect_error(readSIRAP(file), "SIRAP object")
})

test_that("readSIRAP errors when the file is missing", {
  expect_error(readSIRAP(tempfile(fileext = ".rds")), "does not exist")
})
