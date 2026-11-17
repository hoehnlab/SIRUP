square <- function(x0, y0, size = 1) {
  sprintf(
    "[[%s,%s],[%s,%s],[%s,%s],[%s,%s],[%s,%s]]",
    x0, y0, x0 + size, y0, x0 + size, y0 + size, x0, y0 + size, x0, y0
  )
}

qupath_feature <- function(rings, type = "Polygon", properties = "") {
  sprintf(
    '{"type":"Feature","id":"x","geometry":{"type":"%s","coordinates":%s},"properties":{%s}}',
    type, rings, properties
  )
}

qupath_file <- function(...) {
  path <- withr::local_tempfile(fileext = ".geojson", .local_envir = parent.frame())
  writeLines(
    paste0('{"type":"FeatureCollection","features":[', paste(c(...), collapse = ","), "]}"),
    path
  )
  path
}

endo <- '"objectType":"annotation","classification":{"name":"ENDO","color":[255,255,255]}'
peri <- '"objectType":"annotation","classification":{"name":"PERI","color":[0,255,51]}'

test_that("readQuPathPolygons reads polygons, classes, and IDs", {
  skip_if_not_installed("jsonlite")
  path <- qupath_file(
    qupath_feature(paste0("[", square(0, 0), "]"), properties = endo),
    qupath_feature(paste0("[", square(5, 5), "]"), properties = endo),
    qupath_feature(paste0("[", square(9, 9), "]"), properties = peri)
  )
  out <- readQuPathPolygons(path)
  expect_equal(out$polygon_id, c("ENDO_01", "ENDO_02", "PERI_01"))
  expect_equal(out$class, c("ENDO", "ENDO", "PERI"))
  expect_equal(out$vertices[[2]]$x, c(5, 6, 6, 5, 5))
  expect_equal(out$vertices[[2]]$y, c(5, 5, 6, 6, 5))
  expect_false("sample" %in% names(out))

  peri_only <- readQuPathPolygons(path, class = "PERI", sample_id = "s1")
  expect_equal(peri_only$polygon_id, "PERI_01")
  expect_equal(peri_only$sample, "s1")
})

test_that("readQuPathPolygons uses names and skips non-annotations", {
  skip_if_not_installed("jsonlite")
  ring <- paste0("[", square(0, 0), "]")
  path <- qupath_file(
    qupath_feature(ring, properties = paste0(endo, ',"name":"GC1"')),
    qupath_feature(ring, properties = '"objectType":"annotation","name":"box"'),
    qupath_feature(ring, properties = paste0('"objectType":"detection",', '"classification":{"name":"ENDO"}')),
    qupath_feature("[0,0]", type = "Point", properties = endo)
  )
  out <- readQuPathPolygons(path)
  expect_equal(out$polygon_id, c("GC1", "box"))
  expect_true(is.na(out$class[[2]]))
  classified <- readQuPathPolygons(path, keep_unclassified = FALSE)
  expect_equal(classified$polygon_id, "GC1")
  expect_error(readQuPathPolygons(path, class = "none"), "No QuPath annotations")
})

test_that("readQuPathPolygons splits MultiPolygons and warns about holes", {
  skip_if_not_installed("jsonlite")
  holed <- paste0("[", square(0, 0, 10), ",", square(2, 2, 2), "]")
  multi <- paste0("[[", square(20, 20), "],[", square(30, 30), "]]")
  path <- qupath_file(
    qupath_feature(holed, properties = paste0(endo, ',"name":"holed"')),
    qupath_feature(multi, type = "MultiPolygon", properties = paste0(peri, ',"name":"multi"'))
  )
  expect_warning(out <- readQuPathPolygons(path), "holes of 1 annotation")
  expect_equal(out$polygon_id, c("holed", "multi_p1", "multi_p2"))
  # Only the outer ring stays.
  expect_equal(nrow(out$vertices[[1]]), 5L)
  expect_equal(max(out$vertices[[1]]$x), 10)
})

test_that("readQuPathPolygons scales coordinates and checks input", {
  skip_if_not_installed("jsonlite")
  path <- qupath_file(qupath_feature(paste0("[", square(10, 20, 10), "]"), properties = endo))
  out <- readQuPathPolygons(path, scale = 0.5)
  expect_equal(out$vertices[[1]]$x, c(5, 10, 10, 5, 5))
  expect_equal(out$vertices[[1]]$y, c(10, 10, 15, 15, 10))
  expect_error(readQuPathPolygons(path, scale = 0), "positive number")
  expect_error(readQuPathPolygons("missing.geojson"), "does not exist")
})

test_that("readQuPathPolygons output works with addPolygons", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("sf")
  skip_if_not_installed("SeuratObject")
  x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
  path <- qupath_file(qupath_feature(paste0("[", square(2.5, 2.5), "]"), properties = endo))
  x <- addPolygons(x, readQuPathPolygons(path))
  expect_equal(getPolygons(x)$polygon_id, "ENDO_01")
  expect_true(any(!is.na(getSpatial(x)$polygon_id)))
})

test_that("readQuPathPolygons reads older bare lists and malformed properties", {
  skip_if_not_installed("jsonlite")
  ring <- paste0("[", square(0, 0), "]")
  path <- withr::local_tempfile(fileext = ".json")
  writeLines(paste0(
    "[",
    qupath_feature(ring, properties = '"classification":{"name":"Tumor","colorRGB":-3670016}'), ",",
    qupath_feature(ring, properties = '"classification":{"names":["a","b"]}'),
    "]"
  ), path)
  out <- readQuPathPolygons(path)
  expect_equal(out$polygon_id, c("Tumor_01", "annotation_01"))
  expect_equal(out$class, c("Tumor", NA))
})

test_that("addQuPathPolygons stores polygons and assigns spots", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("sf")
  skip_if_not_installed("SeuratObject")
  x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
  path <- qupath_file(qupath_feature(
    paste0("[", square(3.5, 0.5, 2), "]"),
    properties = paste0(endo, ',"name":"GC1"')
  ))
  stored_before <- nrow(getPolygons(x))
  out <- addQuPathPolygons(x, path)
  expect_equal(nrow(getPolygons(out)), stored_before + 1L)
  expect_true("GC1" %in% getPolygons(out)$polygon_id)
  expect_equal(getPolygons(out)$class[getPolygons(out)$polygon_id == "GC1"], "ENDO")
  inside <- rownames(getSpatial(out)@meta.data)[!is.na(getSpatial(out)$polygon_id)]
  expect_true(all(c("spot04", "spot05", "spot09", "spot10") %in% inside))

  # The same ID cannot be added twice, unless the new polygon replaces the old one.
  expect_error(addQuPathPolygons(out, path), "already stored")
  replaced <- addQuPathPolygons(out, path, replace = TRUE)
  expect_equal(getPolygons(replaced)$polygon_id, "GC1")
  expect_error(addQuPathPolygons(x, path, image_scale = "lowres"), "differs")

  # A buffer adds spots whose circles touch the polygon.
  grown <- addQuPathPolygons(x, path, buffer = 1)
  n_inside <- function(object) sum(!is.na(getSpatial(object)$polygon_id))
  expect_gt(n_inside(grown), n_inside(out))
})

test_that("addQuPathPolygons needs an image scale", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("sf")
  skip_if_not_installed("SeuratObject")
  x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
  x$metadata$image_scale <- NULL
  path <- qupath_file(qupath_feature(paste0("[", square(3.5, 0.5, 2), "]"), properties = endo))
  expect_error(addQuPathPolygons(x, path), "image_scale")
  stored <- addQuPathPolygons(x, path, image_scale = "hires", swap_coords = FALSE)
  expect_equal(stored$metadata$image_scale, "hires")
  expect_error(addQuPathPolygons(x, path, image_scale = "hires", buffer = -1), "buffer")
})
