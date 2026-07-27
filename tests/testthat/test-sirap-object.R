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
    v_sequence_start = c(1, 1, 1),
    v_sequence_end = c(101, 111, 121),
    j_sequence_start = c(130, 140, 150),
    j_sequence_end = c(160, 170, 180),
    junction_length = c(30, 33, 36),
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
  expect_identical(names(x), c("spatial", "receptors", "metadata", "analyses"))
  expect_identical(x$analyses, list())
  expect_identical(x$metadata$receptor_fields$barcode, "barcode")
  expect_identical(x$metadata$receptor_fields$clone, "clone_id")
  expect_identical(x$metadata$receptor_fields$locus, "locus")
  expect_identical(receptorTable(x), x$receptors)
  expect_identical(spatialObject(x), x$spatial)
  expect_identical(SIRAP:::.sirap_spatial_input(x), x$spatial)
  expect_identical(SIRAP:::.sirap_receptor_input(x)$receptors, x$receptors)
  expect_invisible(validateSIRAP(x))
})

test_that("analysis accessors store arbitrary named results", {
  x <- make_test_sirap()
  x <- setAnalysis(x, "clone_counts", data.frame(clone_id = "A", n = 2))
  x <- setAnalysis(x, "dowser_trees", list(tree = "(A,B);"))
  x <- setAnalysis(x, "empty_result", NULL)

  expect_identical(
    listAnalyses(x),
    c("clone_counts", "dowser_trees", "empty_result")
  )
  expect_equal(getAnalysis(x, "clone_counts")$n, 2)
  expect_identical(getAnalysis(x, "dowser_trees"), list(tree = "(A,B);"))
  expect_null(getAnalysis(x, "empty_result"))

  x <- removeAnalysis(x, "dowser_trees")
  expect_identical(listAnalyses(x), c("clone_counts", "empty_result"))
  expect_error(getAnalysis(x, "dowser_trees"), "is not stored")
  expect_error(setAnalysis(x, "", 1), "non-empty")
})

test_that("legacy SIRAP objects behave as though analyses are empty", {
  x <- make_test_sirap()
  x$analyses <- NULL

  expect_invisible(validateSIRAP(x))
  expect_identical(listAnalyses(x), character())

  x <- setAnalysis(x, "result", 1)
  expect_identical(getAnalysis(x, "result"), 1)
  expect_identical(names(x), c("spatial", "receptors", "metadata", "analyses"))
})

test_that("createSIRAP validates initial analyses", {
  x <- make_test_sirap()
  spatial <- spatialObject(x)
  receptors <- receptorTable(x)

  with_result <- createSIRAP(
    spatial,
    receptors,
    analyses = list(note = "already analyzed")
  )
  expect_identical(getAnalysis(with_result, "note"), "already analyzed")
  expect_error(
    createSIRAP(spatial, receptors, analyses = list(1)),
    "unique, non-empty names"
  )
  expect_error(
    createSIRAP(spatial, receptors, analyses = "not a list"),
    "named list"
  )
})

test_that("with and update receptor functions operate on the receptor table", {
  x <- make_test_sirap()
  x$metadata$receptor_fields["clone"] <- list(NULL)
  x$receptors$clone_id <- NULL

  n_records <- withReceptors(x, nrow)
  expect_identical(n_records, 3L)
  expect_null(x$metadata$receptor_fields$clone)

  x <- updateReceptors(
    x,
    function(receptors, prefix) {
      receptors$clone_id <- paste0(prefix, c("1", "2", "1"))
      receptors
    },
    prefix = "clone_"
  )
  expect_identical(x$metadata$receptor_fields$clone, "clone_id")
  expect_identical(x$receptors$clone_id, c("clone_1", "clone_2", "clone_1"))

  x <- updateReceptors(
    x,
    transform,
    clonal_group = clone_id,
    receptor_fields = list(clone = "clonal_group")
  )
  expect_identical(x$metadata$receptor_fields$clone, "clonal_group")
  expect_error(updateReceptors(x, function(receptors) 1), "must return a data frame")
})

test_that("with and update spatial functions operate on the Seurat object", {
  x <- make_test_sirap()

  expect_equal(withSpatial(x, ncol), 4)
  x <- updateSpatial(x, function(spatial, value) {
    spatial$new_field <- value
    spatial
  }, value = 1:4)

  expect_identical(unname(x$spatial$new_field), 1:4)
  expect_error(updateSpatial(x, function(spatial) data.frame()), "must return a Seurat")
})

test_that("with and update analysis functions operate on one named result", {
  x <- setAnalysis(make_test_sirap(), "scores", c(3, 1, 2))

  expect_identical(withAnalysis(x, "scores", sort), c(1, 2, 3))
  expect_identical(getAnalysis(x, "scores"), c(3, 1, 2))

  x <- updateAnalysis(x, "scores", sort, decreasing = TRUE)
  expect_identical(getAnalysis(x, "scores"), c(3, 2, 1))
  expect_error(withAnalysis(x, "missing", identity), "is not stored")
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

test_that("receptor summaries and QC plots accept SIRAP objects", {
  x <- make_test_sirap(barcode_col = "spot")

  receptor_summary <- spatialReceptorSummary(x)
  expect_equal(receptor_summary$n_sequences, 3)
  expect_equal(receptor_summary$n_barcodes, 3)
  expect_equal(receptor_summary$n_loci, 1)

  expect_s3_class(plotLengthQc(x), "ggplot")
})

test_that("default optional receptor fields are NULL when absent", {
  skip_if_not_installed("SeuratObject")
  spatial <- SeuratObject::CreateSeuratObject(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )
  x <- createSIRAP(spatial, data.frame(barcode = "cell1"))

  expect_null(x$metadata$receptor_fields$clone)
  expect_null(x$metadata$receptor_fields$locus)
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
  expect_error(
    createSIRAP(spatial, data.frame(barcode = "cell1"), metadata = list(receptor_fields = list())),
    "reserved"
  )
})

test_that("SIRAP metadata supplies receptor inputs to spatial analyses", {
  x <- make_test_sirap(barcode_col = "spot", clone_col = "clonal_group")

  polygon <- data.frame(
    x = c(-1, 4, 4, -1, -1),
    y = c(-1, -1, 4, 4, -1)
  )
  dominance <- spatialCloneDominance(x, polygon = polygon)
  expect_s3_class(dominance, "sirap_spatial_clone_dominance")
  expect_equal(dominance$metadata$receptor_fields$barcode, "spot")
  expect_equal(dominance$metadata$receptor_fields$clone, "clonal_group")
  expect_equal(sort(dominance$dominance$clone_id), c("A", "B"))
  expect_equal(sort(dominance$region_bcrs$.sirap_barcode), c("cell1", "cell2"))
})

test_that("spatialLocalMoran uses the receptor barcode stored in SIRAP metadata", {
  skip_if_not_installed("spdep")
  x <- make_test_sirap(barcode_col = "spot", clone_col = "clonal_group")

  local <- spatialLocalMoran(
    object = x,
    variable = "receptor_present",
    variable_type = "receptor",
    k = 1
  )

  expect_equal(local$value[match(paste0("cell", 1:4), local$cell)], c(1, 1, 0, 0))
})

test_that("SIRAP print reports spatial and receptor dimensions", {
  x <- make_test_sirap()

  expect_output(
    print(x),
    paste0(
      "<SIRAP object>\\n",
      "Spatial data:    Seurat, 4 barcodes, 1 gene, 0 images\\n",
      "Receptor data:   tibble, 3 receptor sequences x 8 columns\\n",
      "Spatial overlap: 2/3 receptor barcodes in spatial data \\(66.7%\\)\\n",
      "Barcode field:   barcode\\n",
      "Clone field:     clone_id\\n",
      "Analyses:        0"
    )
  )
})

test_that("legacy top-level schema objects are rejected", {
  x <- make_test_sirap()
  x$schema <- x$metadata$receptor_fields

  expect_error(validateSIRAP(x), "must contain, in order")
})

test_that("spatial polygon formatting accepts a SIRAP object", {
  x <- make_test_sirap()
  points <- data.frame(
    x = c(0, 3, 3, 0, 0),
    y = c(0, 0, 3, 3, 0),
    polygon = "region_1",
    image = "metadata"
  )

  out <- formatSpatialPolygons(x, points, point_type = "boundary")
  expect_equal(out$polygon, rep("region_1", 5))
  expect_equal(out$x, points$x)
  expect_equal(out$y, points$y)
})
