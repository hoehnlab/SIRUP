make_test_sirap <- function(spatial_id = "spatial_id", clone = "clone_id") {
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
    spatial_id_value = c("cell1", "cell2", "bulk1"),
    clone_value = c("A", "B", "A"),
    locus = c("IGH", "IGH", "IGH"),
    v_sequence_start = c(1, 1, 1),
    v_sequence_end = c(101, 111, 121),
    j_sequence_start = c(130, 140, 150),
    j_sequence_end = c(160, 170, 180),
    junction_length = c(30, 33, 36),
    stringsAsFactors = FALSE
  )
  names(receptors)[names(receptors) == "spatial_id_value"] <- spatial_id
  names(receptors)[names(receptors) == "clone_value"] <- clone
  createSIRAP(
    spatial = spatial,
    receptors = receptors,
    spatial_id = spatial_id,
    clone = clone,
    image_scale = "hires",
    swap_coords = FALSE
  )
}

make_collapse_test_sirap <- function() {
  skip_if_not_installed("SeuratObject")
  coords <- expand.grid(x = seq_len(5), y = seq_len(5))
  cells <- paste0("cell", seq_len(nrow(coords)))
  spatial <- SeuratObject::CreateSeuratObject(
    counts = matrix(
      1,
      nrow = 1,
      ncol = length(cells),
      dimnames = list("gene", cells)
    )
  )
  spatial$x <- coords$x
  spatial$y <- coords$y

  receptors <- data.frame(
    spot = c("cell1", "cell2", "cell25", "cell13"),
    sequence_alignment = c("AAAA", "AAAA", "AAAA", "CCCC"),
    consensus_count = c(1L, 5L, 3L, 2L),
    locus = "IGH",
    stringsAsFactors = FALSE
  )
  createSIRAP(
    spatial = spatial,
    receptors = receptors,
    spatial_id = "spot",
    locus = "locus",
    count = "consensus_count",
    image_scale = "hires",
    swap_coords = FALSE
  )
}

test_that("createSIRAP creates a validated tibble-backed container", {
  x <- make_test_sirap()

  expect_s3_class(x, "sirap")
  expect_true(isSIRAP(x))
  expect_s3_class(x$receptors, "tbl_df")
  expect_identical(names(x), c("spatial", "receptors", "metadata", "polygons", "analyses"))
  expect_identical(x$analyses, list())
  expect_identical(x$metadata$receptor_fields$spatial_id, "spatial_id")
  expect_identical(x$metadata$receptor_fields$clone, "clone_id")
  expect_identical(x$metadata$receptor_fields$locus, "locus")
  expect_identical(x$metadata$image_scale, "hires")
  expect_false(x$metadata$swap_coords)
  expect_equal(nrow(getPolygons(x)), 0)
  expect_identical(getReceptors(x), x$receptors)
  expect_identical(getSpatial(x), x$spatial)
  expect_identical(SIRAP:::.sirap_spatial_input(x), x$spatial)
  expect_identical(SIRAP:::.sirap_receptor_input(x)$receptors, x$receptors)
  expect_invisible(validateSIRAP(x))
})

test_that("SIRAP objects store canonical polygon geometry", {
  x <- make_test_sirap()
  polygon_table <- tibble::tibble(
    polygon_id = c("GC1", "GC2"),
    vertices = list(
      tibble::tibble(spatial_id = c("cell1", "cell2", "cell3", "cell1"), x = c(0, 1, 2, 0), y = c(0, 1, 0, 0)),
      tibble::tibble(spatial_id = c("cell2", "cell3", "cell4", "cell2"), x = c(1, 2, 3, 1), y = c(1, 0, 1, 1))
    )
  )

  with_polygons <- addPolygons(x, polygons = polygon_table)
  expect_equal(getPolygons(with_polygons), polygon_table)
  expect_invisible(validateSIRAP(with_polygons))

  flat <- data.frame(
    polygon = c("GC1", "GC1", "GC1", "GC2", "GC2", "GC2"),
    spatial_id = c("cell1", "cell2", "cell3", "cell2", "cell3", "cell4"),
    x = c(0, 1, 2, 1, 2, 3),
    y = c(0, 1, 0, 1, 0, 1)
  )
  flat_path <- tempfile(fileext = ".csv")
  utils::write.csv(flat, flat_path, row.names = FALSE)
  from_path <- addPolygons(x, polygons_path = flat_path)
  expect_equal(from_path$polygons$polygon_id, c("GC1", "GC2"))
  expect_equal(from_path$metadata$polygons_path, flat_path)
  expect_error(addPolygons(x, polygons = polygon_table, polygons_path = flat_path), "only one")
})

test_that("addPolygons annotates SIRAP spatial metadata and receptors", {
  x <- make_test_sirap()
  polygon <- data.frame(
    polygon = "GC1",
    x = c(-0.5, 1.5, 1.5, -0.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5, -0.5)
  )

  out <- addPolygons(x, polygons = polygon)

  expect_equal(unname(out$spatial$polygon_id), c("GC1", "GC1", NA, NA))
  expect_equal(out$receptors$polygon_id, c("GC1", "GC1", NA))
  expect_invisible(validateSIRAP(out))
})

test_that("createSIRAP accepts a polygon path", {
  x <- make_test_sirap()
  flat <- data.frame(
    polygon_id = rep("GC1", 5),
    x = c(-0.5, 1.5, 1.5, -0.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5, -0.5)
  )
  flat_path <- tempfile(fileext = ".csv")
  utils::write.csv(flat, flat_path, row.names = FALSE)
  out <- createSIRAP(
    spatial = getSpatial(x),
    receptors = getReceptors(x),
    spatial_id = "spatial_id",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE,
    polygons_path = flat_path
  )
  expect_equal(out$polygons$polygon_id, "GC1")
  expect_equal(unname(out$spatial$polygon_id), c("GC1", "GC1", NA, NA))
  expect_equal(out$receptors$polygon_id, c("GC1", "GC1", NA))
  expect_invisible(validateSIRAP(out))
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
  expect_identical(names(x), c("spatial", "receptors", "metadata", "polygons", "analyses"))
})

test_that("createSIRAP validates initial analyses", {
  x <- make_test_sirap()
  spatial <- getSpatial(x)
  receptors <- getReceptors(x)

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

test_that("update receptor functions operate on the receptor table", {
  x <- make_test_sirap()
  x$metadata$receptor_fields["clone"] <- list(NULL)
  x$receptors$clone_id <- NULL

  n_records <- nrow(getReceptors(x))
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

test_that("update spatial functions operate on the Seurat object", {
  x <- make_test_sirap()

  expect_equal(ncol(getSpatial(x)), 4)
  x <- updateSpatial(x, function(spatial, value) {
    spatial$new_field <- value
    spatial
  }, value = 1:4)

  expect_identical(unname(x$spatial$new_field), 1:4)
  expect_error(updateSpatial(x, function(spatial) data.frame()), "must return a Seurat")
})

test_that("update analysis functions operate on one named result", {
  x <- setAnalysis(make_test_sirap(), "scores", c(3, 1, 2))

  expect_identical(sort(getAnalysis(x, "scores")), c(1, 2, 3))
  expect_identical(getAnalysis(x, "scores"), c(3, 1, 2))

  x <- updateAnalysis(x, "scores", sort, decreasing = TRUE)
  expect_identical(getAnalysis(x, "scores"), c(3, 2, 1))
  expect_error(updateAnalysis(x, "missing", identity), "is not stored")
})

test_that("createSIRAP retains unmatched nonspatial receptors", {
  x <- make_test_sirap()
  info <- summary(x)

  expect_equal(nrow(getReceptors(x)), 3)
  expect_equal(info$n_receptor_ids, 3)
  expect_equal(info$n_matched_ids, 2)
  expect_equal(info$n_matched_records, 2)
  expect_equal(info$id_match_rate, 2 / 3)
})

test_that("receptor summaries and QC plots accept SIRAP objects", {
  x <- make_test_sirap(spatial_id = "spot")

  receptor_summary <- spatialReceptorSummary(x)
  expect_equal(receptor_summary$n_sequences, 3)
  expect_equal(receptor_summary$n_spatial_ids, 3)
  expect_equal(receptor_summary$n_loci, 1)

  expect_s3_class(plotLengthQc(x), "ggplot")
})

test_that("collapseIdenticalNeighbors updates SIRAP objects", {
  skip_if_not_installed("kedd")
  skip_if_not_installed("RANN")
  skip_if_not_installed("igraph")

  set.seed(1)
  x <- make_collapse_test_sirap()
  collapsed <- collapseIdenticalNeighbors(
    x,
    sequence = "sequence_alignment",
    consensus_count = "consensus_count",
    num_neighbors_non_boundary = 2L
  )

  expect_s3_class(collapsed, "sirap")
  expect_equal(nrow(getReceptors(collapsed)), 3)
  expect_false("cell1" %in% getReceptors(collapsed)$spot)
  expect_true(all(c("cell2", "cell13", "cell25") %in% getReceptors(collapsed)$spot))
  expect_identical(collapsed$metadata$receptor_fields$spatial_id, "spot")
})

test_that("collapseIdenticalNeighbors accepts receptor tables with a Seurat object", {
  skip_if_not_installed("kedd")
  skip_if_not_installed("RANN")
  skip_if_not_installed("igraph")

  set.seed(1)
  x <- make_collapse_test_sirap()
  collapsed <- collapseIdenticalNeighbors(
    getReceptors(x),
    seurat_object = getSpatial(x),
    spatial_id = "spot",
    sequence = "sequence_alignment",
    consensus_count = "consensus_count",
    num_neighbors_non_boundary = 2L,
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_s3_class(collapsed, "data.frame")
  expect_equal(nrow(collapsed), 3)
  expect_false(isSIRAP(collapsed))
  expect_error(
    collapseIdenticalNeighbors(getReceptors(x), spatial_id = "spot"),
    "`seurat_object` is required"
  )
})

test_that("default optional receptor fields are NULL when absent", {
  skip_if_not_installed("SeuratObject")
  spatial <- SeuratObject::CreateSeuratObject(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )
  x <- createSIRAP(spatial, data.frame(spatial_id = "cell1"))

  expect_null(x$metadata$receptor_fields$clone)
  expect_null(x$metadata$receptor_fields$locus)
  expect_null(x$metadata$image_scale)
  expect_null(x$metadata$swap_coords)
  expect_error(
    createSIRAP(spatial, data.frame(spatial_id = "cell1"), clone = "missing"),
    "missing column"
  )
})

test_that("createSIRAP rejects receptor tables without spatial overlap", {
  skip_if_not_installed("SeuratObject")
  spatial <- SeuratObject::CreateSeuratObject(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )

  expect_error(
    createSIRAP(spatial, data.frame(spatial_id = "bulk1")),
    "No receptor IDs match"
  )
  expect_error(
    createSIRAP(spatial, data.frame(spatial_id = "cell1"), metadata = list("unnamed")),
    "named list"
  )
  expect_error(
    createSIRAP(spatial, data.frame(spatial_id = "cell1"), metadata = list(receptor_fields = list())),
    "reserved"
  )
  expect_error(
    createSIRAP(spatial, data.frame(spatial_id = "cell1"), metadata = list(image_scale = "hires")),
    "reserved"
  )
})

test_that("spatial options require explicit values without SIRAP metadata", {
  skip_if_not_installed("SeuratObject")
  spatial <- SeuratObject::CreateSeuratObject(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )
  x <- createSIRAP(spatial, data.frame(spatial_id = "cell1"))

  expect_error(
    SIRAP:::.sirap_spatial_options(spatial),
    "`image_scale` must be supplied"
  )
  expect_error(
    SIRAP:::.sirap_spatial_options(x),
    "`image_scale` must be supplied"
  )
  expect_identical(
    SIRAP:::.sirap_spatial_options(x, image_scale = "lowres", swap_coords = TRUE),
    list(image_scale = "lowres", swap_coords = TRUE)
  )
})

test_that("SIRAP metadata supplies receptor inputs to spatial analyses", {
  x <- make_test_sirap(spatial_id = "spot", clone = "clonal_group")

  polygon <- data.frame(
    x = c(-1, 4, 4, -1, -1),
    y = c(-1, -1, 4, 4, -1)
  )
  dominance <- spatialCloneDominance(x, polygon = polygon)
  expect_s3_class(dominance, "sirap_spatial_clone_dominance")
  expect_equal(dominance$metadata$receptor_fields$spatial_id, "spot")
  expect_equal(dominance$metadata$receptor_fields$clone, "clonal_group")
  expect_equal(sort(dominance$dominance$clone_id), c("A", "B"))
  expect_equal(sort(dominance$region_receptors$.sirap_spatial_id), c("cell1", "cell2"))
})

test_that("spatialLocalMoran uses the receptor ID stored in SIRAP metadata", {
  skip_if_not_installed("spdep")
  x <- make_test_sirap(spatial_id = "spot", clone = "clonal_group")

  local <- spatialLocalMoran(
    object = x,
    variable = "receptor_present",
    variable_type = "receptor",
    k = 1
  )

  expect_equal(local$value[match(paste0("cell", 1:4), local$cell)], c(1, 1, 0, 0))
})

test_that("spatialCloneJoinCount uses receptor fields stored in SIRAP metadata", {
  skip_if_not_installed("spdep")
  spatial <- SeuratObject::CreateSeuratObject(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 4,
      dimnames = list("gene", paste0("cell", 1:4))
    )
  )
  spatial$x <- c(0, 1, 2, 3)
  spatial$y <- c(0, 0, 0, 0)
  receptors <- data.frame(
    spot = c("cell1", "cell2", "cell3"),
    clonal_group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  x <- createSIRAP(
    spatial = spatial,
    receptors = receptors,
    spatial_id = "spot",
    clone = "clonal_group",
    image_scale = "hires",
    swap_coords = FALSE
  )

  out <- spatialCloneJoinCount(
    object = x,
    clone_val = "A",
    k = 1,
    min_spatial_ids = 2,
    return_data = TRUE
  )

  expect_true("A" %in% out$summary$clone_id)
  expect_equal(out$summary$n_clone_spatial_ids[out$summary$clone_id == "A"], 2)
  expect_true(all(c("join_count_mark", "clone_id", ".sirap_spatial_id") %in% names(out$data)))
  expect_false("style" %in% names(out))
  expect_false("style" %in% names(out$summary))
  expect_equal(sort(unique(out$data$.sirap_spatial_id)), c("cell1", "cell2", "cell3"))
})

test_that("SIRAP print reports spatial and receptor dimensions", {
  x <- make_test_sirap()

  expect_output(
    print(x),
    paste0(
      "<SIRAP object>\\n",
      "Spatial data:    Seurat, 4 spatial IDs, 1 gene, 0 images\\n",
      "Receptor data:   tibble, 3 receptor sequences x 8 columns\\n",
      "Spatial overlap: 2/3 receptor IDs in spatial data \\(66.7%\\)\\n",
      "Spatial ID field: spatial_id\\n",
      "Clone field:     clone_id\\n",
      "Analyses:        0"
    )
  )
})

test_that("SIRAP summaries align values and report stored polygons", {
  x <- make_test_sirap()
  polygon_table <- tibble::tibble(
    polygon_id = c("GC1", "GC2"),
    vertices = list(
      tibble::tibble(x = c(0, 1, 0, 0), y = c(0, 0, 1, 0)),
      tibble::tibble(x = c(1, 2, 1, 1), y = c(1, 1, 2, 1))
    )
  )
  x <- addPolygons(x, polygons = polygon_table)

  expect_equal(summary(x)$n_polygons, 2)
  expect_output(
    print(x),
    "Analyses:        0\\nPolygons:        2"
  )
  expect_output(
    print(summary(x)),
    "(?s)  Spatial IDs:\\s+4.*  Receptor records:\\s+3.*  Polygons:\\s+2",
    perl = TRUE
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
