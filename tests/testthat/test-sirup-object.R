make_test_sirup <- function(spatial_id = "spatial_id", clone = "clone_id") {
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
  createSIRUP(
    spatial = spatial,
    receptors = receptors,
    spatial_id = spatial_id,
    clone = clone,
    image_scale = "hires",
    swap_coords = FALSE
  )
}

make_collapse_test_sirup <- function() {
  coords <- expand.grid(x = seq_len(5), y = seq_len(5))
  cells <- paste0("cell", seq_len(nrow(coords)))
  spatial <- sirup_test_seurat(
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
  createSIRUP(
    spatial = spatial,
    receptors = receptors,
    spatial_id = "spot",
    locus = "locus",
    count = "consensus_count",
    image_scale = "hires",
    swap_coords = FALSE
  )
}

test_that("createSIRUP creates a validated tibble-backed container", {
  x <- make_test_sirup()

  expect_s3_class(x, "sirup")
  expect_true(isSIRUP(x))
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
  expect_identical(SIRUP:::.sirup_spatial_input(x), x$spatial)
  expect_identical(SIRUP:::.sirup_receptor_input(x)$receptors, x$receptors)
  expect_invisible(validateSIRUP(x))
})

test_that("createSIRUP supports named spatial lists keyed by receptor sample", {
  make_spatial <- function(prefix) {
    spatial <- sirup_test_seurat(
      counts = matrix(
        1,
        nrow = 1,
        ncol = 2,
        dimnames = list("gene", c("cell1", "cell2"))
      )
    )
    spatial$x <- c(0, 1)
    spatial$y <- c(0, 1)
    spatial$section <- prefix
    spatial
  }
  spatial <- list(
    section_a = make_spatial("section_a"),
    section_b = make_spatial("section_b")
  )
  receptors <- data.frame(
    barcode = c("cell1", "cell1", "bulk1"),
    sample = c("section_a", "section_b", "blood"),
    clone_id = c("A", "A", "A"),
    locus = "IGH",
    stringsAsFactors = FALSE
  )

  x <- createSIRUP(
    spatial = spatial,
    receptors = receptors,
    spatial_id = "barcode",
    sample_id = "sample",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_s3_class(x, "sirup")
  expect_equal(spatialSamples(x), c("section_a", "section_b"))
  expect_equal(x$metadata$sample_names, c("section_a", "section_b"))
  section_a <- spatial$section_a
  section_a$sample <- "section_a"
  expect_identical(getSpatial(x, sample_id = "section_a"), section_a)
  expect_identical(getReceptors(x), x$receptors)
  expect_equal(summary(x)$n_spatial_ids, 4)
  expect_equal(summary(x)$n_receptor_ids, 3)
  expect_equal(summary(x)$n_matched_ids, 2)
  expect_equal(summary(x)$n_matched_records, 2)

  coords <- spatialCoordinateTable(x)
  expect_true("sample" %in% names(coords))
  expect_equal(sort(unique(coords$sample)), c("section_a", "section_b"))
  expect_equal(sum(coords$cell == "cell1"), 2)
})

test_that("mergeSIRUP combines single-sample objects into a multisample object", {

  make_object <- function(cells, clones) {
    spatial <- sirup_test_seurat(
      counts = matrix(1, nrow = 1, ncol = length(cells), dimnames = list("gene", cells))
    )
    spatial$x <- seq_along(cells) - 1
    spatial$y <- seq_along(cells) - 1
    createSIRUP(
      spatial = spatial,
      receptors = data.frame(
        barcode = cells,
        clone_id = clones,
        locus = "IGH",
        stringsAsFactors = FALSE
      ),
      spatial_id = "barcode",
      clone = "clone_id",
      locus = "locus",
      image_scale = "hires",
      swap_coords = FALSE
    )
  }
  obj1 <- setAnalysis(make_object(c("cell1", "cell2"), c("A", "B")), "score", data.frame(x = 1))
  obj2 <- setAnalysis(make_object(c("cell1", "cell3"), c("C", "D")), "score", data.frame(x = 2))

  merged <- mergeSIRUP(obj1, obj2, samples = c("samp1", "samp2"))

  expect_s3_class(merged, "sirup")
  expect_equal(spatialSamples(merged), c("samp1", "samp2"))
  expect_equal(merged$metadata$sample_names, c("samp1", "samp2"))
  expect_equal(getReceptors(merged)$sample_id, c("samp1", "samp1", "samp2", "samp2"))
  expect_equal(getReceptors(merged, sample_id = "samp2")$clone_id, c("C", "D"))
  expect_identical(getAnalysis(merged, "score", sample_id = "samp1"), data.frame(x = 1))
  expect_identical(getAnalysis(merged, "score", sample_id = "samp2"), data.frame(x = 2))
  expect_equal(merged$metadata$receptor_fields$sample, "sample_id")
  expect_invisible(validateSIRUP(merged))
})

test_that("mergeSIRUP validates existing sample fields", {

  spatial <- sirup_test_seurat(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )
  spatial$x <- 0
  spatial$y <- 0
  make_object <- function(sample_id) {
    createSIRUP(
      spatial = spatial,
      receptors = data.frame(
        barcode = "cell1",
        sample_id = sample_id,
        clone_id = "A",
        stringsAsFactors = FALSE
      ),
      spatial_id = "barcode",
      clone = "clone_id",
      image_scale = "hires",
      swap_coords = FALSE
    )
  }

  expect_error(
    mergeSIRUP(make_object("wrong"), make_object("samp2"), samples = c("samp1", "samp2")),
    "do not match the requested sample"
  )
  merged <- mergeSIRUP(make_object("samp1"), make_object("samp2"), samples = c("samp1", "samp2"))
  expect_equal(getReceptors(merged)$sample_id, c("samp1", "samp2"))
})

test_that("mergeSIRUP preserves sample-tagged polygons when inputs differ", {

  make_object <- function(cell) {
    spatial <- sirup_test_seurat(
      counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", cell))
    )
    spatial$x <- 0
    spatial$y <- 0
    createSIRUP(
      spatial = spatial,
      receptors = data.frame(barcode = cell, clone_id = "A", stringsAsFactors = FALSE),
      spatial_id = "barcode",
      clone = "clone_id",
      image_scale = "hires",
      swap_coords = FALSE
    )
  }
  obj1 <- addPolygons(
    make_object("cell1"),
    polygons = data.frame(polygon = "GC1", x = c(-1, 1, 0, -1), y = c(-1, -1, 1, -1))
  )
  obj2 <- make_object("cell1")

  merged <- mergeSIRUP(obj1, obj2, samples = c("samp1", "samp2"))

  expect_equal(getPolygons(merged)$sample, "samp1")
  expect_equal(getPolygons(merged)$polygon_id, "samp1::GC1")
  expect_equal(getPolygons(merged, sample_id = "samp1")$polygon_id, "GC1")
  expect_equal(getSpatial(merged, "samp1")[[]]["cell1", "polygon_id"], "samp1::GC1")
  expect_true(is.na(getSpatial(merged, "samp2")[[]]["cell1", "polygon_id"]))
  expect_equal(getReceptors(merged)$polygon_id, c("samp1::GC1", NA))
})

test_that("spatial list coordinate options resolve per sample", {

  make_spatial <- function(x, y) {
    spatial <- sirup_test_seurat(
      counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
    )
    spatial$x <- x
    spatial$y <- y
    spatial
  }
  x <- createSIRUP(
    spatial = list(
      section_a = make_spatial(x = 1, y = 10),
      section_b = make_spatial(x = 2, y = 20)
    ),
    receptors = data.frame(
      barcode = c("cell1", "cell1"),
      sample = c("section_a", "section_b"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample_id = "sample",
    image_scale = list(section_a = "hires", section_b = "hires"),
    swap_coords = c(section_a = FALSE, section_b = TRUE)
  )

  coords <- spatialCoordinateTable(x)
  coords <- coords[order(coords$sample), ]

  expect_equal(coords$sample, c("section_a", "section_b"))
  expect_equal(coords$x, c(1, 20))
  expect_equal(coords$y, c(10, 2))
  expect_equal(coords$.sirup_y, c(-10, -2))
})

test_that("spatial list SIRUP objects require a sample receptor field", {
  spatial <- sirup_test_seurat(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )
  receptors <- data.frame(barcode = "cell1", stringsAsFactors = FALSE)

  expect_error(
    createSIRUP(
      spatial = list(section_a = spatial),
      receptors = receptors,
      spatial_id = "barcode",
      image_scale = "hires",
      swap_coords = FALSE
    ),
    "`sample_id` is required"
  )
})

test_that("stored polygons annotate spatial list SIRUP objects by sample", {
  make_spatial <- function() {
    spatial <- sirup_test_seurat(
      counts = matrix(
        1,
        nrow = 1,
        ncol = 2,
        dimnames = list("gene", c("cell1", "cell2"))
      )
    )
    spatial$x <- c(0, 10)
    spatial$y <- c(0, 10)
    spatial
  }
  x <- createSIRUP(
    spatial = list(section_a = make_spatial(), section_b = make_spatial()),
    receptors = data.frame(
      barcode = c("cell1", "cell1", "bulk1"),
      sample = c("section_a", "section_b", "blood"),
      clone_id = c("A", "B", "C"),
      locus = "IGH",
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample_id = "sample",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
  polygons <- data.frame(
    polygon = rep(c("a_gc", "b_gc"), each = 4),
    sample = rep(c("section_a", "section_b"), each = 4),
    x = rep(c(-1, 1, 0, -1), 2),
    y = rep(c(-1, -1, 1, -1), 2),
    stringsAsFactors = FALSE
  )

  x <- addPolygons(x, polygons = polygons)

  expect_equal(getSpatial(x, "section_a")[[]]["cell1", "polygon_id"], "a_gc")
  expect_equal(getSpatial(x, "section_b")[[]]["cell1", "polygon_id"], "b_gc")
  expect_equal(x$receptors$polygon_id, c("a_gc", "b_gc", NA))
  expect_true(x$receptors$polygon_a_gc[[1]])
  expect_false(x$receptors$polygon_b_gc[[1]])
  expect_true(x$receptors$polygon_b_gc[[2]])
  dominance <- spatialCloneDominance(x, locus_val = "IGH")
  expect_equal(sort(unique(dominance$region_receptors$polygon)), c("a_gc", "b_gc"))
  expect_equal(sort(dominance$dominance$clone_id), c("A", "B"))
  expect_invisible(validateSIRUP(x))
})

test_that("direct polygon dominance honors polygon sample columns", {

  make_spatial <- function() {
    spatial <- sirup_test_seurat(
      counts = matrix(
        1,
        nrow = 1,
        ncol = 2,
        dimnames = list("gene", c("cell1", "cell2"))
      )
    )
    spatial$x <- c(0, 10)
    spatial$y <- c(0, 10)
    spatial
  }
  x <- createSIRUP(
    spatial = list(section_a = make_spatial(), section_b = make_spatial()),
    receptors = data.frame(
      barcode = c("cell1", "cell1"),
      sample = c("section_a", "section_b"),
      clone_id = c("A", "B"),
      locus = "IGH",
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample_id = "sample",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
  polygons <- data.frame(
    polygon = rep(c("a_gc", "b_gc"), each = 4),
    sample = rep(c("section_a", "section_b"), each = 4),
    x = rep(c(-1, 1, 0, -1), 2),
    y = rep(c(-1, -1, 1, -1), 2),
    stringsAsFactors = FALSE
  )

  dominance <- spatialCloneDominance(x, polygon = polygons)

  expect_equal(
    dominance$region_receptors[, c("sample", "clone_id", "polygon")],
    data.frame(
      sample = c("section_a", "section_b"),
      clone_id = c("A", "B"),
      polygon = c("a_gc", "b_gc"),
      stringsAsFactors = FALSE
    )
  )
  expect_equal(sort(unique(dominance$polygons$sample)), c("section_a", "section_b"))
})

test_that("addPolygons can tag multisample polygons from named paths", {

  make_spatial <- function() {
    spatial <- sirup_test_seurat(
      counts = matrix(1, nrow = 1, ncol = 2, dimnames = list("gene", c("cell1", "cell2")))
    )
    spatial$x <- c(0, 10)
    spatial$y <- c(0, 10)
    spatial
  }
  x <- createSIRUP(
    spatial = list(section_a = make_spatial(), section_b = make_spatial()),
    receptors = data.frame(
      barcode = c("cell1", "cell1"),
      sample = c("section_a", "section_b"),
      clone_id = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample_id = "sample",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
  poly_a <- data.frame(
    polygon = rep(c("gc", "empty"), each = 4),
    x = c(-1, 1, 0, -1, 100, 101, 100.5, 100),
    y = c(-1, -1, 1, -1, 100, 100, 101, 100),
    stringsAsFactors = FALSE
  )
  poly_b <- data.frame(polygon = "gc", x = c(-1, 1, 0, -1), y = c(-1, -1, 1, -1))
  path_a <- tempfile(fileext = ".csv")
  path_b <- tempfile(fileext = ".csv")
  utils::write.csv(poly_a, path_a, row.names = FALSE)
  utils::write.csv(poly_b, path_b, row.names = FALSE)

  out <- addPolygons(x, polygons_path = c(section_a = path_a, section_b = path_b))

  expect_equal(getSpatial(out, "section_a")[[]]["cell1", "polygon_id"], "section_a::gc")
  expect_equal(getSpatial(out, "section_b")[[]]["cell1", "polygon_id"], "section_b::gc")
  expect_equal(out$receptors$polygon_id, c("section_a::gc", "section_b::gc"))
  one <- filterSample(out, "section_a")
  expect_equal(sort(getPolygons(one)$polygon_id), c("empty", "gc"))
  expect_equal(getSpatial(one)[[]]["cell1", "polygon_id"], "gc")
  expect_equal(getReceptors(one)$polygon_id, "gc")
  dominance <- spatialCloneDominance(one)
  expect_equal(unique(dominance$region_spatial_ids$polygon), "gc")
  expect_equal(unique(dominance$region_receptors$polygon), "gc")
  expect_error(addPolygons(x, polygons = poly_a), "must include `sample` or `spatial_sample`")
})

test_that("SIRUP objects store canonical polygon geometry", {
  x <- make_test_sirup()
  polygon_table <- tibble::tibble(
    polygon_id = c("GC1", "GC2"),
    vertices = list(
      tibble::tibble(spatial_id = c("cell1", "cell2", "cell3", "cell1"), x = c(0, 1, 2, 0), y = c(0, 1, 0, 0)),
      tibble::tibble(spatial_id = c("cell2", "cell3", "cell4", "cell2"), x = c(1, 2, 3, 1), y = c(1, 0, 1, 1))
    )
  )

  with_polygons <- addPolygons(x, polygons = polygon_table)
  expect_equal(getPolygons(with_polygons), polygon_table)
  expect_invisible(validateSIRUP(with_polygons))

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

test_that("addPolygons annotates SIRUP spatial metadata and receptors", {
  x <- make_test_sirup()
  polygon <- data.frame(
    polygon = "GC1",
    x = c(-0.5, 1.5, 1.5, -0.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5, -0.5)
  )

  out <- addPolygons(x, polygons = polygon)

  expect_equal(unname(out$spatial$polygon_id), c("GC1", "GC1", NA, NA))
  expect_equal(out$receptors$polygon_id, c("GC1", "GC1", NA))
  expect_equal(unname(out$spatial$polygon_GC1), c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(out$receptors$polygon_GC1, c(TRUE, TRUE, FALSE))
  expect_invisible(validateSIRUP(out))
})

test_that("addPolygons can exclude boundary points from stored annotations", {
  x <- make_test_sirup()
  polygon <- data.frame(
    polygon = "GC1",
    x = c(0, 2, 2, 0),
    y = c(0, 0, 2, 2)
  )

  out <- addPolygons(x, polygons = polygon, include_boundary = FALSE)

  expect_equal(unname(out$spatial$polygon_id), c(NA, "GC1", NA, NA))
  expect_equal(out$receptors$polygon_id, c(NA, "GC1", NA))
  expect_false(out$metadata$polygons_include_boundary)
})

test_that("polygonReceptorSummary stores one summary row per SIRUP polygon", {
  x <- make_test_sirup()
  x$receptors$mu_freq <- c(0.1, 0.3, 0.5)
  x$receptors$productive <- c(TRUE, FALSE, TRUE)
  polygon <- data.frame(
    polygon = "GC1",
    x = c(-0.5, 1.5, 1.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5)
  )
  x <- addPolygons(x, polygons = polygon)

  out <- polygonReceptorSummary(
    x,
    summary_columns = c("mu_freq", "clone_id", "productive", "productive"),
    summary_operations = c("mean", "cat", "any", "all")
  )
  summary <- getAnalysis(out, "polygon_receptor_summary")

  expect_s3_class(out, "sirup")
  expect_equal(summary$polygon_id, "GC1")
  expect_equal(summary$n_receptors, 2)
  expect_equal(summary$n_spatial_ids, 2)
  expect_equal(summary$n_clones, 2)
  expect_equal(summary$top_clone_fraction, 0.5)
  expect_equal(summary$mean_mu_freq, 0.2)
  expect_equal(summary$cat_clone_id, "A,B")
  expect_true(summary$any_productive)
  expect_false(summary$all_productive)
})

test_that("polygonReceptorSummary can weight receptor summaries", {
  x <- make_test_sirup()
  x$receptors$mu_freq <- c(0.1, 0.3, 0.5)
  x$receptors$consensus_count <- c(10, 1, 100)
  polygon <- data.frame(
    polygon = "GC1",
    x = c(-0.5, 1.5, 1.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5)
  )
  x <- addPolygons(x, polygons = polygon)

  out <- NULL
  expect_warning(
    out <- polygonReceptorSummary(
      x,
      weight = "consensus_count",
      summary_columns = "mu_freq",
      summary_operations = "mean"
    ),
    "Weighted polygon summaries"
  )
  summary <- getAnalysis(out, "polygon_receptor_summary")

  expect_equal(summary$n_receptors, 11)
  expect_equal(summary$top_clone_fraction, 10 / 11)
  expect_equal(summary$mean_mu_freq, stats::weighted.mean(c(0.1, 0.3), c(10, 1)))
})

test_that("polygon analyses use stored SIRUP annotations when polygons are attached", {
  x <- make_test_sirup()
  polygon <- data.frame(
    polygon = "GC1",
    x = c(0, 2, 2, 0),
    y = c(0, 0, 2, 2)
  )
  x <- addPolygons(x, polygons = polygon, include_boundary = FALSE)

  dominance <- spatialCloneDominance(x)
  expect_equal(dominance$dominance$clone_id, "B")
  expect_equal(dominance$region_spatial_ids$spatial_id, "cell2")

  join_data <- SIRUP:::.sirup_joincount_receptor_data_stored_polygons(
    object = x,
    receptors = getReceptors(x),
    clone_col = "clone_id",
    spatial_id_col = "spatial_id",
    polygon_id = NULL,
    images = NULL,
    image_scale = "hires",
    swap_coords = FALSE
  )
  expect_equal(unique(join_data$.sirup_spatial_id), "cell2")
})

test_that("withPolygons applies functions by polygon and across all polygons", {
  x <- make_test_sirup()
  polygon <- data.frame(
    polygon = "GC1",
    x = c(-0.5, 1.5, 1.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5)
  )
  x <- addPolygons(x, polygons = polygon)

  by_polygon <- withPolygons(x, function(dat, polygon_id) {
    data.frame(polygon_id = polygon_id, n = nrow(dat))
  }, pass_polygon_id = TRUE)
  all_polygons <- withPolygons(x, nrow, mode = "all_polygons")

  expect_equal(names(by_polygon), "GC1")
  expect_equal(by_polygon$GC1$n, 2)
  expect_equal(all_polygons, 2)
})

test_that("formatSpatialClones asks for germlines before formatting", {
  x <- make_test_sirup()

  expect_error(
    formatSpatialClones(x),
    "createGermlines"
  )
})

test_that("polygonAbundance and polygonDiversity give an error below threshold", {
  x <- make_test_sirup()
  polygon <- data.frame(
    polygon = "GC1",
    x = c(-0.5, 1.5, 1.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5)
  )
  x <- addPolygons(x, polygons = polygon)

  expect_error(polygonAbundance(x, min_n = 10), "No polygon has at least 10")
  expect_error(polygonDiversity(x, min_n = 10), "No polygon has at least 10")
})

test_that("polygonAbundance runs above threshold with bulk default cell_id", {
  x <- make_test_sirup()
  polygon <- data.frame(
    polygon = "GC1",
    x = c(-0.5, 1.5, 1.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5)
  )
  x <- addPolygons(x, polygons = polygon)

  abundance <- polygonAbundance(x, min_n = 1, nboot = 2)

  expect_s3_class(abundance, "sirup")
  expect_s4_class(getAnalysis(abundance, "polygon_abundance"), "AbundanceCurve")
})

make_polygon_depth_sirup <- function() {
  # Three polygons with 14, 9 and 3 receptors.
  cells <- paste0("cell", 1:9)
  spatial <- sirup_test_seurat(
    counts = matrix(1, nrow = 1, ncol = 9, dimnames = list("gene", cells))
  )
  spatial$x <- c(1, 1, 1, 11, 11, 11, 21, 21, 21)
  spatial$y <- spatial$x
  receptors <- data.frame(
    spatial_id = c(
      rep(cells[1:3], length.out = 14),
      rep(cells[4:6], length.out = 9),
      rep(cells[7:9], length.out = 3)
    ),
    clone_id = c(
      rep(c("A", "B", "C", "D", "E"), length.out = 14),
      rep(c("F", "G", "H"), length.out = 9),
      c("I", "J", "I")
    ),
    locus = "IGH",
    stringsAsFactors = FALSE
  )
  x <- createSIRUP(
    spatial = spatial,
    receptors = receptors,
    spatial_id = "spatial_id",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
  square <- function(id, centre) {
    data.frame(
      polygon = id,
      x = centre + c(-1, 1, 1, -1),
      y = centre + c(-1, -1, 1, 1)
    )
  }
  addPolygons(x, polygons = rbind(square("GC1", 1), square("GC2", 11), square("GC3", 21)))
}

test_that("polygonDiversity and polygonAbundance use min_n = 10 by default", {
  expect_equal(formals(polygonDiversity)$min_n, 10)
  expect_equal(formals(polygonAbundance)$min_n, 10)
})

test_that("polygonDiversity makes one object for all polygons", {
  x <- make_polygon_depth_sirup()
  expect_warning(
    x <- polygonDiversity(x, min_n = 5, max_q = 1, step_q = 1, nboot = 5),
    "Excluded: GC3"
  )
  diversity <- getAnalysis(x, "polygon_diversity")
  expect_s4_class(diversity, "DiversityCurve")
  expect_equal(diversity@group_by, "polygon_id")
  # GC3 has 3 receptors and is left out.
  expect_setequal(unique(as.character(diversity@diversity$polygon_id)), c("GC1", "GC2"))
})

test_that("min_n sets which polygons stay", {
  x <- make_polygon_depth_sirup()
  x <- polygonAbundance(x, min_n = 3, nboot = 5)
  abundance <- getAnalysis(x, "polygon_abundance")
  expect_s4_class(abundance, "AbundanceCurve")
  expect_setequal(unique(as.character(abundance@abundance$polygon_id)), c("GC1", "GC2", "GC3"))
  expect_warning(x2 <- polygonAbundance(x, min_n = 5, nboot = 5), "Excluded: GC3")
  expect_setequal(
    unique(as.character(getAnalysis(x2, "polygon_abundance")@abundance$polygon_id)),
    c("GC1", "GC2")
  )
})

test_that("polygonDiversity does not accept group", {
  x <- make_polygon_depth_sirup()
  expect_error(polygonDiversity(x, min_n = 5, group = "locus"), "`group` is not an argument")
})

test_that("common_size = FALSE analyses each polygon alone and returns a list", {
  x <- make_polygon_depth_sirup()
  x <- polygonDiversity(x, common_size = FALSE, min_n = 5, max_q = 2, step_q = 1, nboot = 5)
  diversity <- getAnalysis(x, "polygon_diversity")
  expect_type(diversity, "list")
  expect_named(diversity, c("GC1", "GC2", "GC3"))
  expect_s4_class(diversity$GC1, "DiversityCurve")
  expect_s4_class(diversity$GC2, "DiversityCurve")
  expect_null(diversity$GC3)
  x <- polygonAbundance(x, common_size = FALSE, min_n = 5, nboot = 5)
  abundance <- getAnalysis(x, "polygon_abundance")
  expect_type(abundance, "list")
  expect_s4_class(abundance$GC1, "AbundanceCurve")
  expect_null(abundance$GC3)
  # no polygon passes: each polygon gives NULL and there is no error
  y <- polygonAbundance(make_polygon_depth_sirup(), common_size = FALSE, min_n = 100, nboot = 5)
  expect_true(all(vapply(getAnalysis(y, "polygon_abundance"), is.null, logical(1))))
})

test_that("summarizePolygonDiversity reads the list from common_size = FALSE", {
  x <- make_polygon_depth_sirup()
  x <- polygonDiversity(x, common_size = FALSE, min_n = 5, max_q = 2, step_q = 1, nboot = 5)
  x <- summarizePolygonDiversity(x, q = 1, analysis_name = "shannon")
  tab <- getAnalysis(x, "shannon")
  expect_equal(tab$polygon_id, c("GC1", "GC2", "GC3"))
  expect_true(is.na(tab$d[tab$polygon_id == "GC3"]))
  expect_false(anyNA(tab$d[tab$polygon_id != "GC3"]))
})

test_that("summarizePolygonDiversity keeps polygons that were left out", {
  x <- make_polygon_depth_sirup()
  expect_warning(
    x <- polygonDiversity(x, min_n = 5, max_q = 2, step_q = 1, nboot = 5),
    "Excluded: GC3"
  )
  x <- summarizePolygonDiversity(x, q = 1, analysis_name = "shannon")
  tab <- getAnalysis(x, "shannon")
  expect_setequal(tab$polygon_id, c("GC1", "GC2", "GC3"))
  expect_true(is.na(tab$d[tab$polygon_id == "GC3"]))
  expect_false(anyNA(tab$d[tab$polygon_id != "GC3"]))
})

test_that("summarizePolygonDiversity stores a q-specific diversity table", {
  x <- make_test_sirup()
  polygon <- data.frame(
    polygon = "GC1",
    x = c(-0.5, 1.5, 1.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5)
  )
  x <- addPolygons(x, polygons = polygon)
  x <- polygonDiversity(
    x,
    min_n = 1,
    max_q = 2,
    step_q = 1,
    nboot = 2
  )
  x <- summarizePolygonDiversity(x, q = 2)
  summary <- getAnalysis(x, "polygon_diversity_q2")

  expect_equal(summary$polygon_id, "GC1")
  expect_equal(summary$q, 2)
  expect_true("d" %in% names(summary))
})

test_that("createSIRUP accepts a polygon path", {
  x <- make_test_sirup()
  flat <- data.frame(
    polygon_id = rep("GC1", 5),
    x = c(-0.5, 1.5, 1.5, -0.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5, -0.5)
  )
  flat_path <- tempfile(fileext = ".csv")
  utils::write.csv(flat, flat_path, row.names = FALSE)
  out <- createSIRUP(
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
  expect_equal(unname(out$spatial$polygon_GC1), c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(out$receptors$polygon_GC1, c(TRUE, TRUE, FALSE))
  expect_invisible(validateSIRUP(out))
})

test_that("analysis accessors store arbitrary named results", {
  x <- make_test_sirup()
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

test_that("legacy SIRUP objects behave as though analyses are empty", {
  x <- make_test_sirup()
  x$analyses <- NULL

  expect_invisible(validateSIRUP(x))
  expect_identical(listAnalyses(x), character())

  x <- setAnalysis(x, "result", 1)
  expect_identical(getAnalysis(x, "result"), 1)
  expect_identical(names(x), c("spatial", "receptors", "metadata", "polygons", "analyses"))
})

test_that("createSIRUP validates initial analyses", {
  x <- make_test_sirup()
  spatial <- getSpatial(x)
  receptors <- getReceptors(x)

  with_result <- createSIRUP(
    spatial,
    receptors,
    analyses = list(note = "already analyzed")
  )
  expect_identical(getAnalysis(with_result, "note"), "already analyzed")
  expect_error(
    createSIRUP(spatial, receptors, analyses = list(1)),
    "unique, non-empty names"
  )
  expect_error(
    createSIRUP(spatial, receptors, analyses = "not a list"),
    "named list"
  )
})

test_that("update receptor functions operate on the receptor table", {
  x <- make_test_sirup()
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
  x <- make_test_sirup()

  expect_equal(ncol(getSpatial(x)), 4)
  x <- updateSpatial(x, function(spatial, value) {
    spatial$new_field <- value
    spatial
  }, value = 1:4)

  expect_identical(unname(x$spatial$new_field), 1:4)
  expect_error(updateSpatial(x, function(spatial) data.frame()), "must return a Seurat")
})

test_that("accessor and update methods can target one spatial sample", {

  make_section <- function() {
    spatial <- sirup_test_seurat(
      counts = matrix(1, nrow = 1, ncol = 2, dimnames = list("gene", c("cell1", "cell2")))
    )
    spatial$x <- c(0, 1)
    spatial$y <- c(0, 1)
    spatial
  }
  x <- createSIRUP(
    spatial = list(section_a = make_section(), section_b = make_section()),
    receptors = data.frame(
      barcode = c("cell1", "cell1", "cell2"),
      sample = c("section_a", "section_b", "section_b"),
      clone_id = c("A", "B", "C"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample_id = "sample",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_equal(getReceptors(x, "section_a")$clone_id, "A")
  expect_equal(getReceptors(x, "section_b")$clone_id, c("B", "C"))

  x <- updateReceptors(x, function(receptors) {
    receptors$sample_only <- "yes"
    receptors$clone_id <- paste0(receptors$clone_id, "_updated")
    receptors
  }, sample_id = "section_b")
  expect_true(all(is.na(getReceptors(x, "section_a")$sample_only)))
  expect_true(all(getReceptors(x, "section_b")$sample_only == "yes"))
  expect_equal(getReceptors(x, "section_b")$clone_id, c("B_updated", "C_updated"))

  x <- updateSpatial(x, function(spatial) {
    spatial$updated <- TRUE
    spatial
  }, sample_id = "section_a")
  expect_true(all(getSpatial(x, "section_a")$updated))
  expect_false("updated" %in% colnames(getSpatial(x, "section_b")[[]]))

  x <- updateSpatial(x, function(spatial) {
    spatial$all_updated <- TRUE
    spatial
  })
  expect_true(all(getSpatial(x, "section_a")$all_updated))
  expect_true(all(getSpatial(x, "section_b")$all_updated))
})

test_that("update analysis functions operate on one named result", {
  x <- setAnalysis(make_test_sirup(), "scores", c(3, 1, 2))

  expect_identical(sort(getAnalysis(x, "scores")), c(1, 2, 3))
  expect_identical(getAnalysis(x, "scores"), c(3, 1, 2))

  x <- updateAnalysis(x, "scores", sort, decreasing = TRUE)
  expect_identical(getAnalysis(x, "scores"), c(3, 2, 1))
  expect_error(updateAnalysis(x, "missing", identity), "is not stored")
})

test_that("createSIRUP retains unmatched nonspatial receptors", {
  x <- make_test_sirup()
  info <- summary(x)

  expect_equal(nrow(getReceptors(x)), 3)
  expect_equal(info$n_receptor_ids, 3)
  expect_equal(info$n_matched_ids, 2)
  expect_equal(info$n_matched_records, 2)
  expect_equal(info$id_match_rate, 2 / 3)
})

test_that("receptor summaries and QC plots accept SIRUP objects", {
  x <- make_test_sirup(spatial_id = "spot")

  receptor_summary <- spatialReceptorSummary(x)
  expect_equal(receptor_summary$n_sequences, 3)
  expect_equal(receptor_summary$n_spatial_ids, 3)
  expect_equal(receptor_summary$n_loci, 1)

  expect_s3_class(plotLengthQc(x), "ggplot")
})

test_that("spatialReceptorSummary counts sample and spatial ID pairs", {

  make_section <- function() {
    spatial <- sirup_test_seurat(
      counts = matrix(1, nrow = 1, ncol = 2, dimnames = list("gene", c("cell1", "cell2")))
    )
    spatial$x <- c(0, 1)
    spatial$y <- c(0, 1)
    spatial
  }
  x <- createSIRUP(
    spatial = list(section_a = make_section(), section_b = make_section()),
    receptors = data.frame(
      barcode = c("cell1", "cell1", "cell2"),
      sample = c("section_a", "section_b", "section_b"),
      clone_id = c("A", "B", "C"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample_id = "sample",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )

  receptor_summary <- spatialReceptorSummary(x)
  expect_equal(receptor_summary$n_sequences, 3)
  expect_equal(receptor_summary$n_spatial_ids, 3)
})

test_that("collapseIdenticalNeighbors updates SIRUP objects", {

  set.seed(1)
  x <- make_collapse_test_sirup()
  x$receptors <- dplyr::bind_rows(
    getReceptors(x),
    data.frame(
      spot = "bulk1",
      sequence_alignment = "GGGG",
      consensus_count = 10L,
      locus = "IGH",
      stringsAsFactors = FALSE
    )
  )
  collapsed <- collapseIdenticalNeighbors(
    x,
    sequence = "sequence_alignment",
    consensus_count = "consensus_count",
    num_neighbors_non_boundary = 2L
  )

  expect_s3_class(collapsed, "sirup")
  expect_equal(nrow(getReceptors(collapsed)), 4)
  expect_false("cell1" %in% getReceptors(collapsed)$spot)
  expect_true(all(c("cell2", "cell13", "cell25", "bulk1") %in% getReceptors(collapsed)$spot))
  expect_identical(collapsed$metadata$receptor_fields$spatial_id, "spot")
})

test_that("collapseIdenticalNeighbors accepts receptor tables with a Seurat object", {

  set.seed(1)
  x <- make_collapse_test_sirup()
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
  expect_false(isSIRUP(collapsed))
  expect_error(
    collapseIdenticalNeighbors(getReceptors(x), spatial_id = "spot"),
    "`seurat_object` is required"
  )
})

test_that("default optional receptor fields are NULL when absent", {
  spatial <- sirup_test_seurat(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )
  x <- createSIRUP(spatial, data.frame(spatial_id = "cell1"))

  expect_null(x$metadata$receptor_fields$clone)
  expect_null(x$metadata$receptor_fields$locus)
  expect_null(x$metadata$image_scale)
  expect_null(x$metadata$swap_coords)
  expect_error(
    createSIRUP(spatial, data.frame(spatial_id = "cell1"), clone = "missing"),
    "missing column"
  )
})

test_that("createSIRUP rejects receptor tables without spatial overlap", {
  spatial <- sirup_test_seurat(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )

  expect_error(
    createSIRUP(spatial, data.frame(spatial_id = "bulk1")),
    "No receptor IDs match"
  )
  expect_error(
    createSIRUP(spatial, data.frame(spatial_id = "cell1"), metadata = list("unnamed")),
    "named list"
  )
  expect_error(
    createSIRUP(spatial, data.frame(spatial_id = "cell1"), metadata = list(receptor_fields = list())),
    "reserved"
  )
  expect_error(
    createSIRUP(spatial, data.frame(spatial_id = "cell1"), metadata = list(image_scale = "hires")),
    "reserved"
  )
})

test_that("spatial options require explicit values without SIRUP metadata", {
  spatial <- sirup_test_seurat(
    counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("gene", "cell1"))
  )
  x <- createSIRUP(spatial, data.frame(spatial_id = "cell1"))

  expect_error(
    SIRUP:::.sirup_spatial_options(spatial),
    "`image_scale` must be supplied"
  )
  expect_error(
    SIRUP:::.sirup_spatial_options(x),
    "`image_scale` must be supplied"
  )
  expect_identical(
    SIRUP:::.sirup_spatial_options(x, image_scale = "lowres", swap_coords = TRUE),
    list(image_scale = "lowres", swap_coords = TRUE)
  )
})

test_that("SIRUP metadata supplies receptor inputs to spatial analyses", {
  x <- make_test_sirup(spatial_id = "spot", clone = "clonal_group")

  polygon <- data.frame(
    x = c(-1, 4, 4, -1, -1),
    y = c(-1, -1, 4, 4, -1)
  )
  dominance <- spatialCloneDominance(x, polygon = polygon)
  expect_s3_class(dominance, "sirup_spatial_clone_dominance")
  expect_equal(dominance$metadata$receptor_fields$spatial_id, "spot")
  expect_equal(dominance$metadata$receptor_fields$clone, "clonal_group")
  expect_equal(sort(dominance$dominance$clone_id), c("A", "B"))
  expect_equal(sort(dominance$region_receptors$.sirup_spatial_id), c("cell1", "cell2"))
})

test_that("spatialLocalMoran uses the receptor ID stored in SIRUP metadata", {
  x <- make_test_sirup(spatial_id = "spot", clone = "clonal_group")

  local <- spatialLocalMoran(
    object = x,
    variable = "receptor_present",
    variable_type = "receptor",
    k = 1
  )

  expect_equal(local$value[match(paste0("cell", 1:4), local$cell)], c(1, 1, 0, 0))
})

test_that("spatialCloneJoinCount uses receptor fields stored in SIRUP metadata", {
  spatial <- sirup_test_seurat(
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
  x <- createSIRUP(
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
  expect_true(all(c("join_count_mark", "clone_id", ".sirup_spatial_id") %in% names(out$data)))
  expect_false("style" %in% names(out))
  expect_false("style" %in% names(out$summary))
  expect_equal(sort(unique(out$data$.sirup_spatial_id)), c("cell1", "cell2", "cell3"))
})

test_that("SIRUP print reports spatial and receptor dimensions", {
  x <- make_test_sirup()

  expect_output(
    print(x),
    paste0(
      "<SIRUP object>\\n",
      "Spatial data:    Seurat, 4 spatial IDs, 1 gene, 0 images\\n",
      "Receptor data:   tibble, 3 receptor sequences x 8 columns\\n",
      "Spatial overlap: 2/3 receptor IDs in spatial data \\(66.7%\\)\\n",
      "Spatial ID field: spatial_id\\n",
      "Clone ID field:  clone_id\\n",
      "Analyses:        0"
    )
  )
})

test_that("SIRUP summaries align values and report stored polygons", {
  x <- make_test_sirup()
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
    "(?s)  Spatial IDs:\\s+4.*  Receptor records:\\s+3.*  Image scale:\\s+hires.*  Swap coords:\\s+FALSE.*  Polygons:\\s+2",
    perl = TRUE
  )
})

test_that("multi-sample SIRUP summaries include per-sample counts", {

  make_section <- function() {
    spatial <- sirup_test_seurat(
      counts = matrix(1, nrow = 1, ncol = 2, dimnames = list("gene", c("cell1", "cell2")))
    )
    spatial$x <- c(0, 1)
    spatial$y <- c(0, 1)
    spatial
  }
  x <- createSIRUP(
    spatial = list(section_a = make_section(), section_b = make_section()),
    receptors = data.frame(
      barcode = c("cell1", "cell1", "cell2"),
      sample = c("section_a", "section_b", "section_b"),
      clone_id = c("A", "B", "C"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample_id = "sample",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )

  out <- summary(x)
  expect_equal(out$samples$sample, c("section_a", "section_b"))
  expect_equal(out$samples$n_receptor_records, c(1, 2))
  expect_equal(out$samples$n_receptor_ids, c(1, 2))
  expect_equal(out$samples$n_matched_ids, c(1, 2))
  expect_equal(out$image_scale, "hires")
  expect_false(out$swap_coords)
  expect_output(
    print(x),
    "(?s)Samples:\\s+2\\nSpatial ID field: barcode\\nSample ID field: sample\\nClone ID field:\\s+clone_id",
    perl = TRUE
  )
  expect_output(print(out), "(?s)Image scale:\\s+hires.*Swap coords:\\s+FALSE.*Per-sample summary", perl = TRUE)

  x_named <- createSIRUP(
    spatial = list(section_a = make_section(), section_b = make_section()),
    receptors = data.frame(
      barcode = c("cell1", "cell1", "cell2"),
      sample = c("section_a", "section_b", "section_b"),
      clone_id = c("A", "B", "C"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample_id = "sample",
    clone = "clone_id",
    image_scale = c(section_a = "hires", section_b = "lowres"),
    swap_coords = c(section_a = FALSE, section_b = TRUE)
  )
  expect_output(
    print(summary(x_named)),
    "(?s)Image scale:\\s+section_a=hires, section_b=lowres.*Swap coords:\\s+section_a=FALSE, section_b=TRUE",
    perl = TRUE
  )
})

test_that("legacy top-level schema objects are rejected", {
  x <- make_test_sirup()
  x$schema <- x$metadata$receptor_fields

  expect_error(validateSIRUP(x), "must contain, in order")
})

test_that("spatial polygon formatting accepts a SIRUP object", {
  x <- make_test_sirup()
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
