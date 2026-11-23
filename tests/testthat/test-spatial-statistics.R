test_that("kcross data frame conversion uses theoretical fallback", {
  k <- data.frame(
    r = c(0, 1, 2),
    border = c(0, 2, 8)
  )
  out <- SIRUP:::.sirup_kcross_df(k, image = "sample")

  expect_equal(out$image, rep("sample", 3))
  expect_equal(out$observed, k$border)
  expect_equal(out$theoretical, pi * k$r^2)
  expect_equal(out$simulated, out$theoretical)
})

test_that("spatial coordinate table falls back to metadata coordinates", {

  object <- sirup_test_seurat(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 3,
      dimnames = list("gene", paste0("cell", 1:3))
    )
  )
  object$x <- c(0, 10, 20)
  object$y <- c(0, 5, 10)

  coords <- SIRUP:::.sirup_spatial_coordinate_table(
    object,
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_equal(coords$cell, paste0("cell", 1:3))
  expect_equal(coords$x, c(0, 10, 20))
  expect_equal(coords$.sirup_y, c(0, -5, -10))
  expect_equal(unique(coords$image), "metadata")
})

test_that("spatialRipleyCrossK keeps all non-missing marks in the point pattern", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")

  object <- sirup_test_seurat(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 6,
      dimnames = list("gene", paste0("cell", 1:6))
    )
  )
  object$x <- c(0, 10, 20, 30, 40, 50)
  object$y <- c(0, 0, 0, 0, 0, 0)
  object$cell_type <- c("A", "B", "Other", "A", "B", "Other")

  out <- spatialRipleyCrossK(
    object,
    mark_col = "cell_type",
    i = "A",
    j = "B",
    nsim = 0,
    max_r = 20,
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_equal(nrow(out$points), 6)
  expect_true("Other" %in% out$points$mark)
  expect_true(all(c("A", "B") %in% out$points$mark))

  sirup <- createSIRUP(
    object,
    data.frame(spatial_id = "cell1"),
    image_scale = "hires",
    swap_coords = FALSE
  )
  sirup_out <- spatialRipleyCrossK(
    sirup,
    mark_col = "cell_type",
    i = "A",
    j = "B",
    nsim = 0,
    max_r = 20
  )
  expect_equal(sirup_out$points, out$points)
})

test_that("spatialMoran can test a particular clone from a SIRUP object", {

  object <- sirup_test_seurat(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 6,
      dimnames = list("gene", paste0("cell", 1:6))
    )
  )
  object$x <- seq(0, 5)
  object$y <- c(0, 0, 0, 0, 0, 0)
  receptors <- data.frame(
    spot = paste0("cell", 1:6),
    clonal_group = c("A", "A", "A", "B", "B", "C"),
    stringsAsFactors = FALSE
  )
  sirup <- createSIRUP(
    object,
    receptors,
    spatial_id = "spot",
    clone = "clonal_group",
    image_scale = "hires",
    swap_coords = FALSE
  )

  out <- spatialMoran(
    sirup,
    variable = "A",
    variable_type = "clone",
    k = 2
  )
  local <- spatialLocalMoran(
    sirup,
    variable = "A",
    variable_type = "clone",
    k = 2
  )

  expect_equal(out$variable, "clone_A")
  expect_equal(out$data$value, c(1, 1, 1, 0, 0, 0))
  expect_true(inherits(out$test, "htest"))
  expect_equal(local$variable, rep("clone_A", 6))
  expect_true(all(c("local_i", "p_value") %in% names(local)))
})

test_that("spatialLee tests bivariate feature association", {

  object <- sirup_test_seurat(
    counts = matrix(
      c(
        1, 2, 3, 4, 5, 6, 7, 8,
        2, 3, 4, 5, 6, 7, 8, 9
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("gene-x", "gene-y"), paste0("cell", 1:8))
    )
  )
  object$x <- seq(0, 7)
  object$y <- rep(0, 8)

  out <- spatialLee(
    object,
    x = "gene-x",
    y = "gene-y",
    k = 2,
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_s3_class(out$test, "htest")
  expect_equal(out$x, "gene-x")
  expect_equal(out$y, "gene-y")
  expect_equal(out$data$value_x, seq(1, 8))
  expect_equal(out$data$value_y, seq(2, 9))
  expect_true("Lee's L statistic" %in% names(out$test$estimate))
})

test_that("spatialCloneJoinCount tests clone marks from receptor IDs", {

  object <- sirup_test_seurat(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 6,
      dimnames = list("gene", paste0("cell", 1:6))
    )
  )
  object$x <- seq(0, 5)
  object$y <- c(0, 0, 0, 0, 0, 0)
  receptors <- data.frame(
    barcode_8um = paste0("cell", 1:6),
    clone_id = c("A", "A", "A", "B", "B", "C"),
    stringsAsFactors = FALSE
  )

  out <- spatialCloneJoinCount(
    object = object,
    receptors = receptors,
    clone = "clone_id",
    spatial_id = "barcode_8um",
    clone_val = "A",
    k = 2,
    min_spatial_ids = 2,
    return_data = TRUE,
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_true("A" %in% out$summary$clone_id)
  expect_equal(out$summary$n_clone_spatial_ids[out$summary$clone_id == "A"], 3)
  expect_true(all(c("join_count_mark", "clone_id", ".sirup_spatial_id") %in% names(out$data)))
})

test_that("plotSpatialCloneJoinCount returns one combined plot", {

  object <- sirup_test_seurat(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 6,
      dimnames = list("gene", paste0("cell", 1:6))
    )
  )
  object$x <- seq(0, 5)
  object$y <- c(0, 0, 0, 0, 0, 0)
  receptors <- data.frame(
    barcode_8um = paste0("cell", 1:6),
    clone_id = c("A", "A", "A", "B", "B", "C"),
    stringsAsFactors = FALSE
  )
  join_count <- spatialCloneJoinCount(
    object = object,
    receptors = receptors,
    clone = "clone_id",
    spatial_id = "barcode_8um",
    clone_val = "A",
    k = 2,
    min_spatial_ids = 2,
    return_data = TRUE,
    image_scale = "hires",
    swap_coords = FALSE
  )

  plot <- plotSpatialCloneJoinCount(
    join_count = join_count,
    object = object,
    clone_val = "A",
    polygon = data.frame(
      x = c(-1, 4, 4, -1, -1),
      y = c(-1, -1, 4, 4, -1),
      polygon = "test"
    ),
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_false(is.list(plot) && !inherits(plot, "ggplot"))
})

test_that("plotSpatialCloneJoinCount requires stored mark data", {

  object <- sirup_test_seurat(
    counts = matrix(1, nrow = 1, ncol = 2, dimnames = list("gene", c("cell1", "cell2")))
  )
  object$x <- c(0, 1)
  object$y <- c(0, 0)
  join_count <- list(summary = data.frame(clone_id = "A"), data = NULL)

  expect_error(
    plotSpatialCloneJoinCount(
      join_count,
      object,
      clone_val = "A",
      image_scale = "hires",
      swap_coords = FALSE
    ),
    "return_data = TRUE",
    fixed = TRUE
  )
})

test_that("spatialLee can use the Monte Carlo method", {

  object <- sirup_test_seurat(
    counts = matrix(
      c(
        1, 2, 3, 4, 5, 6, 7, 8,
        2, 3, 4, 5, 6, 7, 8, 9
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("gene-x", "gene-y"), paste0("cell", 1:8))
    )
  )
  object$x <- seq(0, 7)
  object$y <- rep(0, 8)

  # spdep::lee.test() densifies the weights matrix, so it cannot run on
  # Visium HD. The mc method permutes labels instead and keeps the weights
  # sparse; both must return the same observed L.
  analytical <- spatialLee(
    object, x = "gene-x", y = "gene-y", k = 2,
    image_scale = "hires", swap_coords = FALSE
  )
  mc <- spatialLee(
    object, x = "gene-x", y = "gene-y", k = 2,
    image_scale = "hires", swap_coords = FALSE,
    method = "mc", nsim = 99
  )

  expect_s3_class(mc$test, "htest")
  expect_equal(mc$method, "mc")
  expect_equal(analytical$method, "analytical")
  expect_equal(mc$L, analytical$L, tolerance = 1e-8)
  expect_true(mc$test$p.value >= 0 && mc$test$p.value <= 1)
})

test_that("spatialLee mc matches spdep::lee.mc under the same seed", {

  set.seed(42)
  counts <- matrix(rpois(2 * 40, lambda = 5), nrow = 2,
                   dimnames = list(c("gene-x", "gene-y"), paste0("cell", 1:40)))
  object <- sirup_test_seurat(counts = counts)
  object$x <- rep(0:7, each = 5)
  object$y <- rep(0:4, times = 8)

  # spatialLee runs the permutation loop itself so it can report progress.
  # It must still reproduce spdep::lee.mc exactly: same RNG draws, same L,
  # same rank and p-value.
  set.seed(7)
  ours <- spatialLee(object, x = "gene-x", y = "gene-y", k = 4,
                     image_scale = "hires", swap_coords = FALSE,
                     method = "mc", nsim = 99)

  dat <- ours$data
  listw <- ours$weights
  set.seed(7)
  theirs <- spdep::lee.mc(dat$value_x, dat$value_y, listw, nsim = 99,
                          zero.policy = TRUE, alternative = "greater")

  expect_equal(ours$L, unname(theirs$statistic))
  expect_equal(ours$test$p.value, theirs$p.value)
  expect_equal(unname(ours$test$parameter), unname(theirs$parameter))
  expect_equal(ours$test$res, theirs$res)
})

test_that("clone spatial statistics distinguish duplicate barcodes across spatial samples", {

  make_section <- function() {
    object <- sirup_test_seurat(
      counts = matrix(
        1,
        nrow = 1,
        ncol = 4,
        dimnames = list("gene", paste0("cell", 1:4))
      )
    )
    object$x <- 0:3
    object$y <- c(0, 0, 0, 0)
    object
  }
  sirup <- createSIRUP(
    spatial = list(section_a = make_section(), section_b = make_section()),
    receptors = data.frame(
      barcode = c("cell1", "cell2", "cell3", "cell4", "cell1", "cell2", "cell3", "cell4"),
      sample = c("section_a", "section_a", "section_a", "section_a", "section_b", "section_b", "section_b", "section_b"),
      clone_id = c("A", "A", "B", "C", "B", "B", "A", "C"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample_id = "sample",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )

  values <- SIRUP:::.sirup_spatial_variable_values(
    sirup,
    variable = "A",
    variable_type = "clone",
    value = NULL,
    receptors = NULL,
    spatial_id_col = NULL,
    clone_col = NULL,
    layer = "counts",
    images = NULL,
    image_scale = NULL,
    swap_coords = NULL
  )
  by_key <- stats::setNames(values$value, paste(values$sample, values$cell, sep = ":"))
  expect_equal(unname(by_key[c("section_a:cell1", "section_a:cell2")]), c(1, 1))
  expect_equal(unname(by_key[c("section_b:cell1", "section_b:cell2")]), c(0, 0))

  moran <- spatialMoran(sirup, variable = "A", variable_type = "clone", k = 1)
  expect_true(all(c("section_a::metadata", "section_b::metadata") %in% moran$summary$image))
  expect_equal(
    unname(stats::setNames(moran$data$value, paste(moran$data$sample, moran$data$cell, sep = ":"))[
      c("section_a:cell1", "section_b:cell1")
    ]),
    c(1, 0)
  )
})

test_that("spatialLee builds independent neighbor graphs per spatial sample", {

  make_section <- function(offset = 0) {
    object <- sirup_test_seurat(
      counts = matrix(
        c(1, 2, 3, 4, 4, 3, 2, 1),
        nrow = 2,
        dimnames = list(c("geneA", "geneB"), paste0("cell", 1:4))
      )
    )
    object$x <- offset + 0:3
    object$y <- c(0, 0, 0, 0)
    object
  }
  sirup <- createSIRUP(
    spatial = list(section_a = make_section(), section_b = make_section(offset = 1000)),
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

  out <- spatialLee(sirup, x = "geneA", y = "geneB", k = 1)

  expect_equal(sort(names(out$tests)), c("section_a::metadata", "section_b::metadata"))
  expect_equal(out$summary$image, c("section_a::metadata", "section_b::metadata"))
  expect_true(all(table(out$data$image) == 4))
})

test_that("spatialRipleyCrossK reads marks from multisample spatial metadata", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")

  make_section <- function(mark) {
    object <- sirup_test_seurat(
      counts = matrix(1, nrow = 1, ncol = 4, dimnames = list("gene", paste0("cell", 1:4)))
    )
    object$x <- c(0, 1, 0, 1)
    object$y <- c(0, 0, 1, 1)
    object$mark <- mark
    object
  }
  sirup <- createSIRUP(
    spatial = list(
      section_a = make_section(c("A", "B", "A", "B")),
      section_b = make_section(c("B", "A", "B", "A"))
    ),
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

  out <- spatialRipleyCrossK(sirup, mark_col = "mark", i = "A", j = "B")

  expect_equal(sort(unique(out$points$sample)), c("section_a", "section_b"))
  expect_equal(sort(names(out$observed)), c("section_a::metadata", "section_b::metadata"))
  expect_equal(as.integer(table(out$points$sample)), c(4, 4))
})

test_that("spatialCloneJoinCount uses sample-aware receptor coordinate joins", {

  make_section <- function() {
    object <- sirup_test_seurat(
      counts = matrix(
        1,
        nrow = 1,
        ncol = 4,
        dimnames = list("gene", paste0("cell", 1:4))
      )
    )
    object$x <- 0:3
    object$y <- c(0, 0, 0, 0)
    object
  }
  sirup <- createSIRUP(
    spatial = list(section_a = make_section(), section_b = make_section()),
    receptors = data.frame(
      barcode = c("cell1", "cell2", "cell3", "cell4", "cell1", "cell2", "cell3", "cell4"),
      sample = c("section_a", "section_a", "section_a", "section_a", "section_b", "section_b", "section_b", "section_b"),
      clone_id = c("A", "A", "B", "C", "B", "B", "A", "C"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample_id = "sample",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )

  out <- spatialCloneJoinCount(
    sirup,
    clone_val = c("A", "B"),
    k = 1,
    min_spatial_ids = 2,
    return_data = TRUE
  )

  expect_true(all(c("section_a::metadata", "section_b::metadata") %in% out$summary$image))
  clone_a <- out$data[out$data$clone_id == "A", , drop = FALSE]
  clone_b <- out$data[out$data$clone_id == "B", , drop = FALSE]
  expect_true(all(as.character(clone_a$join_count_mark[clone_a$.sirup_sample == "section_a" & clone_a$cell %in% c("cell1", "cell2")]) == "A"))
  expect_false(any(clone_a$.sirup_sample == "section_b" & clone_a$cell %in% c("cell1", "cell2") & as.character(clone_a$join_count_mark) == "A"))
  expect_true(all(as.character(clone_b$join_count_mark[clone_b$.sirup_sample == "section_b" & clone_b$cell %in% c("cell1", "cell2")]) == "B"))
  expect_false(any(clone_b$.sirup_sample == "section_a" & clone_b$cell %in% c("cell1", "cell2") & as.character(clone_b$join_count_mark) == "B"))
})

test_that("spatialCloneJoinCount honors use_stored_polygons for multisample objects", {

  make_section <- function() {
    object <- sirup_test_seurat(
      counts = matrix(1, nrow = 1, ncol = 4, dimnames = list("gene", paste0("cell", 1:4)))
    )
    object$x <- 0:3
    object$y <- c(0, 0, 0, 0)
    object
  }
  # One GC per section covering cells 1-3 (x 0,1,2); cell 4 is outside.
  gc <- function(samp) data.frame(
    polygon = "gc1", sample = samp,
    x = c(-0.5, 2.5, 2.5, -0.5), y = c(-0.5, -0.5, 0.5, 0.5),
    stringsAsFactors = FALSE
  )
  make_object <- function(with_polygons) {
    obj <- createSIRUP(
      spatial = list(section_a = make_section(), section_b = make_section()),
      receptors = data.frame(
        barcode = rep(paste0("cell", 1:4), 2),
        sample = rep(c("section_a", "section_b"), each = 4),
        clone_id = rep(c("A", "A", "C", "C"), 2),
        stringsAsFactors = FALSE
      ),
      spatial_id = "barcode", sample_id = "sample", clone = "clone_id",
      image_scale = "hires", swap_coords = FALSE
    )
    if (with_polygons) obj <- addPolygons(obj, polygons = rbind(gc("section_a"), gc("section_b")))
    obj
  }
  sirup <- make_object(with_polygons = TRUE)

  # Default (NULL) auto-detects stored polygons -> one graph PER GC (3 spots each).
  per_gc <- spatialCloneJoinCount(sirup, clone_val = "A", k = 1, min_spatial_ids = 2)
  expect_true(all(per_gc$summary$n_spatial_ids == 3))

  # use_stored_polygons = FALSE -> one tissue-wide graph PER SECTION (4 spots each).
  tissue_wide <- spatialCloneJoinCount(
    sirup, clone_val = "A", k = 1, min_spatial_ids = 2, use_stored_polygons = FALSE
  )
  expect_true(all(tissue_wide$summary$n_spatial_ids == 4))
  expect_true(all(c("section_a::metadata", "section_b::metadata") %in% tissue_wide$summary$image))

  # Explicit TRUE requires a SIRUP object with stored polygons.
  expect_error(
    spatialCloneJoinCount(make_object(with_polygons = FALSE), clone_val = "A", k = 1,
                          use_stored_polygons = TRUE),
    "stored polygons"
  )
})

test_that("Moran statistics run with the default variable_type", {
  skip_if_not_installed("spdep")
  skip_if_not_installed("SeuratObject")
  path <- system.file("extdata", "sirup_synthetic.rds", package = "SIRUP")
  # The shipped object has the installed version, so reading it gives no warning.
  expect_no_warning(x <- readSIRUP(path))
  expect_no_error(spatialMoran(x, "gene1", k = 4))
  expect_no_error(spatialLocalMoran(x, "gene1", k = 4))
})
