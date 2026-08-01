test_that("kcross data frame conversion uses theoretical fallback", {
  k <- data.frame(
    r = c(0, 1, 2),
    border = c(0, 2, 8)
  )
  out <- SIRAP:::.sirap_kcross_df(k, image = "sample")

  expect_equal(out$image, rep("sample", 3))
  expect_equal(out$observed, k$border)
  expect_equal(out$theoretical, pi * k$r^2)
  expect_equal(out$simulated, out$theoretical)
})

test_that("spatial coordinate table falls back to metadata coordinates", {
  skip_if_not_installed("SeuratObject")

  object <- SeuratObject::CreateSeuratObject(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 3,
      dimnames = list("gene", paste0("cell", 1:3))
    )
  )
  object$x <- c(0, 10, 20)
  object$y <- c(0, 5, 10)

  coords <- SIRAP:::.sirap_spatial_coordinate_table(
    object,
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_equal(coords$cell, paste0("cell", 1:3))
  expect_equal(coords$x, c(0, 10, 20))
  expect_equal(coords$.sirap_y, c(0, -5, -10))
  expect_equal(unique(coords$image), "metadata")
})

test_that("spatialRipleyCrossK keeps all non-missing marks in the point pattern", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")

  object <- SeuratObject::CreateSeuratObject(
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

  sirap <- createSIRAP(
    object,
    data.frame(barcode = "cell1"),
    barcode = "barcode",
    image_scale = "hires",
    swap_coords = FALSE
  )
  sirap_out <- spatialRipleyCrossK(
    sirap,
    mark_col = "cell_type",
    i = "A",
    j = "B",
    nsim = 0,
    max_r = 20
  )
  expect_equal(sirap_out$points, out$points)
})

test_that("spatialMoran can test a particular clone from a SIRAP object", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("spdep")

  object <- SeuratObject::CreateSeuratObject(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 6,
      dimnames = list("gene", paste0("cell", 1:6))
    )
  )
  object$x <- c(0, 1, 2, 10, 11, 12)
  object$y <- c(0, 0, 0, 0, 0, 0)
  receptors <- data.frame(
    spot = paste0("cell", 1:6),
    clonal_group = c("A", "A", "A", "B", "B", "C"),
    stringsAsFactors = FALSE
  )
  sirap <- createSIRAP(
    object,
    receptors,
    barcode = "spot",
    clone = "clonal_group",
    image_scale = "hires",
    swap_coords = FALSE
  )

  out <- spatialMoran(
    sirap,
    variable = "A",
    variable_type = "clone",
    k = 1
  )
  local <- spatialLocalMoran(
    sirap,
    variable = "A",
    variable_type = "clone",
    k = 1
  )

  expect_equal(out$variable, "clone_A")
  expect_equal(out$data$value, c(1, 1, 1, 0, 0, 0))
  expect_true(inherits(out$test, "htest"))
  expect_equal(local$variable, rep("clone_A", 6))
  expect_true(all(c("local_i", "p_value") %in% names(local)))
})

test_that("spatialCloneJoinCount tests clone marks from receptor barcodes", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("spdep")

  object <- SeuratObject::CreateSeuratObject(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 6,
      dimnames = list("gene", paste0("cell", 1:6))
    )
  )
  object$x <- c(0, 1, 2, 10, 11, 12)
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
    barcode = "barcode_8um",
    clone_val = "A",
    k = 1,
    min_barcodes = 2,
    return_data = TRUE,
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_true("A" %in% out$summary$clone_id)
  expect_equal(out$summary$n_clone_barcodes[out$summary$clone_id == "A"], 3)
  expect_true(all(c("join_count_mark", "clone_id", ".sirap_barcode") %in% names(out$data)))
})

test_that("plotSpatialCloneJoinCount returns one combined plot", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("spdep")
  skip_if_not_installed("cowplot")

  object <- SeuratObject::CreateSeuratObject(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 6,
      dimnames = list("gene", paste0("cell", 1:6))
    )
  )
  object$x <- c(0, 1, 2, 10, 11, 12)
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
    barcode = "barcode_8um",
    clone_val = "A",
    k = 1,
    min_barcodes = 2,
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
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("cowplot")

  object <- SeuratObject::CreateSeuratObject(
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
