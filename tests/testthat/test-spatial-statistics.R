test_that("kcross data frame conversion uses theoretical fallback", {
  k <- data.frame(
    r = c(0, 1, 2),
    border = c(0, 2, 8)
  )
  out <- SIRAP:::.sirap_kcross_df(k, correction = "border", image = "sample")

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

  coords <- SIRAP:::.sirap_spatial_coordinate_table(object)

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
    max_r = 20
  )

  expect_equal(nrow(out$points), 6)
  expect_true("Other" %in% out$points$mark)
  expect_true(all(c("A", "B") %in% out$points$mark))

  sirap <- createSIRAP(
    object,
    data.frame(barcode = "cell1"),
    barcode_col = "barcode"
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
    clone_ids = "A",
    k = 1,
    min_barcodes = 2,
    return_data = TRUE
  )

  expect_true("A" %in% out$summary$clone_id)
  expect_equal(out$summary$n_clone_barcodes[out$summary$clone_id == "A"], 3)
  expect_true(all(c("join_count_mark", "clone_id", ".sirap_barcode") %in% names(out$data)))
})
