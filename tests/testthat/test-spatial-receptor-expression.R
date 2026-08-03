test_that("plotReceptorsExpression layers grouped receptors over expression", {
  testthat::local_mocked_bindings(
    .sirap_resolve_inputs = function(object, receptors = NULL, spatial_id = NULL) {
      list(object = object, receptors = receptors, spatial_id_col = spatial_id)
    },
    .require_namespace = function(package) invisible(TRUE),
    .as_data_frame = as.data.frame,
    .guess_spatial_id_col = function(receptors, spatial_id = NULL) spatial_id,
    .check_required_cols = function(x, cols, label) {
      stopifnot(all(cols %in% colnames(x)))
    },
    .sirap_fetch_data = function(object, vars, layer) {
      out <- data.frame(GENE = c(1, 2, 0), row.names = c("spot1", "spot2", "spot3"))
      colnames(out) <- vars
      out
    },
    .sirap_spatial_images = function(object, images = NULL) "image1",
    .sirap_spatial_data = function(object, image, image_scale, swap_coords) {
      list(
        coords = data.frame(
          cell = c("spot1", "spot2", "spot3"),
          x = c(1, 2, 3),
          .sirap_y = c(1, 2, 3)
        ),
        image = NULL
      )
    },
    .sirap_spatial_base = function(...) ggplot2::ggplot(),
    .sirap_return_spatial_plots = function(plots, combine = TRUE) {
      if (combine) plots[[1L]] else plots
    },
    .package = "SIRAP"
  )

  receptors <- data.frame(
    barcode = c("spot1", "spot2", "spot2"),
    receptor_type = c("BCR", "BCR", "TCR")
  )

  plot <- plotReceptorsExpression(
    object = list(),
    receptors = receptors,
    group = "receptor_type",
    group_val = c("BCR", "TCR"),
    spatial_id = "barcode",
    feature = "GENE",
    colors = c(BCR = "blue", TCR = "orange"),
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_length(plot$layers, 3L)
  expect_identical(plot$layers[[2L]]$aes_params$colour, "blue")
  expect_identical(plot$layers[[3L]]$aes_params$colour, "orange")
  expect_true(any(vapply(
    plot$scales$scales,
    function(scale) "shape" %in% scale$aesthetics,
    logical(1)
  )))
})
