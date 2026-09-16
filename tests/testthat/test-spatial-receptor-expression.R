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
  # spot2 carries both BCR and TCR, so it is drawn as a single combined,
  # colour-averaged "BCR + TCR" ring rather than two overlapping rings.
  expect_identical(plot$layers[[2L]]$aes_params$colour, "blue")
  expect_identical(
    plot$layers[[3L]]$aes_params$colour,
    SIRAP:::.sirap_average_colors(c("blue", "orange"))
  )
  expect_true(any(vapply(
    plot$scales$scales,
    function(scale) "shape" %in% scale$aesthetics,
    logical(1)
  )))
})

test_that(".sirap_combine_receptor_groups merges multi-group spots with a blended color", {
  lookup <- data.frame(
    barcode = c("spot1", "spot2", "spot3", "spot3"),
    .sirap_receptor_group = c("BCR", "TCR", "BCR", "TCR"),
    stringsAsFactors = FALSE
  )
  res <- SIRAP:::.sirap_combine_receptor_groups(
    lookup, "barcode", c("BCR", "TCR"), c(BCR = "blue", TCR = "orange")
  )

  # One row per spot; the multi-group spot becomes the combined group.
  expect_identical(nrow(res$lookup), 3L)
  expect_identical(levels(res$lookup$.sirap_receptor_group), c("BCR", "TCR", "BCR + TCR"))
  expect_identical(
    as.character(res$lookup$.sirap_receptor_group[res$lookup$barcode == "spot3"]),
    "BCR + TCR"
  )
  # Single groups keep their colors; the combination is the average.
  expect_identical(res$colors[["BCR"]], "blue")
  expect_identical(res$colors[["TCR"]], "orange")
  expect_identical(res$colors[["BCR + TCR"]], SIRAP:::.sirap_average_colors(c("blue", "orange")))

  # overlap = a fixed color colors all combined spots with that color.
  res_fixed <- SIRAP:::.sirap_combine_receptor_groups(
    lookup, "barcode", c("BCR", "TCR"), c(BCR = "blue", TCR = "orange"), overlap = "black"
  )
  expect_identical(res_fixed$colors[["BCR + TCR"]], "black")
  expect_identical(res_fixed$colors[["BCR"]], "blue")
  expect_error(
    SIRAP:::.sirap_combine_receptor_groups(
      lookup, "barcode", c("BCR", "TCR"), c(BCR = "blue", TCR = "orange"), overlap = "not_a_color"
    ),
    "overlap"
  )
})

# Helper: values on the continuous receptor layer. Receptors map colour;
# the expression layer maps fill.
expression_receptor_fill_values <- function(plot) {
  layers <- Filter(function(l) any(c("colour", "color") %in% names(l$mapping)), plot$layers)
  sort(as.numeric(layers[[1]]$data$.sirap_count))
}

test_that("plotReceptorsExpression colours receptors by a continuous value", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  counts <- Matrix::Matrix(
    c(
      1, 2, 0, 3,
      0, 1, 2, 4
    ),
    nrow = 2, ncol = 4, sparse = TRUE
  )
  dimnames(counts) <- list(c("AICDA", "geneB"), paste0("cell", 1:4))
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object$x <- c(0, 1, 2, 3)
  object$y <- c(0, 1, 1, 0)

  receptors <- data.frame(
    spatial_id = c("cell1", "cell1", "cell2", "cell4"),
    receptor_type = c("BCR", "BCR", "BCR", "TCR"),
    shm = c(2, 8, 5, 1),
    stringsAsFactors = FALSE
  )
  common <- list(
    object, receptors = receptors, spatial_id = "spatial_id",
    feature = "AICDA", image_scale = "hires", swap_coords = FALSE
  )

  # Default mean: cell1 = mean(2, 8) = 5, cell2 = 5, cell4 = 1
  mean_plot <- do.call(plotReceptorsExpression, c(common, list(value = "shm")))
  expect_s3_class(mean_plot, "ggplot")
  expect_equal(expression_receptor_fill_values(mean_plot), c(1, 5, 5))

  max_plot <- do.call(plotReceptorsExpression, c(common, list(value = "shm", value_combine = "max")))
  expect_equal(expression_receptor_fill_values(max_plot), c(1, 5, 8))

  sum_plot <- do.call(plotReceptorsExpression, c(common, list(value = "shm", value_combine = "sum")))
  expect_equal(expression_receptor_fill_values(sum_plot), c(1, 5, 10))

  count_plot <- do.call(plotReceptorsExpression, c(common, list(value = "count")))
  expect_equal(expression_receptor_fill_values(count_plot), c(1, 1, 2))

  # Expression owns the fill scale; receptors own the colour scale.
  aes_used <- unlist(lapply(mean_plot$scales$scales, function(sc) sc$aesthetics))
  expect_true("fill" %in% aes_used)
  expect_true(any(c("colour", "color") %in% aes_used))

  # Receptors stay hollow rings in continuous mode, exactly as in grouped mode.
  rec_layer <- Filter(function(l) any(c("colour", "color") %in% names(l$mapping)), mean_plot$layers)[[1]]
  expect_identical(rec_layer$aes_params$shape, 21)
  expect_true(is.na(rec_layer$aes_params$fill))

  # The expression layer is the one mapping fill.
  expr_layer <- Filter(function(l) "fill" %in% names(l$mapping), mean_plot$layers)[[1]]
  expect_identical(expr_layer$aes_params$shape, 21)
  # colour = NA would make geom_point drop every row, so the spots must use a
  # transparent border instead. Building the plot must not drop anything.
  expect_false(is.na(expr_layer$aes_params$colour))
  expect_silent(built <- ggplot2::ggplot_build(mean_plot))
  expect_equal(nrow(built$data[[1]]), nrow(expr_layer$data %||% built$data[[1]]))
  expect_true(all(!is.na(built$data[[1]]$fill)))

  # `group` still works alongside `value`, restricting which receptors are drawn.
  grouped <- do.call(plotReceptorsExpression, c(common, list(
    value = "shm", group = "receptor_type", group_val = "TCR"
  )))
  expect_equal(expression_receptor_fill_values(grouped), 1)

  # Without value, group is still required.
  expect_error(do.call(plotReceptorsExpression, common), "`group` is required")

  # And the grouped path still behaves.
  grouped_only <- do.call(plotReceptorsExpression, c(common, list(group = "receptor_type")))
  expect_s3_class(grouped_only, "ggplot")
})
