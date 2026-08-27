make_image_free_spatial <- function() {
  counts <- Matrix::Matrix(
    c(
      1, 2, 0, 3,
      0, 1, 2, 4
    ),
    nrow = 2,
    ncol = 4,
    sparse = TRUE
  )
  dimnames(counts) <- list(c("geneA", "geneB"), paste0("cell", 1:4))
  object <- SeuratObject::CreateSeuratObject(
    counts = counts
  )
  object$x <- c(0, 1, 2, 3)
  object$y <- c(0, 1, 1, 0)
  object$cell_type <- c("B", "B", "T", "T")
  object
}

test_that("shared spatial plotting helpers use metadata coordinates without images", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  object <- make_image_free_spatial()

  expect_equal(SIRAP:::.sirap_spatial_images(object), "metadata")
  spatial <- SIRAP:::.sirap_spatial_data(
    object,
    image = "metadata",
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_null(spatial$image)
  expect_equal(spatial$coords$cell, paste0("cell", 1:4))
  expect_equal(spatial$coords$x, c(0, 1, 2, 3))
  expect_equal(spatial$coords$.sirap_y, c(0, -1, -1, 0))

  plot <- SIRAP:::.sirap_spatial_base(spatial$coords, image = NULL)
  expect_s3_class(plot, "ggplot")
})

test_that("core spatial ggplot functions work without Seurat images", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("patchwork")

  object <- make_image_free_spatial()
  receptors <- data.frame(
    spatial_id = c("cell1", "cell2", "cell4"),
    receptor_type = c("BCR", "BCR", "TCR"),
    stringsAsFactors = FALSE
  )

  expression_plot <- spatialExpressionPlot(
    object,
    features = "geneA",
    image_scale = "hires",
    swap_coords = FALSE
  )
  discrete_plot <- spatialDiscretePlot(
    object,
    meta_col = "cell_type",
    image_scale = "hires",
    swap_coords = FALSE
  )
  receptor_plot <- plotReceptors(
    object,
    receptors = receptors,
    group = "receptor_type",
    spatial_id = "spatial_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
  receptor_expression_plot <- plotReceptorsExpression(
    object,
    receptors = receptors,
    group = "receptor_type",
    spatial_id = "spatial_id",
    feature = "geneA",
    image_scale = "hires",
    swap_coords = FALSE
  )
  blend_plot <- spatialExpressionBlendPlot(
    object,
    features = c("geneA", "geneB"),
    only_blended = TRUE,
    image_scale = "hires",
    swap_coords = FALSE
  )
  multi_plot <- spatialExpressionMultiPlot(
    object,
    features = c("geneA", "geneB"),
    image_scale = "hires",
    swap_coords = FALSE
  )
  color_tip_plot <- plotSpatialColorTips(
    object,
    tip_data = data.frame(
      spatial_id = c("cell1", "cell3"),
      tip_color = c("#111111", "#CC5500"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "spatial_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
  clone_color_plot <- plotSpatialCloneColorTips(
    data.frame(
      clone_id = c("clone1", "clone1", "clone2"),
      spatial_id = c("cell1", "cell2", "cell3"),
      tip_color = c("#111111", "#222222", "#333333"),
      stringsAsFactors = FALSE
    ),
    clone_val = "clone1",
    tissue_positions = object,
    spatial_id = "spatial_id",
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_s3_class(expression_plot, "ggplot")
  expect_s3_class(discrete_plot, "ggplot")
  expect_s3_class(receptor_plot, "ggplot")
  expect_s3_class(receptor_expression_plot, "ggplot")
  expect_true(inherits(blend_plot, "ggplot") || inherits(blend_plot, "patchwork"))
  expect_s3_class(multi_plot, "ggplot")
  expect_s3_class(color_tip_plot, "ggplot")
  expect_s3_class(clone_color_plot, "ggplot")
})

test_that("spatialExpressionPlot uses feature colors or palette", {
  expect_equal(length(.sirap_feature_color_values()), 256)
  expect_identical(
    .sirap_feature_color_values(feature_colors = c("black", "cyan"), feature_color_palette = "plasma"),
    c("black", "cyan")
  )
  expect_error(
    .sirap_feature_color_values(feature_color_palette = "not_a_palette"),
    "`feature_color_palette`"
  )
})

test_that(".sirap_fetch_data tries layer before legacy slot fallback", {
  calls <- character()
  testthat::local_mocked_bindings(
    .sirap_fetch_data_layer = function(object, vars, layer) {
      calls <<- c(calls, paste0("layer:", layer))
      data.frame(geneA = 1, row.names = "cell1")
    },
    .sirap_fetch_data_slot = function(object, vars, slot) {
      calls <<- c(calls, paste0("slot:", slot))
      data.frame(geneA = 2, row.names = "cell1")
    },
    .package = "SIRAP"
  )

  out <- .sirap_fetch_data(list(), vars = "geneA", layer = "counts")

  expect_equal(out$geneA, 1)
  expect_identical(calls, "layer:counts")
})

test_that(".sirap_fetch_data falls back to slot only after layer errors", {
  calls <- character()
  testthat::local_mocked_bindings(
    .sirap_fetch_data_layer = function(object, vars, layer) {
      calls <<- c(calls, paste0("layer:", layer))
      rlang::abort("layer failed")
    },
    .sirap_fetch_data_slot = function(object, vars, slot) {
      calls <<- c(calls, paste0("slot:", slot))
      data.frame(geneA = 2, row.names = "cell1")
    },
    .package = "SIRAP"
  )

  out <- .sirap_fetch_data(list(), vars = "geneA", layer = "counts")

  expect_equal(out$geneA, 2)
  expect_identical(calls, c("layer:counts", "slot:counts"))
})

test_that(".sirap_fetch_data reports ordered layer and slot failures", {
  testthat::local_mocked_bindings(
    .sirap_fetch_data_assay = function(object) "RNA",
    .sirap_fetch_data_layer = function(object, vars, layer) {
      rlang::abort("layer failed")
    },
    .sirap_fetch_data_slot = function(object, vars, slot) {
      rlang::abort("slot failed")
    },
    .package = "SIRAP"
  )

  expect_error(
    .sirap_fetch_data(list(), vars = "geneA", layer = "counts"),
    "SIRAP tried layer = counts in assay = RNA first, then tried slot = counts as a fall back"
  )
})

test_that("spatial statistic and polygon plots work without Seurat images", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  object <- make_image_free_spatial()
  local_moran <- data.frame(
    cell = paste0("cell", 1:4),
    x = c(0, 1, 2, 3),
    y = c(0, 1, 1, 0),
    .sirap_y = c(0, -1, -1, 0),
    image = "metadata",
    variable = "geneA",
    z_i = c(-1, 0, 1, 2),
    stringsAsFactors = FALSE
  )
  polygon <- data.frame(
    polygon = "region",
    x = c(-0.5, 3.5, 3.5, -0.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5, -0.5),
    stringsAsFactors = FALSE
  )

  local_plot <- plotSpatialLocalMoran(
    object,
    local_moran = local_moran,
    image_scale = "hires",
    swap_coords = FALSE
  )
  polygon_plot <- plotSpatialPolygons(
    object,
    polygon = polygon,
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_s3_class(local_plot, "ggplot")
  expect_s3_class(polygon_plot, "ggplot")
})

test_that("3D clone density prepares image-free plot data", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("pals")

  object <- make_image_free_spatial()
  receptors <- data.frame(
    spatial_id = paste0("cell", 1:4),
    clone_id = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  polygon <- data.frame(
    polygon = "region",
    x = c(-0.5, 3.5, 3.5, -0.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5, -0.5),
    stringsAsFactors = FALSE
  )
  dominance <- spatialCloneDominance(
    receptors = receptors,
    object = object,
    polygon = polygon,
    spatial_id = "spatial_id",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )

  plot_data <- plotSpatialCloneDensity3D(
    dominance,
    object = object,
    image_scale = "hires",
    swap_coords = FALSE,
    top_n = 2,
    return_data = TRUE
  )

  expect_null(plot_data$image)
  expect_true(all(c("xlim", "ylim") %in% names(plot_data$image_extent)))
  expect_true(all(is.finite(plot_data$image_extent$xlim)))
  expect_true(all(is.finite(plot_data$image_extent$ylim)))
})

# Helper: count the GeomPolygon layers on a ggplot (works on the resolved plot).
n_polygon_layers <- function(plot) {
  sum(vapply(plot$layers, function(l) inherits(l$geom, "GeomPolygon"), logical(1)))
}

test_that("plotReceptors overlays polygons from a data frame and a CSV path", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  object <- make_image_free_spatial()
  receptors <- data.frame(
    spatial_id = c("cell1", "cell2", "cell4"),
    receptor_type = c("BCR", "BCR", "TCR"),
    stringsAsFactors = FALSE
  )
  polygon <- data.frame(
    polygon = c(rep("left", 5), rep("right", 5)),
    x = c(-0.5, 1.5, 1.5, -0.5, -0.5, 1.5, 3.5, 3.5, 1.5, 1.5),
    y = c(-0.5, -0.5, 1.5, 1.5, -0.5, -0.5, -0.5, 1.5, 1.5, -0.5),
    stringsAsFactors = FALSE
  )

  # A data frame of vertices draws one polygon fill layer plus a label layer.
  df_plot <- plotReceptors(
    object,
    receptors = receptors,
    group = "receptor_type",
    spatial_id = "spatial_id",
    image_scale = "hires",
    swap_coords = FALSE,
    polygons = polygon,
    polygon_label = TRUE
  )
  expect_s3_class(df_plot, "ggplot")
  expect_equal(n_polygon_layers(df_plot), 1L)

  # A CSV path resolves to the same overlay.
  csv_path <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(polygon, csv_path, row.names = FALSE)
  csv_plot <- plotReceptors(
    object,
    receptors = receptors,
    group = "receptor_type",
    spatial_id = "spatial_id",
    image_scale = "hires",
    swap_coords = FALSE,
    polygons = csv_path,
    polygon_val = "left"
  )
  expect_s3_class(csv_plot, "ggplot")
  expect_equal(n_polygon_layers(csv_plot), 1L)

  # Default `polygons = FALSE` draws no polygon layer.
  none_plot <- plotReceptors(
    object,
    receptors = receptors,
    group = "receptor_type",
    spatial_id = "spatial_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
  expect_equal(n_polygon_layers(none_plot), 0L)
})

test_that("plotSpatialColorTips overlays SIRAP object polygons", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  object <- make_image_free_spatial()
  receptors <- data.frame(
    spatial_id = c("cell1", "cell2", "cell3"),
    sequence_id = paste0("seq", 1:3),
    stringsAsFactors = FALSE
  )
  polygon <- data.frame(
    polygon_id = "region",
    x = c(-0.5, 3.5, 3.5, -0.5, -0.5),
    y = c(-0.5, -0.5, 1.5, 1.5, -0.5),
    stringsAsFactors = FALSE
  )
  sirap <- createSIRAP(
    spatial = object,
    receptors = receptors,
    spatial_id = "spatial_id",
    locus = NULL,
    image_scale = "hires",
    swap_coords = FALSE,
    polygons = polygon
  )

  tip_plot <- plotSpatialColorTips(
    sirap,
    tip_data = data.frame(
      spatial_id = c("cell1", "cell3"),
      tip_color = c("#111111", "#CC5500"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "spatial_id",
    polygons = TRUE
  )
  expect_s3_class(tip_plot, "ggplot")
  expect_equal(n_polygon_layers(tip_plot), 1L)
})
