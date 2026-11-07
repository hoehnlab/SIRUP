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
  object <- sirup_test_seurat(
    counts = counts
  )
  object$x <- c(0, 1, 2, 3)
  object$y <- c(0, 1, 1, 0)
  object$cell_type <- c("B", "B", "T", "T")
  object
}

test_that("shared spatial plotting helpers use metadata coordinates without images", {

  object <- make_image_free_spatial()

  expect_equal(SIRUP:::.sirup_spatial_images(object), "metadata")
  spatial <- SIRUP:::.sirup_spatial_data(
    object,
    image = "metadata",
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_null(spatial$image)
  expect_equal(spatial$coords$cell, paste0("cell", 1:4))
  expect_equal(spatial$coords$x, c(0, 1, 2, 3))
  expect_equal(spatial$coords$.sirup_y, c(0, -1, -1, 0))

  plot <- SIRUP:::.sirup_spatial_base(spatial$coords, image = NULL)
  expect_s3_class(plot, "ggplot")
})

test_that("core spatial ggplot functions work without Seurat images", {

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

test_that("spatialDiscretePlot reads multisample metadata without barcode bleed", {

  make_section <- function(cell_type) {
    object <- make_image_free_spatial()
    object$cell_type <- cell_type
    object
  }
  sirup <- createSIRUP(
    spatial = list(
      section_a = make_section(c("B", "B", "T", "T")),
      section_b = make_section(c("T", "T", "B", "B"))
    ),
    receptors = data.frame(
      barcode = c("cell1", "cell1"),
      sample = c("section_a", "section_b"),
      clone_id = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample = "sample",
    clone = "clone_id",
    image_scale = "hires",
    swap_coords = FALSE
  )

  plots <- spatialDiscretePlot(sirup, meta_col = "cell_type", combine = FALSE)

  expect_equal(sort(names(plots)), c("section_a", "section_b"))
  expect_equal(as.character(unique(plots[["section_a"]]$data$cell_type)), c("B", "T"))
  expect_equal(as.character(unique(plots[["section_b"]]$data$cell_type)), c("T", "B"))

  section_b_plot <- spatialDiscretePlot(
    sirup,
    meta_col = "cell_type",
    sample_id = "section_b",
    combine = FALSE
  )
  expect_equal(names(section_b_plot), "section_b")
  expect_equal(as.character(unique(section_b_plot[["section_b"]]$data$cell_type)), c("T", "B"))
})

test_that("spatialExpressionPlot uses feature colors or palette", {
  expect_equal(length(.sirup_feature_color_values()), 256)
  expect_identical(
    .sirup_feature_color_values(feature_colors = c("black", "cyan"), feature_color_palette = "plasma"),
    c("black", "cyan")
  )
  expect_error(
    .sirup_feature_color_values(feature_color_palette = "not_a_palette"),
    "`feature_color_palette`"
  )
})

test_that(".sirup_fetch_data tries layer before legacy slot fallback", {
  calls <- character()
  testthat::local_mocked_bindings(
    .sirup_fetch_data_layer = function(object, vars, layer) {
      calls <<- c(calls, paste0("layer:", layer))
      data.frame(geneA = 1, row.names = "cell1")
    },
    .sirup_fetch_data_slot = function(object, vars, slot) {
      calls <<- c(calls, paste0("slot:", slot))
      data.frame(geneA = 2, row.names = "cell1")
    },
    .package = "SIRUP"
  )

  out <- .sirup_fetch_data(list(), vars = "geneA", layer = "counts")

  expect_equal(out$geneA, 1)
  expect_identical(calls, "layer:counts")
})

test_that(".sirup_fetch_data falls back to slot only after layer errors", {
  calls <- character()
  testthat::local_mocked_bindings(
    .sirup_fetch_data_layer = function(object, vars, layer) {
      calls <<- c(calls, paste0("layer:", layer))
      rlang::abort("layer failed")
    },
    .sirup_fetch_data_slot = function(object, vars, slot) {
      calls <<- c(calls, paste0("slot:", slot))
      data.frame(geneA = 2, row.names = "cell1")
    },
    .package = "SIRUP"
  )

  out <- .sirup_fetch_data(list(), vars = "geneA", layer = "counts")

  expect_equal(out$geneA, 2)
  expect_identical(calls, c("layer:counts", "slot:counts"))
})

test_that(".sirup_fetch_data reports ordered layer and slot failures", {
  testthat::local_mocked_bindings(
    .sirup_fetch_data_assay = function(object) "RNA",
    .sirup_fetch_data_layer = function(object, vars, layer) {
      rlang::abort("layer failed")
    },
    .sirup_fetch_data_slot = function(object, vars, slot) {
      rlang::abort("slot failed")
    },
    .package = "SIRUP"
  )

  expect_error(
    .sirup_fetch_data(list(), vars = "geneA", layer = "counts"),
    "SIRUP tried layer = counts in assay = RNA first, then tried slot = counts as a fall back"
  )
})

test_that("spatial statistic and polygon plots work without Seurat images", {

  object <- make_image_free_spatial()
  local_moran <- data.frame(
    cell = paste0("cell", 1:4),
    x = c(0, 1, 2, 3),
    y = c(0, 1, 1, 0),
    .sirup_y = c(0, -1, -1, 0),
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

polygon_layer_data <- function(plot) {
  polygon_layers <- which(vapply(plot$layers, function(l) inherits(l$geom, "GeomPolygon"), logical(1)))
  plot$layers[[polygon_layers[[1L]]]]$data
}

test_that("plotReceptors overlays polygons from a data frame and a CSV path", {

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

test_that("plotSpatialColorTips overlays SIRUP object polygons", {

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
  sirup <- createSIRUP(
    spatial = object,
    receptors = receptors,
    spatial_id = "spatial_id",
    locus = NULL,
    image_scale = "hires",
    swap_coords = FALSE,
    polygons = polygon
  )

  tip_plot <- plotSpatialColorTips(
    sirup,
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

test_that("stored polygon overlays are filtered per spatial sample", {

  make_section <- function() {
    make_image_free_spatial()
  }
  sirup <- createSIRUP(
    spatial = list(section_a = make_section(), section_b = make_section()),
    receptors = data.frame(
      spatial_id = c("cell1", "cell1"),
      sample = c("section_a", "section_b"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "spatial_id",
    sample = "sample",
    image_scale = "hires",
    swap_coords = FALSE
  )
  polygons <- data.frame(
    polygon = rep(c("a_gc", "b_gc"), each = 5),
    sample = rep(c("section_a", "section_b"), each = 5),
    x = rep(c(-0.5, 1.5, 1.5, -0.5, -0.5), 2),
    y = rep(c(-0.5, -0.5, 1.5, 1.5, -0.5), 2),
    stringsAsFactors = FALSE
  )
  sirup <- addPolygons(sirup, polygons = polygons)

  plots <- plotSpatialColorTips(
    sirup,
    tip_data = data.frame(
      spatial_id = c("cell1", "cell1"),
      sample = c("section_a", "section_b"),
      tip_color = c("#111111", "#CC5500"),
      stringsAsFactors = FALSE
    ),
    spatial_id = "spatial_id",
    sample = "sample",
    polygons = TRUE,
    combine = FALSE
  )

  expect_equal(names(plots), c("section_a", "section_b"))
  expect_equal(unique(as.character(polygon_layer_data(plots[["section_a"]])$polygon)), "a_gc")
  expect_equal(unique(as.character(polygon_layer_data(plots[["section_b"]])$polygon)), "b_gc")

  base_plots <- list(
    section_a = ggplot2::ggplot() + ggplot2::geom_point(data = data.frame(x = 0, y = 0), ggplot2::aes(.data$x, .data$y)),
    section_b = ggplot2::ggplot() + ggplot2::geom_point(data = data.frame(x = 0, y = 0), ggplot2::aes(.data$x, .data$y))
  )
  overlay <- SIRUP:::.sirup_resolve_overlay_polygons(
    polygons = TRUE,
    object = sirup,
    polygon_id = "polygon",
    polygon_val = NULL
  )
  finished <- SIRUP:::.sirup_add_polygon_overlay(base_plots, overlay)
  expect_equal(unique(as.character(polygon_layer_data(finished[["section_a"]])$polygon)), "a_gc")
  expect_equal(unique(as.character(polygon_layer_data(finished[["section_b"]])$polygon)), "b_gc")

  single_finished <- SIRUP:::.sirup_add_polygon_overlay(base_plots[["section_b"]], overlay, image = "section_b")
  expect_equal(unique(as.character(polygon_layer_data(single_finished)$polygon)), "b_gc")
})

test_that("spatialExpressionBlendPlot forwards assay and layer to per-feature panels", {

  object <- make_image_free_spatial()

  # A data-only assay, as produced by ALRA imputation: no counts layer.
  imputed <- Matrix::Matrix(
    c(
      0.5, 1.5, 0.2, 2.5,
      0.1, 0.9, 1.2, 3.5
    ),
    nrow = 2,
    ncol = 4,
    sparse = TRUE
  )
  dimnames(imputed) <- list(c("geneA", "geneB"), paste0("cell", 1:4))
  suppressWarnings(object[["alra"]] <- SeuratObject::CreateAssay5Object(data = imputed))

  # only_blended = FALSE builds the two per-feature panels. Those calls used to
  # drop `assay`/`layer` and fall back to spatialExpressionPlot's default
  # layer = "counts", which does not exist in a data-only assay. The existing
  # blend test only covers only_blended = TRUE, which skips this branch.
  expect_error(
    spatialExpressionBlendPlot(
      object,
      features = c("geneA", "geneB"),
      assay = "alra",
      layer = "data",
      only_blended = FALSE,
      image_scale = "hires",
      swap_coords = FALSE
    ),
    NA
  )
})

# Helper: pull the continuous colour values off the resolved receptor layer.
receptor_color_values <- function(plot) {
  layer <- Filter(function(l) inherits(l$geom, "GeomPoint"), plot$layers)[[1]]
  sort(as.numeric(layer$data$.sirup_count))
}

# Helper: the legend title lives on the colour scale, not in plot$labels.
receptor_color_scale_name <- function(plot) {
  scales <- Filter(function(s) "colour" %in% s$aesthetics, plot$scales$scales)
  scales[[1]]$name
}

test_that("plotReceptors colours spots by a continuous value column", {

  object <- make_image_free_spatial()
  # cell1 carries two receptors, so it must be combined to one number.
  receptors <- data.frame(
    spatial_id = c("cell1", "cell1", "cell2", "cell4"),
    receptor_type = c("BCR", "BCR", "BCR", "TCR"),
    shm = c(2, 8, 5, 1),
    stringsAsFactors = FALSE
  )

  common <- list(
    object, receptors = receptors, spatial_id = "spatial_id",
    image_scale = "hires", swap_coords = FALSE
  )

  mean_plot <- do.call(plotReceptors, c(common, list(value = "shm")))
  expect_s3_class(mean_plot, "ggplot")
  # cell1 = mean(2, 8) = 5, then cell2 = 5 and cell4 = 1
  expect_equal(receptor_color_values(mean_plot), c(1, 5, 5))

  max_plot <- do.call(plotReceptors, c(common, list(value = "shm", value_combine = "max")))
  expect_equal(receptor_color_values(max_plot), c(1, 5, 8))

  min_plot <- do.call(plotReceptors, c(common, list(value = "shm", value_combine = "min")))
  expect_equal(receptor_color_values(min_plot), c(1, 2, 5))

  median_plot <- do.call(plotReceptors, c(common, list(value = "shm", value_combine = "median")))
  expect_equal(receptor_color_values(median_plot), c(1, 5, 5))

  sum_plot <- do.call(plotReceptors, c(common, list(value = "shm", value_combine = "sum")))
  # cell1 = 2 + 8 = 10, cell2 = 5, cell4 = 1
  expect_equal(receptor_color_values(sum_plot), c(1, 5, 10))

  # mode takes the most frequent value, ties resolving to the smallest.
  tied <- data.frame(
    spatial_id = c("cell1", "cell1", "cell1", "cell2"),
    shm = c(3, 3, 9, 4),
    stringsAsFactors = FALSE
  )
  mode_plot <- plotReceptors(
    object, receptors = tied, spatial_id = "spatial_id",
    image_scale = "hires", swap_coords = FALSE,
    value = "shm", value_combine = "mode"
  )
  expect_equal(receptor_color_values(mode_plot), c(3, 4))

  # The legend takes the column name unless a title is given.
  expect_equal(receptor_color_scale_name(mean_plot), "shm")
  titled <- do.call(plotReceptors, c(common, list(value = "shm", receptor_legend_title = "SHM (%)")))
  expect_equal(receptor_color_scale_name(titled), "SHM (%)")

  # A group filter restricts which receptors are summarised.
  grouped <- do.call(plotReceptors, c(common, list(value = "shm", group = "receptor_type", group_val = "TCR")))
  expect_equal(receptor_color_values(grouped), 1)
})

test_that("plotReceptors value = 'count' counts receptors per spot", {

  object <- make_image_free_spatial()
  receptors <- data.frame(
    spatial_id = c("cell1", "cell1", "cell2", "cell4"),
    receptor_type = c("BCR", "BCR", "BCR", "TCR"),
    stringsAsFactors = FALSE
  )
  plot <- plotReceptors(
    object, receptors = receptors, spatial_id = "spatial_id",
    image_scale = "hires", swap_coords = FALSE, value = "count"
  )
  expect_s3_class(plot, "ggplot")
  expect_equal(receptor_color_values(plot), c(1, 1, 2))

  # A real `count` column does not override the receptor tally, but it warns.
  receptors$count <- c(10, 10, 20, 30)
  expect_warning(
    plot2 <- plotReceptors(
      object, receptors = receptors, spatial_id = "spatial_id",
      image_scale = "hires", swap_coords = FALSE, value = "count"
    ),
    "counts receptors per spot"
  )
  expect_equal(receptor_color_values(plot2), c(1, 1, 2))
})

test_that("plotReceptors validates value and accepts multi-stop palettes", {

  object <- make_image_free_spatial()
  receptors <- data.frame(
    spatial_id = c("cell1", "cell2"),
    receptor_type = c("BCR", "BCR"),
    shm = c(1, 2),
    label = c("a", "b"),
    stringsAsFactors = FALSE
  )
  base <- list(object, receptors = receptors, spatial_id = "spatial_id",
               image_scale = "hires", swap_coords = FALSE)

  expect_error(do.call(plotReceptors, c(base, list(value = "label"))), "must be numeric")
  expect_error(
    do.call(plotReceptors, c(base, list(value = "shm", value_combine = "nonsense"))),
    "should be one of"
  )
  expect_error(do.call(plotReceptors, c(base, list(value = "missing_col"))), "receptors")
  expect_error(do.call(plotReceptors, c(base, list(value = c("shm", "shm")))), "single column name")
  expect_error(
    do.call(plotReceptors, c(base, list(value = "shm", value_colors = "red"))),
    "at least two colours"
  )
  # Three or more colours switch to a multi-stop gradient.
  multi <- do.call(plotReceptors, c(base, list(
    value = "shm", value_colors = c("#440154", "#21918C", "#FDE725")
  )))
  expect_s3_class(multi, "ggplot")

  # `group` is still required when no value is given.
  expect_error(do.call(plotReceptors, base), "`group` is required")
})
