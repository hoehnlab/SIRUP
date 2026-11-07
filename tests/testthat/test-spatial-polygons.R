test_that("point-in-polygon includes boundary points when requested", {
  vx <- c(0, 10, 10, 0, 0)
  vy <- c(0, 0, 10, 10, 0)

  inside <- .sirup_points_in_polygon(
    px = c(5, 10, 12),
    py = c(5, 5, 5),
    vx = vx,
    vy = vy,
    include_boundary = TRUE
  )

  expect_equal(inside, c(TRUE, TRUE, FALSE))

  inside_no_boundary <- .sirup_points_in_polygon(
    px = c(5, 10),
    py = c(5, 5),
    vx = vx,
    vy = vy,
    include_boundary = FALSE
  )

  expect_equal(inside_no_boundary, c(TRUE, FALSE))
})

test_that("polygon assignment annotates spatial metadata and receptor tables", {
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
  polygon <- data.frame(
    polygon = c(rep("left", 5), rep("right", 5)),
    x = c(-0.5, 1.5, 1.5, -0.5, -0.5, 1.5, 3.5, 3.5, 1.5, 1.5),
    y = c(-0.5, -0.5, 1.5, 1.5, -0.5, 1.5, 1.5, 3.5, 3.5, 1.5),
    stringsAsFactors = FALSE
  )
  receptors <- data.frame(
    spatial_id = c("cell1", "cell2", "cell3", "missing"),
    sequence_id = paste0("seq", 1:4),
    stringsAsFactors = FALSE
  )

  spatial_out <- assignSpatialPolygons(
    spatial,
    polygon = polygon,
    image_scale = "hires",
    swap_coords = FALSE
  )
  receptors_out <- assignReceptorPolygons(
    receptors,
    object = spatial,
    polygon = polygon,
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_equal(unname(spatial_out$polygon_id), c("left", "left", "right", "right"))
  expect_equal(receptors_out$polygon_id, c("left", "left", "right", NA))
  expect_equal(unname(spatial_out$polygon_left), c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(unname(spatial_out$polygon_right), c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(receptors_out$polygon_left, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(receptors_out$polygon_right, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(ncol(spatial_out), ncol(spatial))
})

test_that("polygon assignment can use receptor coordinate columns directly", {
  polygon <- data.frame(
    polygon = c(rep("left", 5), rep("right", 5)),
    x = c(-0.5, 1.5, 1.5, -0.5, -0.5, 1.5, 3.5, 3.5, 1.5, 1.5),
    y = c(-0.5, -0.5, 1.5, 1.5, -0.5, 1.5, 1.5, 3.5, 3.5, 1.5),
    stringsAsFactors = FALSE
  )
  receptors <- data.frame(
    sequence_id = paste0("seq", 1:4),
    x = c(0, 1, 2, 10),
    y = c(0, 1, 2, 10),
    stringsAsFactors = FALSE
  )

  receptors_out <- assignReceptorPolygons(
    receptors,
    polygon = polygon
  )

  expect_equal(receptors_out$polygon_id, c("left", "left", "right", NA))
  expect_equal(receptors_out$polygon_left, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(receptors_out$polygon_right, c(FALSE, FALSE, TRUE, FALSE))
  expect_false(".sirup_receptor_row" %in% colnames(receptors_out))
})

test_that("overlapping polygon assignment keeps flags and chooses smallest polygon", {
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
  polygon <- data.frame(
    polygon = c(rep("outer", 5), rep("inner", 5)),
    x = c(-1, 4, 4, -1, -1, 0.5, 2.5, 2.5, 0.5, 0.5),
    y = c(-1, -1, 4, 4, -1, 0.5, 0.5, 2.5, 2.5, 0.5),
    stringsAsFactors = FALSE
  )
  receptors <- data.frame(
    spatial_id = c("cell1", "cell2", "cell3", "cell4", "missing"),
    sequence_id = paste0("seq", 1:5),
    stringsAsFactors = FALSE
  )

  spatial_out <- assignSpatialPolygons(
    spatial,
    polygon = polygon,
    image_scale = "hires",
    swap_coords = FALSE
  )
  receptors_out <- assignReceptorPolygons(
    receptors,
    object = spatial,
    polygon = polygon,
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_equal(unname(spatial_out$polygon_id), c("outer", "inner", "inner", "outer"))
  expect_equal(unname(spatial_out$polygon_outer), c(TRUE, TRUE, TRUE, TRUE))
  expect_equal(unname(spatial_out$polygon_inner), c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(receptors_out$polygon_id, c("outer", "inner", "inner", "outer", NA))
  expect_equal(receptors_out$polygon_outer, c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(receptors_out$polygon_inner, c(FALSE, TRUE, TRUE, FALSE, FALSE))
})

test_that("direct receptor coordinate polygon assignment keeps overlap flags", {
  polygon <- data.frame(
    polygon = c(rep("outer", 5), rep("inner", 5)),
    vertex_x = c(-1, 4, 4, -1, -1, 0.5, 2.5, 2.5, 0.5, 0.5),
    vertex_y = c(-1, -1, 4, 4, -1, 0.5, 0.5, 2.5, 2.5, 0.5),
    stringsAsFactors = FALSE
  )
  receptors <- data.frame(
    sequence_id = paste0("seq", 1:5),
    receptor_x = c(0, 1, 2, 3, 10),
    receptor_y = c(0, 1, 2, 3, 10),
    stringsAsFactors = FALSE
  )

  receptors_out <- assignReceptorPolygons(
    receptors,
    polygon = polygon,
    x = "vertex_x",
    y = "vertex_y",
    receptor_x = "receptor_x",
    receptor_y = "receptor_y"
  )

  expect_equal(receptors_out$polygon_id, c("outer", "inner", "inner", "outer", NA))
  expect_equal(receptors_out$polygon_outer, c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(receptors_out$polygon_inner, c(FALSE, TRUE, TRUE, FALSE, FALSE))
})

test_that("polygon vertex selection closes open polygons automatically", {
  polygon <- data.frame(
    x = c(0, 10, 10, 0),
    y = c(0, 0, 10, 10),
    stringsAsFactors = FALSE
  )

  closed <- .sirup_select_polygon_vertices(
    polygon = polygon,
    polygon_id = NULL,
    polygon_col = "polygon",
    x_col = "x",
    y_col = "y"
  )

  expect_equal(nrow(closed), 5)
  expect_equal(unname(unlist(closed[1, ])), unname(unlist(closed[5, ])))
})

test_that("polygon plotting coordinates can be downsampled reproducibly", {
  coords <- data.frame(
    cell = paste0("cell", 1:10),
    x = seq_len(10),
    y = seq_len(10),
    stringsAsFactors = FALSE
  )

  sampled_a <- .sirup_downsample_polygon_plot_coords(coords, downsample_n = 4, seed = 10)
  sampled_b <- .sirup_downsample_polygon_plot_coords(coords, downsample_n = 4, seed = 10)

  expect_equal(nrow(sampled_a), 4)
  expect_equal(sampled_a, sampled_b)
  expect_equal(.sirup_downsample_polygon_plot_coords(coords, downsample_n = 20), coords)
  expect_error(.sirup_downsample_polygon_plot_coords(coords, downsample_n = 0), "positive number")
})

test_that("selected polygon table displays only current polygon without image", {
  points <- data.frame(
    spatial_id = c("cell1", "cell2"),
    x = c(1, 2),
    y = c(1, 2),
    .sirup_plot_y = c(-1, -2),
    num = c(1L, 2L),
    polygon = c("1", "renamed"),
    image = c("slice1", "slice1"),
    stringsAsFactors = FALSE
  )

  current <- points[points$num == 2L, , drop = FALSE]
  shown <- .sirup_public_polygon_vertices(current)
  shown <- shown[, setdiff(colnames(shown), "image"), drop = FALSE]

  expect_equal(shown$spatial_id, "cell2")
  expect_equal(shown$num, 2L)
  expect_equal(shown$polygon, "renamed")
  expect_false("image" %in% names(shown))
  expect_false(".sirup_plot_y" %in% names(shown))
})

test_that("polygon closing uses num while keeping editable polygon names", {
  points <- data.frame(
    spatial_id = paste0("cell", 1:3),
    x = c(0, 1, 0),
    y = c(0, 0, 1),
    .sirup_plot_y = c(0, 0, -1),
    num = c(1L, 1L, 1L),
    polygon = c("custom", "custom", "custom"),
    image = c("slice1", "slice1", "slice1"),
    stringsAsFactors = FALSE
  )

  out <- .sirup_close_current_polygon(points, 1L)

  expect_equal(nrow(out), 4L)
  expect_equal(out$num, c(1L, 1L, 1L, 1L))
  expect_equal(out$polygon, rep("custom", 4L))
  expect_equal(out$spatial_id[[4L]], "cell1")
})

test_that("polygon plotting data preserves polygon groups and flips plotting y", {
  polygon <- rbind(
    data.frame(
      polygon = "left",
      x = c(0, 5, 5, 0, 0),
      y = c(0, 0, 5, 5, 0),
      stringsAsFactors = FALSE
    ),
    data.frame(
      polygon = "right",
      x = c(10, 15, 15, 10, 10),
      y = c(10, 10, 15, 15, 10),
      stringsAsFactors = FALSE
    )
  )

  plot_data <- .sirup_polygon_plot_data(
    polygon = polygon,
    polygon_id = NULL,
    polygon_col = "polygon",
    x_col = "x",
    y_col = "y"
  )

  expect_equal(sort(unique(plot_data$polygon)), c("left", "right"))
  expect_equal(plot_data$.sirup_y, -plot_data$y)
})

test_that("plotSpatialPolygons colors stored polygons by an analysis column", {
  spatial <- sirup_test_seurat(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 4,
      dimnames = list("gene", paste0("cell", 1:4))
    )
  )
  spatial$x <- 0:3
  spatial$y <- 0:3
  x <- createSIRUP(
    spatial = spatial,
    receptors = data.frame(spatial_id = paste0("cell", 1:4)),
    spatial_id = "spatial_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
  x <- addPolygons(x, polygons = tibble::tibble(
    polygon_id = c("left", "right"),
    polygon_score = c(3, 9),
    vertices = list(
      tibble::tibble(x = c(0, 1, 1, 0, 0), y = c(0, 0, 1, 1, 0)),
      tibble::tibble(x = c(2, 3, 3, 2, 2), y = c(2, 2, 3, 3, 2))
    )
  ))
  x <- setAnalysis(x, "gc_metrics", data.frame(
    polygon_id = c("left", "right"),
    effective_shannon = c(2, 8)
  ))

  ns <- environment(plotSpatialPolygons)
  old_images <- get(".sirup_spatial_images", envir = ns)
  old_data <- get(".sirup_spatial_data", envir = ns)
  old_base <- get(".sirup_spatial_base", envir = ns)
  locked <- vapply(c(".sirup_spatial_images", ".sirup_spatial_data", ".sirup_spatial_base"),
                   bindingIsLocked, logical(1), env = ns)
  for (name in names(locked)[locked]) unlockBinding(name, ns)
  assign(".sirup_spatial_images", function(object, images = NULL, sample_id = NULL) "metadata", envir = ns)
  assign(".sirup_spatial_data", function(object, image, image_scale = NULL, swap_coords = NULL) {
    list(coords = data.frame(cell = paste0("cell", 1:4), x = 0:3, y = 0:3, .sirup_y = 0:-3), image = NULL)
  }, envir = ns)
  assign(".sirup_spatial_base", function(plot_data, image, image_alpha = 1, zoom = FALSE, view_data = NULL) {
    ggplot2::ggplot()
  }, envir = ns)
  withr::defer({
    assign(".sirup_spatial_images", old_images, envir = ns)
    assign(".sirup_spatial_data", old_data, envir = ns)
    assign(".sirup_spatial_base", old_base, envir = ns)
    for (name in names(locked)[locked]) lockBinding(name, ns)
  })

  plot <- plotSpatialPolygons(
    x,
    fill_analysis = "gc_metrics",
    fill_column = "effective_shannon"
  )
  expect_s3_class(plot, "ggplot")
  expect_true(any(vapply(plot$layers, function(layer) "fill" %in% names(layer$mapping), logical(1))))
  fill_scale <- plot$scales$get_scales("fill")
  expect_s3_class(fill_scale, "ScaleContinuous")
  expect_identical(unname(fill_scale$palette(c(0, 1))), c("#D3D3D3", "#FFFF00"))

  custom_plot <- plotSpatialPolygons(
    x,
    fill_analysis = "gc_metrics",
    fill_column = "effective_shannon",
    analysis_fill = c("black", "cyan")
  )
  custom_scale <- custom_plot$scales$get_scales("fill")
  expect_identical(unname(custom_scale$palette(c(0, 1))), c("#000000", "#00FFFF"))

  polygon_fill_plot <- plotSpatialPolygons(
    x,
    fill_column = "polygon_score"
  )
  polygon_fill_scale <- polygon_fill_plot$scales$get_scales("fill")
  expect_s3_class(polygon_fill_scale, "ScaleContinuous")
  expect_identical(unname(polygon_fill_scale$palette(c(0, 1))), c("#D3D3D3", "#FFFF00"))

  expect_error(
    plotSpatialPolygons(x, fill_column = "not_in_polygons"),
    "was not found in polygons.*fill_analysis"
  )

  x <- setAnalysis(x, "gc_classes", data.frame(
    polygon_id = c("left", "right"),
    class = c("A", "B")
  ))
  discrete_plot <- plotSpatialPolygons(
    x,
    fill_analysis = "gc_classes",
    fill_column = "class"
  )
  discrete_scale <- discrete_plot$scales$get_scales("fill")
  expect_equal(
    unname(discrete_scale$palette(2)),
    pals::glasbey(2)
  )
})

test_that("plotSpatialPolygons can place labels outside polygons", {
  spatial <- sirup_test_seurat(
    counts = matrix(
      1,
      nrow = 1,
      ncol = 4,
      dimnames = list("gene", paste0("cell", 1:4))
    )
  )
  spatial$x <- 0:3
  spatial$y <- 0:3
  x <- createSIRUP(
    spatial = spatial,
    receptors = data.frame(spatial_id = paste0("cell", 1:4)),
    spatial_id = "spatial_id",
    image_scale = "hires",
    swap_coords = FALSE
  )
  x <- addPolygons(x, polygons = tibble::tibble(
    polygon_id = c("left", "right"),
    vertices = list(
      tibble::tibble(x = c(0, 1, 1, 0, 0), y = c(0, 0, 1, 1, 0)),
      tibble::tibble(x = c(4, 5, 5, 4, 4), y = c(0, 0, 1, 1, 0))
    )
  ))

  ns <- environment(plotSpatialPolygons)
  old_images <- get(".sirup_spatial_images", envir = ns)
  old_data <- get(".sirup_spatial_data", envir = ns)
  old_base <- get(".sirup_spatial_base", envir = ns)
  locked <- vapply(c(".sirup_spatial_images", ".sirup_spatial_data", ".sirup_spatial_base"),
                   bindingIsLocked, logical(1), env = ns)
  for (name in names(locked)[locked]) unlockBinding(name, ns)
  assign(".sirup_spatial_images", function(object, images = NULL, sample_id = NULL) "metadata", envir = ns)
  assign(".sirup_spatial_data", function(object, image, image_scale = NULL, swap_coords = NULL) {
    list(coords = data.frame(cell = paste0("cell", 1:4), x = 0:3, y = 0:3, .sirup_y = 0:-3), image = NULL)
  }, envir = ns)
  assign(".sirup_spatial_base", function(plot_data, image, image_alpha = 1, zoom = FALSE, view_data = NULL) {
    ggplot2::ggplot()
  }, envir = ns)
  withr::defer({
    assign(".sirup_spatial_images", old_images, envir = ns)
    assign(".sirup_spatial_data", old_data, envir = ns)
    assign(".sirup_spatial_base", old_base, envir = ns)
    for (name in names(locked)[locked]) lockBinding(name, ns)
  })

  plot <- plotSpatialPolygons(x, label_outside = TRUE, label_offset = 0.1)

  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$layers[[2]]$geom, "GeomSegment")
  expect_s3_class(plot$layers[[3]]$geom, "GeomText")
  label_data <- plot$layers[[3]]$data
  expect_true(all(c(".sirup_lx", ".sirup_ly", ".sirup_ex", ".sirup_ey") %in% names(label_data)))
  expect_false(isTRUE(all.equal(label_data$x, label_data$.sirup_lx)))
  expect_error(plotSpatialPolygons(x, label_outside = NA), "`label_outside` must be TRUE or FALSE")
  expect_error(plotSpatialPolygons(x, label_leader = NA), "`label_leader` must be TRUE or FALSE")
  expect_error(plotSpatialPolygons(x, label_outside = TRUE, label_offset = -0.1), "`label_offset` must be a non-negative number")
})

test_that("V gene helpers retain alleles and extract families and superfamilies", {
  genes <- .sirup_v_genes(c("IGHV3-23*01", "IGHV3-30*02,IGHV3-23*01"))
  families <- .sirup_v_gene_family(c(genes, "IGKV1D-39*01"))

  expect_equal(genes, c("IGHV3-23*01", "IGHV3-30*02"))
  expect_equal(families, c("IGHV3-23", "IGHV3-30", "IGKV1D-39"))
  expect_equal(.sirup_v_gene_superfamily(families), c("IGHV3", "IGHV3", "IGKV1"))
})

test_that("polygon pie summaries group lower-frequency levels and label top levels", {
  region_receptors <- data.frame(
    polygon = c(rep("poly1", 5), rep("poly2", 3), rep("poly3", 3)),
    clone_id = c("A", "A", "B", "C", "D", "E", "F", "F", "C", "C", "G"),
    stringsAsFactors = FALSE
  )

  pies <- .sirup_polygon_pie_summary(
    region_receptors,
    fill = "clone_id",
    top_n = 2,
    label_top_n = 1,
    other_label = "Other"
  )

  poly1 <- pies[pies$polygon == "poly1", , drop = FALSE]
  expect_true("Other" %in% poly1$level)
  expect_true("C" %in% poly1$level)
  expect_equal(poly1$level, c("A", "B", "C", "Other"))
  expect_equal(poly1$n, c(2, 1, 1, 1))
  expect_equal(sum(poly1$proportion), 1)
  expect_equal(sum(poly1$label != ""), 1)
  expect_equal(sort(unique(pies$polygon_label)), c("poly1", "poly2", "poly3"))
  expect_false(".sirup_order" %in% names(pies))

  pies_without_other <- .sirup_polygon_pie_summary(
    region_receptors,
    fill = "clone_id",
    top_n = 1,
    label_top_n = 1,
    other_label = "Other",
    plot_other = FALSE
  )
  expect_false("Other" %in% pies_without_other$level)
  expect_equal(pies_without_other$level[pies_without_other$polygon == "poly1"], c("A", "C"))
  expect_equal(pies_without_other$level[pies_without_other$polygon == "poly2"], "F")
  expect_equal(pies_without_other$level[pies_without_other$polygon == "poly3"], "C")
})

test_that("discrete palette uses pals glasbey and configured other color", {
  values <- .sirup_discrete_palette_values(
    c("A", "B", "Other"),
    other_label = "Other",
    other_color = "#EEEEEE"
  )

  expect_equal(unname(values[c("A", "B")]), pals::glasbey(2))
  expect_equal(unname(values["Other"]), "#EEEEEE")

  many_levels <- paste0("clone", seq_len(46))
  many_values <- .sirup_discrete_palette_values(
    c(many_levels, "Other"),
    other_label = "Other",
    other_color = "#EEEEEE"
  )
  expect_equal(names(many_values), c(many_levels, "Other"))
  expect_equal(length(many_values), 47)
  expect_equal(unname(many_values["Other"]), "#EEEEEE")
})

test_that("spatial pie sectors start at the visual top", {
  sector <- .sirup_pie_sector(center_x = 10, center_y = 20, radius = 5, ymin = 0, ymax = 0.25)

  expect_equal(sector$x[2], 10, tolerance = 1e-8)
  expect_equal(sector$y[2], 15, tolerance = 1e-8)
})

test_that("3D clone density data groups clones and builds surfaces", {
  ns <- environment(spatialCloneDominance)
  old_image_dim <- get(".sirup_image_dim", envir = ns)
  image_dim_binding_locked <- bindingIsLocked(".sirup_image_dim", ns)
  if (image_dim_binding_locked) {
    unlockBinding(".sirup_image_dim", ns)
  }
  assign(
    ".sirup_image_dim",
    function(image) {
      list(height = 100, width = 100)
    },
    envir = ns
  )
  if (image_dim_binding_locked) {
    lockBinding(".sirup_image_dim", ns)
  }
  withr::defer({
    if (image_dim_binding_locked) {
      unlockBinding(".sirup_image_dim", ns)
    }
    assign(".sirup_image_dim", old_image_dim, envir = ns)
    if (image_dim_binding_locked) {
      lockBinding(".sirup_image_dim", ns)
    }
  })

  region_receptors <- data.frame(
    clone_id = c("A", "A", "B", "C"),
    x = c(20, 22, 70, 72),
    y = c(20, 22, 70, 72),
    .sirup_y = -c(20, 22, 70, 72),
    polygon = "poly1",
    image = "metadata",
    stringsAsFactors = FALSE
  )
  polygons <- data.frame(
    x = c(10, 80, 80, 10, 10),
    y = c(10, 10, 80, 80, 10),
    .sirup_y = -c(10, 10, 80, 80, 10),
    polygon = "poly1",
    stringsAsFactors = FALSE
  )

  out <- .sirup_clone_density_3d_data(
    region_receptors = region_receptors,
    polygons = polygons,
    image = list(),
    spatial_coords = data.frame(
      x = c(0, 100),
      .sirup_y = c(-100, 0)
    ),
    clone = "clone_id",
    weight = NULL,
    top_n = 1,
    other_label = "Other",
    other_color = "#EEEEEE",
    plot_other = TRUE,
    bandwidth = 8,
    grid_n = 20,
    plot_type = "density",
    zoom = TRUE,
    height_scale = 0.2,
    z_offset = 0.05,
    density_floor = 0.03
  )

  expect_equal(names(out$surfaces), c("A", "Other"))
  expect_equal(dim(out$surfaces$A$z), c(20, 20))
  expect_true(max(out$surfaces$A$z, na.rm = TRUE) > out$z_base)
  expect_equal(unname(out$colors["Other"]), "#EEEEEE")
  expect_true(all(out$image_extent$xlim >= 0 & out$image_extent$xlim <= 100))
  expect_true(all(out$image_extent$ylim >= -100 & out$image_extent$ylim <= 0))
})

test_that("polygon plots filter tagged polygons per spatial sample", {

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
      locus = "IGH",
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample = "sample",
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

  plots <- plotSpatialPolygons(x, combine = FALSE)

  expect_equal(unique(as.character(plots[["section_a"]]$data$polygon)), "a_gc")
  expect_equal(unique(as.character(plots[["section_b"]]$data$polygon)), "b_gc")

  section_b_plot <- plotSpatialPolygons(x, sample_id = "section_b", combine = FALSE)
  expect_equal(names(section_b_plot), "section_b")
  expect_equal(unique(as.character(section_b_plot[["section_b"]]$data$polygon)), "b_gc")
})

test_that("polygon pies filter tagged polygons per spatial sample", {

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
      locus = "IGH",
      stringsAsFactors = FALSE
    ),
    spatial_id = "barcode",
    sample = "sample",
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

  side_by_side <- plotSpatialPolygonPies(dominance, x)
  expect_equal(names(side_by_side), c("section_a", "section_b"))
  expect_false(is.null(side_by_side[["section_a"]]))
  expect_false(is.null(side_by_side[["section_b"]]))

  side_by_side_b <- plotSpatialPolygonPies(dominance, x, sample_id = "section_b")
  expect_false(is.null(side_by_side_b))

  plot_a <- plotSpatialPolygonPies(dominance, x, mode = "spatial_clipped", images = "section_a")
  plot_b <- plotSpatialPolygonPies(dominance, x, mode = "spatial_clipped", images = "section_b")

  expect_equal(unique(as.character(plot_a$layers[[1]]$data$polygon)), "a_gc")
  expect_equal(unique(as.character(plot_b$layers[[1]]$data$polygon)), "b_gc")
})

test_that("spatialCloneDominance counts clone proportions inside a polygon", {
  ns <- environment(spatialCloneDominance)
  old_coordinate_table <- get(".sirup_spatial_coordinate_table", envir = ns)
  coordinate_binding_locked <- bindingIsLocked(".sirup_spatial_coordinate_table", ns)
  if (coordinate_binding_locked) {
    unlockBinding(".sirup_spatial_coordinate_table", ns)
  }
  assign(
    ".sirup_spatial_coordinate_table",
    function(object, images = NULL, image_scale = NULL, swap_coords = NULL) {
      data.frame(
        cell = paste0("cell", 1:5),
        x = c(1, 2, 8, 12, 1),
        y = c(1, 2, 8, 5, 12),
        image = "metadata",
        stringsAsFactors = FALSE
      )
    },
    envir = ns
  )
  if (coordinate_binding_locked) {
    lockBinding(".sirup_spatial_coordinate_table", ns)
  }
  withr::defer({
    if (coordinate_binding_locked) {
      unlockBinding(".sirup_spatial_coordinate_table", ns)
    }
    assign(".sirup_spatial_coordinate_table", old_coordinate_table, envir = ns)
    if (coordinate_binding_locked) {
      lockBinding(".sirup_spatial_coordinate_table", ns)
    }
  })

  receptors <- data.frame(
    barcode = c("cell1", "cell1", "cell2", "cell3", "cell4", "cell5"),
    clone_id = c("A", "A", "B", "A", "C", "D"),
    locus = c("IGH", "IGK", "IGH", "IGH", "IGH", "IGH"),
    v_call = c("IGHV3-23*01", "IGKV1-5*01", "IGHV4-39*01", "IGHV3-23*02", "IGHV1-69*01", "IGHV2-5*01"),
    x = seq(101, 106),
    y = seq(201, 206),
    stringsAsFactors = FALSE
  )
  polygon <- data.frame(
    x = c(0, 10, 10, 0, 0),
    y = c(0, 0, 10, 10, 0),
    barcode = c("v1", "v2", "v3", "v4", "v1"),
    stringsAsFactors = FALSE
  )

  out <- spatialCloneDominance(
    receptors = receptors,
    object = list(),
    polygon = polygon,
    spatial_id = "barcode",
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_s3_class(out, "sirup_spatial_clone_dominance")
  expect_equal(out$dominance$clone_id, c("A", "B"))
  expect_equal(out$dominance$v_gene[out$dominance$clone_id == "A"], "IGHV3-23*01,IGHV3-23*02,IGKV1-5*01")
  expect_equal(out$dominance$v_gene_family[out$dominance$clone_id == "A"], "IGHV3-23,IGKV1-5")
  expect_equal(out$dominance$v_gene_superfamily[out$dominance$clone_id == "A"], "IGHV3,IGKV1")
  expect_equal(out$dominance$n_receptors, c(3, 1))
  expect_equal(out$dominance$proportion, c(0.75, 0.25))
  expect_equal(unique(out$dominance$region_receptor_total), 4)
  expect_equal(unique(out$dominance$n_region_spatial_ids), 3)
  expect_equal(sort(unique(out$region_receptors$.sirup_spatial_id)), c("cell1", "cell2", "cell3"))
  expect_true(all(c("x", "y", ".sirup_y", "polygon") %in% names(out$region_receptors)))
  expect_false(any(c("x.x", "x.y", "y.x", "y.y") %in% names(out$region_receptors)))
  expect_equal(sort(unique(out$region_receptors$x)), c(1, 2, 8))
  expect_equal(sort(unique(out$region_receptors$y)), c(1, 2, 8))
  expect_true(all(c("x", "y", ".sirup_y", "polygon") %in% names(out$polygons)))
  expect_false(any(c("polygon_col", "spatial_id_col", "clone_col", "v_call_col") %in% names(out)))
  expect_equal(out$metadata$receptor_fields$spatial_id, "barcode")
  expect_equal(out$metadata$receptor_fields$clone, "clone_id")
  expect_equal(out$metadata$receptor_fields$locus, "locus")
  expect_null(out$metadata$receptor_fields$count)
  expect_equal(out$metadata$receptor_fields$v_call, "v_call")
  expect_equal(names(formals(spatialCloneDominance))[1], "receptors")

  raw_out <- spatialCloneDominance(
    receptors = receptors,
    object = list(),
    polygon = polygon,
    spatial_id = "barcode",
    image_scale = "hires",
    swap_coords = FALSE
  )
  expect_equal(raw_out$dominance, out$dominance)

  igh <- spatialCloneDominance(
    receptors = receptors,
    object = list(),
    polygon = polygon,
    spatial_id = "barcode",
    locus_val = "IGH",
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_equal(igh$dominance$clone_id, c("A", "B"))
  expect_equal(igh$dominance$n_receptors, c(2, 1))
  expect_equal(igh$dominance$proportion, c(2 / 3, 1 / 3))

  no_singletons <- spatialCloneDominance(
    receptors = receptors,
    object = list(),
    polygon = polygon,
    spatial_id = "barcode",
    exclude_singletons = TRUE,
    image_scale = "hires",
    swap_coords = FALSE
  )
  expect_equal(no_singletons$dominance$clone_id, "A")
  expect_equal(no_singletons$dominance$n_receptors, 3)
  expect_equal(unique(no_singletons$region_receptors$clone_id), "A")

  region <- .sirup_region_from_dominance(out)
  expect_equal(region$receptors, out$region_receptors)
  expect_equal(region$polygons, out$polygons)
})

test_that("regional plotting APIs require dominance and spatial objects", {
  region_args <- names(formals(plotSpatialRegionReceptors))
  pie_args <- names(formals(plotSpatialPolygonPies))
  density_args <- names(formals(plotSpatialCloneDensity3D))
  redundant_args <- c(
    "receptors", "bcrs", "polygon", "spatial_id_col", "x_col", "y_col",
    "polygon_col", "require_closed", "color_by", "fill_by"
  )

  expect_equal(region_args[1:2], c("dominance", "object"))
  expect_equal(pie_args[1:2], c("dominance", "object"))
  expect_equal(density_args[1:2], c("dominance", "object"))
  expect_false(any(redundant_args %in% region_args))
  expect_false(any(redundant_args %in% pie_args))
  expect_error(
    plotSpatialRegionReceptors(data.frame(), object = list()),
    "returned by `spatialCloneDominance\\(\\)`"
  )
  expect_error(
    plotSpatialPolygonPies(data.frame(), object = list()),
    "returned by `spatialCloneDominance\\(\\)`"
  )
})

test_that("spatialCloneDominance summarizes multiple polygons independently", {
  ns <- environment(spatialCloneDominance)
  old_coordinate_table <- get(".sirup_spatial_coordinate_table", envir = ns)
  coordinate_binding_locked <- bindingIsLocked(".sirup_spatial_coordinate_table", ns)
  if (coordinate_binding_locked) {
    unlockBinding(".sirup_spatial_coordinate_table", ns)
  }
  assign(
    ".sirup_spatial_coordinate_table",
    function(object, images = NULL, image_scale = NULL, swap_coords = NULL) {
      data.frame(
        cell = paste0("cell", 1:5),
        x = c(1, 2, 8, 12, 13),
        y = c(1, 2, 8, 12, 13),
        image = "metadata",
        stringsAsFactors = FALSE
      )
    },
    envir = ns
  )
  if (coordinate_binding_locked) {
    lockBinding(".sirup_spatial_coordinate_table", ns)
  }
  withr::defer({
    if (coordinate_binding_locked) {
      unlockBinding(".sirup_spatial_coordinate_table", ns)
    }
    assign(".sirup_spatial_coordinate_table", old_coordinate_table, envir = ns)
    if (coordinate_binding_locked) {
      lockBinding(".sirup_spatial_coordinate_table", ns)
    }
  })

  receptors <- data.frame(
    barcode = c("cell1", "cell2", "cell3", "cell4", "cell5"),
    clone_id = c("A", "B", "C", "D", "D"),
    stringsAsFactors = FALSE
  )
  polygon <- rbind(
    data.frame(
      polygon = "left",
      x = c(0, 5, 5, 0, 0),
      y = c(0, 0, 5, 5, 0),
      stringsAsFactors = FALSE
    ),
    data.frame(
      polygon = "right",
      x = c(10, 15, 15, 10, 10),
      y = c(10, 10, 15, 15, 10),
      stringsAsFactors = FALSE
    )
  )

  out <- spatialCloneDominance(
    receptors = receptors,
    object = list(),
    polygon = polygon,
    spatial_id = "barcode",
    image_scale = "hires",
    swap_coords = FALSE
  )

  expect_equal(sort(unique(out$dominance$polygon)), c("left", "right"))
  expect_equal(out$dominance$n_receptors[out$dominance$polygon == "left" & out$dominance$clone_id == "A"], 1)
  expect_equal(out$dominance$n_receptors[out$dominance$polygon == "left" & out$dominance$clone_id == "B"], 1)
  expect_equal(out$dominance$n_receptors[out$dominance$polygon == "right" & out$dominance$clone_id == "D"], 2)
  expect_equal(out$dominance$proportion[out$dominance$polygon == "right" & out$dominance$clone_id == "D"], 1)

  right <- spatialCloneDominance(
    receptors = receptors,
    object = list(),
    polygon = polygon,
    spatial_id = "barcode",
    polygon_val = "right",
    image_scale = "hires",
    swap_coords = FALSE
  )
  expect_equal(unique(right$dominance$polygon), "right")
  expect_equal(right$dominance$clone_id, "D")
})

test_that(".sirup_spatial_base crops to spatial barcodes by default", {
  plot_data <- data.frame(x = c(10, 20), .sirup_y = c(-10, -20))
  spots <- data.frame(x = c(0, 100), .sirup_y = c(0, -100))

  # crop = TRUE with spot coordinates as view_data -> limits follow the spots.
  p_crop <- .sirup_spatial_base(plot_data, image = NULL, crop = TRUE, view_data = spots)
  xlim <- p_crop$coordinates$limits$x
  ylim <- p_crop$coordinates$limits$y
  expect_true(xlim[[1]] <= 0 && xlim[[2]] >= 100)
  expect_true(ylim[[1]] <= -100 && ylim[[2]] >= 0)

  # crop = FALSE with no image -> default unit box, not the spot extent.
  p_full <- .sirup_spatial_base(plot_data, image = NULL, crop = FALSE, view_data = spots)
  expect_equal(p_full$coordinates$limits$x, c(0, 1))

  # An empty view_data is treated as absent and falls back to plot_data.
  p_empty <- .sirup_spatial_base(
    plot_data,
    image = NULL,
    crop = TRUE,
    view_data = data.frame(x = numeric(), .sirup_y = numeric())
  )
  xr <- p_empty$coordinates$limits$x
  expect_true(xr[[1]] <= 10 && xr[[2]] >= 20)
})

test_that("side-by-side pies order polygons by receptor count", {
  # Polygon 1 = 12 receptors, polygon 2 = 70, polygon 3 = 16. The default
  # "receptors" order shows the most-populated polygon first; each pie keeps
  # its map label. "polygon" keeps the numbered map order. The helper returns
  # stable polygon IDs, not display labels.
  pie_data <- data.frame(
    polygon = c("GC1", "GC2", "GC3"),
    polygon_label = c("1", "2", "3"),
    n = c(12L, 70L, 16L),
    stringsAsFactors = FALSE
  )

  expect_equal(
    .sirup_pie_polygon_order(pie_data, "receptors"),
    c("GC2", "GC3", "GC1")
  )
  expect_equal(
    .sirup_pie_polygon_order(pie_data, "polygon"),
    c("GC1", "GC2", "GC3")
  )

  # Several clones per polygon: the totals per polygon drive the order.
  multi <- data.frame(
    polygon = c("GC1", "GC1", "GC2", "GC3", "GC3"),
    polygon_label = c("1", "1", "2", "3", "3"),
    n = c(5L, 7L, 70L, 10L, 6L),
    stringsAsFactors = FALSE
  )
  expect_equal(.sirup_pie_polygon_order(multi, "receptors"), c("GC2", "GC3", "GC1"))

  # A tie in receptor count breaks by the polygon label, ascending.
  tied <- data.frame(
    polygon = c("GC1", "GC2", "GC3"),
    polygon_label = c("1", "2", "3"),
    n = c(20L, 20L, 5L),
    stringsAsFactors = FALSE
  )
  expect_equal(.sirup_pie_polygon_order(tied, "receptors"), c("GC1", "GC2", "GC3"))
})

test_that("stacked polygon dominance bars use clean polygon labels", {

  region_bcrs <- data.frame(
    polygon = c("section_a::1", "section_a::1", "section_a::2"),
    polygon_label = c("1", "1", "2"),
    clone_id = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )

  pie_data <- .sirup_polygon_pie_summary(
    region_bcrs,
    fill = "clone_id",
    top_n = NULL,
    label_top_n = 1L,
    other_label = "Other"
  )
  expect_equal(unique(pie_data$polygon), c("section_a::1", "section_a::2"))
  expect_equal(unique(pie_data$polygon_label), c("1", "2"))

  p <- .sirup_stacked_bar_polygon_dominance(
    pie_data,
    fill_values = c(A = "#111111", B = "#222222"),
    label_size = 3,
    pie_order = "polygon"
  )
  expect_equal(p$scales$get_scales("x")$labels, c("1", "2"))
})

test_that("stacked polygon dominance bars keep same-numbered polygons distinct", {

  region_bcrs <- data.frame(
    polygon = c("section_a::5", "section_b::5"),
    polygon_label = c("5", "5"),
    clone_id = c("A", "B"),
    stringsAsFactors = FALSE
  )

  pie_data <- .sirup_polygon_pie_summary(
    region_bcrs,
    fill = "clone_id",
    top_n = NULL,
    label_top_n = 1L,
    other_label = "Other"
  )
  p <- .sirup_stacked_bar_polygon_dominance(
    pie_data,
    fill_values = c(A = "#111111", B = "#222222"),
    label_size = 3,
    pie_order = "polygon"
  )

  x_scale <- p$scales$get_scales("x")
  expect_equal(x_scale$breaks, 1:2)
  expect_equal(x_scale$labels, c("5", "5"))
})

test_that("expandPolygons moves vertices outward from the centroid", {
  # Closed unit square centered at the origin. Centroid is (0, 0), so each
  # vertex is sqrt(2) from the centroid. With delta = 1, each vertex moves to
  # distance sqrt(2) + 1, i.e. it scales by (sqrt(2) + 1) / sqrt(2).
  square <- data.frame(
    polygon = "A",
    x = c(-1, 1, 1, -1, -1),
    y = c(-1, -1, 1, 1, -1),
    stringsAsFactors = FALSE
  )
  expanded <- expandPolygons(square, delta = 1)

  scale <- (sqrt(2) + 1) / sqrt(2)
  expect_equal(expanded$x, square$x * scale)
  expect_equal(expanded$y, square$y * scale)

  # Shape is kept: same vertex count and the closing vertex still matches.
  last <- nrow(expanded)
  expect_equal(nrow(expanded), nrow(square))
  expect_equal(c(expanded$x[[1]], expanded$y[[1]]), c(expanded$x[[last]], expanded$y[[last]]))

  # Expanding grows the area; a negative delta shrinks it.
  area <- function(p) .sirup_polygon_area(p$x, p$y)
  expect_gt(area(expanded), area(square))
  expect_lt(area(expandPolygons(square, delta = -0.2)), area(square))
})

test_that("expandPolygons uses the area centroid and needs a delta", {
  # An L-free convex quad whose vertex mean differs from its area centroid.
  rect <- data.frame(x = c(0, 4, 4, 0), y = c(0, 0, 2, 2), stringsAsFactors = FALSE)
  centroid <- .sirup_polygon_centroid(rect$x, rect$y)
  expect_equal(centroid$x, 2)
  expect_equal(centroid$y, 1)

  # `delta` is required (no default).
  expect_error(expandPolygons(rect), "delta")
  expect_error(expandPolygons(rect, delta = c(1, 2)), "single finite number")
})

test_that("expandPolygons handles multiple polygons and SIRUP objects", {

  # Each polygon expands independently around its own centroid.
  multi <- data.frame(
    polygon = c("A", "A", "A", "B", "B", "B"),
    x = c(0, 2, 1, 10, 12, 11),
    y = c(0, 0, 2, 0, 0, 2),
    stringsAsFactors = FALSE
  )
  out <- expandPolygons(multi, delta = 0.5)
  expect_equal(nrow(out), nrow(multi))
  expect_equal(as.character(out$polygon), as.character(multi$polygon))

  # SIRUP object input returns a valid, expanded object.
  object <- sirup_test_seurat(
    counts = matrix(1, 1, 4, dimnames = list("gene", paste0("cell", 1:4)))
  )
  object$x <- c(0, 1, 2, 3)
  object$y <- c(0, 1, 2, 3)
  receptors <- data.frame(
    spatial_id = paste0("cell", 1:4),
    sequence_id = paste0("seq", 1:4),
    stringsAsFactors = FALSE
  )
  polygon <- data.frame(
    polygon_id = "r",
    x = c(0.5, 2.5, 2.5, 0.5, 0.5),
    y = c(0.5, 0.5, 2.5, 2.5, 0.5),
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
  grown <- expandPolygons(sirup, delta = 1)
  expect_true(isSIRUP(grown))
  before <- getPolygons(sirup)$vertices[[1]]
  after <- getPolygons(grown)$vertices[[1]]
  expect_gt(diff(range(after$x)), diff(range(before$x)))
})

test_that(".sirup_polygon_feature_values labels source and guards non-numeric", {
  counts <- matrix(c(0, 1, 2, 3, 4, 5, 6, 7),
                   nrow = 2,
                   dimnames = list(c("GeneA", "GeneB"), paste0("cell", 1:4)))
  obj <- sirup_test_seurat(counts = counts)
  gene <- "GeneA"

  # assay hit keeps its source label
  from_assay <- .sirup_polygon_feature_values(obj, gene, layer = "counts")
  expect_false(is.null(from_assay))
  expect_identical(attr(from_assay, "source"), "assay")
  expect_identical(colnames(from_assay), gene)

  # a numeric metadata column (a module score) is reported as metadata, not
  # as an assay hit -- FetchData resolves both, so only the label distinguishes
  # a real gene from a same-named metadata column
  obj$my_module_score <- seq_len(ncol(obj)) / ncol(obj)
  from_meta <- .sirup_polygon_feature_values(obj, "my_module_score", layer = "counts")
  expect_false(is.null(from_meta))
  expect_identical(attr(from_meta, "source"), "metadata")
  expect_identical(rownames(from_meta), colnames(obj))
  expect_equal(from_meta[["my_module_score"]], seq_len(ncol(obj)) / ncol(obj))

  # a factor/character column would break the continuous colour scale, so it
  # is refused here rather than erroring inside the Shiny reactive
  obj$a_label <- rep("x", ncol(obj))
  expect_null(.sirup_polygon_feature_values(obj, "a_label", layer = "counts"))
  expect_null(.sirup_polygon_feature_values(obj, "not_a_real_feature", layer = "counts"))
  expect_null(.sirup_polygon_feature_values(obj, "", layer = "counts"))
})

test_that("selectSpatialPolygons gains a show_zeros argument", {
  expect_true("show_zeros" %in% names(formals(selectSpatialPolygons)))
  expect_true(isTRUE(eval(formals(selectSpatialPolygons)$show_zeros)))
})

test_that(".sirup_nn_index returns the true nearest point on irregular data", {
  # Irregular on purpose: the index must not assume a lattice, since Visium is
  # hexagonal and single-cell data has no grid at all.
  set.seed(42)
  xy <- cbind(runif(500, 0, 100), runif(500, -100, 0))
  idx <- .sirup_nn_index(xy)
  expect_true(is.character(idx$backend) && nzchar(idx$backend))
  expect_true(is.function(idx$query))

  set.seed(7)
  queries <- cbind(runif(25, 0, 100), runif(25, -100, 0))
  for (i in seq_len(nrow(queries))) {
    qx <- queries[i, 1]; qy <- queries[i, 2]
    expected <- which.min((xy[, 1] - qx)^2 + (xy[, 2] - qy)^2)
    expect_identical(as.integer(idx$query(qx, qy)), as.integer(expected))
  }

  # an exact hit on a data point returns that point
  expect_identical(as.integer(idx$query(xy[13, 1], xy[13, 2])), 13L)
})

test_that(".sirup_nn_index scan backend agrees with whichever backend is chosen", {
  set.seed(1)
  xy <- cbind(runif(200), runif(200))
  chosen <- .sirup_nn_index(xy)
  xs <- xy[, 1]; ys <- xy[, 2]
  scan_query <- function(x, y) which.min((xs - x)^2 + (ys - y)^2)
  for (i in seq_len(20)) {
    qx <- runif(1); qy <- runif(1)
    expect_identical(as.integer(chosen$query(qx, qy)), as.integer(scan_query(qx, qy)))
  }
})

test_that(".sirup_region_from_dominance selects by polygon_id or polygon label", {

  make_spatial <- function() {
    s <- sirup_test_seurat(
      counts = matrix(1, nrow = 1, ncol = 2, dimnames = list("gene", c("cell1", "cell2")))
    )
    s$x <- c(0, 10); s$y <- c(0, 10); s
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
    spatial_id = "barcode", sample = "sample", clone = "clone_id",
    image_scale = "hires", swap_coords = FALSE
  )
  # Same polygon label "gc1" in BOTH sections -> composite, section-unique ids
  # (section_a::gc1 / section_b::gc1) with a shared clean `polygon` label "gc1".
  polygons <- data.frame(
    polygon = rep("gc1", 8),
    sample = rep(c("section_a", "section_b"), each = 4),
    x = rep(c(-1, 1, 0, -1), 2),
    y = rep(c(-1, -1, 1, -1), 2),
    stringsAsFactors = FALSE
  )
  x <- addPolygons(x, polygons = polygons)
  dom <- spatialCloneDominance(x, locus_val = "IGH")

  # Composite polygon_id selects ONE section only.
  r_a <- SIRUP:::.sirup_region_from_dominance(dom, polygon_val = "section_a::gc1")
  expect_equal(nrow(r_a$receptors), 1L)
  expect_equal(as.character(r_a$receptors$sample), "section_a")

  # Clean label "gc1" matches both sections (backwards compatible).
  r_all <- SIRUP:::.sirup_region_from_dominance(dom, polygon_val = "gc1")
  expect_equal(nrow(r_all$receptors), 2L)
})
