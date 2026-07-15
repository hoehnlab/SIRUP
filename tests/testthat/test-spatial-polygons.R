test_that("point-in-polygon includes boundary points when requested", {
  vx <- c(0, 10, 10, 0, 0)
  vy <- c(0, 0, 10, 10, 0)

  inside <- .sirap_points_in_polygon(
    px = c(5, 10, 12),
    py = c(5, 5, 5),
    vx = vx,
    vy = vy,
    include_boundary = TRUE
  )

  expect_equal(inside, c(TRUE, TRUE, FALSE))

  inside_no_boundary <- .sirap_points_in_polygon(
    px = c(5, 10),
    py = c(5, 5),
    vx = vx,
    vy = vy,
    include_boundary = FALSE
  )

  expect_equal(inside_no_boundary, c(TRUE, FALSE))
})

test_that("polygon vertex selection requires closed polygons by default", {
  polygon <- data.frame(
    x = c(0, 10, 10, 0),
    y = c(0, 0, 10, 10),
    stringsAsFactors = FALSE
  )

  expect_error(
    .sirap_select_polygon_vertices(
      polygon = polygon,
      polygon_id = NULL,
      polygon_col = "polygon",
      x_col = "x",
      y_col = "y",
      require_closed = TRUE
    ),
    "must be closed"
  )

  closed <- .sirap_select_polygon_vertices(
    polygon = polygon,
    polygon_id = NULL,
    polygon_col = "polygon",
    x_col = "x",
    y_col = "y",
    require_closed = FALSE
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

  sampled_a <- .sirap_downsample_polygon_plot_coords(coords, downsample_n = 4, seed = 10)
  sampled_b <- .sirap_downsample_polygon_plot_coords(coords, downsample_n = 4, seed = 10)

  expect_equal(nrow(sampled_a), 4)
  expect_equal(sampled_a, sampled_b)
  expect_equal(.sirap_downsample_polygon_plot_coords(coords, downsample_n = 20), coords)
  expect_error(.sirap_downsample_polygon_plot_coords(coords, downsample_n = 0), "positive number")
})

test_that("selected polygon table displays only current polygon without image", {
  points <- data.frame(
    barcode = c("cell1", "cell2"),
    x = c(1, 2),
    y = c(1, 2),
    .sirap_plot_y = c(-1, -2),
    polygon = c(1L, 2L),
    image = c("slice1", "slice1"),
    stringsAsFactors = FALSE
  )

  current <- points[points$polygon == 2L, , drop = FALSE]
  shown <- .sirap_public_polygon_vertices(current)
  shown <- shown[, setdiff(colnames(shown), "image"), drop = FALSE]

  expect_equal(shown$barcode, "cell2")
  expect_false("image" %in% names(shown))
  expect_false(".sirap_plot_y" %in% names(shown))
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

  plot_data <- .sirap_polygon_plot_data(
    polygon = polygon,
    polygon_id = NULL,
    polygon_col = "polygon",
    x_col = "x",
    y_col = "y",
    require_closed = TRUE
  )

  expect_equal(sort(unique(plot_data$polygon)), c("left", "right"))
  expect_equal(plot_data$.sirap_y, -plot_data$y)
})

test_that("V gene helpers retain alleles and extract families and superfamilies", {
  genes <- .sirap_v_genes(c("IGHV3-23*01", "IGHV3-30*02,IGHV3-23*01"))
  families <- .sirap_v_gene_family(c(genes, "IGKV1D-39*01"))

  expect_equal(genes, c("IGHV3-23*01", "IGHV3-30*02"))
  expect_equal(families, c("IGHV3-23", "IGHV3-30", "IGKV1D-39"))
  expect_equal(.sirap_v_gene_superfamily(families), c("IGHV3", "IGHV3", "IGKV1"))
})

test_that("polygon pie summaries group lower-frequency levels and label top levels", {
  region_bcrs <- data.frame(
    polygon = c(rep("poly1", 5), rep("poly2", 3), rep("poly3", 3)),
    clone_id = c("A", "A", "B", "C", "D", "E", "F", "F", "C", "C", "G"),
    stringsAsFactors = FALSE
  )

  pies <- .sirap_polygon_pie_summary(
    region_bcrs,
    fill_by = "clone_id",
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
  expect_equal(sort(unique(pies$polygon_label)), c("1", "2", "3"))
  expect_false(".sirap_order" %in% names(pies))

  pies_without_other <- .sirap_polygon_pie_summary(
    region_bcrs,
    fill_by = "clone_id",
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
  values <- .sirap_discrete_palette_values(
    c("A", "B", "Other"),
    other_label = "Other",
    other_color = "#EEEEEE"
  )

  expect_equal(unname(values[c("A", "B")]), pals::glasbey(2))
  expect_equal(unname(values["Other"]), "#EEEEEE")

  many_levels <- paste0("clone", seq_len(46))
  many_values <- .sirap_discrete_palette_values(
    c(many_levels, "Other"),
    other_label = "Other",
    other_color = "#EEEEEE"
  )
  expect_equal(names(many_values), c(many_levels, "Other"))
  expect_equal(length(many_values), 47)
  expect_equal(unname(many_values["Other"]), "#EEEEEE")
})

test_that("spatial pie sectors start at the visual top", {
  sector <- .sirap_pie_sector(center_x = 10, center_y = 20, radius = 5, ymin = 0, ymax = 0.25)

  expect_equal(sector$x[2], 10, tolerance = 1e-8)
  expect_equal(sector$y[2], 15, tolerance = 1e-8)
})

test_that("3D clone density data groups clones and builds surfaces", {
  ns <- environment(spatialCloneDominance)
  old_image_dim <- get(".sirap_image_dim", envir = ns)
  image_dim_binding_locked <- bindingIsLocked(".sirap_image_dim", ns)
  if (image_dim_binding_locked) {
    unlockBinding(".sirap_image_dim", ns)
  }
  assign(
    ".sirap_image_dim",
    function(image) {
      list(height = 100, width = 100)
    },
    envir = ns
  )
  if (image_dim_binding_locked) {
    lockBinding(".sirap_image_dim", ns)
  }
  withr::defer({
    if (image_dim_binding_locked) {
      unlockBinding(".sirap_image_dim", ns)
    }
    assign(".sirap_image_dim", old_image_dim, envir = ns)
    if (image_dim_binding_locked) {
      lockBinding(".sirap_image_dim", ns)
    }
  })

  region_bcrs <- data.frame(
    clone_id = c("A", "A", "B", "C"),
    x = c(20, 22, 70, 72),
    y = c(20, 22, 70, 72),
    .sirap_y = -c(20, 22, 70, 72),
    polygon = "poly1",
    image = "metadata",
    stringsAsFactors = FALSE
  )
  polygons <- data.frame(
    x = c(10, 80, 80, 10, 10),
    y = c(10, 10, 80, 80, 10),
    .sirap_y = -c(10, 10, 80, 80, 10),
    polygon = "poly1",
    stringsAsFactors = FALSE
  )

  out <- .sirap_clone_density_3d_data(
    region_bcrs = region_bcrs,
    polygons = polygons,
    image = list(),
    spatial_coords = data.frame(
      x = c(0, 100),
      .sirap_y = c(-100, 0)
    ),
    clone_col = "clone_id",
    weight_col = NULL,
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

test_that("spatialCloneDominance counts clone proportions inside a polygon", {
  ns <- environment(spatialCloneDominance)
  old_coordinate_table <- get(".sirap_spatial_coordinate_table", envir = ns)
  coordinate_binding_locked <- bindingIsLocked(".sirap_spatial_coordinate_table", ns)
  if (coordinate_binding_locked) {
    unlockBinding(".sirap_spatial_coordinate_table", ns)
  }
  assign(
    ".sirap_spatial_coordinate_table",
    function(object, images = NULL, image.scale = "hires", swap_coords = FALSE) {
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
    lockBinding(".sirap_spatial_coordinate_table", ns)
  }
  withr::defer({
    if (coordinate_binding_locked) {
      unlockBinding(".sirap_spatial_coordinate_table", ns)
    }
    assign(".sirap_spatial_coordinate_table", old_coordinate_table, envir = ns)
    if (coordinate_binding_locked) {
      lockBinding(".sirap_spatial_coordinate_table", ns)
    }
  })

  bcrs <- data.frame(
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
    bcrs = bcrs,
    object = list(),
    polygon = polygon,
    barcode_col = "barcode"
  )

  expect_s3_class(out, "sirap_spatial_clone_dominance")
  expect_equal(out$dominance$clone_id, c("A", "B"))
  expect_equal(out$dominance$v_gene[out$dominance$clone_id == "A"], "IGHV3-23*01,IGHV3-23*02,IGKV1-5*01")
  expect_equal(out$dominance$v_gene_family[out$dominance$clone_id == "A"], "IGHV3-23,IGKV1-5")
  expect_equal(out$dominance$v_gene_superfamily[out$dominance$clone_id == "A"], "IGHV3,IGKV1")
  expect_equal(out$dominance$n_bcrs, c(3, 1))
  expect_equal(out$dominance$proportion, c(0.75, 0.25))
  expect_equal(unique(out$dominance$region_bcr_total), 4)
  expect_equal(unique(out$dominance$n_region_barcodes), 3)
  expect_equal(sort(unique(out$region_bcrs$.sirap_barcode)), c("cell1", "cell2", "cell3"))
  expect_true(all(c("x", "y", ".sirap_y", "polygon") %in% names(out$region_bcrs)))
  expect_false(any(c("x.x", "x.y", "y.x", "y.y") %in% names(out$region_bcrs)))
  expect_equal(sort(unique(out$region_bcrs$x)), c(1, 2, 8))
  expect_equal(sort(unique(out$region_bcrs$y)), c(1, 2, 8))
  expect_true(all(c("x", "y", ".sirap_y", "polygon") %in% names(out$polygons)))
  expect_false(any(c("polygon_col", "barcode_col", "clone_col", "v_call_col") %in% names(out)))

  igh <- spatialCloneDominance(
    bcrs = bcrs,
    object = list(),
    polygon = polygon,
    barcode_col = "barcode",
    locus = "IGH"
  )

  expect_equal(igh$dominance$clone_id, c("A", "B"))
  expect_equal(igh$dominance$n_bcrs, c(2, 1))
  expect_equal(igh$dominance$proportion, c(2 / 3, 1 / 3))

  no_singletons <- spatialCloneDominance(
    bcrs = bcrs,
    object = list(),
    polygon = polygon,
    barcode_col = "barcode",
    exclude_singletons = TRUE
  )
  expect_equal(no_singletons$dominance$clone_id, "A")
  expect_equal(no_singletons$dominance$n_bcrs, 3)
  expect_equal(unique(no_singletons$region_bcrs$clone_id), "A")

  region <- .sirap_region_from_plot_input(
    bcrs = out,
    object = list(),
    polygon = NULL,
    barcode_col = NULL,
    images = NULL,
    image.scale = "hires",
    swap_coords = FALSE,
    polygon_id = NULL,
    x_col = "x",
    y_col = "y",
    polygon_col = "polygon",
    require_closed = TRUE
  )
  expect_equal(region$bcrs, out$region_bcrs)
  expect_equal(region$polygons, out$polygons)
})

test_that("spatialCloneDominance summarizes multiple polygons independently", {
  ns <- environment(spatialCloneDominance)
  old_coordinate_table <- get(".sirap_spatial_coordinate_table", envir = ns)
  coordinate_binding_locked <- bindingIsLocked(".sirap_spatial_coordinate_table", ns)
  if (coordinate_binding_locked) {
    unlockBinding(".sirap_spatial_coordinate_table", ns)
  }
  assign(
    ".sirap_spatial_coordinate_table",
    function(object, images = NULL, image.scale = "hires", swap_coords = FALSE) {
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
    lockBinding(".sirap_spatial_coordinate_table", ns)
  }
  withr::defer({
    if (coordinate_binding_locked) {
      unlockBinding(".sirap_spatial_coordinate_table", ns)
    }
    assign(".sirap_spatial_coordinate_table", old_coordinate_table, envir = ns)
    if (coordinate_binding_locked) {
      lockBinding(".sirap_spatial_coordinate_table", ns)
    }
  })

  bcrs <- data.frame(
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
    bcrs = bcrs,
    object = list(),
    polygon = polygon,
    barcode_col = "barcode"
  )

  expect_equal(sort(unique(out$dominance$polygon)), c("left", "right"))
  expect_equal(out$dominance$n_bcrs[out$dominance$polygon == "left" & out$dominance$clone_id == "A"], 1)
  expect_equal(out$dominance$n_bcrs[out$dominance$polygon == "left" & out$dominance$clone_id == "B"], 1)
  expect_equal(out$dominance$n_bcrs[out$dominance$polygon == "right" & out$dominance$clone_id == "D"], 2)
  expect_equal(out$dominance$proportion[out$dominance$polygon == "right" & out$dominance$clone_id == "D"], 1)

  right <- spatialCloneDominance(
    bcrs = bcrs,
    object = list(),
    polygon = polygon,
    barcode_col = "barcode",
    polygon_id = "right"
  )
  expect_equal(unique(right$dominance$polygon), "right")
  expect_equal(right$dominance$clone_id, "D")
})
