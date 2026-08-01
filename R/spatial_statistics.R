
#' Calculate global Moran's I for a spatial variable
#'
#' @param object Seurat object with spatial images, or a SIRAP object.
#' @param variable Feature or metadata column to analyze.
#' @param variable_type One of `feature`, `metadata`, `receptor`, or `clone`.
#'   Use `clone` to calculate Moran's I for a binary indicator marking spots
#'   that contain a particular clone ID.
#' @param value Optional metadata value to convert to a binary indicator.
#' @param receptors AIRR-style receptor data frame. Required when
#'   `variable_type = "receptor"` or `"clone"` unless `object` is a SIRAP
#'   object.
#' @param barcode Barcode column in `receptors`. Defaults to `NULL`; when
#'   `object` is a SIRAP object, the stored receptor barcode field is used.
#'   Pass this explicitly when using separate Seurat and receptor objects.
#' @param clone Clone column in `receptors`; used only when
#'   `variable_type = "clone"`. Defaults to `NULL`; when `object` is a SIRAP
#'   object, the stored receptor clone field is used when available.
#' @param layer Seurat expression layer to fetch for features.
#' @param images Optional Seurat image names.
#' @param image_scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param radius Optional distance radius defining neighborhoods.
#' @param k Number of nearest neighbors when `radius` is `NULL`.
#' @param alternative Alternative hypothesis passed to `spdep::moran.test()`.
#' @return List containing the Moran test, input data, neighbors, and weights.
#' @export
spatialMoran <- function(object,
                         variable,
                         variable_type = c("feature", "metadata", "receptor", "clone"),
                         value = NULL,
                         receptors = NULL,
                         barcode = NULL,
                         clone = NULL,
                         layer = "counts",
                         images = NULL,
                         image_scale = NULL,
                         swap_coords = NULL,
                         radius = NULL,
                         k = 6,
                         alternative = "greater") {
  .require_namespace("spdep")
  dat <- .sirap_spatial_variable_values(
    object = object,
    variable = variable,
    variable_type = variable_type,
    value = value,
    receptors = receptors,
    barcode_col = barcode,
    clone_col = clone,
    layer = layer,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  dat <- dat[is.finite(dat$value), , drop = FALSE]
  .check_spatial_values(dat$value)

  nb <- .sirap_spatial_neighbors(dat, radius = radius, k = k)
  listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  test <- spdep::moran.test(dat$value, listw, zero.policy = TRUE, alternative = alternative)

  list(
    test = test,
    data = dat,
    neighbors = nb,
    weights = listw,
    radius = radius,
    k = k,
    variable = unique(dat$variable)[[1]]
  )
}

#' Calculate local Moran's I/LISA for a spatial variable
#'
#' @inheritParams spatialMoran
#' @param radius Optional distance radius defining local neighborhoods. When
#'   `NULL`, k-nearest neighbors are used.
#' @param k Number of nearest neighbors when `radius` is `NULL`.
#' @return Data frame with coordinates, values, local Moran's I, z scores, and p values.
#' @export
spatialLocalMoran <- function(object,
                              variable,
                              variable_type = c("feature", "metadata", "receptor", "clone"),
                              value = NULL,
                              receptors = NULL,
                              barcode = NULL,
                              clone = NULL,
                              layer = "counts",
                              images = NULL,
                              image_scale = NULL,
                              swap_coords = NULL,
                              radius = NULL,
                              k = 6) {
  .require_namespace("spdep")

  dat <- .sirap_spatial_variable_values(
    object = object,
    variable = variable,
    variable_type = variable_type,
    value = value,
    receptors = receptors,
    barcode_col = barcode,
    clone_col = clone,
    layer = layer,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  dat <- dat[is.finite(dat$value), , drop = FALSE]
  .check_spatial_values(dat$value)

  nb <- .sirap_spatial_neighbors(dat, radius = radius, k = k)
  listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  lisa <- spdep::localmoran(dat$value, listw, zero.policy = TRUE)
  lisa <- as.data.frame(lisa)

  out <- cbind(dat, lisa)
  names(out)[names(out) == "Ii"] <- "local_i"
  names(out)[names(out) == "E.Ii"] <- "expected_i"
  names(out)[names(out) == "Var.Ii"] <- "variance_i"
  names(out)[names(out) == "Z.Ii"] <- "z_i"
  names(out)[grepl("^Pr", names(out))] <- "p_value"
  out$radius <- radius
  out$k <- if (is.null(radius)) k else NA_integer_
  out
}

#' Test clone localization with join-count statistics
#'
#' Runs `spdep::joincount.test()` for each clone to test whether neighboring
#' receptor-positive spatial barcodes are more likely to share that clone than
#' expected under random labeling. This is useful for asking whether clones
#' localize within regions of interest such as germinal centers.
#'
#' @param object Seurat object with spatial images, or a SIRAP object.
#' @param receptors AIRR-style receptor data frame containing spatial barcodes
#'   and clone calls. May be omitted when `object` is a SIRAP object.
#' @param clone Column in `receptors` containing clone IDs. Defaults to
#'   `NULL`; when `object` is a SIRAP object, the stored receptor clone field
#'   is used when available.
#' @param barcode Column in `receptors` matching Seurat spot/cell barcodes.
#'   Defaults to `NULL`; when `object` is a SIRAP object, the stored receptor
#'   barcode field is used. Pass this explicitly when using separate Seurat and
#'   receptor objects.
#' @param clone_val Optional clone ID values to test. Defaults to all clones
#'   observed in at least `min_barcodes` spatial barcodes.
#' @param polygon Optional polygon vertex data frame used to restrict the test
#'   to one or more regions.
#' @param polygon_val Optional polygon ID value(s) to select when `polygon`
#'   exists.
#' @param images Optional Seurat image names.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param x,y Coordinate columns in `polygon`.
#' @param polygon_id Optional polygon grouping column in `polygon`.
#' @param include_boundary Include Seurat spots on polygon edges.
#' @param require_closed Require the final polygon vertex to equal the first.
#' @param min_barcodes Minimum number of spatial barcodes containing a clone for
#'   that clone to be tested.
#' @param radius Optional distance radius defining neighborhoods.
#' @param k Number of nearest neighbors when `radius` is `NULL`.
#' @param alternative Alternative hypothesis passed to
#'   `spdep::joincount.test()`.
#' @param other_label Label used for non-target clone barcodes in each binary
#'   test.
#' @param return_data Return the barcode-level mark data used for each test.
#' @return List containing a summary table, raw `spdep` tests, barcode-level
#'   data, neighbors, and weights.
#' @export
spatialCloneJoinCount <- function(object,
                                  receptors = NULL,
                                  clone = NULL,
                                  barcode = NULL,
                                  clone_val = NULL,
                                  polygon = NULL,
                                  polygon_val = NULL,
                                  images = NULL,
                                  image_scale = NULL,
                                  swap_coords = NULL,
                                  x = "x",
                                  y = "y",
                                  polygon_id = "polygon",
                                  include_boundary = TRUE,
                                  require_closed = TRUE,
                                  min_barcodes = 2,
                                  radius = NULL,
                                  k = 6,
                                  alternative = "greater",
                                  other_label = "Other",
                                  return_data = FALSE) {
  options <- .sirap_spatial_options(object, image_scale = image_scale, swap_coords = swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  inputs <- .sirap_resolve_inputs(
    object,
    receptors,
    barcode
  )
  object <- inputs$object
  receptors <- inputs$receptors
  barcode <- inputs$barcode_col
  if (is.null(receptors)) {
    rlang::abort("`receptors` is required when `object` is not a SIRAP object.")
  }
  if (is.null(clone) && !is.null(inputs$receptor_fields$clone)) {
    clone <- inputs$receptor_fields$clone
  }
  if (is.null(barcode)) {
    rlang::abort("`barcode` is required when it cannot be resolved from a SIRAP object.")
  }
  if (is.null(clone)) {
    rlang::abort("`clone` is required when it cannot be resolved from a SIRAP object.")
  }
  .require_namespace("spdep")
  if (!is.numeric(min_barcodes) || length(min_barcodes) != 1L || is.na(min_barcodes) || min_barcodes < 1) {
    rlang::abort("`min_barcodes` must be a positive number.")
  }
  min_barcodes <- as.integer(min_barcodes)

  dat <- .sirap_joincount_receptor_data(
    object = object,
    receptors = receptors,
    clone_col = clone,
    barcode_col = barcode,
    polygon = polygon,
    polygon_id = polygon_val,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords,
    x_col = x,
    y_col = y,
    polygon_col = polygon_id,
    include_boundary = include_boundary,
    require_closed = require_closed
  )
  .check_required_cols(dat, c(clone, ".sirap_barcode", "x", ".sirap_y", "image", "polygon"), "receptors")
  dat <- dat[!is.na(dat[[clone]]) & nzchar(as.character(dat[[clone]])), , drop = FALSE]
  if (nrow(dat) == 0L) {
    rlang::abort("No receptor rows with clone IDs were found in the selected region(s).")
  }

  if (is.null(clone_val)) {
    clone_counts <- sort(table(unique(dat[, c(".sirap_barcode", clone), drop = FALSE])[[clone]]), decreasing = TRUE)
    clone_val <- names(clone_counts)[clone_counts >= min_barcodes]
  } else {
    clone_val <- as.character(clone_val)
  }
  clone_val <- unique(clone_val[!is.na(clone_val) & nzchar(clone_val)])
  if (length(clone_val) == 0L) {
    rlang::abort("No clones met `min_barcodes` in the selected region(s).")
  }

  tests <- list()
  summaries <- list()
  mark_data <- list()
  neighbors <- list()
  weights <- list()
  groups <- split(dat, interaction(dat$image, dat$polygon, drop = TRUE, sep = "||"))
  for (group_name in names(groups)) {
    group_dat <- groups[[group_name]]
    group_key <- strsplit(group_name, "\\|\\|", fixed = FALSE)[[1L]]
    image <- group_key[[1L]]
    polygon_value <- group_key[[2L]]
    sites <- group_dat[, c(".sirap_barcode", "x", "y", ".sirap_y", "image", "polygon"), drop = FALSE]
    sites <- sites[!duplicated(sites$.sirap_barcode), , drop = FALSE]
    if (nrow(sites) < 2L) {
      next
    }
    nb <- .sirap_spatial_neighbors(sites, radius = radius, k = k)
    listw <- spdep::nb2listw(nb, style = "B", zero.policy = TRUE)
    neighbors[[group_name]] <- nb
    weights[[group_name]] <- listw

    for (clone_id in clone_val) {
      clone_cells <- unique(group_dat$.sirap_barcode[as.character(group_dat[[clone]]) == clone_id])
      n_clone_barcodes <- length(clone_cells)
      if (n_clone_barcodes < min_barcodes || n_clone_barcodes == nrow(sites)) {
        next
      }
      local_other_label <- if (identical(clone_id, other_label)) paste0(other_label, "_background") else other_label
      sites$join_count_mark <- factor(
        ifelse(sites$.sirap_barcode %in% clone_cells, clone_id, local_other_label),
        levels = c(local_other_label, clone_id)
      )
      test <- spdep::joincount.test(
        sites$join_count_mark,
        listw,
        zero.policy = TRUE,
        alternative = alternative
      )
      test_name <- paste(image, polygon_value, clone_id, sep = "||")
      tests[[test_name]] <- test
      summaries[[test_name]] <- .sirap_joincount_summary_row(
        test = test,
        clone_id = clone_id,
        image = image,
        polygon = polygon_value,
        n_barcodes = nrow(sites),
        n_clone_barcodes = n_clone_barcodes,
        radius = radius,
        k = if (is.null(radius)) k else NA_integer_,
        alternative = alternative
      )
      if (isTRUE(return_data)) {
        mark_data[[test_name]] <- cbind(
          sites,
          data.frame(clone_id = clone_id, stringsAsFactors = FALSE)
        )
      }
    }
  }

  summary <- dplyr::bind_rows(summaries)
  if (nrow(summary) == 0L) {
    rlang::abort("No clones could be tested after applying `min_barcodes` and region grouping.")
  }
  summary <- summary[order(summary$p_value, summary$image, summary$polygon, summary$clone_id), , drop = FALSE]
  rownames(summary) <- NULL

  list(
    summary = summary,
    tests = tests,
    data = if (isTRUE(return_data)) dplyr::bind_rows(mark_data) else NULL,
    neighbors = neighbors,
    weights = weights,
    radius = radius,
    k = k,
    alternative = alternative,
    min_barcodes = min_barcodes
  )
}

#' Plot a spatial clone join-count result
#'
#' Combines a tissue map of one tested clone with its observed and expected
#' same-clone join counts. The map uses the same receptor-positive spots and
#' neighborhood graph that were used by [spatialCloneJoinCount()].
#'
#' @param join_count Output from [spatialCloneJoinCount()] created with
#'   `return_data = TRUE`.
#' @param object Seurat object with spatial images, or a SIRAP object.
#' @param clone_val Clone ID to plot.
#' @param images Optional image names to include.
#' @param polygon_val Optional polygon IDs to include.
#' @param polygon Optional polygon vertex data frame to draw over the tissue.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate axes.
#' @param show_neighbors Draw edges between neighboring spots that both contain
#'   the selected clone.
#' @param clone_color Color for selected-clone spots and neighbor edges.
#' @param other_color Color for other receptor-positive spots.
#' @param point_size,point_alpha Size and opacity of plotted spots.
#' @param edge_width,edge_alpha Width and opacity of neighbor edges.
#' @param polygon_color,polygon_linewidth Color and width of polygon outlines.
#' @param image_alpha Tissue image opacity.
#' @param zoom Zoom the spatial panel to the spots used in the test.
#' @param title Optional title for the combined plot.
#' @return One combined ggplot object.
#' @export
plotSpatialCloneJoinCount <- function(join_count,
                                      object,
                                      clone_val,
                                      images = NULL,
                                      polygon_val = NULL,
                                      polygon = NULL,
                                      image_scale = NULL,
                                      swap_coords = NULL,
                                      show_neighbors = TRUE,
                                      clone_color = "#D55E00",
                                      other_color = "grey75",
                                      point_size = 1.8,
                                      point_alpha = 0.9,
                                      edge_width = 0.6,
                                      edge_alpha = 0.75,
                                      polygon_color = "black",
                                      polygon_linewidth = 0.6,
                                      image_alpha = 1,
                                      zoom = FALSE,
                                      title = NULL) {
  .require_namespace("cowplot")
  context <- .sirap_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (!is.list(join_count) || is.null(join_count$summary)) {
    rlang::abort("`join_count` must be output from `spatialCloneJoinCount()`.")
  }
  if (is.null(join_count$data)) {
    rlang::abort(
      "`join_count$data` is missing. Re-run `spatialCloneJoinCount()` with `return_data = TRUE`."
    )
  }
  if (length(clone_val) != 1L || is.na(clone_val) || !nzchar(as.character(clone_val))) {
    rlang::abort("`clone_val` must be one non-missing clone ID.")
  }
  clone_val <- as.character(clone_val)

  summary_data <- .as_data_frame(join_count$summary)
  mark_data <- .as_data_frame(join_count$data)
  .check_required_cols(
    summary_data,
    c(
      "clone_id", "image", "polygon", "n_barcodes", "n_clone_barcodes",
      "observed_join_count", "expected_join_count", "statistic", "p_value"
    ),
    "join_count$summary"
  )
  .check_required_cols(
    mark_data,
    c("clone_id", "image", "polygon", ".sirap_barcode", "x", ".sirap_y", "join_count_mark"),
    "join_count$data"
  )

  summary_data <- summary_data[as.character(summary_data$clone_id) == clone_val, , drop = FALSE]
  mark_data <- mark_data[as.character(mark_data$clone_id) == clone_val, , drop = FALSE]
  if (!is.null(images)) {
    summary_data <- summary_data[summary_data$image %in% images, , drop = FALSE]
    mark_data <- mark_data[mark_data$image %in% images, , drop = FALSE]
  }
  if (!is.null(polygon_val)) {
    polygon_val <- as.character(polygon_val)
    summary_data <- summary_data[as.character(summary_data$polygon) %in% polygon_val, , drop = FALSE]
    mark_data <- mark_data[as.character(mark_data$polygon) %in% polygon_val, , drop = FALSE]
  }
  if (nrow(summary_data) == 0L || nrow(mark_data) == 0L) {
    rlang::abort("The selected clone was not found in the requested join-count result.")
  }

  summary_data$.sirap_group <- paste(summary_data$image, summary_data$polygon, sep = "||")
  mark_data$.sirap_group <- paste(mark_data$image, mark_data$polygon, sep = "||")
  keep_groups <- intersect(summary_data$.sirap_group, unique(mark_data$.sirap_group))
  summary_data <- summary_data[summary_data$.sirap_group %in% keep_groups, , drop = FALSE]
  mark_data <- mark_data[mark_data$.sirap_group %in% keep_groups, , drop = FALSE]
  summary_data$.sirap_label <- ifelse(
    summary_data$polygon == "all",
    summary_data$image,
    paste(summary_data$image, summary_data$polygon, sep = ": ")
  )

  polygon_data <- NULL
  if (!is.null(polygon)) {
    polygon_data <- .sirap_polygon_plot_data(
      polygon = .as_data_frame(polygon),
      polygon_id = polygon_val,
      polygon_col = "polygon",
      x_col = "x",
      y_col = "y",
      require_closed = TRUE
    )
  }

  spatial_plots <- lapply(keep_groups, function(group_name) {
    group_data <- mark_data[mark_data$.sirap_group == group_name, , drop = FALSE]
    image_name <- as.character(group_data$image[[1L]])
    polygon_name <- as.character(group_data$polygon[[1L]])
    group_data$.sirap_selected <- as.character(group_data$join_count_mark) == clone_val

    image_names <- tryCatch(SeuratObject::Images(object), error = function(e) character())
    if (image_name %in% image_names) {
      spatial <- .sirap_spatial_data(
        object,
        image = image_name,
        image_scale = image_scale,
        swap_coords = swap_coords
      )
      p <- .sirap_spatial_base(
        plot_data = group_data,
        image = spatial$image,
        image_alpha = image_alpha,
        zoom = zoom,
        view_data = spatial$coords
      )
    } else {
      view_data <- .sirap_spatial_coordinate_table(
        object,
        images = if (length(image_names) > 0L) image_name else NULL,
        image_scale = image_scale,
        swap_coords = swap_coords
      )
      p <- ggplot2::ggplot(group_data, ggplot2::aes(x = .data$x, y = .data$.sirap_y)) +
        ggplot2::theme_void()
      if (!isTRUE(zoom) && nrow(view_data) > 0L) {
        p <- p + ggplot2::coord_fixed(
          xlim = range(view_data$x, na.rm = TRUE),
          ylim = range(view_data$.sirap_y, na.rm = TRUE),
          expand = FALSE
        )
      } else {
        p <- p + ggplot2::coord_fixed()
      }
    }

    if (isTRUE(show_neighbors)) {
      nb <- join_count$neighbors[[group_name]]
      if (!is.null(nb)) {
        edge_data <- .sirap_joincount_clone_edges(group_data, nb)
        if (nrow(edge_data) > 0L) {
          p <- p + ggplot2::geom_segment(
            data = edge_data,
            ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
            inherit.aes = FALSE,
            color = clone_color,
            linewidth = edge_width,
            alpha = edge_alpha
          )
        }
      }
    }
    p <- p +
      ggplot2::geom_point(color = other_color, size = point_size, alpha = point_alpha) +
      ggplot2::geom_point(
        data = group_data[group_data$.sirap_selected, , drop = FALSE],
        color = clone_color,
        size = point_size,
        alpha = point_alpha
      )
    if (!is.null(polygon_data)) {
      this_polygon <- polygon_data
      if (polygon_name != "all") {
        this_polygon <- this_polygon[as.character(this_polygon$polygon) == polygon_name, , drop = FALSE]
      }
      if (nrow(this_polygon) > 0L) {
        p <- p + ggplot2::geom_path(
          data = this_polygon,
          ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$polygon),
          inherit.aes = FALSE,
          color = polygon_color,
          linewidth = polygon_linewidth
        )
      }
    }
    p + ggplot2::labs(title = summary_data$.sirap_label[match(group_name, summary_data$.sirap_group)])
  })
  spatial_panel <- cowplot::plot_grid(plotlist = spatial_plots, ncol = 1, align = "v")

  count_data <- dplyr::bind_rows(
    data.frame(
      .sirap_group = summary_data$.sirap_group,
      .sirap_label = summary_data$.sirap_label,
      count_type = "Expected",
      join_count = summary_data$expected_join_count,
      stringsAsFactors = FALSE
    ),
    data.frame(
      .sirap_group = summary_data$.sirap_group,
      .sirap_label = summary_data$.sirap_label,
      count_type = "Observed",
      join_count = summary_data$observed_join_count,
      stringsAsFactors = FALSE
    )
  )
  count_data$count_type <- factor(count_data$count_type, levels = c("Expected", "Observed"))
  annotation_data <- summary_data[, c(".sirap_group", ".sirap_label"), drop = FALSE]
  annotation_data$.sirap_annotation <- paste0(
    "z = ", formatC(summary_data$statistic, digits = 2, format = "f"),
    "\np = ", format.pval(summary_data$p_value, digits = 2, eps = 0.001),
    "\n", summary_data$n_clone_barcodes, "/", summary_data$n_barcodes, " clone-positive spots"
  )
  annotation_data$x <- 1.5
  annotation_data$y <- max(c(count_data$join_count, 0), na.rm = TRUE) * 1.12

  count_panel <- ggplot2::ggplot(
    count_data,
    ggplot2::aes(x = .data$count_type, y = .data$join_count, fill = .data$count_type)
  ) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = signif(.data$join_count, 3)),
      vjust = -0.4,
      size = 3.5
    ) +
    ggplot2::geom_text(
      data = annotation_data,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$.sirap_annotation),
      inherit.aes = FALSE,
      vjust = 0,
      size = 3.3
    ) +
    ggplot2::scale_fill_manual(values = c("Expected" = "grey70", "Observed" = clone_color), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.35))) +
    ggplot2::labs(x = NULL, y = "Same-clone joins") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  if (nrow(summary_data) > 1L) {
    count_panel <- count_panel + ggplot2::facet_wrap(ggplot2::vars(.data$.sirap_label), ncol = 1)
  }

  combined <- cowplot::plot_grid(
    spatial_panel,
    count_panel,
    nrow = 1,
    rel_widths = c(1.7, 1),
    align = "h",
    axis = "tb"
  )
  plot_title <- title %||% paste("Spatial join count for clone", clone_val)
  cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(plot_title, fontface = "bold", x = 0.5, hjust = 0.5),
    combined,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
}

# Convert a join-count neighbor graph into line segments connecting selected
# clone-positive neighboring spots for plotting.
.sirap_joincount_clone_edges <- function(group_data, neighbors) {
  if (length(neighbors) != nrow(group_data)) {
    rlang::abort("The stored neighbor graph does not match `join_count$data`.")
  }
  selected <- group_data$.sirap_selected
  edges <- vector("list", length(neighbors))
  edge_index <- 1L
  for (i in seq_along(neighbors)) {
    if (!isTRUE(selected[[i]])) {
      next
    }
    for (j in neighbors[[i]]) {
      if (j > i && isTRUE(selected[[j]])) {
        edges[[edge_index]] <- data.frame(
          x = group_data$x[[i]],
          y = group_data$.sirap_y[[i]],
          xend = group_data$x[[j]],
          yend = group_data$.sirap_y[[j]]
        )
        edge_index <- edge_index + 1L
      }
    }
  }
  dplyr::bind_rows(edges)
}

#' Plot local Moran's I over a tissue image
#'
#' @param object SIRAP object, or a spatial Seurat object.
#' @param local_moran Output from `spatialLocalMoran()`.
#' @param value_col Column to color. Defaults to `z_i`.
#' @param images Optional Seurat image names.
#' @param image_scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param colors Three colors for low, midpoint, and high values.
#' @param midpoint Color midpoint.
#' @param point_size Point size.
#' @param point_alpha Point opacity.
#' @param order_points Plot lower values first so higher values appear on top.
#' @param image_alpha Tissue image opacity.
#' @param zoom Zoom to plotted points.
#' @param title Optional title.
#' @return A ggplot/patchwork object.
#' @export
plotSpatialLocalMoran <- function(object,
                                  local_moran,
                                  value_col = "z_i",
                                  images = NULL,
                                  image_scale = NULL,
                                  swap_coords = NULL,
                                  colors = c("#2c7bb6", "#f7f7f7", "#d7191c"),
                                  midpoint = 0,
                                  point_size = 0.7,
                                  point_alpha = 1,
                                  order_points = TRUE,
                                  image_alpha = 1,
                                  zoom = FALSE,
                                  title = NULL) {
  context <- .sirap_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  .check_required_cols(local_moran, c("cell", "x", ".sirap_y", value_col), "local_moran")
  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    plot_data <- local_moran[local_moran$image == image, , drop = FALSE]
    if (nrow(plot_data) == 0L) {
      next
    }
    if (isTRUE(order_points)) {
      plot_data <- plot_data[order(plot_data[[value_col]], na.last = TRUE), , drop = FALSE]
    }
    plots[[image]] <- .sirap_spatial_base(
      plot_data,
      spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[value_col]]),
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_gradient2(
        low = colors[[1]],
        mid = colors[[2]],
        high = colors[[3]],
        midpoint = midpoint,
        name = value_col
      ) +
      ggplot2::labs(title = title %||% paste("Local Moran's I z-score:", unique(plot_data$variable)[[1]]))
  }

  .sirap_return_spatial_plots(plots, combine = TRUE)
}

#' Calculate Ripley's cross K for marked spatial points
#'
#' @param object SIRAP object, or a spatial Seurat object.
#' @param mark_col Metadata column containing point types.
#' @param i,j Mark levels to compare.
#' @param nsim Number of random-label simulations for an envelope. Use `0` for
#'   observed K only.
#' @param max_r Optional maximum radius to evaluate. Ignored if `r` is supplied.
#' @param r Optional distance vector passed to `spatstat.explore::Kcross()`.
#' @param image_scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @return List with observed K, optional envelope, plotting data, and point data.
#' @export
spatialRipleyCrossK <- function(object,
                                mark_col,
                                i,
                                j,
                                nsim = 0,
                                max_r = NULL,
                                r = NULL,
                                image_scale = NULL,
                                swap_coords = NULL) {
  context <- .sirap_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  .require_namespace("spatstat.geom")
  .require_namespace("spatstat.explore")
  if (nsim > 0L) {
    .require_namespace("spatstat.random")
  }

  # Get the spatial coordinates
  dat <- .sirap_spatial_coordinate_table(object, image_scale = image_scale, swap_coords = swap_coords)
  .check_required_cols(object@meta.data, mark_col, "object@meta.data")
  dat$mark <- object@meta.data[dat$cell, mark_col]
  dat <- dat[!is.na(dat$mark), , drop = FALSE]
  if (!all(c(i, j) %in% dat$mark)) {
    rlang::abort("Both `i` and `j` must be present in the marked spatial data.")
  }
  dat$mark <- factor(dat$mark)

  observed <- list()
  envelopes <- list()
  plot_data <- list()
  for (image in unique(dat$image)) {
    image_data <- dat[dat$image == image, , drop = FALSE]
    window <- spatstat.geom::owin(
      xrange = range(image_data$x, na.rm = TRUE),
      yrange = range(image_data$.sirap_y, na.rm = TRUE)
    )
    pp <- spatstat.geom::ppp(
      x = image_data$x,
      y = image_data$.sirap_y,
      marks = image_data$mark,
      window = window,
      checkdup = FALSE
    )
    r_eval <- r
    if (is.null(r_eval) && !is.null(max_r)) {
      r_eval <- seq(0, max_r, length.out = 513L)
    }
    k_obs <- spatstat.explore::Kcross(pp, i = i, j = j, correction = "border", r = r_eval)
    observed[[image]] <- k_obs

    if (nsim > 0L) {
      envelopes[[image]] <- spatstat.explore::envelope(
        pp,
        fun = spatstat.explore::Kcross,
        nsim = nsim,
        simulate = expression(spatstat.random::rlabel(pp)),
        i = i,
        j = j,
        correction = "border",
        r = r_eval
      )
      plot_data[[image]] <- .sirap_kcross_df(envelopes[[image]], image = image)
    } else {
      plot_data[[image]] <- .sirap_kcross_df(k_obs, image = image)
    }
  }

  list(
    observed = observed,
    envelope = envelopes,
    plot_data = dplyr::bind_rows(plot_data),
    points = dat,
    i = i,
    j = j,
    correction = "border",
    nsim = nsim,
    max_r = max_r,
    r = r
  )
}

#' Plot Ripley's cross K versus distance
#'
#' @param kcross Output from `spatialRipleyCrossK()`.
#' @param title Optional title.
#' @param observed_color Color for the observed cross-K curve.
#' @param simulated_color Color for the dotted simulated curve.
#' @param envelope_fill Fill color for the simulation envelope.
#' @param envelope_alpha Alpha for the simulation envelope.
#' @param facet Facet by image when multiple spatial images are present.
#' @param log_y Use a log10 y scale.
#' @return A ggplot object.
#' @export
plotRipleyCrossK <- function(kcross,
                             title = NULL,
                             observed_color = "black",
                             simulated_color = "grey30",
                             envelope_fill = "grey75",
                             envelope_alpha = 0.45,
                             facet = TRUE,
                             log_y = FALSE) {
  dat <- kcross$plot_data
  .check_required_cols(dat, c("r", "observed", "simulated"), "kcross$plot_data")
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$r))
  if (all(c("lo", "hi") %in% colnames(dat))) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$lo,
          ymax = .data$hi,
          group = .data$image,
          fill = "Simulation envelope"
        ),
        alpha = envelope_alpha,
        inherit.aes = TRUE
      )
  }
  p <- p +
    ggplot2::geom_line(
      ggplot2::aes(
        y = .data$simulated,
        color = "Simulated",
        linetype = "Simulated"
      ),
      size = 0.8
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = .data$observed,
        color = "Observed",
        linetype = "Observed"
      ),
      size = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = c("Observed" = observed_color, "Simulated" = simulated_color),
      breaks = c("Observed", "Simulated"),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Observed" = "solid", "Simulated" = "dotted"),
      breaks = c("Observed", "Simulated"),
      name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c("Simulation envelope" = envelope_fill),
      breaks = "Simulation envelope",
      name = NULL
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% paste("Ripley's cross K:", kcross$i, "to", kcross$j),
      x = "Distance",
      y = "Kcross"
    )
  if (isTRUE(facet) && length(unique(dat$image)) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$image))
  }
  if (isTRUE(log_y)) {
    p <- p + ggplot2::scale_y_log10()
  }
  p
}

# Resolve feature, metadata, receptor, or clone values onto spatial coordinates.
.sirap_spatial_variable_values <- function(object,
                                           variable,
                                           variable_type,
                                           value,
                                           receptors,
                                           barcode_col,
                                           clone_col,
                                           layer,
                                           images,
                                           image_scale,
                                           swap_coords) {
  options <- .sirap_spatial_options(object, image_scale = image_scale, swap_coords = swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  inputs <- .sirap_resolve_inputs(object, receptors, barcode_col)
  object <- inputs$object
  receptors <- inputs$receptors
  barcode_col <- inputs$barcode_col
  if (is.null(clone_col) && !is.null(inputs$receptor_fields$clone)) {
    clone_col <- inputs$receptor_fields$clone
  }
  variable_type <- match.arg(variable_type, c("feature", "metadata", "receptor", "clone"))
  .require_namespace("SeuratObject")

  # Start from one row per spatial coordinate so all statistic modes share the
  # same neighbor graph construction.
  dat <- .sirap_spatial_coordinate_table(object, images = images, image_scale = image_scale, swap_coords = swap_coords)

  if (variable_type == "feature") {
    # Features are pulled from the requested expression layer and analyzed as
    # quantitative values.
    expr <- .sirap_fetch_data(object, vars = variable, layer = layer)
    dat <- merge(dat, expr, by.x = "cell", by.y = 0, all.x = TRUE)
    dat$value <- dat[[variable]]
    dat$variable <- variable
  } else if (variable_type == "metadata") {
    # Numeric/logical metadata are analyzed directly. Categorical metadata need
    # an explicit level so the statistic is a binary indicator.
    .check_required_cols(object@meta.data, variable, "object@meta.data")
    meta <- object@meta.data[, variable, drop = FALSE]
    meta$cell <- rownames(meta)
    dat <- merge(dat, meta, by = "cell", all.x = TRUE)
    if (!is.null(value)) {
      dat$value <- as.integer(dat[[variable]] %in% value)
      dat$variable <- paste0(variable, "_", paste(value, collapse = "_"))
    } else if (is.numeric(dat[[variable]]) || is.integer(dat[[variable]]) || is.logical(dat[[variable]])) {
      dat$value <- as.numeric(dat[[variable]])
      dat$variable <- variable
    } else {
      rlang::abort("For non-numeric metadata, provide `value` to analyze a binary indicator.")
    }
  } else if (variable_type == "receptor") {
    # Receptor mode asks whether receptor-positive spots are spatially
    # autocorrelated, independent of clone identity.
    if (is.null(receptors)) {
      rlang::abort("`receptors` is required when `variable_type = 'receptor'`.")
    }
    if (is.null(barcode_col)) {
      rlang::abort("`barcode` is required when it cannot be resolved from a SIRAP object.")
    }
    receptors <- .as_data_frame(receptors)
    .check_required_cols(receptors, barcode_col, "receptors")
    receptor_barcodes <- unique(as.character(receptors[[barcode_col]]))
    dat$value <- as.integer(dat$cell %in% receptor_barcodes)
    dat$variable <- variable
  } else {
    # Clone mode is a clone-specific receptor indicator: spots containing any
    # row from the requested clone(s) are 1 and all other spatial spots are 0.
    if (is.null(receptors)) {
      rlang::abort("`receptors` is required when `variable_type = 'clone'`.")
    }
    if (is.null(barcode_col)) {
      rlang::abort("`barcode` is required when it cannot be resolved from a SIRAP object.")
    }
    if (is.null(clone_col)) {
      rlang::abort("`clone` is required when it cannot be resolved from a SIRAP object.")
    }
    receptors <- .as_data_frame(receptors)
    .check_required_cols(receptors, c(barcode_col, clone_col), "receptors")
    clone_ids <- as.character(value %||% variable)
    clone_ids <- unique(clone_ids[!is.na(clone_ids) & nzchar(clone_ids)])
    if (length(clone_ids) == 0L) {
      rlang::abort("`variable` or `value` must name at least one clone ID when `variable_type = 'clone'`.")
    }
    receptor_barcodes <- unique(as.character(receptors[[barcode_col]][as.character(receptors[[clone_col]]) %in% clone_ids]))
    dat$value <- as.integer(dat$cell %in% receptor_barcodes)
    dat$variable <- paste0("clone_", paste(clone_ids, collapse = "_"))
  }

  dat <- dat[!is.na(dat$value), , drop = FALSE]
  dat
}

# Collect spatial coordinates from Seurat images when present; otherwise use
# spatial coordinate columns stored in object metadata.
.sirap_spatial_coordinate_table <- function(object, images = NULL, image_scale = NULL, swap_coords = NULL) {
  options <- .sirap_require_spatial_options(image_scale, swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  .require_namespace("SeuratObject")
  image_names <- tryCatch(
    SeuratObject::Images(object),
    error = function(e) character()
  )
  if (is.null(images) && length(image_names) == 0L) {
    return(.sirap_metadata_coordinate_table(object))
  }

  coords <- lapply(.sirap_spatial_images(object, images = images), function(image) {
    spatial <- .sirap_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    out <- spatial$coords
    out$image <- image
    out
  })
  dplyr::bind_rows(coords)
}

# Fallback for image-free objects used in tests or non-Visium spatial datasets.
.sirap_metadata_coordinate_table <- function(object) {
  meta <- object@meta.data
  x_col <- intersect(c("x", "imagecol", "pxl_col_in_fullres", "col", "array_col", "spatial_x"), colnames(meta))
  y_col <- intersect(c("y", "imagerow", "pxl_row_in_fullres", "row", "array_row", "spatial_y"), colnames(meta))
  if (length(x_col) == 0L || length(y_col) == 0L) {
    rlang::abort(
      paste0(
        "Could not find spatial images or coordinate columns in `object@meta.data`. ",
        "For image-free spatial statistics, add coordinate columns such as `x` and `y`."
      )
    )
  }
  coords <- data.frame(
    cell = rownames(meta),
    x = meta[[x_col[[1]]]],
    y = meta[[y_col[[1]]]],
    image = "metadata",
    stringsAsFactors = FALSE
  )
  coords <- coords[coords$cell %in% colnames(object), , drop = FALSE]
  coords$.sirap_y <- -coords$y
  coords
}

# Build an spdep neighbor graph either by distance radius or k-nearest
# neighbors over SIRAP's plotting coordinate convention.
.sirap_spatial_neighbors <- function(dat, radius = NULL, k = 6) {
  if (nrow(dat) < 2L) {
    rlang::abort("At least two spatial points are required to build neighborhoods.")
  }
  coords <- as.matrix(dat[, c("x", ".sirap_y")])
  if (!is.null(radius)) {
    spdep::dnearneigh(coords, d1 = 0, d2 = radius, longlat = FALSE)
  } else {
    k <- min(k, nrow(dat) - 1L)
    knn <- spdep::knearneigh(coords, k = k, longlat = FALSE)
    spdep::knn2nb(knn)
  }
}

# Guard against degenerate values before asking spdep for a statistic.
.check_spatial_values <- function(x) {
  if (length(unique(x)) < 2L) {
    rlang::abort("Spatial statistic requires at least two distinct values.")
  }
  invisible(TRUE)
}

# Join receptor clone calls to spatial coordinates and, when supplied, restrict
# them to polygon regions before join-count testing.
.sirap_joincount_receptor_data <- function(object,
                                           receptors,
                                           clone_col,
                                           barcode_col,
                                           polygon,
                                           polygon_id,
                                           images,
                                           image_scale,
                                           swap_coords,
                                           x_col,
                                           y_col,
                                           polygon_col,
                                           include_boundary,
                                           require_closed) {
  receptors <- .as_data_frame(receptors)
  .check_required_cols(receptors, c(barcode_col, clone_col), "receptors")
  coords <- .sirap_spatial_coordinate_table(object, images = images, image_scale = image_scale, swap_coords = swap_coords)
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  if (nrow(coords) == 0L) {
    rlang::abort("No finite spatial coordinates were found in `object`.")
  }

  receptor_cols <- setdiff(colnames(receptors), setdiff(c("cell", "x", "y", ".sirap_y", "image", "polygon"), c(barcode_col, clone_col)))
  receptors <- receptors[, receptor_cols, drop = FALSE]
  receptors$.sirap_barcode <- as.character(receptors[[barcode_col]])

  if (is.null(polygon)) {
    region_coords <- coords
    region_coords$polygon <- "all"
    out <- merge(
      receptors,
      region_coords[, c("cell", "x", "y", ".sirap_y", "image", "polygon"), drop = FALSE],
      by.x = ".sirap_barcode",
      by.y = "cell",
      all.x = FALSE,
      all.y = FALSE,
      sort = FALSE
    )
    return(out)
  }

  polygon <- .as_data_frame(polygon)
  .check_required_cols(polygon, c(x_col, y_col), "polygon")
  polygon_specs <- .sirap_polygon_specs(polygon, polygon_id = polygon_id, polygon_col = polygon_col)
  regions <- lapply(polygon_specs, function(spec) {
    vertices <- .sirap_select_polygon_vertices(
      polygon = polygon,
      polygon_id = spec$id,
      polygon_col = polygon_col,
      x_col = x_col,
      y_col = y_col,
      require_closed = require_closed
    )
    inside <- .sirap_points_in_polygon(
      px = coords$x,
      py = coords$y,
      vx = vertices$x,
      vy = vertices$y,
      include_boundary = include_boundary
    )
    region_coords <- coords[inside, , drop = FALSE]
    if (nrow(region_coords) == 0L) {
      return(data.frame())
    }
    region_coords$polygon <- as.character(if (isTRUE(spec$add_polygon_col)) spec$id else "polygon_1")
    merge(
      receptors,
      region_coords[, c("cell", "x", "y", ".sirap_y", "image", "polygon"), drop = FALSE],
      by.x = ".sirap_barcode",
      by.y = "cell",
      all.x = FALSE,
      all.y = FALSE,
      sort = FALSE
    )
  })
  dplyr::bind_rows(regions)
}

# Convert spdep's join-count htest object into one stable summary row.
.sirap_joincount_summary_row <- function(test,
                                         clone_id,
                                         image,
                                         polygon,
                                         n_barcodes,
                                         n_clone_barcodes,
                                         radius,
                                         k,
                                         alternative) {
  htest <- .sirap_joincount_clone_htest(test, clone_id)
  estimate <- htest$estimate %||% numeric()
  estimate_names <- names(estimate)
  get_estimate <- function(pattern, fallback_index) {
    idx <- grep(pattern, estimate_names, ignore.case = TRUE)
    if (length(idx) > 0L) {
      return(unname(estimate[[idx[[1L]]]]))
    }
    if (length(estimate) >= fallback_index) {
      return(unname(estimate[[fallback_index]]))
    }
    NA_real_
  }
  data.frame(
    clone_id = as.character(clone_id),
    image = as.character(image),
    polygon = as.character(polygon),
    n_barcodes = as.integer(n_barcodes),
    n_clone_barcodes = as.integer(n_clone_barcodes),
    observed_join_count = get_estimate("join|observed", 1L),
    expected_join_count = get_estimate("expect", 2L),
    variance = get_estimate("var", 3L),
    statistic = if (length(htest$statistic) > 0L) unname(htest$statistic[[1L]]) else NA_real_,
    p_value = htest$p.value %||% NA_real_,
    method = htest$method %||% NA_character_,
    radius = radius %||% NA_real_,
    k = k,
    alternative = alternative,
    stringsAsFactors = FALSE
  )
}

# spdep may return either an htest or a list of htests depending on factor
# levels; this extracts the htest for the requested clone.
.sirap_joincount_clone_htest <- function(test, clone_id) {
  if (inherits(test, "htest")) {
    return(test)
  }
  if (is.list(test)) {
    named_match <- which(names(test) == as.character(clone_id))
    if (length(named_match) > 0L && inherits(test[[named_match[[1L]]]], "htest")) {
      return(test[[named_match[[1L]]]])
    }
    htest_idx <- which(vapply(test, inherits, logical(1), what = "htest"))
    if (length(htest_idx) > 0L) {
      method_text <- vapply(test[htest_idx], function(x) x$method %||% "", character(1))
      method_match <- grep(as.character(clone_id), method_text, fixed = TRUE)
      if (length(method_match) > 0L) {
        return(test[[htest_idx[[method_match[[1L]]]]]])
      }
      return(test[[utils::tail(htest_idx, 1L)]])
    }
  }
  rlang::abort("Could not extract an htest result from `spdep::joincount.test()`.")
}

# Convert a spatstat Kcross object to a consistent plotting data frame.
.sirap_kcross_df <- function(x, image) {
  dat <- as.data.frame(x)
  theoretical_col <- intersect(c("theo", "theoretical"), colnames(dat))
  simulated_col <- intersect(c("mmean", "mean", "simulated", "simmean", theoretical_col), colnames(dat))
  observed_col <- if ("border" %in% colnames(dat)) {
    "border"
  } else if ("obs" %in% colnames(dat)) {
    "obs"
  } else {
    setdiff(colnames(dat), c("r", theoretical_col, simulated_col, "lo", "hi"))[[1]]
  }
  theoretical <- if (length(theoretical_col) > 0L) dat[[theoretical_col[[1]]]] else pi * dat$r^2
  simulated <- if (length(simulated_col) > 0L) dat[[simulated_col[[1]]]] else theoretical
  out <- data.frame(
    image = image,
    r = dat$r,
    observed = dat[[observed_col]],
    theoretical = theoretical,
    simulated = simulated,
    stringsAsFactors = FALSE
  )
  if (all(c("lo", "hi") %in% colnames(dat))) {
    out$lo <- dat$lo
    out$hi <- dat$hi
  }
  out
}
