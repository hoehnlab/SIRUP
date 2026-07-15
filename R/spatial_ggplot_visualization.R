#' Plot V/J/junction length QC histograms
#'
#' @param airr_df AIRR-style data frame.
#' @param title Optional plot title prefix.
#' @return A ggplot object.
#' @export
plotLengthQc <- function(airr_df, title = NULL) {
  .check_required_cols(airr_df, c("v_sequence_end", "v_sequence_start", "j_sequence_end", "j_sequence_start", "junction_length"), "airr_df")
  out <- .as_data_frame(airr_df)
  out$vlength <- out$v_sequence_end - out$v_sequence_start
  out$jlength <- out$j_sequence_end - out$j_sequence_start
  long <- tidyr::pivot_longer(
    out,
    cols = dplyr::all_of(c("vlength", "jlength", "junction_length")),
    names_to = "metric",
    values_to = "length"
  ) %>%
    dplyr::mutate(metric = factor(.data$metric, levels = c("vlength", "jlength", "junction_length"), labels = c("V Length", "J Length", "Junction Length")))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$length)) +
    ggplot2::geom_histogram(binwidth = 5, na.rm = TRUE) +
    ggplot2::facet_wrap(ggplot2::vars(.data$metric), scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Length", y = "Count")
}

#' Plot spatial features after excluding zero-expression spots
#'
#' @param object Seurat object with spatial images.
#' @param features Feature names to plot.
#' @param assay Optional assay to set before plotting.
#' @param layer Seurat expression layer to fetch.
#' @param images Optional Seurat image names.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param feature_colors Colors used, in order, for the continuous expression
#'   scale. The first color is the low value and the last color is the high
#'   value.
#' @param zero_tol Values at or below this threshold are treated as zero.
#' @param point_size Size of plotted spots.
#' @param point_alpha Opacity for plotted spots.
#' @param order_points Plot lower values first so higher values appear on top.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to non-zero spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param combine Combine plots with patchwork when possible.
#' @return A ggplot/patchwork object, list, or `NULL` when all spots are zero.
#' @export
spatialExpressionPlot <- function(object,
                                      features,
                                      assay = NULL,
                                      layer = "counts",
                                      images = NULL,
                                      image.scale = "hires",
                                      swap_coords = FALSE,
                                      feature_colors = c("#D3D3D3", "#FFFF00"),
                                      zero_tol = 0,
                                      point_size = 0.6,
                                      point_alpha = 1,
                                      order_points = TRUE,
                                      image.alpha = 1,
                                      zoom = FALSE,
                                      combine = TRUE) {
  .require_namespace("SeuratObject")
  .require_namespace("patchwork")
  if (!is.null(assay)) {
    SeuratObject::DefaultAssay(object) <- assay
  }
  if (length(feature_colors) < 2L) {
    rlang::abort("`feature_colors` must contain at least two colors.")
  }

  expr <- .sirap_fetch_data(object, vars = features, layer = layer)
  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    dat <- merge(spatial$coords, expr, by.x = "cell", by.y = 0, all.x = TRUE)
    for (feature in features) {
      keep <- !is.na(dat[[feature]]) & dat[[feature]] > zero_tol
      if (!any(keep)) {
        next
      }
      plot_data <- dat[keep, , drop = FALSE]
      if (isTRUE(order_points)) {
        plot_data <- plot_data[order(plot_data[[feature]], na.last = TRUE), , drop = FALSE]
      }
      plots[[paste(image, feature, sep = "_")]] <- .sirap_spatial_base(
        plot_data,
        spatial$image,
        image.alpha = image.alpha,
        zoom = zoom,
        view_data = spatial$coords
      ) +
        ggplot2::geom_point(
          ggplot2::aes(color = .data[[feature]]),
          size = point_size,
          alpha = point_alpha,
          na.rm = TRUE
        ) +
        ggplot2::scale_color_gradientn(colors = feature_colors, name = feature) +
        ggplot2::labs(title = feature)
    }
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

#' Plot spatial metadata groups after excluding zero/blank spots
#'
#' @param object Seurat object with spatial images.
#' @param meta_col Metadata column to plot.
#' @param zero_value Value treated as not present.
#' @param images Optional Seurat image names.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param cols Optional named color vector.
#' @param point_size Size of plotted spots.
#' @param point_alpha Opacity for plotted spots.
#' @param order_points Plot lower numeric values or earlier factor levels first
#'   so higher/later values appear on top.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param combine Combine plots with patchwork when possible.
#' @return A ggplot/patchwork object, list, or `NULL` when all spots are blank.
#' @export
spatialDiscretePlot <- function(object,
                                  meta_col,
                                  zero_value = 0,
                                  images = NULL,
                                  image.scale = "hires",
                                  swap_coords = FALSE,
                                  cols = NULL,
                                  point_size = 1.6,
                                  point_alpha = 1,
                                  order_points = TRUE,
                                  image.alpha = 1,
                                  zoom = FALSE,
                                  combine = TRUE) {
  .check_required_cols(object@meta.data, meta_col, "object@meta.data")

  values <- object@meta.data[[meta_col]]
  names(values) <- rownames(object@meta.data)
  keep_cells <- names(values)[!is.na(values) & as.character(values) != as.character(zero_value)]
  if (length(keep_cells) == 0L) {
    message("All spots are zero/blank for `", meta_col, "`.")
    return(NULL)
  }

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    plot_data <- spatial$coords[spatial$coords$cell %in% keep_cells, , drop = FALSE]
    if (nrow(plot_data) == 0L) {
      next
    }
    plot_data[[meta_col]] <- values[plot_data$cell]
    plot_data[[meta_col]] <- factor(plot_data[[meta_col]], levels = unique(plot_data[[meta_col]]))
    if (isTRUE(order_points)) {
      plot_order <- suppressWarnings(as.numeric(as.character(plot_data[[meta_col]])))
      if (all(is.na(plot_order))) {
        plot_order <- as.integer(plot_data[[meta_col]])
      }
      plot_data <- plot_data[order(plot_order, na.last = TRUE), , drop = FALSE]
    }
    color_values <- cols %||% scales::hue_pal()(length(levels(plot_data[[meta_col]])))
    if (is.null(names(color_values))) {
      names(color_values) <- levels(plot_data[[meta_col]])
    }

    plots[[image]] <- .sirap_spatial_base(
      plot_data,
      spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[meta_col]]),
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = color_values, name = meta_col) +
      ggplot2::labs(title = meta_col)
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

#' Blend two spatial features after excluding double-zero spots
#'
#' @param object Seurat object with spatial images.
#' @param features Two feature names.
#' @param combine Combine plots with patchwork.
#' @param feature_1_alt_name,feature_2_alt_name Optional axis labels for the
#'   blend legend.
#' @param assay Optional assay to set before plotting.
#' @param layer Seurat expression layer to fetch.
#' @param only_blended Return only the blended spatial panel and legend.
#' @param bottom_left,bottom_right,top_left,top_right Corner colors for the
#'   two-feature blend grid.
#' @param zero_tol Values at or below this threshold are treated as zero.
#' @param images Optional Seurat image names.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param point_size Size of plotted spots.
#' @param point_alpha Opacity for plotted spots.
#' @param order_points Plot lower blended expression intensity first so higher
#'   intensity spots appear on top.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @return A ggplot/patchwork object or a list when `combine = FALSE`.
#' @export
spatialExpressionBlendPlot <- function(object,
                                                features,
                                                combine = TRUE,
                                                feature_1_alt_name = NULL,
                                                feature_2_alt_name = NULL,
                                                assay = NULL,
                                                layer = "counts",
                                                only_blended = FALSE,
                                                bottom_left = "#000000",
                                                bottom_right = "#FF0000",
                                                top_left = "#00FF00",
                                                top_right = "#FFFF00",
                                                zero_tol = 0,
                                                images = NULL,
                                                image.scale = "hires",
                                                swap_coords = FALSE,
                                                point_size = 0.6,
                                                point_alpha = 1,
                                                order_points = TRUE,
                                                image.alpha = 1,
                                                zoom = FALSE) {
  .require_namespace("SeuratObject")
  .require_namespace("patchwork")
  if (length(features) != 2L) {
    rlang::abort(paste0("Requires two features, received ", length(features), "."))
  }
  if (!is.null(assay)) {
    SeuratObject::DefaultAssay(object) <- assay
  }

  expr <- .sirap_fetch_data(object, vars = features, layer = layer)
  color_grid <- .sirap_blend_color_grid(
    bottom_left = bottom_left,
    bottom_right = bottom_right,
    top_left = top_left,
    top_right = top_right
  )
  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    dat <- merge(spatial$coords, expr, by.x = "cell", by.y = 0, all.x = TRUE)
    keep <- rowSums(dat[, features, drop = FALSE] > zero_tol, na.rm = TRUE) > 0
    if (!any(keep)) {
      next
    }
    dat <- dat[keep, , drop = FALSE]
    if (isTRUE(order_points)) {
      dat$.sirap_blend_order <- rowSums(dat[, features, drop = FALSE], na.rm = TRUE)
      dat <- dat[order(dat$.sirap_blend_order, na.last = TRUE), , drop = FALSE]
    }
    dat$.sirap_blend_color <- .sirap_blend_colors(dat[, features, drop = FALSE], color_grid)
    color_values <- stats::setNames(unique(dat$.sirap_blend_color), unique(dat$.sirap_blend_color))

    blended <- .sirap_spatial_base(
      dat,
      spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$.sirap_blend_color),
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::labs(title = paste(features, collapse = " + "))

    legend <- .sirap_blend_legend(
      color_grid,
      xlab = feature_1_alt_name %||% features[[1]],
      ylab = feature_2_alt_name %||% features[[2]]
    )

    if (only_blended) {
      plots[[image]] <- patchwork::wrap_plots(list(blended, legend), nrow = 1, widths = c(0.82, 0.18))
    } else {
      f1 <- spatialExpressionPlot(
        object,
        features[[1]],
        images = image,
        image.scale = image.scale,
        swap_coords = swap_coords,
        feature_colors = c(bottom_left, bottom_right),
        zero_tol = zero_tol,
        point_size = point_size,
        point_alpha = point_alpha,
        order_points = order_points,
        image.alpha = image.alpha,
        zoom = zoom,
        combine = TRUE
      )
      f2 <- spatialExpressionPlot(
        object,
        features[[2]],
        images = image,
        image.scale = image.scale,
        swap_coords = swap_coords,
        feature_colors = c(bottom_left, top_left),
        zero_tol = zero_tol,
        point_size = point_size,
        point_alpha = point_alpha,
        order_points = order_points,
        image.alpha = image.alpha,
        zoom = zoom,
        combine = TRUE
      )
      plots[[image]] <- patchwork::wrap_plots(list(f1, f2, blended, legend), nrow = 1, widths = c(0.28, 0.28, 0.28, 0.16))
    }
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

#' Plot multiple spatial features as ordered color layers
#'
#' @param object Seurat object with spatial images.
#' @param features Feature names.
#' @param feat_order Optional plotting priority, broad/low priority first and
#'   rare/high priority last.
#' @param combine Combine spatial panels and legends with patchwork.
#' @param assay Optional assay to set before plotting.
#' @param layer Seurat expression layer to fetch.
#' @param colors One color per feature.
#' @param breadth_metric Metric used to choose default feature order.
#' @param images Optional Seurat image names.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param point_size Size of plotted spots.
#' @param point_alpha Opacity for plotted spots.
#' @param order_points Plot lower-priority feature assignments first so later
#'   feature levels appear on top.
#' @param zero_tol Values at or below this threshold are treated as zero.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @return A ggplot/patchwork object or a list when `combine = FALSE`.
#' @export
spatialExpressionMultiPlot <- function(object,
                                    features,
                                    feat_order = NULL,
                                    combine = TRUE,
                                    assay = NULL,
                                    layer = "counts",
                                    colors = NULL,
                                    breadth_metric = c("n_expressing"),
                                    images = NULL,
                                    image.scale = "hires",
                                    swap_coords = FALSE,
                                    point_size = 0.6,
                                    point_alpha = 1,
                                    order_points = TRUE,
                                    zero_tol = 0,
                                    image.alpha = 1,
                                    zoom = FALSE) {
  .require_namespace("SeuratObject")
  breadth_metric <- match.arg(breadth_metric)
  if (length(features) < 2L) {
    rlang::abort("Requires at least 2 features.")
  }
  if (!is.null(assay)) {
    SeuratObject::DefaultAssay(object) <- assay
  }
  colors <- colors %||% c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999")
  if (length(colors) < length(features)) {
    rlang::abort("`colors` must be at least as long as `features`.")
  }
  names(colors) <- features

  expr <- .sirap_fetch_data(object, vars = features, layer = layer)
  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    dat <- merge(spatial$coords, expr, by.x = "cell", by.y = 0, all.x = TRUE)
    keep <- rowSums(dat[, features, drop = FALSE] > zero_tol, na.rm = TRUE) > 0
    if (!any(keep)) {
      next
    }
    dat <- dat[keep, , drop = FALSE]
    breadths <- colSums(dat[, features, drop = FALSE] > zero_tol, na.rm = TRUE)
    feature_order <- feat_order %||% names(sort(breadths, decreasing = TRUE))
    missing_features <- setdiff(feature_order, features)
    if (length(missing_features) > 0L) {
      rlang::abort(paste0("`feat_order` includes feature(s) not found in `features`: ", paste(missing_features, collapse = ", ")))
    }

    dat$.sirap_feature <- NA_character_
    for (feature in feature_order) {
      dat$.sirap_feature[dat[[feature]] > zero_tol] <- feature
    }
    dat$.sirap_feature <- factor(dat$.sirap_feature, levels = feature_order)
    if (isTRUE(order_points)) {
      dat <- dat[order(as.integer(dat$.sirap_feature), na.last = TRUE), , drop = FALSE]
    }
    ordered_colors <- colors[feature_order]

    spatial_plot <- .sirap_spatial_base(
      dat,
      spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$.sirap_feature),
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = ordered_colors, name = "Feature") +
      ggplot2::labs(title = paste(features, collapse = " / "))

    plots[[image]] <- spatial_plot
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

#' Plot receptor-positive spatial spots over a tissue image
#'
#' @param object Seurat object, or a SIRAP object containing both spatial and
#'   receptor data.
#' @param receptors AIRR-style receptor table containing a spatial barcode.
#'   May be omitted when `object` is a SIRAP object.
#' @param barcode_col Column in `receptors` matching Seurat spot names. If
#'   `NULL`, the first available barcode-like column is used.
#' @param feature Optional gene/feature to plot underneath receptor-positive
#'   spots.
#' @param gene Alias for `feature`, retained for readability in receptor plots.
#' @param receptor_color Color for receptor-positive spots.
#' @param receptor_size Size for receptor-positive spots.
#' @param receptor_stroke Stroke width for hollow receptor circles.
#' @param receptor_alpha Alpha for receptor circles.
#' @param feature_colors Colors used, in order, to build the continuous gene
#'   expression ramp when `feature` is set.
#' @param feature_size Size for gene expression points.
#' @param feature_alpha Alpha for gene expression points.
#' @param zero_tol Values at or below this threshold are treated as zero.
#' @param title Optional title.
#' @param order_points Plot lower feature expression first so higher expression
#'   points appear on top.
#' @param assay,layer,image.scale Additional plotting controls.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate axes.
#' @param images Optional Seurat image names.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param receptor_label Legend label for receptor-positive spots.
#' @param show_receptor_legend Show the receptor legend. Defaults to `TRUE`
#'   when a feature is shown and `FALSE` for receptor-only plots.
#' @return A ggplot/patchwork object.
#' @export
plotReceptorExpression <- function(object,
                                      receptors = NULL,
                                      barcode_col = NULL,
                                      feature = NULL,
                                      gene = NULL,
                                      receptor_color = "black",
                                      receptor_size = 1,
                                      receptor_stroke = 0.7,
                                      receptor_alpha = 1,
                                      feature_colors = c("#939393", "#FFFF00"),
                                      feature_size = 1,
                                      feature_alpha = 1,
                                      zero_tol = 0,
                                      title = NULL,
                                      order_points = TRUE,
                                      assay = NULL,
                                      layer = "counts",
                                      image.scale = "hires",
                                      swap_coords = FALSE,
                                      images = NULL,
                                      image.alpha = 1,
                                      zoom = FALSE,
                                      receptor_label = "Receptor",
                                      show_receptor_legend = NULL) {
  inputs <- .sirap_resolve_inputs(object, receptors, barcode_col)
  object <- inputs$object
  receptors <- inputs$receptors
  barcode_col <- inputs$barcode_col
  if (is.null(receptors)) {
    rlang::abort("`receptors` is required when `object` is not a SIRAP object.")
  }
  .require_namespace("SeuratObject")
  feature <- feature %||% gene
  show_receptor_legend <- show_receptor_legend %||% !is.null(feature)
  receptors <- .as_data_frame(receptors)
  barcode_col <- .guess_barcode_col(receptors, barcode_col)
  receptor_barcodes <- unique(stats::na.omit(as.character(receptors[[barcode_col]])))
  if (!is.null(assay)) {
    SeuratObject::DefaultAssay(object) <- assay
  }
  if (!is.null(feature) && length(feature_colors) < 2L) {
    rlang::abort("`feature_colors` must contain at least two colors.")
  }

  expr <- if (!is.null(feature)) {
    .sirap_fetch_data(object, vars = feature, layer = layer)
  } else {
    NULL
  }

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    dat <- spatial$coords
    dat$.sirap_receptor <- dat$cell %in% receptor_barcodes
    if (!is.null(feature)) {
      dat <- merge(dat, expr, by.x = "cell", by.y = 0, all.x = TRUE)
      dat$.sirap_expressed <- !is.na(dat[[feature]]) & dat[[feature]] > zero_tol
      keep <- dat$.sirap_expressed | dat$.sirap_receptor
    } else {
      keep <- dat$.sirap_receptor
    }
    if (!any(keep)) {
      next
    }
    dat <- dat[keep, , drop = FALSE]

    p <- .sirap_spatial_base(
      dat,
      spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      view_data = spatial$coords
    )

    if (!is.null(feature)) {
      gene_data <- dat[dat$.sirap_expressed, , drop = FALSE]
      if (nrow(gene_data) > 0L) {
        if (isTRUE(order_points)) {
          gene_data <- gene_data[order(gene_data[[feature]], na.last = TRUE), , drop = FALSE]
        }
        p <- p +
          ggplot2::geom_point(
            data = gene_data,
            ggplot2::aes(x = .data$x, y = .data$.sirap_y, color = .data[[feature]]),
            size = feature_size,
            alpha = feature_alpha,
            inherit.aes = FALSE,
            na.rm = TRUE
          ) +
          ggplot2::scale_color_gradientn(colors = feature_colors, name = feature)
      }
    }

    receptor_data <- dat[dat$.sirap_receptor, , drop = FALSE]
    if (nrow(receptor_data) > 0L) {
      if (isTRUE(show_receptor_legend)) {
        receptor_data$.sirap_receptor_label <- receptor_label
        p <- p +
          ggplot2::geom_point(
            data = receptor_data,
            ggplot2::aes(x = .data$x, y = .data$.sirap_y, shape = .data$.sirap_receptor_label),
            color = receptor_color,
            fill = NA,
            size = receptor_size,
            stroke = receptor_stroke,
            alpha = receptor_alpha,
            inherit.aes = FALSE,
            na.rm = TRUE
          ) +
          ggplot2::scale_shape_manual(
            values = stats::setNames(21, receptor_label),
            name = "Receptor"
          ) +
          ggplot2::guides(
            shape = ggplot2::guide_legend(
              override.aes = list(color = receptor_color, fill = NA, stroke = receptor_stroke)
            )
          )
      } else {
        p <- p +
          ggplot2::geom_point(
            data = receptor_data,
            ggplot2::aes(x = .data$x, y = .data$.sirap_y),
            color = receptor_color,
            fill = NA,
            size = receptor_size,
            stroke = receptor_stroke,
            shape = 21,
            alpha = receptor_alpha,
            inherit.aes = FALSE,
            na.rm = TRUE,
            show.legend = FALSE
          )
      }
    }
    plots[[image]] <- p + ggplot2::labs(title = title %||% feature %||% receptor_label)
  }

  .sirap_return_spatial_plots(plots, combine = TRUE)
}

#' Plot one or more receptor groups over a tissue image
#'
#' @param object Seurat object, or a SIRAP object containing both spatial and
#'   receptor data.
#' @param receptors AIRR-style receptor table containing a spatial barcode and
#'   grouping column. May be omitted when `object` is a SIRAP object.
#' @param group_col Column in `receptors` used to select and color receptor
#'   groups, for example `"clone_id"`.
#' @param values Optional values from `group_col` to plot. Defaults to all
#'   non-missing values in `group_col`.
#' @param barcode_col Column in `receptors` matching Seurat spot names. If
#'   `NULL`, the first available barcode-like column is used.
#' @param colors Optional colors for receptor groups. Named colors are matched
#'   to `values`; otherwise colors are used in value order.
#' @param receptor_size Size for receptor points.
#' @param receptor_alpha Alpha for receptor points.
#' @param title Optional title.
#' @param order_points Plot earlier receptor groups first so later groups appear
#'   on top.
#' @param image.scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param images Optional Seurat image names.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param receptor_legend_title Legend title for receptor groups. Defaults to
#'   `group_col`.
#' @param show_receptor_legend Show the receptor group legend.
#' @param combine Combine image panels with patchwork.
#' @return A ggplot/patchwork object or a list when `combine = FALSE`.
#' @export
plotReceptors <- function(object,
                          receptors = NULL,
                          group_col,
                          values = NULL,
                          barcode_col = NULL,
                          colors = NULL,
                          receptor_size = 1,
                          receptor_alpha = 1,
                          title = NULL,
                          order_points = TRUE,
                          image.scale = "hires",
                          swap_coords = FALSE,
                          images = NULL,
                          image.alpha = 1,
                          zoom = FALSE,
                          receptor_legend_title = NULL,
                          show_receptor_legend = TRUE,
                          combine = TRUE) {
  inputs <- .sirap_resolve_inputs(object, receptors, barcode_col)
  object <- inputs$object
  receptors <- inputs$receptors
  barcode_col <- inputs$barcode_col
  if (is.null(receptors)) {
    rlang::abort("`receptors` is required when `object` is not a SIRAP object.")
  }
  .require_namespace("SeuratObject")
  receptors <- .as_data_frame(receptors)
  barcode_col <- .guess_barcode_col(receptors, barcode_col)
  .check_required_cols(receptors, c(barcode_col, group_col), "receptors")

  receptors$.sirap_receptor_group <- as.character(receptors[[group_col]])
  if (is.null(values)) {
    values <- unique(stats::na.omit(receptors$.sirap_receptor_group))
  }
  values <- as.character(values)
  values <- values[!is.na(values)]
  if (length(values) == 0L) {
    rlang::abort("`values` must contain at least one non-missing receptor group.")
  }

  receptors <- receptors[
    !is.na(receptors$.sirap_receptor_group) &
      receptors$.sirap_receptor_group %in% values &
      !is.na(receptors[[barcode_col]]),
    ,
    drop = FALSE
  ]
  if (nrow(receptors) == 0L) {
    rlang::abort("No receptor rows matched `values` in `group_col`.")
  }
  receptors[[barcode_col]] <- as.character(receptors[[barcode_col]])
  receptors$.sirap_receptor_group <- factor(receptors$.sirap_receptor_group, levels = values)
  colors <- .sirap_receptor_group_colors(values, colors)
  receptor_legend_title <- receptor_legend_title %||% group_col

  receptor_lookup <- receptors[, c(barcode_col, ".sirap_receptor_group"), drop = FALSE]
  receptor_lookup <- receptor_lookup[!duplicated(receptor_lookup), , drop = FALSE]
  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    receptor_data <- merge(
      spatial$coords,
      receptor_lookup,
      by.x = "cell",
      by.y = barcode_col,
      all = FALSE
    )
    receptor_data$.sirap_receptor_group <- factor(
      receptor_data$.sirap_receptor_group,
      levels = values
    )
    if (isTRUE(order_points)) {
      receptor_data <- receptor_data[
        order(as.integer(receptor_data$.sirap_receptor_group), na.last = TRUE),
        ,
        drop = FALSE
      ]
    }

    if (nrow(receptor_data) == 0L) {
      next
    }

    p <- .sirap_spatial_base(
      receptor_data,
      spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      view_data = spatial$coords
    )

    p <- p +
      ggplot2::geom_point(
        data = receptor_data,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, color = .data$.sirap_receptor_group),
        size = receptor_size,
        alpha = receptor_alpha,
        inherit.aes = FALSE,
        na.rm = TRUE,
        show.legend = isTRUE(show_receptor_legend)
      ) +
      ggplot2::scale_color_manual(
        values = colors,
        name = receptor_legend_title,
        drop = FALSE
      )

    plots[[image]] <- p + ggplot2::labs(title = title %||% receptor_legend_title)
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

.sirap_receptor_group_colors <- function(values, colors = NULL) {
  n <- length(values)
  if (is.null(colors)) {
    if (requireNamespace("pals", quietly = TRUE)) {
      colors <- .sirap_extend_palette(pals::glasbey(32), n)
    } else {
      colors <- .sirap_extend_palette(grDevices::hcl.colors(max(3L, n), palette = "Dark 3"), n)
    }
    names(colors) <- values
    return(colors)
  }

  colors <- as.character(colors)
  if (!is.null(names(colors)) && all(values %in% names(colors))) {
    return(colors[values])
  }
  if (length(colors) < n) {
    colors <- .sirap_extend_palette(colors, n)
  }
  colors <- colors[seq_len(n)]
  names(colors) <- values
  colors
}

.guess_barcode_col <- function(data, barcode_col = NULL) {
  if (!is.null(barcode_col)) {
    .check_required_cols(data, barcode_col, "receptors")
    return(barcode_col)
  }

  candidates <- c(
    "barcode",
    "barcode_8um",
    "barcode_2um",
    "spatial_barcode",
    "cell",
    "cell_barcode",
    "sequence_id"
  )
  barcode_col <- candidates[candidates %in% colnames(data)][1]
  if (is.na(barcode_col)) {
    rlang::abort("Could not guess `barcode_col`; pass the receptor barcode column explicitly.")
  }
  barcode_col
}

#' Plot tip colors over a tissue image
#'
#' @param object Seurat object.
#' @param tip_data Data frame with spatial barcode and color columns.
#' @param barcode_col Column in `tip_data` matching Seurat spot names.
#' @param color_col Column containing colors.
#' @param label_col Optional column containing legend labels for colors.
#' @param legend_title Title for the color legend.
#' @param legend_colors Optional named color vector used for the legend.
#' @param title Optional title.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param point_size Size of plotted tips.
#' @param point_alpha Opacity for plotted tips.
#' @param order_points Plot earlier tip-color groups first so later groups
#'   appear on top.
#' @param images Optional Seurat image names.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param show.legend Show color legend.
#' @return A ggplot/patchwork object.
#' @export
plotSpatialColorTips <- function(object,
                                tip_data,
                                barcode_col = "barcode_8um",
                                color_col = "tip_color",
                                label_col = NULL,
                                legend_title = "Tip color",
                                legend_colors = NULL,
                                title = NULL,
                                image.scale = "hires",
                                swap_coords = FALSE,
                                point_size = 1.6,
                                point_alpha = 1,
                                order_points = TRUE,
                                images = NULL,
                                image.alpha = 1,
                                zoom = FALSE,
                                show.legend = FALSE) {
  .require_namespace("SeuratObject")
  tip_data <- .as_data_frame(tip_data)
  .check_required_cols(tip_data, c(barcode_col, color_col), "tip_data")
  if (!is.null(label_col)) {
    .check_required_cols(tip_data, label_col, "tip_data")
  }
  tip_data <- tip_data[!is.na(tip_data[[barcode_col]]) & !is.na(tip_data[[color_col]]), , drop = FALSE]
  spot_color <- stats::setNames(as.character(tip_data[[color_col]]), as.character(tip_data[[barcode_col]]))
  spot_color <- spot_color[!duplicated(names(spot_color))]
  if (is.null(label_col)) {
    colors <- unique(stats::na.omit(as.character(spot_color)))
    names(colors) <- colors
  } else {
    tip_data <- tip_data[!is.na(tip_data[[label_col]]), , drop = FALSE]
    spot_label <- stats::setNames(as.character(tip_data[[label_col]]), as.character(tip_data[[barcode_col]]))
    spot_label <- spot_label[!duplicated(names(spot_label))]
    if (is.null(legend_colors)) {
      legend_df <- tip_data[!duplicated(tip_data[[label_col]]), c(label_col, color_col), drop = FALSE]
      legend_colors <- stats::setNames(as.character(legend_df[[color_col]]), as.character(legend_df[[label_col]]))
    }
    colors <- legend_colors
  }

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale, swap_coords = swap_coords)
    dat <- spatial$coords[spatial$coords$cell %in% names(spot_color), , drop = FALSE]
    if (nrow(dat) == 0L) {
      next
    }
    dat$.sirap_tip_color <- spot_color[dat$cell]
    if (!is.null(label_col)) {
      dat$.sirap_tip_label <- factor(spot_label[dat$cell], levels = names(colors))
      dat$.sirap_tip_legend <- dat$.sirap_tip_label
    } else {
      dat$.sirap_tip_legend <- factor(dat$.sirap_tip_color, levels = names(colors))
    }
    if (isTRUE(order_points)) {
      dat <- dat[order(as.integer(dat$.sirap_tip_legend), na.last = TRUE), , drop = FALSE]
    }
    p <- .sirap_spatial_base(
      dat,
      spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$.sirap_tip_legend),
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = colors,
        name = legend_title,
        breaks = names(colors),
        limits = names(colors),
        drop = FALSE
      ) +
      ggplot2::labs(title = title)
    if (!isTRUE(show.legend)) {
      p <- p + ggplot2::theme(legend.position = "none")
    }
    plots[[image]] <- p
  }

  .sirap_return_spatial_plots(plots, combine = TRUE)
}

#' Plot a clone with tip colors over a spatial tissue image
#'
#' @param colored_clones Output from `assignTipColors()`.
#' @param clone_id Clone ID to plot.
#' @param tissue_positions Seurat object with spatial images.
#' @param clone_col Clone ID column.
#' @param barcode_col Spatial barcode column in `colored_clones`.
#' @param ... Additional arguments passed to `plotSpatialColorTips()`.
#' @return A ggplot object.
#' @export
plotSpatialCloneColorTips <- function(colored_clones,
                                    clone_id,
                                    tissue_positions = NULL,
                                    clone_col = "clone_id",
                                    barcode_col = "barcode_8um",
                                    ...) {
  .check_required_cols(colored_clones, c(clone_col, barcode_col, "tip_color"), "colored_clones")
  clone_df <- .as_data_frame(colored_clones)
  clone_df <- clone_df[clone_df[[clone_col]] == clone_id, , drop = FALSE]

  if (!inherits(tissue_positions, "Seurat")) {
    rlang::abort("plotSpatialCloneColorTips() expects a Seurat object. Use ggplot directly for coordinate data frames.")
  }
  plotSpatialColorTips(
    tissue_positions,
    tip_data = clone_df,
    barcode_col = barcode_col,
    color_col = "tip_color",
    title = paste("Clone", clone_id),
    ...
  )
}

#' Plot paired BCR lineage trees and matching spatial clone colors
#'
#' Tips without spatial coordinates or spatial barcodes can still appear in the
#' tree panel. The spatial panel only plots tips whose `barcode_col` matches a
#' spot in `object`.
#'
#' @param trees Dowser `getTrees()` output, or a list containing a `trees`
#'   element with that output.
#' @param object Seurat object with spatial images.
#' @param barcode_col Column in the dowser clone data matching Seurat spot names.
#' @param clone_ids Optional clone IDs to plot.
#' @param id_col Column in the dowser clone data matching tree tip labels.
#' @param clone_col Clone ID column.
#' @param color_by Column used to assign matching tree and spatial colors.
#'   Use `"identity"` to give each tip its own color. Use another column name,
#'   such as `"tree_cluster"` or a spatial barcode column, to color tips by
#'   groups in that column.
#' @param palette Optional colors. If fewer colors are supplied than needed,
#'   additional colors are generated by averaging supplied colors.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. Defaults to `FALSE`,
#'   which uses Seurat image coordinates natively.
#' @param images Optional Seurat image names.
#' @param zoom Zoom spatial plots to the clone. Defaults to `FALSE`, which
#'   shows the full tissue area covered by spatial coordinates.
#' @param tree_layout Layout passed to `ggtree::ggtree()`.
#' @param tree_tip_size Size of tree tip points.
#' @param treescale Mutational distance shown in the tree scale bar.
#' @param spatial_tip_size Size of spatial tip points.
#' @param spatial_tip_alpha Opacity for spatial tip points.
#' @param order_points Plot earlier spatial tip-color groups first so later
#'   groups appear on top.
#' @param show_spatial_legend Show a color legend based on the tree tip colors
#'   below the spatial plot.
#' @param nproc Number of processes used to build clone tree/spatial plot pairs.
#' @param ... Additional arguments passed to `plotSpatialColorTips()`.
#' @return Named list. Each clone element contains `tree` and `spatial_clone`.
#' @export
plotSpatialCloneTreePairs <- function(trees,
                                      object,
                                      barcode_col,
                                      clone_ids = NULL,
                                      id_col = "sequence_id",
                                      clone_col = "clone_id",
                                      color_by = "identity",
                                      palette = NULL,
                                      image.scale = "hires",
                                      swap_coords = FALSE,
                                      images = NULL,
                                      zoom = FALSE,
                                      tree_layout = "rectangular",
                                      tree_tip_size = 2.5,
                                      treescale = 0.01,
                                      spatial_tip_size = 1.6,
                                      spatial_tip_alpha = 1,
                                      order_points = TRUE,
                                      show_spatial_legend = FALSE,
                                      nproc = 1L,
                                      ...) {
  .require_namespace("ggtree")
  nproc <- .sirap_plot_check_nproc(nproc)
  records <- .sirap_dowser_tree_records(trees)
  if (!is.null(clone_ids)) {
    records <- records[names(records) %in% as.character(clone_ids)]
  }
  if (length(records) == 0L) {
    rlang::abort("No dowser trees matched the requested clone IDs.")
  }

  plot_one_clone <- function(clone_id) {
    record <- records[[clone_id]]
    tree <- record$tree
    clone_data <- .sirap_airr_clone_data(record$data, clone_id = clone_id, clone_col = clone_col)
    .check_required_cols(clone_data, c(id_col, barcode_col, clone_col), "dowser clone data")

    tip_ids <- setdiff(tree$tip.label, "Germline")
    tip_data <- clone_data[clone_data[[id_col]] %in% tip_ids, , drop = FALSE]
    if (nrow(tip_data) == 0L) {
      return(NULL)
    }
    shared_color_col <- NULL
    if (!identical(color_by, "identity")) {
      .check_required_cols(tip_data, color_by, "dowser clone data")
      shared_color_col <- color_by
    }
    colored_tips <- assignTipColors(
      tip_data,
      clone_ids = clone_id,
      clone_col = clone_col,
      id_col = id_col,
      shared_color_col = shared_color_col,
      palette = palette
    )
    legend_info <- .sirap_spatial_tip_legend_info(colored_tips, color_by = color_by)
    colored_tips <- legend_info$tip_data

    tree_plot <- .sirap_plot_colored_tree(
      tree,
      colored_tips,
      id_col = id_col,
      color_col = "tip_color",
      title = paste("Clone", clone_id),
      layout = tree_layout,
      tip_size = tree_tip_size,
      treescale = treescale,
      label_col = ".sirap_tip_legend_label",
      legend_title = "Tip cluster",
      legend_colors = legend_info$colors,
      show.legend = show_spatial_legend,
      legend.position = "bottom",
      legend.ncol = min(3L, length(legend_info$colors))
    )
    tree_legend <- if (isTRUE(show_spatial_legend)) {
      cowplot::get_legend(tree_plot)
    } else {
      NULL
    }
    tree_plot <- tree_plot + ggplot2::theme(legend.position = "none")

    spatial_plot <- plotSpatialColorTips(
      object,
      tip_data = colored_tips,
      barcode_col = barcode_col,
      color_col = "tip_color",
      title = NULL,
      image.scale = image.scale,
      swap_coords = swap_coords,
      images = images,
      point_size = spatial_tip_size,
      point_alpha = spatial_tip_alpha,
      order_points = order_points,
      zoom = zoom,
      show.legend = FALSE,
      ...
    )
    if (!is.null(tree_legend)) {
      spatial_plot <- cowplot::plot_grid(
        spatial_plot,
        tree_legend,
        ncol = 1,
        rel_heights = c(1, max(0.18, 0.08 * ceiling(length(legend_info$colors) / 3)))
      )
    }

    list(
      tree = tree_plot,
      spatial_clone = spatial_plot
    )
  }

  out <- .sirap_plot_lapply(names(records), plot_one_clone, nproc = nproc)
  names(out) <- names(records)
  out <- out[!vapply(out, is.null, logical(1))]

  out
}

.sirap_spatial_tip_legend_info <- function(tip_data, color_by) {
  tip_data$.sirap_tip_legend_label <- .sirap_spatial_tip_labels(tip_data, color_by = color_by)
  legend_df <- tip_data[
    !is.na(tip_data$.sirap_tip_legend_label) &
      !is.na(tip_data$tip_color),
    c(".sirap_tip_legend_label", "tip_color", "tip_color_id"),
    drop = FALSE
  ]
  legend_df <- legend_df[order(as.integer(legend_df$tip_color_id)), , drop = FALSE]
  legend_df <- legend_df[!duplicated(legend_df$.sirap_tip_legend_label), , drop = FALSE]

  spatial <- grepl("^spatial_tip_cluster_[0-9]+$", legend_df$.sirap_tip_legend_label)
  spatial_df <- legend_df[spatial, , drop = FALSE]
  if (nrow(spatial_df) > 0L) {
    spatial_id <- as.integer(sub("^spatial_tip_cluster_", "", spatial_df$.sirap_tip_legend_label))
    spatial_df <- spatial_df[order(spatial_id), , drop = FALSE]
  }
  other_df <- legend_df[!spatial, , drop = FALSE]
  legend_df <- rbind(spatial_df, other_df)

  list(
    tip_data = tip_data,
    colors = stats::setNames(
      as.character(legend_df$tip_color),
      as.character(legend_df$.sirap_tip_legend_label)
    )
  )
}

.sirap_spatial_tip_labels <- function(tip_data, color_by) {
  if (!identical(color_by, "identity") && color_by %in% colnames(tip_data)) {
    values <- as.character(tip_data[[color_by]])
    spatial_cluster <- grepl("^tree_cluster_[0-9]+$", values)
    values[spatial_cluster] <- paste0(
      "spatial_tip_cluster_",
      sub("^tree_cluster_", "", values[spatial_cluster])
    )
    values
  } else {
    paste0("spatial_tip_cluster_", as.integer(tip_data$tip_color_id))
  }
}

.sirap_plot_check_nproc <- function(nproc) {
  nproc <- as.integer(nproc)
  if (is.na(nproc) || nproc < 1L) {
    rlang::abort("`nproc` must be a positive integer.")
  }
  nproc
}

.sirap_plot_lapply <- function(x, fun, nproc) {
  if (nproc == 1L || length(x) <= 1L) {
    return(lapply(x, fun))
  }
  if (.Platform$OS.type == "windows") {
    warning("Parallel plotting with `nproc > 1` is not supported on Windows; using `nproc = 1`.", call. = FALSE)
    return(lapply(x, fun))
  }
  parallel::mclapply(x, fun, mc.cores = nproc)
}

#' Save paired clone tree and spatial plots to a PDF
#'
#' @param pairs Output from `plotSpatialCloneTreePairs()`.
#' @param path PDF path.
#' @param width,height PDF page size in inches.
#' @param rel_widths Relative widths for tree and spatial panels.
#' @param nproc Number of processes used to render PDF pages. Parallel PDF
#'   rendering requires the optional `qpdf` package to combine pages; otherwise
#'   the function falls back to serial rendering.
#' @return Invisibly returns `path`.
#' @export
spatialCloneTreePairsToPDF <- function(pairs,
                                       path,
                                       width = 12,
                                       height = 6,
                                       rel_widths = c(1.2, 1),
                                       nproc = 1L) {
  if (length(pairs) == 0L) {
    rlang::abort("`pairs` does not contain any clone tree/spatial plots.")
  }
  nproc <- .sirap_plot_check_nproc(nproc)
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (nproc > 1L && .Platform$OS.type != "windows" && requireNamespace("qpdf", quietly = TRUE)) {
    tmp_dir <- tempfile("sirap_clone_tree_pairs_pdf_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

    page_paths <- file.path(tmp_dir, sprintf("page_%05d.pdf", seq_along(pairs)))
    .sirap_plot_lapply(seq_along(pairs), function(i) {
      .sirap_write_clone_tree_pair_pdf(
        pair = pairs[[i]],
        path = page_paths[[i]],
        width = width,
        height = height,
        rel_widths = rel_widths
      )
      page_paths[[i]]
    }, nproc = nproc)

    qpdf::pdf_combine(input = page_paths, output = path)
    return(invisible(path))
  }

  if (nproc > 1L) {
    warning(
      "Parallel PDF rendering requires a non-Windows platform and the optional `qpdf` package; using `nproc = 1`.",
      call. = FALSE
    )
  }

  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (pair in pairs) {
    print(.sirap_clone_tree_pair_grid(pair, rel_widths = rel_widths))
  }

  invisible(path)
}

.sirap_write_clone_tree_pair_pdf <- function(pair,
                                             path,
                                             width,
                                             height,
                                             rel_widths) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(.sirap_clone_tree_pair_grid(pair, rel_widths = rel_widths))
  invisible(path)
}

.sirap_clone_tree_pair_grid <- function(pair, rel_widths) {
  cowplot::plot_grid(
    pair$tree,
    pair$spatial_clone,
    ncol = 2,
    rel_widths = rel_widths
  )
}

.sirap_dowser_tree_records <- function(trees) {
  if (is.list(trees) && !is.data.frame(trees) && "trees" %in% names(trees)) {
    trees <- trees$trees
  }
  if (is.data.frame(trees) && all(c("trees", "data") %in% colnames(trees))) {
    records <- lapply(seq_len(nrow(trees)), function(i) {
      clone_id <- if ("clone_id" %in% colnames(trees)) as.character(trees$clone_id[[i]]) else NULL
      tree <- trees$trees[[i]]
      clone_data <- trees$data[[i]]
      clone_id <- clone_id %||% tree$name %||% .sirap_airr_clone_id(clone_data)
      list(tree = tree, data = clone_data)
    })
    if ("clone_id" %in% colnames(trees)) {
      names(records) <- as.character(trees$clone_id)
    } else {
      names(records) <- vapply(seq_along(records), function(i) {
        records[[i]]$tree$name %||% .sirap_airr_clone_id(records[[i]]$data) %||% paste0("clone_", i)
      }, character(1))
    }
    return(records)
  }
  if (inherits(trees, "phylo")) {
    return(stats::setNames(list(list(tree = trees, data = NULL)), trees$name %||% "clone_1"))
  }
  if (is.list(trees)) {
    records <- lapply(seq_along(trees), function(i) {
      x <- trees[[i]]
      tree <- if (inherits(x, "phylo")) x else x$tree %||% x$trees
      if (is.list(tree) && length(tree) == 1L && inherits(tree[[1]], "phylo")) {
        tree <- tree[[1]]
      }
      clone_data <- x$data %||% x$clone_data %||% NULL
      list(tree = tree, data = clone_data)
    })
    record_names <- names(trees)
    if (is.null(record_names) || any(!nzchar(record_names))) {
      record_names <- vapply(seq_along(records), function(i) {
        records[[i]]$tree$name %||% paste0("clone_", i)
      }, character(1))
    }
    names(records) <- record_names
    return(records)
  }
  rlang::abort("Could not interpret `trees` as dowser getTrees() output.")
}

.sirap_airr_clone_id <- function(x) {
  if (inherits(x, "airrClone")) {
    return(x@clone)
  }
  NULL
}

.sirap_airr_clone_data <- function(x, clone_id, clone_col) {
  if (inherits(x, "airrClone")) {
    out <- x@data
    clone_id <- x@clone %||% clone_id
  } else {
    out <- .as_data_frame(x)
  }
  if (!clone_col %in% colnames(out)) {
    out[[clone_col]] <- clone_id
  }
  out
}

.sirap_plot_colored_tree <- function(tree,
                                     tip_data,
                                     id_col,
                                     color_col,
                                     title,
                                     layout,
                                     tip_size,
                                     treescale,
                                     label_col = NULL,
                                     legend_title = "Tip cluster",
                                     legend_colors = NULL,
                                     show.legend = FALSE,
                                     legend.position = "right",
                                     legend.ncol = NULL) {
  p <- ggtree::ggtree(tree, layout = layout)
  join_cols <- c(id_col, color_col)
  if (!is.null(label_col)) {
    join_cols <- c(join_cols, label_col)
  }
  join_data <- tip_data[, join_cols, drop = FALSE]
  join_data <- join_data[!duplicated(join_data[[id_col]]), , drop = FALSE]
  idx <- match(p$data$label, join_data[[id_col]])
  p$data[[color_col]] <- join_data[[color_col]][idx]
  if (!is.null(label_col)) {
    p$data[[label_col]] <- factor(join_data[[label_col]][idx], levels = names(legend_colors))
  }
  color_aes_col <- if (is.null(label_col)) color_col else label_col

  p <- p +
    ggtree::geom_tippoint(
      ggplot2::aes(color = .data[[color_aes_col]]),
      size = tip_size,
      na.rm = TRUE
    ) +
    ggtree::geom_treescale(width = treescale) +
    ggplot2::ggtitle(title)

  if (is.null(label_col)) {
    p <- p + ggplot2::scale_color_identity()
  } else {
    p <- p +
      ggplot2::scale_color_manual(
        values = legend_colors,
        name = legend_title,
        breaks = names(legend_colors),
        limits = names(legend_colors),
        drop = FALSE
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          ncol = legend.ncol,
          byrow = TRUE,
          override.aes = list(size = tip_size)
        )
      )
  }

  p +
    ggplot2::theme(
      legend.position = if (isTRUE(show.legend)) legend.position else "none",
      plot.title = ggplot2::element_text(hjust = 0)
    )
}

.sirap_spatial_images <- function(object, images = NULL) {
  .require_namespace("SeuratObject")
  images <- images %||% SeuratObject::Images(object)
  if (length(images) == 0L) {
    rlang::abort("Could not find any spatial image information in `object`.")
  }
  missing_images <- setdiff(images, SeuratObject::Images(object))
  if (length(missing_images) > 0L) {
    rlang::abort(paste0("Could not find image(s): ", paste(missing_images, collapse = ", ")))
  }
  images
}

.sirap_fetch_data <- function(object, vars, layer = "counts") {
  out <- tryCatch(
    SeuratObject::FetchData(object, vars = vars, layer = layer),
    error = function(e) {
      SeuratObject::FetchData(object, vars = vars, slot = layer)
    }
  )
  as.data.frame(out, stringsAsFactors = FALSE)
}

.sirap_spatial_data <- function(object, image, image.scale = "hires", swap_coords = FALSE) {
  .require_namespace("SeuratObject")
  image_obj <- object[[image]]
  if (isTRUE(swap_coords)) {
    coords <- .sirap_swapped_spatial_coordinates(image_obj, image.scale = image.scale)
  } else {
    coords <- .sirap_tissue_spatial_coordinates(image_obj, image.scale = image.scale)
  }
  if (is.null(coords)) {
    rlang::abort("Could not identify x/y coordinate columns in Seurat tissue coordinates.")
  }
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  coords <- coords[coords$cell %in% colnames(object), , drop = FALSE]
  coords$.sirap_y <- -coords$y

  list(coords = coords, image = image_obj)
}

.sirap_tissue_spatial_coordinates <- function(image, image.scale = "hires") {
  coords <- SeuratObject::GetTissueCoordinates(image, scale = image.scale)
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  if ("cell" %in% colnames(coords)) {
    coords$cell <- as.character(coords$cell)
  } else {
    coords$cell <- rownames(coords)
  }

  x_col <- intersect(c("imagecol", "x", "pxl_col_in_fullres", "col"), colnames(coords))
  y_col <- intersect(c("imagerow", "y", "pxl_row_in_fullres", "row"), colnames(coords))
  if (length(x_col) == 0L || length(y_col) == 0L) {
    return(NULL)
  }
  data.frame(
    cell = coords$cell,
    x = coords[[x_col[[1L]]]],
    y = coords[[y_col[[1L]]]],
    stringsAsFactors = FALSE
  )
}

.sirap_swapped_spatial_coordinates <- function(image, image.scale = "hires") {
  coords <- SeuratObject::GetTissueCoordinates(image, scale = image.scale)
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  if ("cell" %in% colnames(coords)) {
    coords$cell <- as.character(coords$cell)
  } else {
    coords$cell <- rownames(coords)
  }
  if (all(c("x", "y") %in% colnames(coords))) {
    x <- coords$x
    y <- coords$y
  } else if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
    x <- coords$imagecol
    y <- coords$imagerow
  } else {
    rlang::abort(
      "`swap_coords = TRUE` requires `SeuratObject::GetTissueCoordinates()` to return `x`/`y` or `imagecol`/`imagerow` coordinates."
    )
  }
  data.frame(
    cell = coords$cell,
    x = suppressWarnings(as.numeric(y)),
    y = suppressWarnings(as.numeric(x)),
    stringsAsFactors = FALSE
  )
}

.sirap_spatial_coordinate_candidates <- function(coords,
                                                image,
                                                image.scale = "hires",
                                                include_all_scales = FALSE) {
  raster <- .sirap_image_raster(image, alpha = 1)
  image_dim <- dim(raster)
  height <- image_dim[[1L]]
  width <- image_dim[[2L]]
  scale_factors <- .sirap_image_scale_factors(image)
  requested_scale_factor <- .sirap_image_scale_factor(image, image.scale)
  scale_transforms <- c(raw = 1)
  if (is.finite(requested_scale_factor) && requested_scale_factor > 0) {
    scale_transforms <- c(requested = requested_scale_factor, scale_transforms)
  }
  if (isTRUE(include_all_scales) && length(scale_factors) > 0L) {
    scale_transforms <- c(scale_transforms, scale_factors)
  }
  scale_transforms <- scale_transforms[
    is.finite(scale_transforms) & scale_transforms > 0 & !duplicated(scale_transforms)
  ]

  pairs <- list(
    c("x", "y"),
    c("imagecol", "imagerow"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    c("col", "row")
  )
  pairs <- Filter(function(pair) all(pair %in% colnames(coords)), pairs)
  candidates <- list()
  for (pair_idx in seq_along(pairs)) {
    pair <- pairs[[pair_idx]]
    raw_x <- suppressWarnings(as.numeric(coords[[pair[[1L]]]]))
    raw_y <- suppressWarnings(as.numeric(coords[[pair[[2L]]]]))
    for (transform_idx in seq_along(scale_transforms)) {
      scale_transform <- scale_transforms[[transform_idx]]
      candidates[[length(candidates) + 1L]] <- .sirap_spatial_coordinate_candidate(
        x = raw_x * scale_transform,
        y = raw_y * scale_transform,
        pair_idx = pair_idx,
        transform_idx = transform_idx,
        width = width,
        height = height
      )
    }
  }

  candidates <- Filter(Negate(is.null), candidates)
  if (length(candidates) == 0L) {
    return(list())
  }
  scores <- vapply(candidates, function(candidate) {
    candidate$score
  }, numeric(4L))
  candidates[order(scores[1L, ], scores[2L, ], scores[3L, ], scores[4L, ])]
}

.sirap_spatial_coordinate_candidate <- function(x,
                                                y,
                                                pair_idx,
                                                transform_idx,
                                                width,
                                                height) {
  finite <- is.finite(x) & is.finite(y)
  if (!any(finite)) {
    return(NULL)
  }
  overflow_x <- pmax(0, -x[finite]) + pmax(0, x[finite] - width)
  overflow_y <- pmax(0, -y[finite]) + pmax(0, y[finite] - height)
  overflow <- mean(overflow_x / max(width, 1) + overflow_y / max(height, 1))
  inside <- mean(x[finite] >= 0 & x[finite] <= width & y[finite] >= 0 & y[finite] <= height)
  span <- max(
    diff(range(x[finite])) / max(width, 1),
    diff(range(y[finite])) / max(height, 1)
  )
  list(
    x = x,
    y = y,
    score = c(-inside, overflow, -span, pair_idx + transform_idx / 10)
  )
}

.sirap_image_scale_factor <- function(image, image.scale = "hires") {
  if (is.null(image.scale) || length(image.scale) == 0L || is.na(image.scale[[1L]])) {
    return(NA_real_)
  }
  scale_factors <- tryCatch(methods::slot(image, "scale.factors"), error = function(e) NULL)
  if (is.null(scale_factors)) {
    return(NA_real_)
  }
  scale_value <- tryCatch(scale_factors[[image.scale[[1L]]]], error = function(e) NULL)
  if (is.null(scale_value)) {
    scale_value <- tryCatch(scale_factors[[paste0(image.scale[[1L]], ".scale")]], error = function(e) NULL)
  }
  scale_value <- suppressWarnings(as.numeric(scale_value))
  if (length(scale_value) == 0L || !is.finite(scale_value[[1L]])) {
    return(NA_real_)
  }
  scale_value[[1L]]
}

.sirap_image_scale_factors <- function(image) {
  scale_factors <- tryCatch(methods::slot(image, "scale.factors"), error = function(e) NULL)
  if (is.null(scale_factors)) {
    return(numeric())
  }
  scale_values <- tryCatch(unlist(as.list(scale_factors), use.names = TRUE), error = function(e) NULL)
  if (is.null(scale_values)) {
    return(numeric())
  }
  scale_names <- names(scale_values)
  scale_values <- suppressWarnings(as.numeric(scale_values))
  names(scale_values) <- scale_names
  scale_values[is.finite(scale_values) & scale_values > 0]
}

.sirap_spatial_base <- function(plot_data,
                                image,
                                image.alpha = 1,
                                zoom = FALSE,
                                view_data = NULL) {
  raster <- .sirap_image_raster(image, alpha = image.alpha)
  image_dim <- dim(raster)
  height <- image_dim[[1]]
  width <- image_dim[[2]]

  xlim <- c(0, width)
  ylim <- c(-height, 0)
  limit_data <- if (isTRUE(zoom)) plot_data else view_data %||% plot_data
  if (!is.null(limit_data) && nrow(limit_data) > 0L) {
    xlim <- .sirap_pad_range(limit_data$x, fallback = xlim)
    ylim <- .sirap_pad_range(limit_data$.sirap_y, fallback = ylim)
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x, y = .data$.sirap_y)) +
    ggplot2::annotation_raster(raster, xmin = 0, xmax = width, ymin = -height, ymax = 0, interpolate = FALSE) +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "right"
    )
}

.sirap_image_raster <- function(image, alpha = 1) {
  .require_namespace("SeuratObject")
  img <- tryCatch(methods::slot(image, "image"), error = function(e) NULL)
  if (is.null(img) || length(dim(img)) < 2L) {
    img <- tryCatch(SeuratObject::GetImage(image, mode = "raw"), error = function(e) NULL)
  }
  if (is.null(img) || length(dim(img)) < 2L) {
    rlang::abort("Could not extract image raster from Seurat image object.")
  }
  if (alpha < 1) {
    rgba <- array(NA_real_, dim = c(dim(img)[1:2], 4))
    rgba[, , 1:3] <- img[, , 1:3]
    rgba[, , 4] <- alpha
    img <- rgba
  }
  grDevices::as.raster(img)
}

.sirap_pad_range <- function(x, fallback) {
  out <- range(x, na.rm = TRUE)
  if (!all(is.finite(out))) {
    return(fallback)
  }
  if (out[[1]] == out[[2]]) {
    return(out + c(-1, 1))
  }
  pad <- diff(out) * 0.04
  out + c(-pad, pad)
}

.sirap_return_spatial_plots <- function(plots, combine = TRUE) {
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0L) {
    return(NULL)
  }
  if (!isTRUE(combine) || length(plots) == 1L) {
    return(if (length(plots) == 1L && isTRUE(combine)) plots[[1L]] else plots)
  }
  .require_namespace("patchwork")
  patchwork::wrap_plots(plots, ncol = 1)
}

.sirap_blend_color_grid <- function(bottom_left,
                                    bottom_right,
                                    top_left,
                                    top_right,
                                    side_length = 100L) {
  grad <- function(start, end) grDevices::colorRampPalette(c(start, end))(side_length)
  bottom <- grad(bottom_left, bottom_right)
  top <- grad(top_left, top_right)
  grid <- lapply(seq_along(bottom), function(i) grad(bottom[[i]], top[[i]]))
  t(matrix(unlist(grid), ncol = side_length, nrow = side_length))
}

.sirap_blend_colors <- function(dat, color_grid) {
  side_length <- nrow(color_grid)
  dat_norm <- apply(dat, 2, function(x) {
    if (max(x, na.rm = TRUE) <= 0) {
      return(rep(1L, length(x)))
    }
    idx <- round((side_length - 1L) * x / max(x, na.rm = TRUE)) + 1L
    pmax(1L, pmin(side_length, idx))
  })
  vapply(seq_len(nrow(dat_norm)), function(i) {
    color_grid[dat_norm[i, 1], dat_norm[i, 2]]
  }, character(1))
}

.sirap_blend_legend <- function(color_grid, xlab, ylab) {
  side_length <- nrow(color_grid)
  legend_grid <- expand.grid(
    x = seq(0, 1, length.out = side_length),
    y = seq(0, 1, length.out = side_length)
  )
  legend_grid$color <- c(color_grid)
  legend_colors <- stats::setNames(legend_grid$color, legend_grid$color)
  ggplot2::ggplot(
    legend_grid,
    ggplot2::aes(x = .data$x, y = .data$y, color = .data$color)
  ) +
    ggplot2::geom_point(shape = 15, size = 1.9) +
    ggplot2::scale_color_manual(values = legend_colors) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      aspect.ratio = 1,
      axis.title = ggplot2::element_text(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
}
