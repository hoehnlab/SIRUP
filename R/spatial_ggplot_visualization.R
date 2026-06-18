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
#' @param slot Deprecated Seurat slot alias for `layer`.
#' @param images Optional Seurat image names.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param feature_colors Colors used, in order, for the continuous expression
#'   scale. The first color is the low value and the last color is the high
#'   value.
#' @param zero_tol Values at or below this threshold are treated as zero.
#' @param point_size Size of plotted spots.
#' @param pt.size.factor Alias for `point_size` for compatibility.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to non-zero spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param crop Deprecated alias for `zoom`.
#' @param combine Combine plots with patchwork when possible.
#' @param ... Ignored; accepted for compatibility with older plotting
#'   calls.
#' @return A ggplot/patchwork object, list, or `NULL` when all spots are zero.
#' @export
spatialExpressionPlot <- function(object,
                                      features,
                                      assay = NULL,
                                      layer = "counts",
                                      slot = NULL,
                                      images = NULL,
                                      image.scale = "hires",
                                      feature_colors = c("#D3D3D3", "#FFFF00"),
                                      zero_tol = 0,
                                      point_size = 0.6,
                                      pt.size.factor = NULL,
                                      image.alpha = 1,
                                      zoom = FALSE,
                                      crop = NULL,
                                      combine = TRUE,
                                      ...) {
  .require_namespace("SeuratObject")
  .require_namespace("patchwork")
  point_size <- pt.size.factor %||% point_size
  if (!is.null(assay)) {
    SeuratObject::DefaultAssay(object) <- assay
  }
  if (length(feature_colors) < 2L) {
    rlang::abort("`feature_colors` must contain at least two colors.")
  }

  expr <- .sirap_fetch_data(object, vars = features, layer = layer, slot = slot)
  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
    dat <- merge(spatial$coords, expr, by.x = "cell", by.y = 0, all.x = TRUE)
    for (feature in features) {
      keep <- !is.na(dat[[feature]]) & dat[[feature]] > zero_tol
      if (!any(keep)) {
        next
      }
      plot_data <- dat[keep, , drop = FALSE]
      plots[[paste(image, feature, sep = "_")]] <- .sirap_spatial_base(
        plot_data,
        spatial$image,
        image.alpha = image.alpha,
        zoom = zoom,
        crop = crop,
        view_data = spatial$coords
      ) +
        ggplot2::geom_point(
          ggplot2::aes(color = .data[[feature]]),
          size = point_size,
          alpha = 1,
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
#' @param cols Optional named color vector.
#' @param point_size Size of plotted spots.
#' @param pt.size.factor Alias for `point_size` for compatibility.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param crop Deprecated alias for `zoom`.
#' @param combine Combine plots with patchwork when possible.
#' @param ... Ignored; accepted for compatibility with older plotting
#'   calls.
#' @return A ggplot/patchwork object, list, or `NULL` when all spots are blank.
#' @export
spatialDiscretePlot <- function(object,
                                  meta_col,
                                  zero_value = 0,
                                  images = NULL,
                                  image.scale = "hires",
                                  cols = NULL,
                                  point_size = 1.6,
                                  pt.size.factor = NULL,
                                  image.alpha = 1,
                                  zoom = FALSE,
                                  crop = NULL,
                                  combine = TRUE,
                                  ...) {
  .check_required_cols(object@meta.data, meta_col, "object@meta.data")
  point_size <- pt.size.factor %||% point_size

  values <- object@meta.data[[meta_col]]
  names(values) <- rownames(object@meta.data)
  keep_cells <- names(values)[!is.na(values) & as.character(values) != as.character(zero_value)]
  if (length(keep_cells) == 0L) {
    message("All spots are zero/blank for `", meta_col, "`.")
    return(NULL)
  }

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
    plot_data <- spatial$coords[spatial$coords$cell %in% keep_cells, , drop = FALSE]
    if (nrow(plot_data) == 0L) {
      next
    }
    plot_data[[meta_col]] <- values[plot_data$cell]
    plot_data[[meta_col]] <- factor(plot_data[[meta_col]], levels = unique(plot_data[[meta_col]]))
    color_values <- cols %||% scales::hue_pal()(length(levels(plot_data[[meta_col]])))
    if (is.null(names(color_values))) {
      names(color_values) <- levels(plot_data[[meta_col]])
    }

    plots[[image]] <- .sirap_spatial_base(
      plot_data,
      spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      crop = crop,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[meta_col]]),
        size = point_size,
        alpha = 1,
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
#' @param slot Deprecated Seurat slot alias for `layer`.
#' @param only_blended Return only the blended spatial panel and legend.
#' @param bottom_left,bottom_right,top_left,top_right Corner colors for the
#'   two-feature blend grid.
#' @param zero_tol Values at or below this threshold are treated as zero.
#' @param images Optional Seurat image names.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param point_size Size of plotted spots.
#' @param pt.size.factor Alias for `point_size` for compatibility.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param crop Deprecated alias for `zoom`.
#' @param ... Ignored; accepted for compatibility with older plotting
#'   calls.
#' @return A ggplot/patchwork object or a list when `combine = FALSE`.
#' @export
spatialExpressionBlendPlot <- function(object,
                                                features,
                                                combine = TRUE,
                                                feature_1_alt_name = NULL,
                                                feature_2_alt_name = NULL,
                                                assay = NULL,
                                                layer = "counts",
                                                slot = NULL,
                                                only_blended = FALSE,
                                                bottom_left = "#000000",
                                                bottom_right = "#FF0000",
                                                top_left = "#00FF00",
                                                top_right = "#FFFF00",
                                                zero_tol = 0,
                                                images = NULL,
                                                image.scale = "hires",
                                                point_size = 0.6,
                                                pt.size.factor = NULL,
                                                image.alpha = 1,
                                                zoom = FALSE,
                                                crop = NULL,
                                                ...) {
  .require_namespace("SeuratObject")
  .require_namespace("patchwork")
  point_size <- pt.size.factor %||% point_size
  if (length(features) != 2L) {
    rlang::abort(paste0("Requires two features, received ", length(features), "."))
  }
  if (!is.null(assay)) {
    SeuratObject::DefaultAssay(object) <- assay
  }

  expr <- .sirap_fetch_data(object, vars = features, layer = layer, slot = slot)
  color_grid <- .sirap_blend_color_grid(
    bottom_left = bottom_left,
    bottom_right = bottom_right,
    top_left = top_left,
    top_right = top_right
  )
  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
    dat <- merge(spatial$coords, expr, by.x = "cell", by.y = 0, all.x = TRUE)
    keep <- rowSums(dat[, features, drop = FALSE] > zero_tol, na.rm = TRUE) > 0
    if (!any(keep)) {
      next
    }
    dat <- dat[keep, , drop = FALSE]
    dat$.sirap_blend_color <- .sirap_blend_colors(dat[, features, drop = FALSE], color_grid)
    color_values <- stats::setNames(unique(dat$.sirap_blend_color), unique(dat$.sirap_blend_color))

    blended <- .sirap_spatial_base(
      dat,
      spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      crop = crop,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$.sirap_blend_color),
        size = point_size,
        alpha = 1,
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
        feature_colors = c(bottom_left, bottom_right),
        zero_tol = zero_tol,
        point_size = point_size,
        image.alpha = image.alpha,
        zoom = zoom,
        crop = crop,
        combine = TRUE
      )
      f2 <- spatialExpressionPlot(
        object,
        features[[2]],
        images = image,
        image.scale = image.scale,
        feature_colors = c(bottom_left, top_left),
        zero_tol = zero_tol,
        point_size = point_size,
        image.alpha = image.alpha,
        zoom = zoom,
        crop = crop,
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
#' @param slot Deprecated Seurat slot alias for `layer`.
#' @param colors One color per feature.
#' @param breadth_metric Metric used to choose default feature order.
#' @param images Optional Seurat image names.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param point_size Size of plotted spots.
#' @param pt.size.factor Alias for `point_size` for compatibility.
#' @param zero_tol Values at or below this threshold are treated as zero.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param crop Deprecated alias for `zoom`.
#' @param ... Ignored; accepted for compatibility with older plotting
#'   calls.
#' @return A ggplot/patchwork object or a list when `combine = FALSE`.
#' @export
spatialExpressionMultiPlot <- function(object,
                                    features,
                                    feat_order = NULL,
                                    combine = TRUE,
                                    assay = NULL,
                                    layer = "counts",
                                    slot = NULL,
                                    colors = NULL,
                                    breadth_metric = c("n_expressing"),
                                    images = NULL,
                                    image.scale = "hires",
                                    point_size = 0.6,
                                    pt.size.factor = NULL,
                                    zero_tol = 0,
                                    image.alpha = 1,
                                    zoom = FALSE,
                                    crop = NULL,
                                    ...) {
  .require_namespace("SeuratObject")
  point_size <- pt.size.factor %||% point_size
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

  expr <- .sirap_fetch_data(object, vars = features, layer = layer, slot = slot)
  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
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
    ordered_colors <- colors[feature_order]

    spatial_plot <- .sirap_spatial_base(
      dat,
      spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      crop = crop,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$.sirap_feature),
        size = point_size,
        alpha = 1,
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
#' @param object Seurat object.
#' @param receptors AIRR-style receptor table containing a spatial barcode.
#' @param barcode_col Column in `receptors` matching Seurat spot names. If
#'   `NULL`, the first available barcode-like column is used.
#' @param feature Optional gene/feature to plot underneath receptor-positive
#'   spots.
#' @param gene Alias for `feature`, retained for readability in receptor plots.
#' @param overlay_col Ignored; retained for compatibility with the old
#'   plotting implementation.
#' @param receptor_color Color for receptor-positive spots.
#' @param receptor_size Size for receptor-positive spots.
#' @param receptor_stroke Stroke width for hollow receptor circles.
#' @param receptor_alpha Alpha for receptor circles.
#' @param feature_colors Colors used, in order, to build the continuous gene
#'   expression ramp when `feature` is set.
#' @param feature_size Size for gene expression points.
#' @param feature_alpha Alpha for gene expression points.
#' @param zero_tol Values at or below this threshold are treated as zero.
#' @param gene_color_steps Ignored; retained for compatibility.
#' @param title Optional title.
#' @param assay,layer,image.scale Additional plotting controls.
#' @param images Optional Seurat image names.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param crop Deprecated alias for `zoom`.
#' @param receptor_label Legend label for receptor-positive spots.
#' @param show_receptor_legend Show the receptor legend. Defaults to `TRUE`
#'   when a feature is shown and `FALSE` for receptor-only plots.
#' @param ... Ignored; accepted for compatibility with older plotting
#'   calls.
#' @return A ggplot/patchwork object.
#' @export
plotReceptorExpression <- function(object,
                                      receptors,
                                      barcode_col = NULL,
                                      feature = NULL,
                                      gene = NULL,
                                      overlay_col = ".sirap_receptor_overlay",
                                      receptor_color = "black",
                                      receptor_size = 1,
                                      receptor_stroke = 0.7,
                                      receptor_alpha = 1,
                                      feature_colors = c("#939393", "#FFFF00"),
                                      feature_size = 1,
                                      feature_alpha = 1,
                                      zero_tol = 0,
                                      gene_color_steps = 100L,
                                      title = NULL,
                                      assay = NULL,
                                      layer = "counts",
                                      image.scale = "hires",
                                      images = NULL,
                                      image.alpha = 1,
                                      zoom = FALSE,
                                      crop = NULL,
                                      receptor_label = "Receptor",
                                      show_receptor_legend = NULL,
                                      ...) {
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
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
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
      crop = crop,
      view_data = spatial$coords
    )

    if (!is.null(feature)) {
      gene_data <- dat[dat$.sirap_expressed, , drop = FALSE]
      if (nrow(gene_data) > 0L) {
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
            alpha = 1,
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

.guess_barcode_col <- function(data, barcode_col = NULL) {
  if (!is.null(barcode_col)) {
    .check_required_cols(data, barcode_col, "receptors")
    return(barcode_col)
  }

  candidates <- c(
    "barcode",
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
#' @param overlay_col Ignored; retained for compatibility.
#' @param title Optional title.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param pt.size.factor Size of plotted tips.
#' @param images Optional Seurat image names.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted spots. Defaults to `FALSE`, which shows the
#'   full tissue area covered by the spatial coordinates.
#' @param crop Deprecated alias for `zoom`.
#' @param show.legend Show color legend.
#' @param ... Ignored; accepted for compatibility with older plotting
#'   calls.
#' @return A ggplot/patchwork object.
#' @export
plotSpatialColorTips <- function(object,
                                tip_data,
                                barcode_col = "barcode_8um",
                                color_col = "tip_color",
                                overlay_col = ".sirap_tip_color",
                                title = NULL,
                                image.scale = "hires",
                                pt.size.factor = 1.6,
                                images = NULL,
                                image.alpha = 1,
                                zoom = FALSE,
                                crop = NULL,
                                show.legend = FALSE,
                                ...) {
  .require_namespace("SeuratObject")
  tip_data <- .as_data_frame(tip_data)
  .check_required_cols(tip_data, c(barcode_col, color_col), "tip_data")
  tip_data <- tip_data[!is.na(tip_data[[barcode_col]]) & !is.na(tip_data[[color_col]]), , drop = FALSE]
  spot_color <- stats::setNames(as.character(tip_data[[color_col]]), as.character(tip_data[[barcode_col]]))
  spot_color <- spot_color[!duplicated(names(spot_color))]
  colors <- unique(stats::na.omit(as.character(spot_color)))
  names(colors) <- colors

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
    dat <- spatial$coords[spatial$coords$cell %in% names(spot_color), , drop = FALSE]
    if (nrow(dat) == 0L) {
      next
    }
    dat$.sirap_tip_color <- spot_color[dat$cell]
    p <- .sirap_spatial_base(
      dat,
      spatial$image,
      image.alpha = image.alpha,
      zoom = zoom,
      crop = crop,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$.sirap_tip_color),
        size = pt.size.factor,
        alpha = 1,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = colors, name = "Tip color") +
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
#' @param images Optional Seurat image names.
#' @param zoom Zoom spatial plots to the clone. Defaults to `FALSE`, which
#'   shows the full tissue area covered by spatial coordinates.
#' @param tree_layout Layout passed to `ggtree::ggtree()`.
#' @param tree_tip_size Size of tree tip points.
#' @param treescale Mutational distance shown in the tree scale bar.
#' @param spatial_tip_size Size of spatial tip points.
#' @param show_spatial_legend Show the spatial color legend.
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
                                      images = NULL,
                                      zoom = FALSE,
                                      tree_layout = "rectangular",
                                      tree_tip_size = 2.5,
                                      treescale = 0.01,
                                      spatial_tip_size = 1.6,
                                      show_spatial_legend = FALSE,
                                      ...) {
  .require_namespace("ggtree")
  records <- .sirap_dowser_tree_records(trees)
  if (!is.null(clone_ids)) {
    records <- records[names(records) %in% as.character(clone_ids)]
  }
  if (length(records) == 0L) {
    rlang::abort("No dowser trees matched the requested clone IDs.")
  }

  out <- list()
  for (clone_id in names(records)) {
    record <- records[[clone_id]]
    tree <- record$tree
    clone_data <- .sirap_airr_clone_data(record$data, clone_id = clone_id, clone_col = clone_col)
    .check_required_cols(clone_data, c(id_col, barcode_col, clone_col), "dowser clone data")

    tip_ids <- setdiff(tree$tip.label, "Germline")
    tip_data <- clone_data[clone_data[[id_col]] %in% tip_ids, , drop = FALSE]
    if (nrow(tip_data) == 0L) {
      next
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

    tree_plot <- .sirap_plot_colored_tree(
      tree,
      colored_tips,
      id_col = id_col,
      color_col = "tip_color",
      title = paste("Clone", clone_id),
      layout = tree_layout,
      tip_size = tree_tip_size,
      treescale = treescale
    )

    spatial_plot <- plotSpatialColorTips(
      object,
      tip_data = colored_tips,
      barcode_col = barcode_col,
      color_col = "tip_color",
      title = NULL,
      image.scale = image.scale,
      images = images,
      pt.size.factor = spatial_tip_size,
      zoom = zoom,
      show.legend = show_spatial_legend,
      ...
    )

    out[[as.character(clone_id)]] <- list(
      tree = tree_plot,
      spatial_clone = spatial_plot
    )
  }

  out
}

#' Save paired clone tree and spatial plots to a PDF
#'
#' @param pairs Output from `plotSpatialCloneTreePairs()`.
#' @param path PDF path.
#' @param width,height PDF page size in inches.
#' @param rel_widths Relative widths for tree and spatial panels.
#' @return Invisibly returns `path`.
#' @export
spatialCloneTreePairsToPDF <- function(pairs,
                                       path,
                                       width = 12,
                                       height = 6,
                                       rel_widths = c(1.2, 1)) {
  if (length(pairs) == 0L) {
    rlang::abort("`pairs` does not contain any clone tree/spatial plots.")
  }
  .ensure_dir(path)

  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (pair in pairs) {
    print(
      cowplot::plot_grid(
        pair$tree,
        pair$spatial_clone,
        ncol = 2,
        rel_widths = rel_widths
      )
    )
  }

  invisible(path)
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
                                     treescale) {
  p <- ggtree::ggtree(tree, layout = layout)
  join_data <- tip_data[, c(id_col, color_col), drop = FALSE]
  join_data <- join_data[!duplicated(join_data[[id_col]]), , drop = FALSE]
  idx <- match(p$data$label, join_data[[id_col]])
  p$data[[color_col]] <- join_data[[color_col]][idx]
  p +
    ggtree::geom_tippoint(
      ggplot2::aes(color = .data[[color_col]]),
      size = tip_size,
      na.rm = TRUE
    ) +
    ggtree::geom_treescale(width = treescale) +
    ggplot2::scale_color_identity() +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      legend.position = "none",
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

.sirap_fetch_data <- function(object, vars, layer = "counts", slot = NULL) {
  layer <- slot %||% layer
  out <- tryCatch(
    SeuratObject::FetchData(object, vars = vars, layer = layer),
    error = function(e) {
      SeuratObject::FetchData(object, vars = vars, slot = layer)
    }
  )
  as.data.frame(out, stringsAsFactors = FALSE)
}

.sirap_spatial_data <- function(object, image, image.scale = "hires") {
  .require_namespace("SeuratObject")
  image_obj <- object[[image]]
  coords <- SeuratObject::GetTissueCoordinates(image_obj, scale = image.scale)
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  if ("cell" %in% colnames(coords)) {
    coords$cell <- as.character(coords$cell)
  } else {
    coords$cell <- rownames(coords)
  }

  x_col <- intersect(c("imagecol", "x", "pxl_col_in_fullres", "col"), colnames(coords))
  y_col <- intersect(c("imagerow", "y", "pxl_row_in_fullres", "row"), colnames(coords))
  if (length(x_col) == 0L || length(y_col) == 0L) {
    rlang::abort("Could not identify x/y coordinate columns in Seurat tissue coordinates.")
  }
  coords <- data.frame(
    cell = coords$cell,
    x = coords[[x_col[[1]]]],
    y = coords[[y_col[[1]]]],
    stringsAsFactors = FALSE
  )
  coords <- coords[coords$cell %in% colnames(object), , drop = FALSE]
  coords$.sirap_y <- -coords$y

  list(coords = coords, image = image_obj)
}

.sirap_spatial_base <- function(plot_data,
                                image,
                                image.alpha = 1,
                                zoom = FALSE,
                                crop = NULL,
                                view_data = NULL) {
  raster <- .sirap_image_raster(image, alpha = image.alpha)
  image_dim <- dim(raster)
  height <- image_dim[[1]]
  width <- image_dim[[2]]

  xlim <- c(0, width)
  ylim <- c(-height, 0)
  if (!is.null(crop)) {
    zoom <- crop
  }
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
  img <- SeuratObject::GetImage(image, mode = "raw")
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
