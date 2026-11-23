#' Plot length histograms of V, J, and junction regions
#'
#' `plotLengthQc` plots one histogram each for the V length, the J length, and
#' the junction length of the receptors. Used to choose length thresholds for
#' [filterAirr()].
#'
#' @param airr_df AIRR-style receptor table, or a SIRUP object. The table needs
#'   the columns `v_sequence_start`, `v_sequence_end`, `j_sequence_start`,
#'   `j_sequence_end`, and `junction_length`.
#' @param bin_width Width of the histogram bins (default 5).
#' @param title Plot title.
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#'
#' @return A ggplot object with one panel for each length.
#' @seealso [filterAirr()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' plotLengthQc(x)
#' @export
plotLengthQc <- function(airr_df, bin_width = 5, title = NULL, slot = NULL) {
  airr_df <- .sirup_receptor_input(airr_df, slot = slot)$receptors
  .check_required_cols(airr_df, c("v_sequence_end", "v_sequence_start", "j_sequence_end", "j_sequence_start", "junction_length"), "airr_df")
  out <- .as_data_frame(airr_df)
  out$v_sequence_length <- out$v_sequence_end - out$v_sequence_start + 1L
  out$j_sequence_length <- out$j_sequence_end - out$j_sequence_start + 1L
  long <- tidyr::pivot_longer(
    out,
    cols = dplyr::all_of(c("v_sequence_length", "j_sequence_length", "junction_length")),
    names_to = "metric",
    values_to = "length"
  ) %>%
    dplyr::mutate(metric = factor(.data$metric, levels = c("v_sequence_length", "j_sequence_length", "junction_length"), labels = c("V Length", "J Length", "Junction Length")))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$length)) +
    ggplot2::geom_histogram(binwidth = bin_width, na.rm = TRUE) +
    ggplot2::facet_wrap(ggplot2::vars(.data$metric), scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Length", y = "Count")
}

#' Plot gene expression on a tissue image
#'
#' `spatialExpressionPlot` plots the expression of one or more features. By default, it
#' draws only the spots with expression above zero.
#' 
#' If the object has multiple images or samples, the function draws one plot
#' for each image or sample. If given several features, the function draws one
#' panel for each feature. 
#'
#' @param object Seurat object or SIRUP object.
#' @param features Names of the features to plot.
#' @param assay Assay to use.
#' @param layer Seurat expression layer to use (default `"counts"`).
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param feature_colors Colors for the continuous expression scale, from the
#'   low value to the high value. If you give them, they replace
#'   `feature_color_palette`.
#' @param feature_color_palette Color palette for the expression scale, if
#'   `feature_colors` is `NULL`: `"viridis"` (default), `"inferno"`, `"plasma"`,
#'   `"cividis"`, `"rocket"`, or `"mako"`.
#' @param zero_tol Values at or below this threshold count as zero (default 0).
#'    Set to -Inf to plot all spots.
#' @param point_size Size of the spots. The default is 0.6.
#' @param point_alpha Opacity of the spots.
#' @param order_points Draw the low values first, so that the high values show
#'   on top (default `TRUE`)?
#' @param title Plot title. For one panel, it replaces the feature name. For
#'   several panels that the function combines, it is the title of the whole
#'   patchwork, and each panel keeps its feature name.
#'   Give one title for each panel, in the order of the panels, to title every
#'   panel. This works for a multi-sample object with one feature, where the
#'   panels are the samples.
#' @param image_alpha Opacity of the tissue image.
#' @param zoom Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot
#'   shows the full tissue area that the spatial coordinates cover.
#' @param combine Combine the panels into one figure with patchwork, if there
#'   are several images or features (default `TRUE`)?
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#'
#' @return A ggplot object, a patchwork object, or a list of plots. The function
#'   returns `NULL` if all spots have zero expression.
#' @seealso [spatialExpressionBlendPlot()], [spatialExpressionMultiPlot()],
#'   [plotReceptorsExpression()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' spatialExpressionPlot(x, "gene1")
#'
#' # Plot two features. The function combines the panels with patchwork.
#' spatialExpressionPlot(x, c("gene1", "gene2"))
#' @export
spatialExpressionPlot <- function(object,
                                      features,
                                      assay = NULL,
                                      layer = "counts",
                                      images = NULL,
                                      sample_id = NULL,
                                      image_scale = NULL,
                                      swap_coords = NULL,
                                      feature_colors = NULL,
                                      feature_color_palette = "viridis",
                                      zero_tol = 0,
                                      point_size = 0.6,
                                      point_alpha = 1,
                                      order_points = TRUE,
                                      title = NULL,
                                      image_alpha = 1,
                                      zoom = FALSE,
                                      combine = TRUE,
                                      combine_sample_plots = FALSE) {
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  panel_titles <- .sirup_panel_titles(title)
  if (!is.null(panel_titles)) title <- NULL
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (!is.null(assay)) {
    object <- .sirup_set_default_assay(object, assay)
  }
  feature_colors <- .sirup_feature_color_values(
    feature_colors = feature_colors,
    feature_color_palette = feature_color_palette
  )

  expr <- .sirup_spatial_fetch_data(object, vars = features, layer = layer, assay = assay)
  plots <- list()
  spatial_images <- .sirup_spatial_images(object, images = images, sample_id = sample_id)
  single_panel <- length(features) == 1L && length(spatial_images) == 1L
  for (image in spatial_images) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    dat <- .sirup_merge_spatial_values(spatial$coords, expr)
    for (feature in features) {
      keep <- !is.na(dat[[feature]]) & dat[[feature]] > zero_tol
      if (!any(keep)) {
        next
      }
      plot_data <- dat[keep, , drop = FALSE]
      if (isTRUE(order_points)) {
        plot_data <- plot_data[order(plot_data[[feature]], na.last = TRUE), , drop = FALSE]
      }
      plots[[paste(image, feature, sep = "_")]] <- .sirup_spatial_base(
        plot_data,
        spatial$image,
        image_alpha = image_alpha,
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
        ggplot2::labs(title = if (single_panel) title %||% feature else feature)
    }
  }

  out <- .sirup_return_spatial_plots(plots, combine = combine, combine_sample_plots = combine_sample_plots, panel_titles = panel_titles)
  if (!is.null(out) && !is.list(out) && !single_panel && isTRUE(combine) && !is.null(title)) {
    out <- out + patchwork::plot_annotation(title = title)
  }
  out
}

.sirup_feature_color_values <- function(feature_colors = NULL,
                                        feature_color_palette = "viridis",
                                        n = 256L) {
  if (!is.null(feature_colors)) {
    if (length(feature_colors) < 2L) {
      rlang::abort("`feature_colors` must contain at least two colors.")
    }
    return(feature_colors)
  }
  palettes <- c(
    viridis = "Viridis",
    inferno = "Inferno",
    plasma = "Plasma",
    cividis = "Cividis",
    rocket = "Rocket",
    mako = "Mako"
  )
  feature_color_palette <- rlang::arg_match(
    feature_color_palette,
    values = names(palettes),
    error_arg = "feature_color_palette"
  )
  tryCatch(
    grDevices::hcl.colors(n, palette = palettes[[feature_color_palette]]),
    error = function(e) grDevices::hcl.colors(n, palette = "Viridis")
  )
}

#' Plot a metadata column on a tissue image
#'
#' `spatialDiscretePlot` colors each spot by the value of a metadata column. It
#' leaves out the spots that have the "zero_value" value.
#' 
#' #' If the object has multiple images or samples, the function draws one plot
#' for each image or sample.
#'
#' @param object Seurat object or SIRUP object.
#' @param meta_col Name of the metadata column to plot.
#' @param zero_value Value that means "not present" (default 0). The function
#'   does not draw these spots.
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param cols Named vector of colors, one for each value.
#' @param point_size Size of the spots. The default is 1.6.
#' @param point_alpha Opacity of the spots.
#' @param order_points Draw the low numeric values, or the early factor levels,
#'   first, so that the high or late values show on top (default `TRUE`)?
#' @param image_alpha Opacity of the tissue image.
#' @param zoom Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot
#'   shows the full tissue area that the spatial coordinates cover.
#' @param combine Combine the panels with patchwork, if possible (default
#'   `TRUE`)?
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#'
#' @return A ggplot object, a patchwork object, or a list of plots. The function
#'   returns `NULL` if all spots have the blank value.
#' @seealso [spatialExpressionPlot()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' spatialDiscretePlot(x, "cell_type")
#' @export
spatialDiscretePlot <- function(object,
                                  meta_col,
                                  zero_value = 0,
                                  images = NULL,
                                  sample_id = NULL,
                                  image_scale = NULL,
                                  swap_coords = NULL,
                                  cols = NULL,
                                  point_size = 1.6,
                                  point_alpha = 1,
                                  order_points = TRUE,
                                  image_alpha = 1,
                                  zoom = FALSE,
                                  combine = TRUE,
                                  combine_sample_plots = FALSE) {
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  values <- .sirup_spatial_metadata_values(object, meta_col)
  values <- values[!is.na(values[[meta_col]]) & as.character(values[[meta_col]]) != as.character(zero_value), , drop = FALSE]
  if (nrow(values) == 0L) {
    message("All spots are zero/blank for `", meta_col, "`.")
    return(NULL)
  }

  plots <- list()
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    plot_data <- .sirup_merge_spatial_values(spatial$coords, values)
    plot_data <- plot_data[!is.na(plot_data[[meta_col]]) & as.character(plot_data[[meta_col]]) != as.character(zero_value), , drop = FALSE]
    if (nrow(plot_data) == 0L) {
      next
    }
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

    plots[[image]] <- .sirup_spatial_base(
      plot_data,
      spatial$image,
      image_alpha = image_alpha,
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

  .sirup_return_spatial_plots(plots, combine = combine, combine_sample_plots = combine_sample_plots)
}

#' Blend the expression of two features on a tissue image
#'
#' `spatialExpressionBlendPlot` colors each spot by the expression of two
#' features at once. By default, it leaves out the spots where both features are zero.
#' 
#' #' If the object has multiple images or samples, the function draws one plot
#' for each image or sample. 
#'
#' @param object Seurat object or SIRUP object.
#' @param features Names of two features.
#' @param combine Combine the panels and the legend with patchwork (default
#'   `TRUE`)?
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#' @param feature_1_alt_name Label of the first feature in the blend legend. If
#'   `NULL`, the function uses the feature name.
#' @param feature_2_alt_name Label of the second feature in the blend legend. If
#'   `NULL`, the function uses the feature name.
#' @param assay Assay to use.
#' @param layer Seurat expression layer to use (default `"counts"`).
#' @param only_blended Return only the blended panel and the legend (default
#'   `FALSE`)?
#' @param bottom_left Color of the bottom-left corner of the blend grid (low in
#'   both features). The default is `"#FFFFFF"`.
#' @param bottom_right Color of the bottom-right corner of the blend grid (high
#'   in the first feature). The default is `"#0000FF"`.
#' @param top_left Color of the top-left corner of the blend grid (high in the
#'   second feature). The default is `"#FFFF00"`.
#' @param top_right Color of the top-right corner of the blend grid (high in
#'   both features). The default is `"#00FF00"`.
#' @param zero_tol Values at or below this threshold count as zero (default 0).
#'   Set to -Inf to plot all spots.
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param point_size Size of the spots. The default is 0.6.
#' @param point_alpha Opacity of the spots.
#' @param order_points Draw the spots with low blended expression first, so that
#'   the spots with high expression show on top (default `TRUE`)?
#' @param image_alpha Opacity of the tissue image.
#' @param zoom Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot
#'   shows the full tissue area that the spatial coordinates cover.
#'
#' @return A ggplot object or a patchwork object. If `combine = FALSE`, a list.
#' @seealso [spatialExpressionPlot()], [spatialExpressionMultiPlot()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' spatialExpressionBlendPlot(x, c("gene1", "gene2"))
#' @export
spatialExpressionBlendPlot <- function(object,
                                                features,
                                                combine = TRUE,
                                                combine_sample_plots = FALSE,
                                                feature_1_alt_name = NULL,
                                                feature_2_alt_name = NULL,
                                                assay = NULL,
                                                layer = "counts",
                                                only_blended = FALSE,
                                                bottom_left = "#FFFFFF",
                                                bottom_right = "#0000FF",
                                                top_left = "#FFFF00",
                                                top_right = "#00FF00",
                                                zero_tol = 0,
                                                images = NULL,
                                                sample_id = NULL,
                                                image_scale = NULL,
                                                swap_coords = NULL,
                                                point_size = 0.6,
                                                point_alpha = 1,
                                                order_points = TRUE,
                                                image_alpha = 1,
                                                zoom = FALSE) {
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (length(features) != 2L) {
    rlang::abort(paste0("Requires two features, received ", length(features), "."))
  }
  if (!is.null(assay)) {
    object <- .sirup_set_default_assay(object, assay)
  }

  expr <- .sirup_spatial_fetch_data(object, vars = features, layer = layer, assay = assay)
  color_grid <- .sirup_blend_color_grid(
    bottom_left = bottom_left,
    bottom_right = bottom_right,
    top_left = top_left,
    top_right = top_right
  )
  plots <- list()
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    dat <- .sirup_merge_spatial_values(spatial$coords, expr)
    keep <- rowSums(dat[, features, drop = FALSE] > zero_tol, na.rm = TRUE) > 0
    if (!any(keep)) {
      next
    }
    dat <- dat[keep, , drop = FALSE]
    if (isTRUE(order_points)) {
      dat$.sirup_blend_order <- rowSums(dat[, features, drop = FALSE], na.rm = TRUE)
      dat <- dat[order(dat$.sirup_blend_order, na.last = TRUE), , drop = FALSE]
    }
    dat$.sirup_blend_color <- .sirup_blend_colors(dat[, features, drop = FALSE], color_grid)
    color_values <- stats::setNames(unique(dat$.sirup_blend_color), unique(dat$.sirup_blend_color))

    blended <- .sirup_spatial_base(
      dat,
      spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$.sirup_blend_color),
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::labs(title = paste(features, collapse = " + "))

    legend <- .sirup_blend_legend(
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
        assay = assay,
        layer = layer,
        images = image,
        image_scale = image_scale,
        swap_coords = swap_coords,
        feature_colors = c(bottom_left, bottom_right),
        zero_tol = zero_tol,
        point_size = point_size,
        point_alpha = point_alpha,
        order_points = order_points,
        image_alpha = image_alpha,
        zoom = zoom,
        combine = TRUE
      )
      f2 <- spatialExpressionPlot(
        object,
        features[[2]],
        assay = assay,
        layer = layer,
        images = image,
        image_scale = image_scale,
        swap_coords = swap_coords,
        feature_colors = c(bottom_left, top_left),
        zero_tol = zero_tol,
        point_size = point_size,
        point_alpha = point_alpha,
        order_points = order_points,
        image_alpha = image_alpha,
        zoom = zoom,
        combine = TRUE
      )
      plots[[image]] <- patchwork::wrap_plots(list(f1, f2, blended, legend), nrow = 1, widths = c(0.28, 0.28, 0.28, 0.16))
    }
  }

  .sirup_return_spatial_plots(plots, combine = combine, combine_sample_plots = combine_sample_plots)
}

#' Plot several features as ordered color layers
#'
#' `spatialExpressionMultiPlot` draws each feature in its own color on one
#' tissue image. The features are layers. A feature that the function draws
#' later shows on top.
#' 
#' If the object has multiple images or samples, the function draws one plot
#' for each image or sample.
#'
#' @param object Seurat object or SIRUP object.
#' @param features Names of the features to plot.
#' @param feat_order Order of the layers. Put broad features (low priority)
#'   first and rare features (high priority) last. If `NULL`, the function
#'   orders the features with `breadth_metric`.
#' @param combine Combine the spatial panels and the legends with patchwork
#'   (default `TRUE`)?
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#' @param assay Assay to use.
#' @param layer Seurat expression layer to use (default `"counts"`).
#' @param colors Colors, one for each feature.
#' @param breadth_metric Measure that sets the default order of the features
#'   (default `"n_expressing"`, the number of spots that express the feature).
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param point_size Size of the spots. The default is 0.6.
#' @param point_alpha Opacity of the spots.
#' @param order_points Draw the low-priority features first, so that the
#'   high-priority features show on top (default `TRUE`)?
#' @param zero_tol Values at or below this threshold count as zero (default 0).
#'  Set to -Inf to plot all spots.
#' @param image_alpha Opacity of the tissue image.
#' @param zoom Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot
#'   shows the full tissue area that the spatial coordinates cover.
#'
#' @return A ggplot object or a patchwork object. If `combine = FALSE`, a list.
#' @seealso [spatialExpressionPlot()], [spatialExpressionBlendPlot()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' spatialExpressionMultiPlot(x, c("gene1", "gene2", "gene3"))
#' @export
spatialExpressionMultiPlot <- function(object,
                                    features,
                                    feat_order = NULL,
                                    combine = TRUE,
                                    combine_sample_plots = FALSE,
                                    assay = NULL,
                                    layer = "counts",
                                    colors = NULL,
                                    breadth_metric = c("n_expressing"),
                                    images = NULL,
                                    sample_id = NULL,
                                    image_scale = NULL,
                                    swap_coords = NULL,
                                    point_size = 0.6,
                                    point_alpha = 1,
                                    order_points = TRUE,
                                    zero_tol = 0,
                                    image_alpha = 1,
                                    zoom = FALSE) {
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  breadth_metric <- match.arg(breadth_metric)
  if (length(features) < 2L) {
    rlang::abort("Requires at least 2 features.")
  }
  if (!is.null(assay)) {
    object <- .sirup_set_default_assay(object, assay)
  }
  colors <- colors %||% c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999")
  if (length(colors) < length(features)) {
    rlang::abort("`colors` must be at least as long as `features`.")
  }
  names(colors) <- features

  expr <- .sirup_spatial_fetch_data(object, vars = features, layer = layer, assay = assay)
  plots <- list()
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    dat <- .sirup_merge_spatial_values(spatial$coords, expr)
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

    dat$.sirup_feature <- NA_character_
    for (feature in feature_order) {
      dat$.sirup_feature[dat[[feature]] > zero_tol] <- feature
    }
    dat$.sirup_feature <- factor(dat$.sirup_feature, levels = feature_order)
    if (isTRUE(order_points)) {
      dat <- dat[order(as.integer(dat$.sirup_feature), na.last = TRUE), , drop = FALSE]
    }
    ordered_colors <- colors[feature_order]

    spatial_plot <- .sirup_spatial_base(
      dat,
      spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = spatial$coords
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$.sirup_feature),
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = ordered_colors, name = "Feature") +
      ggplot2::labs(title = paste(features, collapse = " / "))

    plots[[image]] <- spatial_plot
  }

  .sirup_return_spatial_plots(plots, combine = combine, combine_sample_plots = combine_sample_plots)
}

#' Plot receptor groups over gene expression
#'
#' `plotReceptorsExpression` plots one feature together with one or more groups
#' of receptors on a tissue image. For example, it shows the BCR-positive and
#' TCR-positive spots on one expression map.
#'
#' @details
#' The function draws the spots with expression above zero with a continuous
#' fill scale. It draws the spots with receptors on top as hollow circles.
#' 
#' If the object has multiple images or samples, the function draws one plot
#' for each image or sample. 
#'
#' @param object Seurat object, or a SIRUP object with spatial data and
#'   receptors.
#' @param receptors AIRR-style receptor table with a spatial ID column and a
#'   group column. Omit it if `object` is a SIRUP object.
#' @param group Name of the `receptors` column that defines the receptor groups,
#'   for example `"receptor_type"` with the values `"BCR"` and `"TCR"`. You can
#'   omit it if you give `value`. If you give both, `group` and `group_val` only
#'   select the receptors to draw.
#' @param group_val Group values to plot, in the order of the plot. If `NULL`,
#'   the function uses all non-missing values of `group`.
#' @param spatial_id Name of the `receptors` column that holds the Seurat spot
#'   names. If `NULL`, the function uses the stored field of a SIRUP object, or
#'   looks for a standard spatial ID column.
#' @param feature Feature to plot under the spots with receptors.
#' @param gene Alias for `feature`.
#' @param colors Colors of the receptor groups. Named colors match the group
#'   values. Unnamed colors follow the order of the values.
#' @param overlap Color of a spot that holds more than one receptor group. The
#'   function draws such a spot once, as a combined group with the member names
#'   joined by `" + "`. With `"average"` (default), the spot gets the mean color
#'   of its members. With one color, all combined spots get that color.
#' @param value Color the receptor spots on a continuous gradient and not by
#'   group. Give the name of a numeric `receptors` column, or `"count"` for the
#'   number of receptors in the spot.
#' @param value_colors Colors of the `value` gradient (default white to navy).
#'   Two colors give a low-to-high gradient. More colors give a gradient with
#'   several stops.
#' @param value_limits Limits of the `value` scale (default `c(0, NA)`).
#' @param value_combine How to combine `value` if a spot holds more than one
#'   receptor: `"mean"` (default), `"median"`, `"mode"`, `"max"`, `"min"`, or
#'   `"sum"`. The function ignores it if `value = "count"`.
#' @param receptor_size Size of the receptor circles (default 1).
#' @param receptor_stroke Line width of the receptor circles (default 0.7).
#' @param receptor_alpha Opacity of the receptor circles (default 1).
#' @param feature_colors Colors that make the expression scale (default
#'   `c("#939393", "#FFFF00")`).
#' @param feature_size Size of the expression points (default 1).
#' @param feature_alpha Opacity of the expression points (default 1).
#' @param zero_tol Values at or below this threshold count as zero (default 0).
#' @param title Plot title. Give one title for each plot, in the order of the
#'   plots, to title every plot of a multi-sample object.
#' @param order_points Draw the low expression values first, so that the high
#'   values show on top (default `TRUE`)? The function draws the receptor groups
#'   in the order of `group_val`.
#' @param assay Assay to use.
#' @param layer Seurat expression layer to use (default `"counts"`).
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Swap the two Seurat coordinate axes?
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_alpha Opacity of the tissue image.
#' @param zoom Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot
#'   shows the full tissue area that the spatial coordinates cover.
#' @param receptor_legend_title Title of the receptor legend. If `NULL`, the
#'   function uses `group`.
#' @param show_receptor_legend Show the receptor legend (default `TRUE`)?
#' @param combine Combine the image panels with patchwork (default `TRUE`)?
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#'
#' @return A ggplot object or a patchwork object. If `combine = FALSE`, a list.
#' @seealso [plotReceptors()], [spatialExpressionPlot()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' plotReceptorsExpression(x, group = "clone_id", feature = "gene1")
#' @export
plotReceptorsExpression <- function(object,
                                    receptors = NULL,
                                    group = NULL,
                                    group_val = NULL,
                                    spatial_id = NULL,
                                    feature = NULL,
                                    gene = NULL,
                                    colors = NULL,
                                    overlap = "average",
                                    value = NULL,
                                    value_colors = c("#FFFFFF", "#000080"),
                                    value_limits = c(0, NA),
                                    value_combine = c("mean", "median", "mode", "max", "min", "sum"),
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
                                    image_scale = NULL,
                                    swap_coords = NULL,
                                    images = NULL,
                                    sample_id = NULL,
                                    image_alpha = 1,
                                    zoom = FALSE,
                                    receptor_legend_title = NULL,
                                    show_receptor_legend = TRUE,
                                    combine = TRUE,
                                    combine_sample_plots = FALSE,
                                    slot = NULL) {
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  panel_titles <- .sirup_panel_titles(title)
  if (!is.null(panel_titles)) title <- NULL
  options <- .sirup_spatial_options(object, image_scale = image_scale, swap_coords = swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  inputs <- .sirup_resolve_inputs(object, receptors, spatial_id, slot)
  object <- inputs$object
  receptors <- inputs$receptors
  spatial_id <- inputs$spatial_id_col
  sample_col <- .sirup_receptor_sample_col(receptors %||% data.frame(), inputs$receptor_fields, object)
  if (is.null(receptors)) {
    rlang::abort("`receptors` is required when `object` is not a SIRUP object.")
  }

  feature <- feature %||% gene
  if (is.null(feature) || length(feature) != 1L || is.na(feature) || !nzchar(feature)) {
    rlang::abort("Provide one gene or feature using `feature` or `gene`.")
  }
  if (length(feature_colors) < 2L) {
    rlang::abort("`feature_colors` must contain at least two colors.")
  }

  receptors <- .as_data_frame(receptors)
  spatial_id <- .guess_spatial_id_col(receptors, spatial_id)
  if (is.null(group) && is.null(value)) {
    rlang::abort("`group` is required unless `value` is set.")
  }
  .check_required_cols(receptors, c(spatial_id, group), "receptors")

  # Continuous mode: one fill value per spot, so `group` only filters.
  use_value <- !is.null(value)
  if (use_value) {
    value_combine <- match.arg(value_combine)
    receptors <- receptors[!is.na(receptors[[spatial_id]]), , drop = FALSE]
    if (!is.null(group)) {
      gv <- if (is.null(group_val)) unique(stats::na.omit(as.character(receptors[[group]]))) else as.character(group_val)
      receptors <- receptors[as.character(receptors[[group]]) %in% gv, , drop = FALSE]
    }
    if (nrow(receptors) == 0L) {
      rlang::abort("No receptor rows with a spatial ID to summarise.")
    }
    summarised <- .sirup_receptor_spot_values(receptors, spatial_id, value, value_combine, sample_col = sample_col)
    receptor_lookup <- summarised$data
    receptor_legend_title <- receptor_legend_title %||% summarised$legend
  } else {

  receptors$.sirup_receptor_group <- as.character(receptors[[group]])
  if (is.null(group_val)) {
    group_val <- unique(stats::na.omit(receptors$.sirup_receptor_group))
  }
  group_val <- as.character(group_val)
  group_val <- group_val[!is.na(group_val) & nzchar(group_val)]
  if (length(group_val) == 0L) {
    rlang::abort("`group_val` must contain at least one non-missing receptor group.")
  }

  receptors <- receptors[
    !is.na(receptors$.sirup_receptor_group) &
      receptors$.sirup_receptor_group %in% group_val &
      !is.na(receptors[[spatial_id]]),
    ,
    drop = FALSE
  ]
  if (nrow(receptors) == 0L) {
    rlang::abort("No receptor rows matched `group_val` in `group`.")
  }
  receptors[[spatial_id]] <- as.character(receptors[[spatial_id]])
  receptors$.sirup_receptor_group <- factor(
    receptors$.sirup_receptor_group,
    levels = group_val
  )
  single_colors <- .sirup_receptor_group_colors(group_val, colors)
  receptor_legend_title <- receptor_legend_title %||% group

  receptor_lookup <- receptors[, unique(c(sample_col, spatial_id, ".sirup_receptor_group")), drop = FALSE]
  receptor_lookup <- receptor_lookup[!duplicated(receptor_lookup), , drop = FALSE]

  # Spots with more than one receptor group become a combined, colour-averaged
  # group (for example "BCR + TCR"), so each spot is drawn as a single ring.
  combined <- .sirup_combine_receptor_groups(receptor_lookup, spatial_id, group_val, single_colors, overlap = overlap, sample_col = sample_col)
  receptor_lookup <- combined$lookup
  plot_groups <- combined$levels
  colors <- combined$colors
  }

  if (!is.null(assay)) {
    object <- .sirup_set_default_assay(object, assay)
  }
  expr <- .sirup_spatial_fetch_data(object, vars = feature, layer = layer, assay = assay)

  plots <- list()
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(
      object,
      image = image,
      image_scale = image_scale,
      swap_coords = swap_coords
    )
    expression_data <- .sirup_merge_spatial_values(spatial$coords, expr)
    expression_data <- expression_data[
      !is.na(expression_data[[feature]]) & expression_data[[feature]] > zero_tol,
      ,
      drop = FALSE
    ]
    if (isTRUE(order_points) && nrow(expression_data) > 0L) {
      expression_data <- expression_data[
        order(expression_data[[feature]], na.last = TRUE),
        ,
        drop = FALSE
      ]
    }

    receptor_data <- .sirup_merge_receptor_lookup(spatial$coords, receptor_lookup, spatial_id, sample_col)
    if (!use_value) {
      receptor_data$.sirup_receptor_group <- factor(
        receptor_data$.sirup_receptor_group,
        levels = plot_groups
      )
    }

    if (nrow(expression_data) == 0L && nrow(receptor_data) == 0L) {
      next
    }

    p <- .sirup_spatial_base(
      spatial$coords,
      spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = spatial$coords
    )

    if (nrow(expression_data) > 0L) {
      # Expression always uses the fill scale and receptors always use the
      # colour scale, in both grouped and continuous mode. One plot carries only
      # one scale per aesthetic, so this keeps the receptor marks identical
      # whichever mode is in use. Shape 21 with no border still reads as a spot.
      p <- p +
        ggplot2::geom_point(
          data = expression_data,
          ggplot2::aes(x = .data$x, y = .data$.sirup_y, fill = .data[[feature]]),
          shape = 21,
          # "transparent", not NA: colour = NA makes geom_point drop the rows,
          # so every expression spot silently disappears.
          color = "transparent",
          stroke = 0,
          size = feature_size,
          alpha = feature_alpha,
          inherit.aes = FALSE,
          na.rm = TRUE
        ) +
        ggplot2::scale_fill_gradientn(colors = feature_colors, name = feature)
    }

    if (use_value) {
      if (nrow(receptor_data) > 0L) {
        if (isTRUE(order_points)) {
          receptor_data <- receptor_data[order(receptor_data$.sirup_count), , drop = FALSE]
        }
        # Hollow rings, same as grouped mode, with the ring colour carrying the
        # value. The expression layer uses the fill scale, so colour is free.
        p <- p +
          ggplot2::geom_point(
            data = receptor_data,
            ggplot2::aes(x = .data$x, y = .data$.sirup_y, color = .data$.sirup_count),
            shape = 21,
            fill = NA,
            size = receptor_size,
            stroke = receptor_stroke,
            alpha = receptor_alpha,
            inherit.aes = FALSE,
            na.rm = TRUE,
            show.legend = isTRUE(show_receptor_legend)
          ) +
          .sirup_value_color_scale(value_colors, value_limits, receptor_legend_title,
                                   aesthetic = "colour")
      }
    } else {
    for (grp in plot_groups) {
      group_data <- receptor_data[
        receptor_data$.sirup_receptor_group == grp,
        ,
        drop = FALSE
      ]
      if (nrow(group_data) == 0L) {
        next
      }
      group_data$.sirup_receptor_label <- factor(grp, levels = plot_groups)
      p <- p + ggplot2::geom_point(
        data = group_data,
        ggplot2::aes(
          x = .data$x,
          y = .data$.sirup_y,
          shape = .data$.sirup_receptor_label
        ),
        color = unname(colors[[grp]]),
        fill = NA,
        size = receptor_size,
        stroke = receptor_stroke,
        alpha = receptor_alpha,
        inherit.aes = FALSE,
        na.rm = TRUE,
        show.legend = isTRUE(show_receptor_legend)
      )
    }

    if (nrow(receptor_data) > 0L && isTRUE(show_receptor_legend)) {
      p <- p +
        ggplot2::scale_shape_manual(
          values = stats::setNames(rep(21, length(plot_groups)), plot_groups),
          name = receptor_legend_title,
          drop = FALSE
        ) +
        ggplot2::guides(
          shape = ggplot2::guide_legend(
            override.aes = list(
              color = unname(colors),
              fill = NA,
              stroke = receptor_stroke,
              alpha = receptor_alpha
            )
          )
        )
    }
    }

    plots[[image]] <- p + ggplot2::labs(title = title %||% feature)
  }

  .sirup_return_spatial_plots(plots, combine = combine, combine_sample_plots = combine_sample_plots, panel_titles = panel_titles)
}

#' Plot receptor groups on a tissue image
#'
#' `plotReceptors` plots the spots that hold receptors. It colors the spots by
#' receptor group, or by a continuous value.
#' 
#' #' If the object has multiple images or samples, the function draws one plot
#' for each image or sample. 
#'
#' @param object Seurat object, or a SIRUP object with spatial data and
#'   receptors.
#' @param receptors AIRR-style receptor table with a spatial ID column and a
#'   group column. Omit it if `object` is a SIRUP object.
#' @param group Name of the `receptors` column that selects and colors the
#'   receptor groups, for example `"clone_id"`.
#' @param group_val Values of `group` to plot. If `NULL`, the function uses all
#'   non-missing values.
#' @param spatial_id Name of the `receptors` column that holds the Seurat spot
#'   names. If `NULL`, the function uses the stored field of a SIRUP object, or
#'   uses the first standard spatial ID column that it finds.
#' @param colors Colors of the receptor groups. Named colors match the group
#'   values. Unnamed colors follow the order of the values.
#' @param overlap Color of a spot that holds more than one receptor group. The
#'   function draws such a spot once, as a combined group with the member names
#'   joined by `" + "`. With `"average"` (default), the spot gets the mean color
#'   of its members. With one color, all combined spots get that color.
#' @param value Color each spot on a continuous gradient and not by group. Give
#'   the name of a numeric `receptors` column, or `"count"` for the number of
#'   receptors in the spot. You can omit `group` in this mode. If you give it
#'   (with `group_val`), the function uses only the matching receptors.
#' @param value_colors Colors of the `value` gradient (default white to navy).
#'   Two colors give a low-to-high gradient. More colors give a gradient with
#'   several stops, so you can give a palette such as `viridisLite::viridis(9)`.
#' @param value_limits Limits of the `value` scale (default `c(0, NA)`, which
#'   starts the scale at zero).
#' @param value_combine How to combine `value` if a spot holds more than one
#'   receptor: `"mean"` (default), `"median"`, `"mode"`, `"max"`, `"min"`, or
#'   `"sum"`. The function ignores it if `value = "count"`.
#' @param receptor_size Size of the receptor points (default 1).
#' @param receptor_alpha Opacity of the receptor points (default 1).
#' @param title Plot title. Give one title for each plot, in the order of the
#'   plots, to title every plot of a multi-sample object.
#' @param order_points Draw the early receptor groups first, so that the later
#'   groups show on top (default `TRUE`)?
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_alpha Opacity of the tissue image.
#' @param zoom Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot
#'   shows the full tissue area that the spatial coordinates cover.
#' @param receptor_legend_title Title of the receptor legend. If `NULL`, the
#'   function uses `group`.
#' @param show_receptor_legend Show the receptor legend (default `TRUE`)?
#' @param polygons Draw polygons under the points on each spatial panel. With
#'   `FALSE` (default), the function draws no polygons. With `TRUE`, it uses the
#'   polygons stored in the SIRUP object. You can also give a data frame of
#'   polygon vertices, or the path of a CSV file of polygon vertices. A vertex
#'   table needs the columns set by `polygon_x` and `polygon_y`. It needs the
#'   column set by `polygon_id` if it holds more than one polygon.
#' @param polygon_x Name of the x-coordinate column in the polygon table
#'   (default `"x"`).
#' @param polygon_y Name of the y-coordinate column in the polygon table
#'   (default `"y"`).
#' @param polygon_id Name of the column that separates the polygons in the
#'   polygon table (default `"polygon"`). If the column does not exist, the
#'   function uses a `polygon` or `polygon_id` column. If neither exists, all
#'   rows form one polygon.
#' @param polygon_val IDs of the polygons to draw. If `NULL`, the function draws
#'   all polygons.
#' @param polygon_alpha Opacity of the polygon fill (default 0.2).
#' @param polygon_fill Fill color of the polygons (default `"white"`).
#' @param polygon_color Color of the polygon outlines (default `"white"`).
#' @param polygon_linewidth Width of the polygon outlines (default 0.6).
#' @param polygon_label Write a label at the center of each polygon (default
#'   `FALSE`)?
#' @param polygon_label_by Name of the column that holds the labels. If `NULL`,
#'   the function uses the polygon ID.
#' @param polygon_label_size Text size of the polygon labels (default 3).
#' @param polygon_label_color Text color of the polygon labels (default
#'   `"white"`).
#' @param combine Combine the image panels with patchwork (default `TRUE`)?
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#'
#' @return A ggplot object or a patchwork object. If `combine = FALSE`, a list.
#' @seealso [plotReceptorsExpression()], [plotSpatialColorTips()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' plotReceptors(x, group = "clone_id")
#'
#' # Color the spots by the number of receptors in each spot.
#' plotReceptors(x, value = "count")
#' @export
plotReceptors <- function(object,
                          receptors = NULL,
                          group = NULL,
                          group_val = NULL,
                          spatial_id = NULL,
                          colors = NULL,
                          overlap = "average",
                          value = NULL,
                          value_colors = c("#FFFFFF", "#000080"),
                          value_limits = c(0, NA),
                          value_combine = c("mean", "median", "mode", "max", "min", "sum"),
                          receptor_size = 1,
                          receptor_alpha = 1,
                          title = NULL,
                          order_points = TRUE,
                          image_scale = NULL,
                          swap_coords = NULL,
                          images = NULL,
                          sample_id = NULL,
                          image_alpha = 1,
                          zoom = FALSE,
                          receptor_legend_title = NULL,
                          show_receptor_legend = TRUE,
                          polygons = FALSE,
                          polygon_x = "x",
                          polygon_y = "y",
                          polygon_id = "polygon",
                          polygon_val = NULL,
                          polygon_alpha = 0.2,
                          polygon_fill = "white",
                          polygon_color = "white",
                          polygon_linewidth = 0.6,
                          polygon_label = FALSE,
                          polygon_label_by = NULL,
                          polygon_label_size = 3,
                          polygon_label_color = "white",
                          combine = TRUE,
                          combine_sample_plots = FALSE,
                          slot = NULL) {
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  panel_titles <- .sirup_panel_titles(title)
  if (!is.null(panel_titles)) title <- NULL
  options <- .sirup_spatial_options(object, image_scale = image_scale, swap_coords = swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  # Resolve any polygon overlay while `object` may still be the SIRUP object;
  # `.sirup_resolve_inputs()` replaces it with the plain Seurat object below.
  polygon_overlay <- .sirup_resolve_overlay_polygons(
    polygons = polygons,
    object = object,
    polygon_x = polygon_x,
    polygon_y = polygon_y,
    polygon_id = polygon_id,
    polygon_val = polygon_val,
    label_by = polygon_label_by
  )
  inputs <- .sirup_resolve_inputs(object, receptors, spatial_id, slot)
  object <- inputs$object
  receptors <- inputs$receptors
  spatial_id <- inputs$spatial_id_col
  sample_col <- .sirup_receptor_sample_col(receptors %||% data.frame(), inputs$receptor_fields, object)
  if (is.null(receptors)) {
    rlang::abort("`receptors` is required when `object` is not a SIRUP object.")
  }
  receptors <- .as_data_frame(receptors)
  spatial_id <- .guess_spatial_id_col(receptors, spatial_id)

  # Continuous mode: colour each spot on a gradient, either by the number of
  # receptors it carries (value = "count") or by a numeric receptor column.
  if (!is.null(value)) {
    value_combine <- match.arg(value_combine)
    receptors <- receptors[!is.na(receptors[[spatial_id]]), , drop = FALSE]
    if (!is.null(group)) {
      .check_required_cols(receptors, group, "receptors")
      gv <- if (is.null(group_val)) unique(stats::na.omit(as.character(receptors[[group]]))) else as.character(group_val)
      receptors <- receptors[as.character(receptors[[group]]) %in% gv, , drop = FALSE]
    }
    if (nrow(receptors) == 0L) {
      rlang::abort("No receptor rows with a spatial ID to summarise.")
    }
    # A spot is drawn once, so several receptors on one spot are reduced to a
    # single number first.
    summarised <- .sirup_receptor_spot_values(receptors, spatial_id, value, value_combine, sample_col = sample_col)
    counts <- summarised$data
    legend_name <- receptor_legend_title %||% summarised$legend
    plots <- list()
    for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
      spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
      image_polygon_overlay <- .sirup_filter_polygon_plot_sample(polygon_overlay, image)
      receptor_data <- .sirup_merge_receptor_lookup(spatial$coords, counts, spatial_id, sample_col)
      if (nrow(receptor_data) == 0L) {
        next
      }
      # Draw high-count spots last so they sit on top.
      if (isTRUE(order_points)) {
        receptor_data <- receptor_data[order(receptor_data$.sirup_count), , drop = FALSE]
      }
      p <- .sirup_spatial_base(
        receptor_data, spatial$image,
        image_alpha = image_alpha, zoom = zoom, view_data = spatial$coords
      )
      if (!is.null(image_polygon_overlay) && nrow(image_polygon_overlay) > 0L) {
        p <- p + .sirup_polygon_fill_layer(
          image_polygon_overlay, polygon_alpha = polygon_alpha,
          polygon_fill = polygon_fill, polygon_color = polygon_color,
          polygon_linewidth = polygon_linewidth
        )
      }
      p <- p +
        ggplot2::geom_point(
          data = receptor_data,
          ggplot2::aes(x = .data$x, y = .data$.sirup_y, color = .data$.sirup_count),
          size = receptor_size, alpha = receptor_alpha,
          inherit.aes = FALSE, na.rm = TRUE,
          show.legend = isTRUE(show_receptor_legend)
        ) +
        .sirup_value_color_scale(value_colors, value_limits, legend_name)
      if (!is.null(image_polygon_overlay) && nrow(image_polygon_overlay) > 0L) {
        label_layer <- .sirup_polygon_label_layer(
          image_polygon_overlay, polygon_label = polygon_label,
          polygon_label_by = polygon_label_by,
          polygon_label_size = polygon_label_size,
          polygon_label_color = polygon_label_color
        )
        if (!is.null(label_layer)) {
          p <- p + label_layer
        }
      }
      plots[[image]] <- p + ggplot2::labs(title = title %||% legend_name)
    }
    return(.sirup_return_spatial_plots(plots, combine = combine, combine_sample_plots = combine_sample_plots, panel_titles = panel_titles))
  }

  if (is.null(group)) {
    rlang::abort("`group` is required unless `value` is set.")
  }
  .check_required_cols(receptors, c(spatial_id, group), "receptors")

  receptors$.sirup_receptor_group <- as.character(receptors[[group]])
  if (is.null(group_val)) {
    group_val <- unique(stats::na.omit(receptors$.sirup_receptor_group))
  }
  group_val <- as.character(group_val)
  group_val <- group_val[!is.na(group_val)]
  if (length(group_val) == 0L) {
    rlang::abort("`group_val` must contain at least one non-missing receptor group.")
  }

  receptors <- receptors[
    !is.na(receptors$.sirup_receptor_group) &
      receptors$.sirup_receptor_group %in% group_val &
      !is.na(receptors[[spatial_id]]),
    ,
    drop = FALSE
  ]
  if (nrow(receptors) == 0L) {
    rlang::abort("No receptor rows matched `group_val` in `group`.")
  }
  receptors[[spatial_id]] <- as.character(receptors[[spatial_id]])
  receptors$.sirup_receptor_group <- factor(receptors$.sirup_receptor_group, levels = group_val)
  single_colors <- .sirup_receptor_group_colors(group_val, colors)
  receptor_legend_title <- receptor_legend_title %||% group

  receptor_lookup <- receptors[, unique(c(sample_col, spatial_id, ".sirup_receptor_group")), drop = FALSE]
  receptor_lookup <- receptor_lookup[!duplicated(receptor_lookup), , drop = FALSE]

  # Spots with more than one receptor group become a combined, colour-averaged
  # group (for example "BCR + TCR"), so each spot is drawn as a single point.
  combined <- .sirup_combine_receptor_groups(receptor_lookup, spatial_id, group_val, single_colors, overlap = overlap, sample_col = sample_col)
  receptor_lookup <- combined$lookup
  plot_groups <- combined$levels
  colors <- combined$colors
  plots <- list()
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    image_polygon_overlay <- .sirup_filter_polygon_plot_sample(polygon_overlay, image)
    receptor_data <- .sirup_merge_receptor_lookup(spatial$coords, receptor_lookup, spatial_id, sample_col)
    receptor_data$.sirup_receptor_group <- factor(
      receptor_data$.sirup_receptor_group,
      levels = plot_groups
    )
    if (isTRUE(order_points)) {
      receptor_data <- receptor_data[
        order(as.integer(receptor_data$.sirup_receptor_group), na.last = TRUE),
        ,
        drop = FALSE
      ]
    }

    if (nrow(receptor_data) == 0L) {
      next
    }

    p <- .sirup_spatial_base(
      receptor_data,
      spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = spatial$coords
    )

    # Draw polygons under the receptor points so both stay visible.
    if (!is.null(image_polygon_overlay) && nrow(image_polygon_overlay) > 0L) {
      p <- p + .sirup_polygon_fill_layer(
        image_polygon_overlay,
        polygon_alpha = polygon_alpha,
        polygon_fill = polygon_fill,
        polygon_color = polygon_color,
        polygon_linewidth = polygon_linewidth
      )
    }

    p <- p +
      ggplot2::geom_point(
        data = receptor_data,
        ggplot2::aes(x = .data$x, y = .data$.sirup_y, color = .data$.sirup_receptor_group),
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

    # Polygon labels sit on top of the points.
    if (!is.null(image_polygon_overlay) && nrow(image_polygon_overlay) > 0L) {
      label_layer <- .sirup_polygon_label_layer(
        image_polygon_overlay,
        polygon_label = polygon_label,
        polygon_label_by = polygon_label_by,
        polygon_label_size = polygon_label_size,
        polygon_label_color = polygon_label_color
      )
      if (!is.null(label_layer)) {
        p <- p + label_layer
      }
    }

    plots[[image]] <- p + ggplot2::labs(title = title %||% receptor_legend_title)
  }

  .sirup_return_spatial_plots(plots, combine = combine, combine_sample_plots = combine_sample_plots, panel_titles = panel_titles)
}

#' Plot tip colors on a tissue image
#'
#' `plotSpatialColorTips` plots the spots of tree tips. Each spot gets the color
#' of its tip.
#' 
#' If the object has multiple images or samples, the function draws one plot
#' for each image or sample. 
#'
#' @param object Seurat object or SIRUP object.
#' @param tip_data Data frame with a spatial ID column and a color column. For
#'   example, use the result of [assignTipColors()].
#' @param spatial_id Name of the column in `tip_data` that holds the Seurat spot
#'   names (default `"spatial_id_8um"`).
#' @param sample_col Name of the column in `tip_data` that matches the names of a
#'   list of spatial objects.
#' @param color Name of the column that holds the colors (default
#'   `"tip_color"`).
#' @param label Name of the column that holds the legend labels of the colors.
#' @param legend_title Title of the color legend (default `"Tip color"`).
#' @param legend_colors Named vector of colors for the legend.
#' @param title Plot title. Give one title for each plot, in the order of the
#'   plots, to title every plot of a multi-sample object.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param point_size Size of the spots. The default is 1.6.
#' @param point_alpha Opacity of the spots.
#' @param order_points Draw the early color groups first, so that the later
#'   groups show on top (default `TRUE`)?
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_alpha Opacity of the tissue image.
#' @param zoom Zoom to the plotted spots (default `FALSE`)? If `FALSE`, the plot
#'   shows the full tissue area that the spatial coordinates cover.
#' @param show_legend Show the color legend (default `FALSE`)?
#' @param polygons Draw polygons under the points on each spatial panel. With
#'   `FALSE` (default), the function draws no polygons. With `TRUE`, it uses the
#'   polygons stored in the SIRUP object. You can also give a data frame of
#'   polygon vertices, or the path of a CSV file of polygon vertices. A vertex
#'   table needs the columns set by `polygon_x` and `polygon_y`. It needs the
#'   column set by `polygon_id` if it holds more than one polygon.
#' @param polygon_x Name of the x-coordinate column in the polygon table
#'   (default `"x"`).
#' @param polygon_y Name of the y-coordinate column in the polygon table
#'   (default `"y"`).
#' @param polygon_id Name of the column that separates the polygons in the
#'   polygon table (default `"polygon"`). If the column does not exist, the
#'   function uses a `polygon` or `polygon_id` column. If neither exists, all
#'   rows form one polygon.
#' @param polygon_val IDs of the polygons to draw. If `NULL`, the function draws
#'   all polygons.
#' @param polygon_alpha Opacity of the polygon fill (default 0.2).
#' @param polygon_fill Fill color of the polygons (default `"white"`).
#' @param polygon_color Color of the polygon outlines (default `"white"`).
#' @param polygon_linewidth Width of the polygon outlines (default 0.6).
#' @param polygon_label Write a label at the center of each polygon (default
#'   `FALSE`)?
#' @param polygon_label_by Name of the column that holds the labels. If `NULL`,
#'   the function uses the polygon ID.
#' @param polygon_label_size Text size of the polygon labels (default 3).
#' @param polygon_label_color Text color of the polygon labels (default
#'   `"white"`).
#' @param combine Combine the image panels with patchwork (default `TRUE`)?
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#'
#' @return A ggplot object or a patchwork object.
#' @seealso [plotSpatialCloneColorTips()], [plotSpatialCloneTreePairs()],
#'   [assignTipColors()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' colors <- assignTipColors(getReceptors(x), shared_color = "spatial_id")
#' plotSpatialColorTips(x, tip_data = colors, spatial_id = "spatial_id")
#' @export
plotSpatialColorTips <- function(object,
                                 tip_data,
                                spatial_id = "spatial_id_8um",
                                sample_col = NULL,
                                color = "tip_color",
                                label = NULL,
                                legend_title = "Tip color",
                                legend_colors = NULL,
                                title = NULL,
                                image_scale = NULL,
                                swap_coords = NULL,
                                point_size = 1.6,
                                point_alpha = 1,
                                order_points = TRUE,
                                images = NULL,
                                sample_id = NULL,
                                image_alpha = 1,
                                 zoom = FALSE,
                                 show_legend = FALSE,
                                polygons = FALSE,
                                polygon_x = "x",
                                polygon_y = "y",
                                polygon_id = "polygon",
                                polygon_val = NULL,
                                polygon_alpha = 0.2,
                                polygon_fill = "white",
                                polygon_color = "white",
                                polygon_linewidth = 0.6,
                                polygon_label = FALSE,
                                polygon_label_by = NULL,
                                polygon_label_size = 3,
                                polygon_label_color = "white",
                                combine = TRUE,
                                combine_sample_plots = FALSE) {
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  panel_titles <- .sirup_panel_titles(title)
  if (!is.null(panel_titles)) title <- NULL
  spatial_id_missing <- missing(spatial_id)
  if (isSIRUP(object) && spatial_id_missing) {
    spatial_id <- .sirup_receptor_fields(object)$spatial_id
  }
  # Resolve any polygon overlay while `object` may still be the SIRUP object;
  # `.sirup_spatial_context()` replaces it with the plain Seurat object below.
  polygon_overlay <- .sirup_resolve_overlay_polygons(
    polygons = polygons,
    object = object,
    polygon_x = polygon_x,
    polygon_y = polygon_y,
    polygon_id = polygon_id,
    polygon_val = polygon_val,
    label_by = polygon_label_by
  )
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  tip_data <- .as_data_frame(tip_data)
  .check_required_cols(tip_data, c(spatial_id, color), "tip_data")
  sample_col <- sample_col %||% if (.sirup_is_spatial_list(object) && "sample" %in% colnames(tip_data)) "sample" else NULL
  if (!is.null(sample_col)) {
    .check_required_cols(tip_data, sample_col, "tip_data")
  }
  if (!is.null(label)) {
    .check_required_cols(tip_data, label, "tip_data")
  }
  tip_data <- tip_data[!is.na(tip_data[[spatial_id]]) & !is.na(tip_data[[color]]), , drop = FALSE]
  tip_key <- if (!is.null(sample_col)) {
    paste(as.character(tip_data[[sample_col]]), as.character(tip_data[[spatial_id]]), sep = "\r")
  } else {
    as.character(tip_data[[spatial_id]])
  }
  spot_color <- stats::setNames(as.character(tip_data[[color]]), tip_key)
  spot_color <- spot_color[!duplicated(names(spot_color))]
  if (is.null(label)) {
    colors <- unique(stats::na.omit(as.character(spot_color)))
    names(colors) <- colors
  } else {
    tip_data <- tip_data[!is.na(tip_data[[label]]), , drop = FALSE]
    tip_label_key <- if (!is.null(sample_col)) {
      paste(as.character(tip_data[[sample_col]]), as.character(tip_data[[spatial_id]]), sep = "\r")
    } else {
      as.character(tip_data[[spatial_id]])
    }
    spot_label <- stats::setNames(as.character(tip_data[[label]]), tip_label_key)
    spot_label <- spot_label[!duplicated(names(spot_label))]
    if (is.null(legend_colors)) {
      legend_df <- tip_data[!duplicated(tip_data[[label]]), c(label, color), drop = FALSE]
      legend_colors <- stats::setNames(as.character(legend_df[[color]]), as.character(legend_df[[label]]))
    }
    colors <- legend_colors
  }

  plots <- list()
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    image_polygon_overlay <- .sirup_filter_polygon_plot_sample(polygon_overlay, image)
    coord_key <- if (!is.null(sample_col) && "sample" %in% colnames(spatial$coords)) {
      paste(as.character(spatial$coords$sample), as.character(spatial$coords$cell), sep = "\r")
    } else {
      as.character(spatial$coords$cell)
    }
    dat <- spatial$coords[coord_key %in% names(spot_color), , drop = FALSE]
    if (nrow(dat) == 0L) {
      next
    }
    dat_key <- if (!is.null(sample_col) && "sample" %in% colnames(dat)) {
      paste(as.character(dat$sample), as.character(dat$cell), sep = "\r")
    } else {
      as.character(dat$cell)
    }
    dat$.sirup_tip_color <- spot_color[dat_key]
    if (!is.null(label)) {
      dat$.sirup_tip_label <- factor(spot_label[dat_key], levels = names(colors))
      dat$.sirup_tip_legend <- dat$.sirup_tip_label
    } else {
      dat$.sirup_tip_legend <- factor(dat$.sirup_tip_color, levels = names(colors))
    }
    if (isTRUE(order_points)) {
      dat <- dat[order(as.integer(dat$.sirup_tip_legend), na.last = TRUE), , drop = FALSE]
    }
    p <- .sirup_spatial_base(
      dat,
      spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = spatial$coords
    )
    # Draw polygons under the tip points so both stay visible.
    if (!is.null(image_polygon_overlay) && nrow(image_polygon_overlay) > 0L) {
      p <- p + .sirup_polygon_fill_layer(
        image_polygon_overlay,
        polygon_alpha = polygon_alpha,
        polygon_fill = polygon_fill,
        polygon_color = polygon_color,
        polygon_linewidth = polygon_linewidth
      )
    }
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$.sirup_tip_legend),
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
    # Polygon labels sit on top of the points.
    if (!is.null(image_polygon_overlay) && nrow(image_polygon_overlay) > 0L) {
      label_layer <- .sirup_polygon_label_layer(
        image_polygon_overlay,
        polygon_label = polygon_label,
        polygon_label_by = polygon_label_by,
        polygon_label_size = polygon_label_size,
        polygon_label_color = polygon_label_color
      )
      if (!is.null(label_layer)) {
        p <- p + label_layer
      }
    }
    if (!isTRUE(show_legend)) {
      p <- p + ggplot2::theme(legend.position = "none")
    }
    plots[[image]] <- p
  }

  .sirup_return_spatial_plots(plots, combine = combine, combine_sample_plots = combine_sample_plots, panel_titles = panel_titles)
}

#' Plot the tip colors of one clone on a tissue image
#'
#' `plotSpatialCloneColorTips` plots the spots of one clone. Each spot gets the
#' color of its tip. The function calls [plotSpatialColorTips()].
#'
#' @param colored_clones Result of [assignTipColors()], or a SIRUP object whose
#'   receptor table has the tip colors.
#' @param clone_val ID of the clone to plot.
#' @param tissue_positions Seurat object or SIRUP object. Omit it if
#'   `colored_clones` is a SIRUP object.
#' @param clone Name of the column that holds the clone IDs (default
#'   `"clone_id"`). If `colored_clones` is a SIRUP object and you omit this
#'   argument, the function uses the stored clone field.
#' @param spatial_id Name of the column in `colored_clones` that holds the
#'   spatial IDs (default `"spatial_id_8um"`). If `colored_clones` is a SIRUP
#'   object and you omit this argument, the function uses the stored spatial ID
#'   field.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`. If
#'   `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param ... Other arguments for [plotSpatialColorTips()].
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#'
#' @return A ggplot object.
#' @seealso [plotSpatialColorTips()], [assignTipColors()], [plotSpatialCloneTreePairs()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- updateReceptors(x, assignTipColors, shared_color = "spatial_id")
#' plotSpatialCloneColorTips(x, clone_val = "1")
#' @export
plotSpatialCloneColorTips <- function(colored_clones,
                                    clone_val,
                                    tissue_positions = NULL,
                                    clone = "clone_id",
                                    spatial_id = "spatial_id_8um",
                                    image_scale = NULL,
                                    swap_coords = NULL,
                                    sample_id = NULL,
                                    slot = NULL,
                                    ...) {
  clone_missing <- missing(clone)
  spatial_id_missing <- missing(spatial_id)
  spatial_options_object <- tissue_positions
  plot_object <- NULL
  if (isSIRUP(colored_clones)) {
    container <- colored_clones
    receptor_fields <- .sirup_receptor_fields(container)
    spatial_options_object <- container
    plot_object <- tissue_positions %||% container
    tissue_positions <- tissue_positions %||% container$spatial
    colored_clones <- .sirup_get_slot(container, slot)
    receptor_fields <- .sirup_slot_fields(receptor_fields, colored_clones)
    if (clone_missing && !is.null(receptor_fields$clone)) {
      clone <- receptor_fields$clone
    }
    if (spatial_id_missing) {
      spatial_id <- receptor_fields$spatial_id
    }
  }
  if (isSIRUP(tissue_positions) && spatial_id_missing) {
    spatial_id <- .sirup_receptor_fields(tissue_positions)$spatial_id
  }
  spatial_options_object <- spatial_options_object %||% tissue_positions
  options <- .sirup_spatial_options(
    spatial_options_object,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  tissue_positions <- .sirup_spatial_input(tissue_positions)
  clone_df <- .as_data_frame(colored_clones)
  .check_required_cols(clone_df, c(clone, spatial_id, "tip_color"), "colored_clones")
  clone_df <- clone_df[clone_df[[clone]] == clone_val, , drop = FALSE]

  if (!inherits(tissue_positions, "Seurat") && !.sirup_is_spatial_list(tissue_positions)) {
    rlang::abort("plotSpatialCloneColorTips() expects a Seurat object or named list of Seurat objects. Use ggplot directly for coordinate data frames.")
  }
  plot_object <- plot_object %||% tissue_positions
  plotSpatialColorTips(
    plot_object,
    tip_data = clone_df,
    spatial_id = spatial_id,
    sample_col = if ("sample" %in% colnames(clone_df)) "sample" else NULL,
    sample_id = sample_id,
    color = "tip_color",
    title = paste("Clone", clone_val),
    image_scale = image_scale,
    swap_coords = swap_coords,
    ...
  )
}

#' Plot each clone tree next to its tissue map
#'
#' `plotSpatialCloneTreePairs` makes two plots for each clone: the lineage tree,
#' and a tissue map of the tips. Each tip has one color in both plots, so you
#' can find a tip of the tree on the tissue.
#'
#' @details
#' A tip with no spatial ID, or with a spatial ID that is not in `object`, shows
#' in the tree plot only. The map shows only the tips whose spatial ID matches a
#' spot in `object`. In the tree, by default tips with a spot are circles, and tips without
#' a spot are squares.
#' 
#' If the object has multiple images or samples, the function draws one paired plot
#' for each image or sample.
#'
#' To match the tips to receptors, see [formatSpatialClones()].
#'
#' @param trees Result of [dowser::getTrees()], or a list with an element
#'   `trees` that holds this result.
#' @param object Seurat object or SIRUP object.
#' @param spatial_id Name of the column in the dowser clone data that holds the
#'   Seurat spot names. Omit it if `object` is a SIRUP object.
#' @param clone_val IDs of the clones to plot. If `NULL`, the function plots all
#'   clones.
#' @param id Name of the column in the dowser clone data that matches the tree
#'   tip labels (default `"sequence_id"`).
#' @param sequence_id Name of the column in the receptor table of `object` that
#'   decides which tree tips belong to the current spatial object (default: the
#'   value of `id`). The function uses it only if `object` is a SIRUP object. A
#'   tip that is not in this column stays in the tree plot, but the map does not
#'   show it.
#' @param clone Name of the column that holds the clone IDs (default
#'   `"clone_id"`).
#' @param color_by Column that sets the colors of the tips in both plots. With
#'   `"identity"` (default), each tip gets its own color. With another column,
#'   such as `"tree_cluster"` or a spatial ID column, tips that share a value
#'   get one color.
#' @param palette Colors to use. If you give fewer colors than the function
#'   needs, it makes more colors by averaging the colors that you give.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param zoom Zoom the maps to the clone (default `FALSE`)? If `FALSE`, the
#'   maps show the full tissue area that the spatial coordinates cover.
#' @param tipsize Size of the tip points of the tree.
#' @param scale Distance that the scale bar of the tree shows, in mutations for
#'   each site (default 0.01).
#' @param layout Layout of the tree, for `ggtree::ggtree()` (default
#'   `"rectangular"`).
#' @param title Title of each tree plot. With `TRUE` (default), the title is the
#'   clone ID. Give a character string to use it as the title.
#' @param shape_by Name of the column in the dowser clone data that sets the
#'   shape of the tree tips. The default marks a tip as `"spatial"` if it has a
#'   spatial ID, and as `"not_spatial"` if it has none.
#' @param shape_values Named vector of shapes, for
#'   [ggplot2::scale_shape_manual()] (default `c(spatial = 21, not_spatial =
#'   22)`). Shapes 21 to 25 are filled shapes: the tip color fills the shape,
#'   and the outline is black. Other shapes get the tip color on the whole mark.
#' @param tree_layout Alias for `layout`.
#' @param tree_tip_size Alias for `tipsize`.
#' @param tree_linewidth Line width of the tree branches, for `ggtree::ggtree()`
#'   (default 0.5).
#' @param tip_stroke Outline width of the filled tree tips (default 0.3).
#' @param treescale_fontsize Font size of the label of the scale bar (default
#'   3.88).
#' @param treescale Alias for `scale`.
#' @param treescale_position Position of the scale bar: `"bottomleft"`
#'   (default), `"topleft"` (near the root, at the top of the tips), or
#'   `"midleft"` (near the root, halfway down the tips).
#' @param spatial_tip_size Size of the tip points on the map (default 1.6).
#' @param spatial_tip_alpha Opacity of the tip points on the map (default 1).
#' @param order_points Draw the early color groups first on the map, so that the
#'   later groups show on top (default `FALSE`)?
#' @param show_spatial_legend Show a color legend below the map, based on the
#'   tip colors of the tree (default `FALSE`)?
#' @param legend_theme A ggplot2 theme, made with `ggplot2::theme()`, that
#'   changes the look of the shared legend (default `NULL`). Use it, for example,
#'   to make the legend text and keys smaller, so that a long legend fits on one
#'   line. The function uses it only if `show_spatial_legend` is `TRUE`.
#' @param combine_sample_plots Combine the maps of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns one tree and map pair for each sample.
#' @param polygons Draw polygons under the points on each spatial panel. With
#'   `FALSE` (default), the function draws no polygons. With `TRUE`, it uses the
#'   polygons stored in the SIRUP object. You can also give a data frame of
#'   polygon vertices, or the path of a CSV file of polygon vertices. A vertex
#'   table needs the columns set by `polygon_x` and `polygon_y`. It needs the
#'   column set by `polygon_id` if it holds more than one polygon.
#' @param polygon_x Name of the x-coordinate column in the polygon table
#'   (default `"x"`).
#' @param polygon_y Name of the y-coordinate column in the polygon table
#'   (default `"y"`).
#' @param polygon_id Name of the column that separates the polygons in the
#'   polygon table (default `"polygon"`). If the column does not exist, the
#'   function uses a `polygon` or `polygon_id` column. If neither exists, all
#'   rows form one polygon.
#' @param polygon_val IDs of the polygons to draw. If `NULL`, the function draws
#'   all polygons.
#' @param polygon_alpha Opacity of the polygon fill (default 0.2).
#' @param polygon_fill Fill color of the polygons (default `"white"`).
#' @param polygon_color Color of the polygon outlines (default `"white"`).
#' @param polygon_linewidth Width of the polygon outlines (default 0.6).
#' @param polygon_label Write a label at the center of each polygon (default
#'   `FALSE`)?
#' @param polygon_label_by Name of the column that holds the labels. If `NULL`,
#'   the function uses the polygon ID.
#' @param polygon_label_size Text size of the polygon labels (default 3).
#' @param polygon_label_color Text color of the polygon labels (default
#'   `"white"`).
#' @param nproc Number of processes for building the plot pairs (default 1).
#' @param ... Other arguments for [plotSpatialColorTips()].
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#'
#' @return A named list with one element for each clone. Each element is a list
#'   with the items `tree` and `spatial_clone` (the two plots), and `legend` and
#'   `legend_rel_height` (for the optional legend). Print an element to see the
#'   tree and the map side by side, with the legend beneath them. Use the `tree`
#'   and `spatial_clone` items to make your own layout.
#' @seealso [spatialCloneTreePairsToPDF()], [clusterTreesTips()],
#'   [getSpatialTrees()], [plotSpatialColorTips()], [plotSpatialCloneColorTips()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- formatSpatialClones(
#'   x,
#'   id = "sequence_id",
#'   germ = "germline_alignment_d_mask",
#'   traits = "spatial_id"
#' )
#' x <- getSpatialTrees(x)
#' trees <- getAnalysis(x, "trees")
#'
#' pairs <- plotSpatialCloneTreePairs(trees, x, clone_val = "4")
#' names(pairs)
#'
#' # Show the tree and the map side by side.
#' pairs[["4"]]
#' @export
plotSpatialCloneTreePairs <- function(trees,
                                      object,
                                      spatial_id = NULL,
                                      clone_val = NULL,
                                      id = "sequence_id",
                                      sequence_id = id,
                                      clone = "clone_id",
                                      color_by = "identity",
                                      palette = NULL,
                                      image_scale = NULL,
                                      swap_coords = NULL,
                                      images = NULL,
                                      sample_id = NULL,
                                      zoom = FALSE,
                                      tipsize = NULL,
                                      scale = 0.01,
                                      layout = "rectangular",
                                      title = TRUE,
                                      shape_by = ".sirup_spatial_status",
                                      shape_values = c(spatial = 21, not_spatial = 22),
                                      tree_layout = NULL,
                                      tree_tip_size = NULL,
                                      tree_linewidth = 0.5,
                                      tip_stroke = 0.3, treescale_fontsize = 3.88,
                                      treescale = NULL,
                                      treescale_position = c("bottomleft", "topleft", "midleft"),
                                      spatial_tip_size = 1.6,
                                      spatial_tip_alpha = 1,
                                      order_points = FALSE,
                                      show_spatial_legend = FALSE,
                                      legend_theme = NULL,
                                      combine_sample_plots = FALSE,
                                      polygons = FALSE,
                                      polygon_x = "x",
                                      polygon_y = "y",
                                      polygon_id = "polygon",
                                      polygon_val = NULL,
                                      polygon_alpha = 0.2,
                                      polygon_fill = "white",
                                      polygon_color = "white",
                                      polygon_linewidth = 0.6,
                                      polygon_label = FALSE,
                                      polygon_label_by = NULL,
                                      polygon_label_size = 3,
                                      polygon_label_color = "white",
                                      nproc = 1L,
                                      slot = NULL,
                                      ...) {
  if (!is.null(legend_theme) && !inherits(legend_theme, "theme")) {
    rlang::abort("`legend_theme` must be NULL or a ggplot2 theme, made with `ggplot2::theme()`.")
  }
  object_tip_ids <- NULL
  object_spatial_ids <- NULL
  object_samples <- NULL
  sample_field <- NULL
  uid_col <- ".sirup_uid"
  if (!is.null(sample_id)) {
    sample_id <- .sirup_plot_sample_id(object, sample_id)
    images <- unique(c(images %||% character(), as.character(sample_id)))
  }
  if (isSIRUP(object)) {
    if (is.null(spatial_id)) {
      spatial_id <- .sirup_receptor_fields(object)$spatial_id
    }
    receptor_data <- getReceptors(object, slot = slot)
    .check_required_cols(receptor_data, c(sequence_id, spatial_id), "object receptor data")
    object_tip_ids <- unique(stats::na.omit(as.character(receptor_data[[sequence_id]])))
    object_spatial_ids <- unique(stats::na.omit(as.character(receptor_data[[spatial_id]])))
    sample_field <- .sirup_receptor_fields(object)$sample
    if (!is.null(sample_field) && sample_field %in% colnames(receptor_data)) {
      object_samples <- unique(stats::na.omit(as.character(receptor_data[[sample_field]])))
    }
  }
  # Resolve the optional polygon overlay while `object` still holds the SIRUP
  # object, because `polygons = TRUE` reads the stored polygon table from it.
  # `.sirup_spatial_context()` replaces `object` with the plain Seurat object.
  polygon_overlay <- .sirup_resolve_overlay_polygons(
    polygons = polygons,
    object = object,
    polygon_x = polygon_x,
    polygon_y = polygon_y,
    polygon_id = polygon_id,
    polygon_val = polygon_val,
    label_by = polygon_label_by
  )
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (is.null(spatial_id)) {
    rlang::abort("`spatial_id` is required when `object` is not a SIRUP object.")
  }
  nproc <- .sirup_plot_check_nproc(nproc)
  records <- .sirup_dowser_tree_records(trees)
  if (!is.null(clone_val)) {
    records <- records[names(records) %in% as.character(clone_val)]
  }
  if (length(records) == 0L) {
    rlang::abort("No dowser trees matched the requested clone IDs.")
  }
  layout <- tree_layout %||% layout
  tipsize <- tree_tip_size %||% tipsize %||% 2.5
  scale <- treescale %||% scale
  treescale_position <- match.arg(treescale_position)

  plot_one_clone <- function(clone_id) {
    record <- records[[clone_id]]
    tree <- record$tree
    clone_data <- .sirup_airr_clone_data(record$data, clone_id = clone_id, clone_col = clone)
    .check_required_cols(clone_data, c(id, spatial_id, clone), "dowser clone data")

    tip_ids <- setdiff(tree$tip.label, "Germline")
    tip_data <- clone_data[clone_data[[id]] %in% tip_ids, , drop = FALSE]
    if (nrow(tip_data) == 0L) {
      return(NULL)
    }
    if (!is.null(object_tip_ids)) {
      # Match tree tips back to this object's receptors. Prefer the verbatim id
      # copy (`.sirup_uid`) carried by `formatSpatialClones()`: dowser sanitizes
      # tip labels for Newick, so the tree `id` column no longer equals the raw
      # receptor id. dowser also comma-joins text fields when it collapses
      # identical sequences, so a tip may carry several ids; it belongs to this
      # object if ANY of them matches. Trees built before `.sirup_uid` existed
      # fall back to matching on the (verbatim) spatial id.
      if (uid_col %in% colnames(tip_data)) {
        tip_id_lists <- strsplit(as.character(tip_data[[uid_col]]), ",", fixed = TRUE)
        in_current_object <- vapply(
          tip_id_lists,
          function(ids) any(trimws(ids) %in% object_tip_ids),
          logical(1)
        )
      } else {
        in_current_object <- as.character(tip_data[[spatial_id]]) %in% object_spatial_ids
      }
      # When the SIRUP `sample` field is set and carried into the tree data,
      # also require the tip's sample to belong to this object. This keeps spots
      # that share a barcode name across sections from being drawn on the wrong
      # section.
      if (!is.null(object_samples) && sample_field %in% colnames(tip_data)) {
        in_current_object <- in_current_object &
          (as.character(tip_data[[sample_field]]) %in% object_samples)
      }
      tip_data[[spatial_id]] <- as.character(tip_data[[spatial_id]])
      tip_data[[spatial_id]][!in_current_object] <- NA_character_
    }
    shared_color_col <- NULL
    if (!identical(color_by, "identity")) {
      .check_required_cols(tip_data, color_by, "dowser clone data")
      shared_color_col <- color_by
    }
    if (identical(shape_by, ".sirup_spatial_status")) {
      spatial_vals <- as.character(tip_data[[spatial_id]])
      tip_data[[shape_by]] <- ifelse(!is.na(spatial_vals) & nzchar(spatial_vals), "spatial", "not_spatial")
    } else if (!is.null(shape_by)) {
      .check_required_cols(tip_data, shape_by, "dowser clone data")
    }
    colored_tips <- assignTipColors(
      tip_data,
      clone_val = clone_id,
      clone = clone,
      id = id,
      shared_color = shared_color_col,
      palette = palette
    )
    legend_info <- .sirup_spatial_tip_legend_info(colored_tips, color_by = color_by)
    colored_tips <- legend_info$tip_data
    tree_title <- if (isTRUE(title)) {
      paste("Clone", clone_id)
    } else if (isFALSE(title) || is.null(title)) {
      NULL
    } else {
      as.character(title)
    }

    tree_plot <- .sirup_plot_colored_tree(
      tree,
      colored_tips,
      id_col = id,
      color_col = "tip_color",
      title = tree_title,
      layout = layout,
      tip_size = tipsize,
      linewidth = tree_linewidth,
      stroke = tip_stroke, treescale_fontsize = treescale_fontsize,
      treescale = scale, treescale_position = treescale_position,
      shape_col = shape_by,
      shape_values = shape_values,
      label_col = ".sirup_tip_legend_label",
      legend_title = "Tip cluster",
      legend_colors = legend_info$colors,
      show_legend = show_spatial_legend,
      legend_position = "bottom",
      legend_ncol = min(6L, length(legend_info$colors))
    )
    if (!is.null(legend_theme) && isTRUE(show_spatial_legend)) {
      tree_plot <- tree_plot + legend_theme
    }
    tree_legend <- if (isTRUE(show_spatial_legend)) {
      cowplot::get_legend(tree_plot)
    } else {
      NULL
    }
    tree_plot <- tree_plot + ggplot2::theme(legend.position = "none")

    spatial_plot <- plotSpatialColorTips(
      object,
      tip_data = colored_tips,
      spatial_id = spatial_id,
      sample_col = sample_field,
      color = "tip_color",
      title = NULL,
      image_scale = image_scale,
      swap_coords = swap_coords,
      images = images,
      sample_id = sample_id,
      point_size = spatial_tip_size,
      point_alpha = spatial_tip_alpha,
      order_points = order_points,
      zoom = zoom,
      show_legend = FALSE,
      combine_sample_plots = combine_sample_plots,
      ...
    )
    overlay_image <- .sirup_spatial_tip_overlay_image(
      colored_tips,
      spatial_id = spatial_id,
      sample_field = sample_field,
      sample_id = sample_id
    )
    # Draw region polygons beneath the tip points so both stay visible.
    spatial_plot <- .sirup_add_polygon_overlay(
      spatial_plot,
      polygon_overlay,
      image = overlay_image,
      polygon_alpha = polygon_alpha,
      polygon_fill = polygon_fill,
      polygon_color = polygon_color,
      polygon_linewidth = polygon_linewidth,
      polygon_label = polygon_label,
      polygon_label_by = polygon_label_by,
      polygon_label_size = polygon_label_size,
      polygon_label_color = polygon_label_color
    )
    if (is.null(spatial_plot)) {
      return(NULL)
    }
    build_pair <- function(spatial_panel) {
      # Return the shared legend separately so the pair assembler can place it
      # full-width beneath BOTH panels (avoids clipping when there are color and
      # shape guides). legend_rel_height scales with the number of colour rows.
      structure(
        list(
          tree = tree_plot,
          spatial_clone = spatial_panel,
          legend = tree_legend,
          legend_rel_height = max(0.12, 0.07 * ceiling(length(legend_info$colors) / 6))
        ),
        class = c("sirup_clone_tree_pair", "list")
      )
    }
    if (is.list(spatial_plot) && !inherits(spatial_plot, "ggplot")) {
      return(lapply(spatial_plot, build_pair))
    }
    build_pair(spatial_plot)
  }

  out <- .sirup_plot_lapply(names(records), plot_one_clone, nproc = nproc)
  names(out) <- names(records)
  out <- out[!vapply(out, is.null, logical(1))]
  nested <- vapply(out, function(x) is.list(x) && length(x) > 0L && is.null(x$tree), logical(1))
  if (any(nested)) {
    flat <- list()
    for (clone_id in names(out)) {
      value <- out[[clone_id]]
      if (isTRUE(nested[[clone_id]])) {
        parsed_samples <- vapply(names(value), function(sample_name) {
          parsed <- .sirup_parse_spatial_image_label(sample_name)
          if (is.null(parsed)) NA_character_ else parsed$sample
        }, character(1))
        sample_counts <- table(stats::na.omit(parsed_samples))
        for (sample_name in names(value)) {
          parsed <- .sirup_parse_spatial_image_label(sample_name)
          flat_sample_name <- if (!is.null(parsed) &&
              !is.na(parsed$sample) &&
              sample_counts[[parsed$sample]] == 1L) {
            parsed$sample
          } else {
            sample_name
          }
          flat[[paste(clone_id, flat_sample_name, sep = "::")]] <- value[[sample_name]]
        }
      } else {
        flat[[clone_id]] <- value
      }
    }
    out <- flat
  }

  out
}

# Resolve the `polygons` argument shared by the spatial plotting functions into
# a plot-ready vertex table (columns `x`, `.sirup_y`, `polygon`, plus
# `label_by`). `FALSE` returns NULL; `TRUE` reads the SIRUP object polygons; a
# data frame or CSV path is used directly. Reuses the same helpers as
# `plotSpatialPolygons()`.
.sirup_resolve_overlay_polygons <- function(polygons,
                                            object,
                                            polygon_x = "x",
                                            polygon_y = "y",
                                            polygon_id = "polygon",
                                            polygon_val = NULL,
                                            label_by = NULL) {
  if (is.null(polygons) || isFALSE(polygons)) {
    return(NULL)
  }
  if (isTRUE(polygons)) {
    if (!isSIRUP(object)) {
      rlang::abort("`polygons = TRUE` needs a SIRUP object with stored polygons. Pass a data frame or CSV path instead.")
    }
    polygon <- .sirup_polygon_vertices_table(.sirup_polygons(object))
    if (nrow(polygon) == 0L) {
      rlang::abort("The SIRUP object does not contain stored polygons.")
    }
  } else if (is.character(polygons)) {
    polygon <- .sirup_read_polygons_path(polygons)
  } else if (is.data.frame(polygons)) {
    polygon <- polygons
  } else {
    rlang::abort("`polygons` must be FALSE, TRUE, a data frame of polygon vertices, or a path to a CSV file.")
  }
  polygon <- .as_data_frame(polygon)
  # Fall back to the common polygon-id column names when the requested one is
  # absent, so a stored or exported polygon table works without extra arguments.
  if (!polygon_id %in% names(polygon)) {
    alt <- intersect(c("polygon", "polygon_id"), names(polygon))
    if (length(alt) > 0L) {
      polygon_id <- alt[[1L]]
    }
  }
  .check_required_cols(polygon, c(polygon_x, polygon_y), "polygons")
  keep_cols <- unique(c(
    if (is.null(label_by)) character() else label_by,
    .sirup_polygon_sample_col(polygon)
  ))
  .sirup_polygon_plot_data(
    polygon = polygon,
    polygon_id = polygon_val,
    polygon_col = polygon_id,
    x_col = polygon_x,
    y_col = polygon_y,
    keep_cols = keep_cols
  )
}

# The polygon fill/outline layer. Add this to a `.sirup_spatial_base()` plot
# BEFORE the receptor/tip points so the points stay visible on top.
.sirup_polygon_fill_layer <- function(polygon_data,
                                      polygon_alpha = 0.2,
                                      polygon_fill = "white",
                                      polygon_color = "white",
                                      polygon_linewidth = 0.6) {
  ggplot2::geom_polygon(
    data = polygon_data,
    ggplot2::aes(x = .data$x, y = .data$.sirup_y, group = .data$polygon),
    inherit.aes = FALSE,
    fill = polygon_fill,
    color = polygon_color,
    alpha = polygon_alpha,
    linewidth = polygon_linewidth
  )
}

# The optional polygon center-label layer, or NULL when labels are off.
.sirup_polygon_label_layer <- function(polygon_data,
                                       polygon_label = FALSE,
                                       polygon_label_by = NULL,
                                       polygon_label_size = 3,
                                       polygon_label_color = "white") {
  if (!isTRUE(polygon_label)) {
    return(NULL)
  }
  centers <- .sirup_polygon_centers(polygon_data, unique(polygon_data$polygon))
  if (!is.null(polygon_label_by) && polygon_label_by %in% names(polygon_data)) {
    lookup <- polygon_data[!duplicated(polygon_data$polygon), c("polygon", polygon_label_by), drop = FALSE]
    centers$.sirup_label <- as.character(lookup[[polygon_label_by]][match(centers$polygon, lookup$polygon)])
  } else {
    centers$.sirup_label <- as.character(centers$polygon)
  }
  ggplot2::geom_text(
    data = centers,
    ggplot2::aes(x = .data$x, y = .data$.sirup_y, label = .data$.sirup_label),
    inherit.aes = FALSE,
    size = polygon_label_size,
    color = polygon_label_color,
    fontface = "bold"
  )
}

.sirup_spatial_tip_overlay_image <- function(tip_data,
                                             spatial_id,
                                             sample_field = NULL,
                                             sample_id = NULL) {
  if (!is.null(sample_id) && length(sample_id) == 1L) {
    return(as.character(sample_id))
  }
  if (is.null(sample_field) || !sample_field %in% colnames(tip_data)) {
    return(NULL)
  }
  spatial_values <- as.character(tip_data[[spatial_id]])
  keep <- !is.na(spatial_values) & nzchar(spatial_values)
  samples <- unique(as.character(stats::na.omit(tip_data[[sample_field]][keep])))
  samples <- samples[nzchar(samples)]
  if (length(samples) == 1L) samples[[1L]] else NULL
}

# Add a polygon overlay to a FINISHED spatial panel (used by
# `plotSpatialCloneTreePairs()`, which consumes `plotSpatialColorTips()`
# output). The fill is inserted beneath the tip points (after the image raster
# layer) so both stay visible; labels are drawn on top. Handles a single ggplot
# and a patchwork of panels.
.sirup_add_polygon_overlay <- function(plot,
                                       polygon_data,
                                       image = NULL,
                                       polygon_alpha = 0.2,
                                       polygon_fill = "white",
                                       polygon_color = "white",
                                       polygon_linewidth = 0.6,
                                       polygon_label = FALSE,
                                       polygon_label_by = NULL,
                                       polygon_label_size = 3,
                                       polygon_label_color = "white") {
  if (is.null(plot) || is.null(polygon_data) || nrow(polygon_data) == 0L) {
    return(plot)
  }
  if (!is.null(image) && length(image) == 1L && !is.na(image) && nzchar(as.character(image))) {
    polygon_data <- .sirup_filter_polygon_plot_sample(polygon_data, as.character(image))
    if (nrow(polygon_data) == 0L) {
      return(plot)
    }
  }
  if (is.list(plot) && !inherits(plot, "ggplot")) {
    plot_names <- names(plot)
    out <- lapply(seq_along(plot), function(i) {
      name <- if (is.null(plot_names)) NA_character_ else plot_names[[i]]
      .sirup_add_polygon_overlay(
        plot[[i]],
        polygon_data = .sirup_filter_polygon_plot_sample(polygon_data, name),
        image = NULL,
        polygon_alpha = polygon_alpha,
        polygon_fill = polygon_fill,
        polygon_color = polygon_color,
        polygon_linewidth = polygon_linewidth,
        polygon_label = polygon_label,
        polygon_label_by = polygon_label_by,
        polygon_label_size = polygon_label_size,
        polygon_label_color = polygon_label_color
      )
    })
    names(out) <- plot_names
    return(out)
  }
  polygon_geom <- .sirup_polygon_fill_layer(
    polygon_data,
    polygon_alpha = polygon_alpha,
    polygon_fill = polygon_fill,
    polygon_color = polygon_color,
    polygon_linewidth = polygon_linewidth
  )
  label_geom <- .sirup_polygon_label_layer(
    polygon_data,
    polygon_label = polygon_label,
    polygon_label_by = polygon_label_by,
    polygon_label_size = polygon_label_size,
    polygon_label_color = polygon_label_color
  )
  if (inherits(plot, "ggplot")) {
    # Insert the polygon fill just after the tissue-image raster so it sits
    # below the tip points that `plotSpatialColorTips()` drew.
    is_raster <- vapply(plot$layers, function(l) inherits(l$geom, "GeomRasterAnn"), logical(1))
    after <- if (any(is_raster)) max(which(is_raster)) else 0L
    plot$layers <- append(plot$layers, list(polygon_geom), after = after)
    if (!is.null(label_geom)) {
      plot <- plot + label_geom
    }
    return(plot)
  }
  # Patchwork of several image panels: add to every panel (drawn on top).
  plot <- plot & polygon_geom
  if (!is.null(label_geom)) {
    plot <- plot & label_geom
  }
  plot
}

#' Save the tree and map pairs to a PDF file
#'
#' `spatialCloneTreePairsToPDF` writes one page for each clone to a PDF file.
#' Each page shows the tree and the tissue map of the clone side by side.
#'
#' @param pairs Result of [plotSpatialCloneTreePairs()].
#' @param path Path of the PDF file.
#' @param width Width of a page, in inches (default 12).
#' @param height Height of a page, in inches (default 6).
#' @param rel_widths Relative widths of the tree panel and the map panel
#'   (default `c(1.2, 1)`).
#' @param nproc Number of processes for the pages (default 1). The function
#'   needs the `qpdf` package to join the pages that it makes in parallel. If
#'   `qpdf` is not installed, it makes the pages one after the other.
#'
#' @return `path`, invisibly.
#' @seealso [plotSpatialCloneTreePairs()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- formatSpatialClones(
#'   x,
#'   id = "sequence_id",
#'   germ = "germline_alignment_d_mask",
#'   traits = "spatial_id"
#' )
#' x <- getSpatialTrees(x)
#' trees <- getAnalysis(x, "trees")
#'
#' pairs <- plotSpatialCloneTreePairs(trees, x, clone_val = c("1", "4"))
#' spatialCloneTreePairsToPDF(pairs, tempfile(fileext = ".pdf"))
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
  nproc <- .sirup_plot_check_nproc(nproc)
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (nproc > 1L && .Platform$OS.type != "windows" && requireNamespace("qpdf", quietly = TRUE)) {
    tmp_dir <- tempfile("sirup_clone_tree_pairs_pdf_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

    page_paths <- file.path(tmp_dir, sprintf("page_%05d.pdf", seq_along(pairs)))
    .sirup_plot_lapply(seq_along(pairs), function(i) {
      .sirup_write_clone_tree_pair_pdf(
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
    print(.sirup_clone_tree_pair_grid(pair, rel_widths = rel_widths))
  }

  invisible(path)
}

# Average a set of colors in sRGB space and return one hex color.
.sirup_average_colors <- function(colors) {
  rgb_mat <- grDevices::col2rgb(colors)
  avg <- round(rowMeans(rgb_mat))
  grDevices::rgb(avg[[1]], avg[[2]], avg[[3]], maxColorValue = 255)
}

# Collapse per-spot receptor groups into one group per spot. A spot that carries
# more than one distinct group becomes a combined group (members joined by
# " + ", ordered as in `group_val`), coloured by the average of its members'
# colors. Single-group spots are unchanged. Returns the per-spot lookup, the
# ordered plot levels (single groups first, then combinations), and the colors.
.sirup_combine_receptor_groups <- function(lookup, spatial_id, group_val, single_colors, sep = " + ", overlap = "average", sample_col = NULL) {
  if (!identical(overlap, "average")) {
    if (length(overlap) != 1L || is.na(overlap) ||
        inherits(try(grDevices::col2rgb(overlap), silent = TRUE), "try-error")) {
      rlang::abort('`overlap` must be "average" or a single color.')
    }
  }
  lookup[[spatial_id]] <- as.character(lookup[[spatial_id]])
  groups <- as.character(lookup$.sirup_receptor_group)
  split_key <- if (!is.null(sample_col) && sample_col %in% colnames(lookup)) {
    paste(as.character(lookup[[sample_col]]), lookup[[spatial_id]], sep = "\r")
  } else {
    lookup[[spatial_id]]
  }
  by_spot <- split(groups, split_key)
  labels <- vapply(by_spot, function(v) {
    v <- group_val[group_val %in% unique(v)]
    paste(v, collapse = sep)
  }, character(1))

  if (!is.null(sample_col) && sample_col %in% colnames(lookup)) {
    key_parts <- strsplit(names(labels), "\r", fixed = TRUE)
    out <- data.frame(
      sample = vapply(key_parts, `[[`, character(1), 1L),
      spatial_id = vapply(key_parts, `[[`, character(1), 2L),
      .sirup_receptor_group = unname(labels),
      stringsAsFactors = FALSE
    )
    names(out)[1:2] <- c(sample_col, spatial_id)
  } else {
    out <- data.frame(names(labels), unname(labels), stringsAsFactors = FALSE)
    names(out) <- c(spatial_id, ".sirup_receptor_group")
  }

  observed <- unique(out$.sirup_receptor_group)
  singles <- group_val[group_val %in% observed]
  combos <- sort(setdiff(observed, singles))
  levels_out <- c(singles, combos)
  out$.sirup_receptor_group <- factor(out$.sirup_receptor_group, levels = levels_out)

  colors <- stats::setNames(character(length(levels_out)), levels_out)
  for (lv in levels_out) {
    members <- strsplit(lv, sep, fixed = TRUE)[[1]]
    member_colors <- single_colors[members]
    member_colors <- member_colors[!is.na(member_colors)]
    colors[[lv]] <- if (length(member_colors) == 0L) {
      "#808080"
    } else if (length(member_colors) == 1L) {
      unname(member_colors[[1L]])          # keep single-group colors unchanged
    } else if (identical(overlap, "average")) {
      .sirup_average_colors(member_colors)   # combinations get the averaged color
    } else {
      overlap                                # or a single fixed overlap color
    }
  }

  list(lookup = out, levels = levels_out, colors = colors)
}

.sirup_receptor_group_colors <- function(values, colors = NULL) {
  n <- length(values)
  if (is.null(colors)) {
    if (requireNamespace("pals", quietly = TRUE)) {
      colors <- .sirup_extend_palette(pals::glasbey(32), n)
    } else {
      colors <- .sirup_extend_palette(grDevices::hcl.colors(max(3L, n), palette = "Dark 3"), n)
    }
    names(colors) <- values
    return(colors)
  }

  colors <- as.character(colors)
  if (!is.null(names(colors)) && all(values %in% names(colors))) {
    return(colors[values])
  }
  if (length(colors) < n) {
    colors <- .sirup_extend_palette(colors, n)
  }
  colors <- colors[seq_len(n)]
  names(colors) <- values
  colors
}

.guess_spatial_id_col <- function(data, spatial_id_col = NULL) {
  if (!is.null(spatial_id_col)) {
    .check_required_cols(data, spatial_id_col, "receptors")
    return(spatial_id_col)
  }

  candidates <- c(
    "spatial_id",
    "spatial_id_8um",
    "spatial_id_2um",
    "spatial_id_16um",
    "spatial_id",
    "barcode_8um",
    "barcode_2um",
    "cell",
    "cell_barcode"
  )
  spatial_id_col <- candidates[candidates %in% colnames(data)][1]
  if (is.na(spatial_id_col)) {
    rlang::abort("Could not guess `spatial_id`; pass the receptor spatial ID column explicitly.")
  }
  spatial_id_col
}

.sirup_spatial_tip_legend_info <- function(tip_data, color_by) {
  tip_data$.sirup_tip_legend_label <- .sirup_spatial_tip_labels(tip_data, color_by = color_by)
  legend_df <- tip_data[
    !is.na(tip_data$.sirup_tip_legend_label) &
      !is.na(tip_data$tip_color),
    c(".sirup_tip_legend_label", "tip_color", "tip_color_id"),
    drop = FALSE
  ]
  legend_df <- legend_df[order(as.integer(legend_df$tip_color_id)), , drop = FALSE]
  legend_df <- legend_df[!duplicated(legend_df$.sirup_tip_legend_label), , drop = FALSE]

  spatial <- grepl("^spatial_tip_cluster_[0-9]+$", legend_df$.sirup_tip_legend_label)
  spatial_df <- legend_df[spatial, , drop = FALSE]
  if (nrow(spatial_df) > 0L) {
    spatial_id <- as.integer(sub("^spatial_tip_cluster_", "", spatial_df$.sirup_tip_legend_label))
    spatial_df <- spatial_df[order(spatial_id), , drop = FALSE]
  }
  other_df <- legend_df[!spatial, , drop = FALSE]
  legend_df <- rbind(spatial_df, other_df)

  list(
    tip_data = tip_data,
    colors = stats::setNames(
      as.character(legend_df$tip_color),
      as.character(legend_df$.sirup_tip_legend_label)
    )
  )
}

.sirup_spatial_tip_labels <- function(tip_data, color_by) {
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

.sirup_plot_check_nproc <- function(nproc) {
  nproc <- as.integer(nproc)
  if (is.na(nproc) || nproc < 1L) {
    rlang::abort("`nproc` must be a positive integer.")
  }
  nproc
}

.sirup_plot_lapply <- function(x, fun, nproc) {
  if (nproc == 1L || length(x) <= 1L) {
    return(lapply(x, fun))
  }
  if (.Platform$OS.type == "windows") {
    warning("Parallel plotting with `nproc > 1` is not supported on Windows; using `nproc = 1`.", call. = FALSE)
    return(lapply(x, fun))
  }
  parallel::mclapply(x, fun, mc.cores = nproc)
}

.sirup_write_clone_tree_pair_pdf <- function(pair,
                                             path,
                                             width,
                                             height,
                                             rel_widths) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(.sirup_clone_tree_pair_grid(pair, rel_widths = rel_widths))
  invisible(path)
}

# Printing a tree and map pair draws the tree and the map side by side, with the
# shared legend beneath them. The `tree`, `spatial_clone`, and `legend` items
# stay available for custom layouts.
#' @export
#' @noRd
print.sirup_clone_tree_pair <- function(x, ...) {
  print(.sirup_clone_tree_pair_grid(x, rel_widths = c(1.2, 1)))
  invisible(x)
}

.sirup_clone_tree_pair_grid <- function(pair, rel_widths) {
  top_row <- cowplot::plot_grid(
    pair$tree,
    pair$spatial_clone,
    ncol = 2,
    rel_widths = rel_widths
  )
  if (is.null(pair$legend)) {
    return(top_row)
  }
  # Shared legend spans the full width beneath both panels so it is not clipped.
  rel_h <- if (is.null(pair$legend_rel_height)) 0.14 else pair$legend_rel_height
  cowplot::plot_grid(
    top_row,
    pair$legend,
    ncol = 1,
    rel_heights = c(1, rel_h)
  )
}

.sirup_dowser_tree_records <- function(trees) {
  if (is.list(trees) && !is.data.frame(trees) && "trees" %in% names(trees)) {
    trees <- trees$trees
  }
  if (is.data.frame(trees) && all(c("trees", "data") %in% colnames(trees))) {
    records <- lapply(seq_len(nrow(trees)), function(i) {
      clone_id <- if ("clone_id" %in% colnames(trees)) as.character(trees$clone_id[[i]]) else NULL
      tree <- trees$trees[[i]]
      clone_data <- trees$data[[i]]
      clone_id <- clone_id %||% tree$name %||% .sirup_airr_clone_id(clone_data)
      list(tree = tree, data = clone_data)
    })
    if ("clone_id" %in% colnames(trees)) {
      names(records) <- as.character(trees$clone_id)
    } else {
      names(records) <- vapply(seq_along(records), function(i) {
        records[[i]]$tree$name %||% .sirup_airr_clone_id(records[[i]]$data) %||% paste0("clone_", i)
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

.sirup_airr_clone_id <- function(x) {
  if (inherits(x, "airrClone")) {
    return(x@clone)
  }
  NULL
}

.sirup_airr_clone_data <- function(x, clone_id, clone_col) {
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

.sirup_plot_colored_tree <- function(tree,
                                     tip_data,
                                     id_col,
                                     color_col,
                                     title,
                                     layout,
                                     tip_size,
                                     treescale,
                                     shape_col = NULL,
                                     shape_values = NULL,
                                     label_col = NULL,
                                     legend_title = "Tip cluster",
                                     legend_colors = NULL,
                                     show_legend = FALSE,
                                     legend_position = "right",
                                     legend_ncol = NULL,
                                     linewidth = 0.5,
                                     stroke = 0.3, treescale_fontsize = 3.88,
                                     treescale_position = "bottomleft") {
  p <- ggtree::ggtree(tree, layout = layout, linewidth = linewidth)
  join_cols <- c(id_col, color_col)
  if (!is.null(label_col)) {
    join_cols <- c(join_cols, label_col)
  }
  if (!is.null(shape_col)) {
    join_cols <- c(join_cols, shape_col)
  }
  join_data <- tip_data[, join_cols, drop = FALSE]
  join_data <- join_data[!duplicated(join_data[[id_col]]), , drop = FALSE]
  idx <- match(p$data$label, join_data[[id_col]])
  p$data[[color_col]] <- join_data[[color_col]][idx]
  if (!is.null(label_col)) {
    p$data[[label_col]] <- factor(join_data[[label_col]][idx], levels = names(legend_colors))
  }
  if (!is.null(shape_col)) {
    p$data[[shape_col]] <- join_data[[shape_col]][idx]
    if (!is.null(shape_values)) {
      p$data[[shape_col]] <- factor(p$data[[shape_col]], levels = names(shape_values))
    }
  }
  color_aes_col <- if (is.null(label_col)) color_col else label_col
  # Filled shapes (pch 21-25) carry the colour on `fill` and draw a black
  # outline; solid shapes keep the colour on `colour`. This lets the default
  # tip shapes be outlined, filled markers.
  fillable <- !is.null(shape_values) && all(shape_values %in% 21:25)
  tip_aes <- if (is.null(shape_col)) {
    if (fillable) ggplot2::aes(fill = .data[[color_aes_col]]) else ggplot2::aes(color = .data[[color_aes_col]])
  } else if (fillable) {
    ggplot2::aes(fill = .data[[color_aes_col]], shape = .data[[shape_col]])
  } else {
    ggplot2::aes(color = .data[[color_aes_col]], shape = .data[[shape_col]])
  }

  tip_layer <- if (is.null(shape_col)) {
    ggtree::geom_tippoint(tip_aes, size = tip_size, shape = 16, na.rm = TRUE)
  } else if (fillable) {
    ggtree::geom_tippoint(tip_aes, size = tip_size, colour = "black", stroke = stroke, na.rm = TRUE)
  } else {
    ggtree::geom_tippoint(tip_aes, size = tip_size, na.rm = TRUE)
  }

  # Tree scale bar: default lower-left (ggtree default); "topleft"/"midleft" place it
  # at the root (x = 0) near the top / halfway down the tip range, respectively.
  treescale_layer <- if (treescale_position %in% c("topleft", "midleft")) {
    n_tip <- length(tree$tip.label)
    yfrac <- if (identical(treescale_position, "midleft")) 0.5 else 0.92
    ggtree::geom_treescale(x = 0, y = n_tip * yfrac, width = treescale, fontsize = treescale_fontsize)
  } else {
    ggtree::geom_treescale(width = treescale, fontsize = treescale_fontsize)
  }
  p <- p +
    tip_layer +
    treescale_layer +
    ggplot2::ggtitle(title) +
    # Do not clip tip points/labels at the panel edge (tree tops were cut off).
    ggplot2::coord_cartesian(clip = "off")

  if (is.null(label_col)) {
    p <- p + if (fillable) ggplot2::scale_fill_identity() else ggplot2::scale_color_identity()
  } else {
    color_scale <- if (fillable) {
      ggplot2::scale_fill_manual(
        values = legend_colors, name = legend_title, breaks = names(legend_colors),
        limits = names(legend_colors), drop = FALSE, na.translate = FALSE
      )
    } else {
      ggplot2::scale_color_manual(
        values = legend_colors, name = legend_title, breaks = names(legend_colors),
        limits = names(legend_colors), drop = FALSE, na.translate = FALSE
      )
    }
    legend_guide <- if (fillable) {
      ggplot2::guides(fill = ggplot2::guide_legend(
        ncol = legend_ncol, byrow = TRUE, override.aes = list(size = tip_size, shape = 21)))
    } else {
      ggplot2::guides(color = ggplot2::guide_legend(
        ncol = legend_ncol, byrow = TRUE, override.aes = list(size = tip_size)))
    }
    p <- p + color_scale + legend_guide
  }
  if (!is.null(shape_col) && !is.null(shape_values)) {
    p <- p + ggplot2::scale_shape_manual(values = shape_values, drop = FALSE, na.translate = FALSE)
  }

  p +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) legend_position else "none",
      legend.box = "horizontal",
      plot.title = ggplot2::element_text(hjust = 0)
    )
}

.sirup_spatial_images <- function(object, images = NULL, sample_id = NULL) {
  if (.sirup_is_spatial_list(object)) {
    image_index <- .sirup_spatial_image_index(object)
    if (!is.null(sample_id)) {
      missing_samples <- setdiff(as.character(sample_id), image_index$sample)
      if (length(missing_samples) > 0L) {
        rlang::abort(paste0("Could not find spatial sample(s): ", paste(missing_samples, collapse = ", ")))
      }
      image_index <- image_index[image_index$sample %in% as.character(sample_id), , drop = FALSE]
    }
    if (is.null(images)) {
      return(image_index$label)
    }
    keep <- image_index$label %in% images |
      image_index$image %in% images |
      image_index$sample %in% images
    missing_images <- images[!vapply(images, function(image) {
      any(image_index$label == image | image_index$image == image | image_index$sample == image)
    }, logical(1))]
    if (length(missing_images) > 0L) {
      rlang::abort(paste0("Could not find image(s) or spatial sample(s): ", paste(missing_images, collapse = ", ")))
    }
    return(image_index$label[keep])
  }
  if (!is.null(sample_id)) {
    rlang::abort("`sample_id` can only select among samples in a multi-sample SIRUP object.")
  }
  object_images <- SeuratObject::Images(object)
  if (is.null(images)) {
    if (length(object_images) == 0L) {
      return("metadata")
    }
    return(object_images)
  }
  missing_images <- setdiff(images, c(object_images, "metadata"))
  if (length(missing_images) > 0L) {
    rlang::abort(paste0("Could not find image(s): ", paste(missing_images, collapse = ", ")))
  }
  images
}

.sirup_spatial_image_index <- function(object) {
  do.call(rbind, lapply(names(object), function(sample) {
    object_images <- SeuratObject::Images(object[[sample]])
    if (length(object_images) == 0L) {
      object_images <- "metadata"
    }
    data.frame(
      sample = sample,
      image = object_images,
      label = paste(sample, object_images, sep = "::"),
      stringsAsFactors = FALSE
    )
  }))
}

.sirup_parse_spatial_image_label <- function(image) {
  parts <- strsplit(image, "::", fixed = TRUE)[[1L]]
  if (length(parts) < 2L) {
    return(NULL)
  }
  list(sample = parts[[1L]], image = paste(parts[-1L], collapse = "::"))
}

#' Get the sample name from an image label
#'
#' `sirupSampleFromImage` returns the sample name in an image label.
#'
#' @details
#' A multi-sample SIRUP object labels each spatial image as
#' `"<sample>::<image>"`, so that image names stay unique across sections. The
#' function returns the sample part of such a label. For a label without a
#' sample part, it returns `NA`.
#'
#' @param image Character vector of image labels.
#'
#' @return A character vector of sample names, with `NA` for each label that has
#'   no sample part.
#'
#' @examples
#' sirupSampleFromImage(c("sample1::image1", "image2"))
#' @export
sirupSampleFromImage <- function(image) {
  vapply(as.character(image), function(label) {
    parsed <- .sirup_parse_spatial_image_label(label)
    if (is.null(parsed)) NA_character_ else parsed$sample
  }, character(1), USE.NAMES = FALSE)
}

.sirup_fetch_data <- function(object, vars, layer = "counts", assay = NULL) {
  if (!is.null(assay)) {
    old_assay <- .sirup_fetch_data_assay(object)
    SeuratObject::DefaultAssay(object) <- assay
    on.exit({
      if (!is.na(old_assay) && nzchar(old_assay)) {
        try(SeuratObject::DefaultAssay(object) <- old_assay, silent = TRUE)
      }
    }, add = TRUE)
  }
  assay <- .sirup_fetch_data_assay(object)
  out <- tryCatch(
    .sirup_fetch_data_layer(object = object, vars = vars, layer = layer),
    error = function(layer_error) {
      slot_out <- tryCatch(
        .sirup_fetch_data_slot(object = object, vars = vars, slot = layer),
        error = function(slot_error) {
          # Report the LAYER error as the parent, not the slot error. The slot
          # fall back only exists for SeuratObject < 5.0; from 5.0 on,
          # FetchData(slot =) is defunct and always fails the same way, so
          # attaching the slot error hid the real cause (for example a layer
          # that is genuinely absent from the assay) behind a generic
          # "deprecated slot argument" message.
          rlang::abort(
            paste0(
              "Could not fetch feature data from Seurat. SIRUP tried layer = ",
              layer, " in assay = ", assay, " first, then tried slot = ",
              layer, " as a fall back."
            ),
            parent = layer_error
          )
        }
      )
      slot_out
    }
  )
  as.data.frame(out, stringsAsFactors = FALSE)
}

.sirup_set_default_assay <- function(object, assay) {
  if (.sirup_is_spatial_list(object)) {
    for (sample in names(object)) {
      SeuratObject::DefaultAssay(object[[sample]]) <- assay
    }
    return(object)
  }
  SeuratObject::DefaultAssay(object) <- assay
  object
}

.sirup_fetch_data_assay <- function(object) {
  tryCatch(
    as.character(SeuratObject::DefaultAssay(object))[[1L]],
    error = function(e) NA_character_
  )
}

.sirup_fetch_data_layer <- function(object, vars, layer) {
  SeuratObject::FetchData(object, vars = vars, layer = layer)
}

.sirup_fetch_data_slot <- function(object, vars, slot) {
  fetch_formals <- .sirup_fetch_data_method_formals(object)
  if (!"slot" %in% fetch_formals) {
    rlang::abort("This SeuratObject `FetchData()` method does not support the legacy `slot` argument.")
  }
  SeuratObject::FetchData(object, vars = vars, slot = slot)
}

.sirup_fetch_data_method_formals <- function(object) {
  class_name <- class(object)[[1L]]
  method <- utils::getS3method("FetchData", class_name, optional = TRUE)
  if (is.null(method)) {
    method <- SeuratObject::FetchData
  }
  names(formals(method))
}

.sirup_spatial_data <- function(object, image, image_scale = NULL, swap_coords = NULL) {
  if (.sirup_is_spatial_list(object)) {
    parsed <- .sirup_parse_spatial_image_label(image)
    if (is.null(parsed) || !parsed$sample %in% names(object)) {
      rlang::abort("Multi-spatial image labels must use '<sample>::<image>'.")
    }
    options <- .sirup_require_spatial_options(
      image_scale,
      swap_coords,
      sample_id = parsed$sample
    )
    out <- .sirup_spatial_data(
      object = object[[parsed$sample]],
      image = parsed$image,
      image_scale = options$image_scale,
      swap_coords = options$swap_coords
    )
    out$coords$sample <- parsed$sample
    out$coords$image <- image
    return(out)
  }
  options <- .sirup_require_spatial_options(image_scale, swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  object_images <- SeuratObject::Images(object)
  if (identical(image, "metadata") && !"metadata" %in% object_images) {
    coords <- .sirup_metadata_coordinate_table(object)
    if (isTRUE(swap_coords)) {
      coords <- data.frame(
        cell = coords$cell,
        x = coords$y,
        y = coords$x,
        image = coords$image,
        stringsAsFactors = FALSE
      )
      coords$.sirup_y <- -coords$y
    }
    return(list(coords = coords, image = NULL))
  }
  image_obj <- object[[image]]
  if (isTRUE(swap_coords)) {
    coords <- .sirup_swapped_spatial_coordinates(image_obj, image_scale = image_scale)
  } else {
    coords <- .sirup_tissue_spatial_coordinates(image_obj, image_scale = image_scale)
  }
  if (is.null(coords)) {
    rlang::abort("Could not identify x/y coordinate columns in Seurat tissue coordinates.")
  }
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  coords <- coords[coords$cell %in% colnames(object), , drop = FALSE]
  coords$.sirup_y <- -coords$y

  list(coords = coords, image = image_obj)
}

.sirup_tissue_spatial_coordinates <- function(image, image_scale = NULL) {
  image_scale <- .sirup_resolve_image_scale(image_scale)
  if (is.null(image_scale)) {
    rlang::abort("`image_scale` must be supplied explicitly or stored in a SIRUP object.")
  }
  coords <- SeuratObject::GetTissueCoordinates(image, scale = image_scale)
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

.sirup_swapped_spatial_coordinates <- function(image, image_scale = NULL) {
  image_scale <- .sirup_resolve_image_scale(image_scale)
  if (is.null(image_scale)) {
    rlang::abort("`image_scale` must be supplied explicitly or stored in a SIRUP object.")
  }
  coords <- SeuratObject::GetTissueCoordinates(image, scale = image_scale)
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

.sirup_spatial_coordinate_candidates <- function(coords,
                                                image,
                                                image_scale = NULL,
                                                include_all_scales = FALSE) {
  image_scale <- .sirup_resolve_image_scale(image_scale)
  if (is.null(image_scale)) {
    rlang::abort("`image_scale` must be supplied explicitly or stored in a SIRUP object.")
  }
  raster <- .sirup_image_raster(image, alpha = 1)
  image_dim <- dim(raster)
  height <- image_dim[[1L]]
  width <- image_dim[[2L]]
  scale_factors <- .sirup_image_scale_factors(image)
  requested_scale_factor <- .sirup_image_scale_factor(image, image_scale)
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
      candidates[[length(candidates) + 1L]] <- .sirup_spatial_coordinate_candidate(
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

.sirup_spatial_coordinate_candidate <- function(x,
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

.sirup_image_scale_factor <- function(image, image_scale = NULL) {
  if (is.null(image_scale) || length(image_scale) == 0L || is.na(image_scale[[1L]])) {
    return(NA_real_)
  }
  scale_factors <- tryCatch(methods::slot(image, "scale.factors"), error = function(e) NULL)
  if (is.null(scale_factors)) {
    return(NA_real_)
  }
  scale_value <- tryCatch(scale_factors[[image_scale[[1L]]]], error = function(e) NULL)
  if (is.null(scale_value)) {
    scale_value <- tryCatch(scale_factors[[paste0(image_scale[[1L]], ".scale")]], error = function(e) NULL)
  }
  scale_value <- suppressWarnings(as.numeric(scale_value))
  if (length(scale_value) == 0L || !is.finite(scale_value[[1L]])) {
    return(NA_real_)
  }
  scale_value[[1L]]
}

.sirup_image_scale_factors <- function(image) {
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

.sirup_spatial_base <- function(plot_data,
                                image,
                                image_alpha = 1,
                                zoom = FALSE,
                                crop = TRUE,
                                view_data = NULL) {
  raster <- NULL
  if (!is.null(image)) {
    raster <- .sirup_image_raster(image, alpha = image_alpha)
  }
  xlim <- c(0, 1)
  ylim <- c(-1, 0)
  if (!is.null(raster)) {
    image_dim <- dim(raster)
    height <- image_dim[[1]]
    width <- image_dim[[2]]
    xlim <- c(0, width)
    ylim <- c(-height, 0)
  }
  # Choose the extent to display:
  #   zoom = TRUE  -> crop tightly to the plotted geometry (`plot_data`);
  #   crop = TRUE  -> crop to the spatial barcodes (`view_data` when supplied,
  #                   otherwise `plot_data`), matching Seurat's crop behaviour;
  #   crop = FALSE -> show the full tissue image.
  # An empty `view_data` is treated as absent so callers can fall back to
  # `plot_data` cleanly.
  if (!is.null(view_data) && nrow(view_data) == 0L) {
    view_data <- NULL
  }
  limit_data <- if (isTRUE(zoom)) {
    plot_data
  } else if (isTRUE(crop)) {
    view_data %||% plot_data
  } else {
    NULL
  }
  if (!is.null(limit_data) && nrow(limit_data) > 0L) {
    xlim <- .sirup_pad_range(limit_data$x, fallback = xlim)
    ylim <- .sirup_pad_range(limit_data$.sirup_y, fallback = ylim)
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x, y = .data$.sirup_y))
  if (!is.null(raster)) {
    p <- p + ggplot2::annotation_raster(raster, xmin = 0, xmax = width, ymin = -height, ymax = 0, interpolate = FALSE)
  }
  p +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "right"
    )
}

.sirup_image_raster <- function(image, alpha = 1) {
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

.sirup_pad_range <- function(x, fallback) {
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

# A title vector with one title for each plot sets the titles of the panels.
# One title stays a single title, as before.
.sirup_panel_titles <- function(title) {
  if (length(title) > 1L) as.character(title) else NULL
}

.sirup_set_panel_titles <- function(plots, panel_titles) {
  if (is.null(panel_titles)) {
    return(plots)
  }
  if (length(panel_titles) != length(plots)) {
    rlang::abort(paste0(
      "`title` must be one title, or one title for each plot. The function makes ",
      length(plots), " plot(s) and `title` has ", length(panel_titles), "."
    ))
  }
  for (i in seq_along(plots)) {
    plots[[i]] <- plots[[i]] + ggplot2::labs(title = panel_titles[[i]])
  }
  plots
}

.sirup_return_spatial_plots <- function(plots, combine = TRUE, combine_sample_plots = FALSE,
                                       panel_titles = NULL) {
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0L) {
    return(NULL)
  }
  plots <- .sirup_set_panel_titles(plots, panel_titles)
  combine <- .sirup_combine_spatial_plots(plots, combine, combine_sample_plots)
  if (!isTRUE(combine) || length(plots) == 1L) {
    plots <- .sirup_name_spatial_plot_list(plots)
    return(if (length(plots) == 1L && isTRUE(combine)) plots[[1L]] else plots)
  }
  patchwork::wrap_plots(plots, ncol = 1)
}

.sirup_combine_spatial_plots <- function(plots, combine = TRUE, combine_sample_plots = FALSE) {
  if (!isTRUE(combine)) {
    return(FALSE)
  }
  samples <- vapply(names(plots), function(name) {
    parsed <- .sirup_parse_spatial_image_label(name)
    if (is.null(parsed)) NA_character_ else parsed$sample
  }, character(1))
  samples <- unique(stats::na.omit(samples))
  if (length(samples) > 1L && !isTRUE(combine_sample_plots)) {
    return(FALSE)
  }
  TRUE
}

.sirup_name_spatial_plot_list <- function(plots) {
  plot_names <- names(plots)
  if (is.null(plot_names) || length(plot_names) == 0L) {
    return(plots)
  }
  parsed <- lapply(plot_names, .sirup_parse_spatial_image_label)
  samples <- vapply(parsed, function(x) {
    if (is.null(x)) NA_character_ else x$sample
  }, character(1))
  if (any(is.na(samples)) || anyDuplicated(samples)) {
    return(plots)
  }
  names(plots) <- samples
  plots
}

.sirup_plot_sample_id <- function(object, sample_id = NULL) {
  if (is.null(sample_id)) {
    return(NULL)
  }
  sample_id <- unique(as.character(sample_id))
  if (any(is.na(sample_id) | !nzchar(sample_id))) {
    rlang::abort("`sample_id` must contain non-missing spatial sample IDs.")
  }
  if (!isSIRUP(object)) {
    rlang::abort("`sample_id` can only be used with SIRUP objects.")
  }
  sample_names <- spatialSamples(object)
  missing_samples <- setdiff(sample_id, sample_names)
  if (length(missing_samples) > 0L) {
    rlang::abort(paste0("Could not find spatial sample(s): ", paste(missing_samples, collapse = ", ")))
  }
  if (length(sample_names) <= 1L) {
    return(NULL)
  }
  sample_id
}

.sirup_blend_color_grid <- function(bottom_left,
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

.sirup_blend_colors <- function(dat, color_grid) {
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

.sirup_blend_legend <- function(color_grid, xlab, ylab) {
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

# Reduce several receptor values on one spot to a single number.
.sirup_value_combine_fun <- function(value_combine) {
  switch(value_combine,
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE),
    max = function(x) max(x, na.rm = TRUE),
    min = function(x) min(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE),
    mode = .sirup_mode,
    rlang::abort(paste0("Unknown `value_combine`: ", value_combine))
  )
}

# Most frequent value. Ties take the smallest value, so the result is stable.
.sirup_mode <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  tab <- table(x)
  as.numeric(names(tab)[tab == max(tab)][[1]])
}

# Two colours give a low-to-high gradient; more give a multi-stop gradient.
# `aesthetic` is "colour" or "fill": plotReceptorsExpression needs the fill
# scale, because its expression layer already uses the colour scale.
.sirup_value_color_scale <- function(value_colors, value_limits, legend_name,
                                     aesthetic = c("colour", "fill")) {
  aesthetic <- match.arg(aesthetic)
  if (length(value_colors) < 2L) {
    rlang::abort("`value_colors` needs at least two colours.")
  }
  if (identical(aesthetic, "fill")) {
    if (length(value_colors) == 2L) {
      ggplot2::scale_fill_gradient(
        low = value_colors[[1]], high = value_colors[[2]],
        limits = value_limits, name = legend_name
      )
    } else {
      ggplot2::scale_fill_gradientn(
        colours = value_colors, limits = value_limits, name = legend_name
      )
    }
  } else if (length(value_colors) == 2L) {
    ggplot2::scale_color_gradient(
      low = value_colors[[1]], high = value_colors[[2]],
      limits = value_limits, name = legend_name
    )
  } else {
    ggplot2::scale_color_gradientn(
      colours = value_colors, limits = value_limits, name = legend_name
    )
  }
}

# One row per spot with the receptor value reduced to a single number. Shared by
# plotReceptors and plotReceptorsExpression.
.sirup_receptor_spot_values <- function(receptors, spatial_id, value, value_combine, sample_col = NULL) {
  if (length(value) != 1L || !is.character(value)) {
    rlang::abort("`value` must be a single column name, or \"count\".")
  }
  is_count <- identical(value, "count")
  if (is_count && "count" %in% names(receptors)) {
    rlang::warn('value = "count" counts receptors per spot; the `count` column is ignored.')
  }
  if (is_count) {
    group_cols <- unique(c(sample_col, spatial_id))
    out <- stats::aggregate(
      rep(1L, nrow(receptors)),
      by = receptors[, group_cols, drop = FALSE],
      FUN = sum
    )
    names(out)[ncol(out)] <- ".sirup_count"
    out$.sirup_count <- as.integer(out$.sirup_count)
    return(list(data = out, legend = "Receptors/spot"))
  }
  .check_required_cols(receptors, value, "receptors")
  if (!is.numeric(receptors[[value]])) {
    rlang::abort(paste0("`value` column \"", value, "\" must be numeric."))
  }
  keep <- !is.na(receptors[[value]])
  if (!any(keep)) {
    rlang::abort(paste0("`value` column \"", value, "\" has no non-missing values."))
  }
  receptors <- receptors[keep, , drop = FALSE]
  out <- stats::aggregate(
    receptors[[value]],
    by = receptors[, unique(c(sample_col, spatial_id)), drop = FALSE],
    FUN = .sirup_value_combine_fun(value_combine)
  )
  names(out)[ncol(out)] <- ".sirup_count"
  list(data = out, legend = value)
}

.sirup_merge_receptor_lookup <- function(coords, lookup, spatial_id, sample_col = NULL) {
  if (!is.null(sample_col) && sample_col %in% colnames(lookup) && "sample" %in% colnames(coords)) {
    return(merge(
      coords,
      lookup,
      by.x = c("sample", "cell"),
      by.y = c(sample_col, spatial_id),
      all = FALSE,
      sort = FALSE
    ))
  }
  merge(coords, lookup, by.x = "cell", by.y = spatial_id, all = FALSE, sort = FALSE)
}
