#' Build a spatial analysis table from a Seurat object
#'
#' @param object Seurat object with spatial images.
#' @param variable Feature or metadata column to analyze.
#' @param variable_type One of `feature`, `metadata`, or `receptor`.
#' @param value Optional metadata value to convert to a binary indicator.
#' @param receptors AIRR-style receptor data frame. Required when
#'   `variable_type = "receptor"`.
#' @param barcode_col Barcode column in `receptors`.
#' @param positive Convert feature values to a positive-expression indicator.
#' @param layer Seurat expression layer to fetch for features.
#' @param images Optional Seurat image names.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param zero_tol Values at or below this threshold are treated as zero.
#' @param as_sf Return an `sf` object with point geometry.
#' @return Data frame, or `sf` object, with spatial coordinates and `value`.
#' @export
spatialVariableData <- function(object,
                                variable,
                                variable_type = c("feature", "metadata", "receptor"),
                                value = NULL,
                                receptors = NULL,
                                barcode_col = "barcode_8um",
                                positive = FALSE,
                                layer = "counts",
                                images = NULL,
                                image.scale = "hires",
                                zero_tol = 0,
                                as_sf = FALSE) {
  variable_type <- match.arg(variable_type)
  .require_namespace("SeuratObject")

  dat <- .sirap_spatial_coordinate_table(object, images = images, image.scale = image.scale)

  if (variable_type == "feature") {
    expr <- .sirap_fetch_data(object, vars = variable, layer = layer)
    dat <- merge(dat, expr, by.x = "cell", by.y = 0, all.x = TRUE)
    dat$value <- dat[[variable]]
    if (isTRUE(positive)) {
      dat$value <- as.integer(!is.na(dat$value) & dat$value > zero_tol)
      dat$variable <- paste0(variable, "_positive")
    } else {
      dat$variable <- variable
    }
  } else if (variable_type == "metadata") {
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
  } else {
    if (is.null(receptors)) {
      rlang::abort("`receptors` is required when `variable_type = 'receptor'`.")
    }
    receptors <- .as_data_frame(receptors)
    .check_required_cols(receptors, barcode_col, "receptors")
    receptor_barcodes <- unique(as.character(receptors[[barcode_col]]))
    dat$value <- as.integer(dat$cell %in% receptor_barcodes)
    dat$variable <- variable
  }

  dat <- dat[!is.na(dat$value), , drop = FALSE]
  if (isTRUE(as_sf)) {
    .require_namespace("sf")
    return(sf::st_as_sf(dat, coords = c("x", ".sirap_y"), crs = NA, remove = FALSE))
  }
  dat
}

#' Calculate global Moran's I for a spatial variable
#'
#' @inheritParams spatialVariableData
#' @param radius Optional distance radius defining neighborhoods.
#' @param k Number of nearest neighbors when `radius` is `NULL`.
#' @param style Spatial weights style passed to `spdep::nb2listw()`.
#' @param alternative Alternative hypothesis passed to `spdep::moran.test()`.
#' @return List containing the Moran test, input data, neighbors, and weights.
#' @export
spatialMoran <- function(object,
                         variable,
                         variable_type = c("feature", "metadata", "receptor"),
                         value = NULL,
                         receptors = NULL,
                         barcode_col = "barcode_8um",
                         positive = FALSE,
                         layer = "counts",
                         images = NULL,
                         image.scale = "hires",
                         zero_tol = 0,
                         radius = NULL,
                         k = 6,
                         style = "W",
                         alternative = "greater") {
  .require_namespace("spdep")
  dat <- spatialVariableData(
    object = object,
    variable = variable,
    variable_type = variable_type,
    value = value,
    receptors = receptors,
    barcode_col = barcode_col,
    positive = positive,
    layer = layer,
    images = images,
    image.scale = image.scale,
    zero_tol = zero_tol
  )
  dat <- dat[is.finite(dat$value), , drop = FALSE]
  .check_spatial_values(dat$value)

  nb <- .sirap_spatial_neighbors(dat, radius = radius, k = k)
  listw <- spdep::nb2listw(nb, style = style, zero.policy = TRUE)
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
                              variable_type = c("feature", "metadata", "receptor"),
                              value = NULL,
                              receptors = NULL,
                              barcode_col = "barcode_8um",
                              positive = FALSE,
                              layer = "counts",
                              images = NULL,
                              image.scale = "hires",
                              zero_tol = 0,
                              radius = NULL,
                              k = 6,
                              style = "W") {
  .require_namespace("spdep")

  dat <- spatialVariableData(
    object = object,
    variable = variable,
    variable_type = variable_type,
    value = value,
    receptors = receptors,
    barcode_col = barcode_col,
    positive = positive,
    layer = layer,
    images = images,
    image.scale = image.scale,
    zero_tol = zero_tol
  )
  dat <- dat[is.finite(dat$value), , drop = FALSE]
  .check_spatial_values(dat$value)

  nb <- .sirap_spatial_neighbors(dat, radius = radius, k = k)
  listw <- spdep::nb2listw(nb, style = style, zero.policy = TRUE)
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

#' Plot local Moran's I over a tissue image
#'
#' @param object Seurat object with spatial images.
#' @param local_moran Output from `spatialLocalMoran()`.
#' @param value_col Column to color. Defaults to `z_i`.
#' @param images Optional Seurat image names.
#' @param image.scale Image scale passed to `SeuratObject::GetTissueCoordinates()`.
#' @param colors Three colors for low, midpoint, and high values.
#' @param midpoint Color midpoint.
#' @param point_size Point size.
#' @param point_alpha Point opacity.
#' @param order_points Plot lower values first so higher values appear on top.
#' @param image.alpha Tissue image opacity.
#' @param zoom Zoom to plotted points.
#' @param title Optional title.
#' @return A ggplot/patchwork object.
#' @export
plotSpatialLocalMoran <- function(object,
                                  local_moran,
                                  value_col = "z_i",
                                  images = NULL,
                                  image.scale = "hires",
                                  colors = c("#2c7bb6", "#f7f7f7", "#d7191c"),
                                  midpoint = 0,
                                  point_size = 0.7,
                                  point_alpha = 1,
                                  order_points = TRUE,
                                  image.alpha = 1,
                                  zoom = FALSE,
                                  title = NULL) {
  .check_required_cols(local_moran, c("cell", "x", ".sirap_y", value_col), "local_moran")
  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
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
      image.alpha = image.alpha,
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
#' @param object Seurat object with spatial images.
#' @param mark_col Metadata column containing point types.
#' @param i,j Mark levels to compare.
#' @param correction Edge correction to plot/use, usually `border`.
#' @param nsim Number of random-label simulations for an envelope. Use `0` for
#'   observed K only.
#' @param max_r Optional maximum radius to evaluate. Ignored if `r` is supplied.
#' @param r Optional distance vector passed to `spatstat.explore::Kcross()`.
#' @return List with observed K, optional envelope, plotting data, and point data.
#' @export
spatialRipleyCrossK <- function(object,
                                mark_col,
                                i,
                                j,
                                correction = "border",
                                nsim = 0,
                                max_r = NULL,
                                r = NULL) {
  .require_namespace("spatstat.geom")
  .require_namespace("spatstat.explore")
  if (nsim > 0L) {
    .require_namespace("spatstat.random")
  }

  dat <- .sirap_spatial_coordinate_table(object)
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
    k_obs <- spatstat.explore::Kcross(pp, i = i, j = j, correction = correction, r = r_eval)
    observed[[image]] <- k_obs

    if (nsim > 0L) {
      envelopes[[image]] <- spatstat.explore::envelope(
        pp,
        fun = spatstat.explore::Kcross,
        nsim = nsim,
        simulate = expression(spatstat.random::rlabel(pp)),
        i = i,
        j = j,
        correction = correction,
        r = r_eval
      )
      plot_data[[image]] <- .sirap_kcross_df(envelopes[[image]], correction = correction, image = image)
    } else {
      plot_data[[image]] <- .sirap_kcross_df(k_obs, correction = correction, image = image)
    }
  }

  list(
    observed = observed,
    envelope = envelopes,
    plot_data = dplyr::bind_rows(plot_data),
    points = dat,
    i = i,
    j = j,
    correction = correction,
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

.sirap_spatial_coordinate_table <- function(object, images = NULL, image.scale = "hires") {
  .require_namespace("SeuratObject")
  image_names <- tryCatch(
    SeuratObject::Images(object),
    error = function(e) character()
  )
  if (is.null(images) && length(image_names) == 0L) {
    return(.sirap_metadata_coordinate_table(object))
  }

  coords <- lapply(.sirap_spatial_images(object, images = images), function(image) {
    spatial <- .sirap_spatial_data(object, image = image, image.scale = image.scale)
    out <- spatial$coords
    out$image <- image
    out
  })
  dplyr::bind_rows(coords)
}

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

.check_spatial_values <- function(x) {
  if (length(unique(x)) < 2L) {
    rlang::abort("Spatial statistic requires at least two distinct values.")
  }
  invisible(TRUE)
}

.sirap_kcross_df <- function(x, correction, image) {
  dat <- as.data.frame(x)
  theoretical_col <- intersect(c("theo", "theoretical"), colnames(dat))
  simulated_col <- intersect(c("mmean", "mean", "simulated", "simmean", theoretical_col), colnames(dat))
  observed_col <- if (correction %in% colnames(dat)) {
    correction
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
