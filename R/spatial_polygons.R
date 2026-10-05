#' Calculate receptor clone dominance inside spatial polygons
#'
#' Finds Seurat spatial IDs whose coordinates fall inside polygons,
#' subsets receptor rows to those spatial IDs, and returns the fraction of regional
#' receptors represented by each clone. Points on the polygon boundary are
#' included. Open polygons are closed automatically by repeating the first
#' vertex at the end.
#'
#' @param receptors AIRR-style receptor data frame containing clone assignments
#'   and a spatial ID column. This may instead be a SIRAP object when
#'   `object` is omitted.
#' @param object Seurat object with spatial coordinates or metadata coordinate
#'   columns, or a SIRAP object. When this is a SIRAP object, `receptors` may
#'   be omitted.
#' @param polygon Data frame of polygon vertices in order. Must contain `x` and
#'   `y` columns by default. If a `polygon` column is present, each polygon is
#'   summarized independently unless `polygon_val` selects one or more polygons.
#' @param spatial_id Column in `receptors` matching Seurat spot/cell spatial IDs.
#'   If `NULL`, a spatial ID-like column is guessed.
#' @param clone Column in `receptors` containing clone IDs.
#' @param locus_val Optional locus value(s) used to filter `receptors`, for example
#'   `"IGH"`. Use `NULL` to keep all loci.
#' @param locus Column in `receptors` containing locus calls.
#' @param count Optional numeric count column. When `NULL`, each receptor
#'   row counts once.
#' @param v_call Column in `receptors` containing V-gene calls. When
#'   present, distinct allele-level V genes, V genes without alleles, and V-gene
#'   superfamilies are returned for each clone.
#' @param exclude_singletons Exclude clones represented by only one receptor row
#'   after locus filtering.
#' @param images Optional Seurat image names to use for spatial coordinates.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param polygon_val Optional polygon ID value(s) to select when `polygon_id`
#'   exists. When `NULL`, all polygons in `polygon_id` are analyzed.
#' @param x Coordinate column in `polygon`.
#' @param y Coordinate column in `polygon`.
#' @param polygon_id Optional polygon grouping column in `polygon`.
#' @param include_boundary Include Seurat spots on polygon edges when polygon
#'   membership must be calculated from `polygon`. If a SIRAP object with stored
#'   polygon annotations is supplied and `polygon = NULL`, stored annotations
#'   are used and this argument is ignored.
#' @return A spatial clone dominance object containing the dominance table,
#'   regional spatial IDs, regional coordinate rows, regional receptor rows, polygon
#'   vertices prepared for plotting, and analysis metadata containing the
#'   receptor field mappings used for the analysis.
#' @export
spatialCloneDominance <- function(receptors = NULL,
                                  object = NULL,
                                  polygon = NULL,
                                  spatial_id = NULL,
                                  clone = "clone_id",
                                  locus_val = NULL,
                                  locus = "locus",
                                  count = NULL,
                                  v_call = "v_call",
                                  exclude_singletons = FALSE,
                                  images = NULL,
                                  image_scale = NULL,
                                  swap_coords = NULL,
                                  polygon_val = NULL,
                                  x = "x",
                                  y = "y",
                                  polygon_id = "polygon",
                                  include_boundary = TRUE) {
  clone_missing <- missing(clone)
  locus_missing <- missing(locus)
  count_missing <- missing(count)
  if (isSIRAP(receptors) && is.null(object)) {
    object <- receptors
    receptors <- NULL
  }
  options <- .sirap_spatial_options(object, image_scale = image_scale, swap_coords = swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  is_sirap <- isSIRAP(object)
  sirap_container <- if (is_sirap) object else NULL
  use_stored_polygons <- is_sirap && is.null(polygon)
  if (is.null(polygon)) {
    if (!is_sirap) {
      rlang::abort("`polygon` must be supplied when `object` is not a SIRAP object.")
    }
    polygon <- .sirap_polygon_vertices_table(.sirap_polygons(object))
    if (nrow(polygon) == 0L) {
      rlang::abort("The SIRAP object does not contain stored polygons.")
    }
  }
  inputs <- .sirap_resolve_inputs(object, receptors, spatial_id)
  object <- inputs$object
  receptors <- inputs$receptors
  spatial_id <- inputs$spatial_id_col
  if (is.null(receptors)) {
    rlang::abort("`receptors` is required when `object` is not a SIRAP object.")
  }
  if (is_sirap && clone_missing && !is.null(inputs$receptor_fields$clone)) {
    clone <- inputs$receptor_fields$clone
  }
  if (is_sirap && locus_missing && !is.null(inputs$receptor_fields$locus)) {
    locus <- inputs$receptor_fields$locus
  }
  if (is_sirap && count_missing && !is.null(inputs$receptor_fields$count)) {
    count <- inputs$receptor_fields$count
  }
  receptors <- .as_data_frame(receptors)
  polygon <- .as_data_frame(polygon)
  spatial_id <- .guess_spatial_id_col(receptors, spatial_id)
  .check_required_cols(receptors, c(spatial_id, clone), "receptors")
  .check_required_cols(polygon, c(x, y), "polygon")
  if (!is.null(locus_val)) {
    .check_required_cols(receptors, locus, "receptors")
    missing_locus <- setdiff(as.character(locus_val), unique(as.character(receptors[[locus]])))
    if (length(missing_locus) == length(locus_val)) {
      rlang::warn(paste0(
        "No rows in `receptors` match `locus_val = ",
        paste(as.character(locus_val), collapse = ","),
        "` using `locus = \"",
        locus,
        "\"`."
      ))
    }
  }
  if (!is.null(count)) {
    .check_required_cols(receptors, count, "receptors")
  }
  receptors <- receptors[
    !is.na(receptors[[clone]]) &
      as.character(receptors[[clone]]) != "",
    ,
    drop = FALSE
  ]
  if (!is.null(locus_val)) {
    receptors <- receptors[receptors[[locus]] %in% locus_val, , drop = FALSE]
  }
  if (isTRUE(exclude_singletons)) {
    clone_ids <- as.character(receptors[[clone]])
    clone_counts <- table(clone_ids, useNA = "no")
    receptors <- receptors[clone_ids %in% names(clone_counts)[clone_counts > 1L], , drop = FALSE]
  }

  coords <- .sirap_spatial_coordinate_table(object, images = images, image_scale = image_scale, swap_coords = swap_coords)
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  polygon_specs <- .sirap_polygon_specs(polygon, polygon_id = polygon_val, polygon_col = polygon_id)
  polygon_plot_data <- .sirap_polygon_plot_data(
    polygon = polygon,
    polygon_id = polygon_val,
    polygon_col = polygon_id,
    x_col = x,
    y_col = y
  )

  results <- lapply(polygon_specs, function(spec) {
    polygon_value <- if (isTRUE(spec$add_polygon_col)) spec$id else "polygon_1"
    if (isTRUE(use_stored_polygons)) {
      region_coords <- .sirap_stored_polygon_coords(
        coords = coords,
        spatial = getSpatial(sirap_container),
        polygon_value = polygon_value
      )
    } else {
      vertices <- .sirap_select_polygon_vertices(
        polygon = polygon,
        polygon_id = spec$id,
        polygon_col = polygon_id,
        x_col = x,
        y_col = y
      )
      inside <- .sirap_points_in_polygon(
        px = coords$x,
        py = coords$y,
        vx = vertices$x,
        vy = vertices$y,
        include_boundary = include_boundary
      )
      region_coords <- coords[inside, , drop = FALSE]
    }
    region_spatial_ids <- unique(as.character(region_coords$cell))

    regional_receptors <- if (isTRUE(use_stored_polygons)) {
      .sirap_receptors_in_polygon(receptors, polygon_value, polygon_id = "polygon_id")
    } else {
      receptors[as.character(receptors[[spatial_id]]) %in% region_spatial_ids, , drop = FALSE]
    }

    if (nrow(regional_receptors) == 0L) {
      out <- data.frame(
        clone_id = character(),
        v_gene = character(),
        v_gene_family = character(),
        v_gene_superfamily = character(),
        n_receptors = numeric(),
        proportion = numeric(),
        region_receptor_total = numeric(),
        n_region_spatial_ids = integer(),
        stringsAsFactors = FALSE
      )
    } else {
      clone_ids <- as.character(regional_receptors[[clone]])
      weights <- if (is.null(count)) {
        rep(1, nrow(regional_receptors))
      } else {
        as.numeric(regional_receptors[[count]])
      }
      weights[is.na(weights)] <- 0
      out <- stats::aggregate(
        n_receptors ~ clone_id,
        data = data.frame(clone_id = clone_ids, n_receptors = weights, stringsAsFactors = FALSE),
        FUN = sum
      )
      out <- out[out$n_receptors > 0, , drop = FALSE]

      clone_ids <- unique(as.character(regional_receptors[[clone]]))
      if (v_call %in% colnames(regional_receptors)) {
        v_calls <- as.character(regional_receptors[[v_call]])
        names(v_calls) <- as.character(regional_receptors[[clone]])
        v_summary <- lapply(clone_ids, function(clone_id) {
          calls <- v_calls[names(v_calls) == clone_id]
          genes <- .sirap_v_genes(calls)
          families <- unique(stats::na.omit(.sirap_v_gene_family(genes)))
          families <- sort(families)
          superfamilies <- unique(stats::na.omit(.sirap_v_gene_superfamily(.sirap_v_gene_family(genes))))
          superfamilies <- sort(superfamilies)
          data.frame(
            clone_id = clone_id,
            v_gene = if (length(genes) == 0L) NA_character_ else paste(genes, collapse = ","),
            v_gene_family = if (length(families) == 0L) NA_character_ else paste(families, collapse = ","),
            v_gene_superfamily = if (length(superfamilies) == 0L) NA_character_ else paste(superfamilies, collapse = ","),
            stringsAsFactors = FALSE
          )
        })
        v_summary <- do.call(rbind, v_summary)
      } else {
        v_summary <- data.frame(
          clone_id = clone_ids,
          v_gene = rep(NA_character_, length(clone_ids)),
          v_gene_family = rep(NA_character_, length(clone_ids)),
          v_gene_superfamily = rep(NA_character_, length(clone_ids)),
          stringsAsFactors = FALSE
        )
      }

      out <- merge(out, v_summary, by = "clone_id", all.x = TRUE, sort = FALSE)
      out$v_gene[is.na(out$v_gene)] <- NA_character_
      out$v_gene_family[is.na(out$v_gene_family)] <- NA_character_
      out$v_gene_superfamily[is.na(out$v_gene_superfamily)] <- NA_character_
      total <- sum(out$n_receptors)
      out$proportion <- out$n_receptors / total
      out$region_receptor_total <- total
      out$n_region_spatial_ids <- length(region_spatial_ids)
      out <- out[order(out$proportion, out$n_receptors, decreasing = TRUE), , drop = FALSE]
      out <- out[, c("clone_id", "v_gene", "v_gene_family", "v_gene_superfamily", "n_receptors", "proportion", "region_receptor_total", "n_region_spatial_ids"), drop = FALSE]
    }
    if (!is.null(locus_val)) {
      out$locus_filter <- rep(paste(as.character(locus_val), collapse = ","), nrow(out))
    }
    out$polygon <- rep(as.character(polygon_value), nrow(out))

    region_spatial_id_table <- data.frame(
      spatial_id = region_spatial_ids,
      polygon = rep(as.character(polygon_value), length(region_spatial_ids)),
      stringsAsFactors = FALSE
    )
    region_coords$polygon <- rep(as.character(polygon_value), nrow(region_coords))
    if (nrow(regional_receptors) == 0L || nrow(region_coords) == 0L) {
      region_receptor_table <- regional_receptors[0, , drop = FALSE]
      region_receptor_table$.sirap_spatial_id <- character()
      region_receptor_table$x <- numeric()
      region_receptor_table$y <- numeric()
      region_receptor_table$.sirap_y <- numeric()
      region_receptor_table$image <- character()
      region_receptor_table$polygon <- character()
      region_receptor_table <- .sirap_add_v_gene_columns(region_receptor_table, v_call_col = v_call)
    } else {
      region_coord_table <- region_coords[, c("cell", "x", "y", "image"), drop = FALSE]
      region_coord_table$.sirap_y <- -region_coord_table$y
      region_receptor_table <- .sirap_drop_plot_coordinate_cols(regional_receptors)
      region_receptor_table$.sirap_spatial_id <- as.character(region_receptor_table[[spatial_id]])
      region_receptor_table <- merge(
        region_receptor_table,
        region_coord_table,
        by.x = ".sirap_spatial_id",
        by.y = "cell",
        all.x = FALSE,
        all.y = FALSE,
        sort = FALSE
      )
      region_receptor_table$polygon <- as.character(polygon_value)
      region_receptor_table <- .sirap_add_v_gene_columns(region_receptor_table, v_call_col = v_call)
    }
    list(
      dominance = out,
      region_spatial_ids = region_spatial_id_table,
      region_coordinates = region_coords,
      region_receptors = region_receptor_table
    )
  })

  out <- .sirap_rbind(results, "dominance")

  out <- list(
    dominance = out,
    region_spatial_ids = .sirap_rbind(results, "region_spatial_ids"),
    region_coordinates = .sirap_rbind(results, "region_coordinates"),
    region_receptors = .sirap_rbind(results, "region_receptors"),
    polygons = polygon_plot_data,
    metadata = list(
      receptor_fields = list(
        spatial_id = spatial_id,
        clone = clone,
        locus = locus,
        count = count,
        v_call = v_call
      )
    )
  )
  class(out) <- c("sirap_spatial_clone_dominance", class(out))
  out
}

# Fetch one feature for the polygon selector and say where it came from.
#
# SeuratObject::FetchData() already resolves object metadata as well as assay
# features, so a numeric module score (AddModuleScore, UCell, or any numeric
# per-cell column) does plot. Two things still needed fixing:
#
#   * the app reported no provenance, so a typo'd gene name that happened to
#     match a metadata column looked identical to a real assay hit;
#   * a FACTOR or CHARACTER metadata column came back as a data.frame and was
#     handed to a continuous colour scale, which errors inside the reactive
#     and leaves the app blank with no message.
#
# Returns a one-column numeric data.frame indexed by cell, with a "source"
# attribute of "assay" or "metadata", or NULL when the feature is missing or
# not numeric.
.sirap_polygon_feature_values <- function(object, feature, assay = NULL, layer = "counts") {
  if (is.null(feature) || !nzchar(feature)) {
    return(NULL)
  }
  meta_names <- tryCatch(colnames(object[[]]), error = function(e) character())
  source <- if (feature %in% meta_names) "metadata" else "assay"

  fetch_object <- object
  if (!is.null(assay) && nzchar(assay)) {
    SeuratObject::DefaultAssay(fetch_object) <- assay
  }
  expr <- tryCatch(
    .sirap_fetch_data(fetch_object, vars = feature, layer = layer),
    error = function(e) NULL
  )
  if (is.null(expr) || !feature %in% colnames(expr)) {
    return(NULL)
  }
  values <- expr[[feature]]
  if (!is.numeric(values) && !is.logical(values)) {
    return(NULL)
  }
  out <- data.frame(as.numeric(values), row.names = rownames(expr))
  names(out) <- feature
  attr(out, "source") <- source
  out
}

# Build a nearest-neighbour index ONCE for the polygon selector.
#
# The click handler used to call RANN::nn2() on the full coordinate table for
# every click. nn2 has no persistent tree, so each click paid an O(n log n)
# BUILD for a single k = 1 query -- at 180k-420k Visium HD bins that is the
# bulk of the lag between clicking and the app accepting the next click.
# The index below is built once when the app starts and queried per click.
#
# Deliberately not a grid lookup: Visium is hexagonal and single-cell data
# (Xenium, MERFISH) has no lattice at all, so nothing here may assume regular
# spacing.
#
# Backends, in order. Each is probed with a real query before being accepted,
# because an installed package is not the same as a working one -- nabor
# 0.5.0 built under R 4.5.1 constructs a WKNNF object whose Rcpp module
# exposes no methods, so $new() succeeds and $query() then fails.
#
#   nabor::WKNNF     persistent k-d tree, exact
#   BiocNeighbors    persistent KMKNN index, exact
#   vectorised scan  no build at all; O(n) per query but only a few
#                    arithmetic passes, which is competitive for single
#                    queries and always available
#
# Returns a list with `backend` and a `query(x, y)` closure giving the row
# index of the nearest point.
.sirap_nn_index <- function(xy) {
  probe <- function(fun) {
    ok <- tryCatch({
      hit <- fun(xy[1L, 1L], xy[1L, 2L])
      is.numeric(hit) && length(hit) == 1L && !is.na(hit)
    }, error = function(e) FALSE)
    isTRUE(ok)
  }

  if (requireNamespace("nabor", quietly = TRUE)) {
    tree <- tryCatch(nabor::WKNNF$new(xy), error = function(e) NULL)
    if (!is.null(tree)) {
      query <- function(x, y) {
        tree$query(matrix(c(x, y), nrow = 1L), k = 1L, eps = 0)$nn.idx[[1L]]
      }
      if (probe(query)) {
        return(list(backend = "nabor::WKNNF", query = query))
      }
    }
  }

  if (requireNamespace("BiocNeighbors", quietly = TRUE)) {
    built <- tryCatch(BiocNeighbors::buildIndex(xy), error = function(e) NULL)
    if (!is.null(built)) {
      query <- function(x, y) {
        BiocNeighbors::queryKNN(built, query = matrix(c(x, y), nrow = 1L),
                                k = 1L)$index[[1L]]
      }
      if (probe(query)) {
        return(list(backend = "BiocNeighbors", query = query))
      }
    }
  }

  # Columns are pulled out as plain vectors so the per-click path does no
  # data.frame subsetting or matrix copying.
  xs <- xy[, 1L]
  ys <- xy[, 2L]
  list(
    backend = "scan",
    query = function(x, y) which.min((xs - x)^2 + (ys - y)^2)
  )
}

#' Select spatial polygon vertices interactively
#'
#' Opens a Shiny app for clicking polygon vertices on a Seurat spatial
#' coordinate plot. Each click snaps to the nearest spatial ID. Use
#' "Close polygon" or "New polygon" before finishing to repeat the first vertex
#' at the end of each polygon.
#' The returned vertices may describe overlapping or nested polygons. Polygon
#' membership is resolved when these vertices are passed to
#' [assignSpatialPolygons()], [assignReceptorPolygons()], [createSIRAP()], or
#' [addPolygons()]: the primary `polygon_id` uses the smallest matching polygon
#' by area, while logical `polygon_<ID>` columns retain membership in each
#' individual polygon.
#'
#' @param object Seurat object with spatial coordinates, or a SIRAP object.
#' @param feature Optional feature shown as a color scale under the clicked
#'   vertices. The app also lets you change this interactively. Looked up in
#'   the assay first; if it is not there, the object metadata is checked, so
#'   module scores (`AddModuleScore()`, `UCell`, or any numeric per-cell
#'   column) can be displayed without adding them to an assay.
#' @param assay Optional assay to set before fetching feature values.
#' @param layer Expression layer/slot used for `feature`.
#' @param images Optional Seurat image names to show.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param save_path Optional CSV path for saved vertices.
#' @param plot_image Show the Seurat tissue image under the spatial points.
#'   Defaults to `FALSE`.
#' @param downsample_n Optional number of spatial points to show in the Shiny
#'   plot. Use `NULL` to show all points. Clicks still snap to the nearest
#'   spatial_id in the full coordinate table.
#' @param downsample_seed Random seed used when `downsample_n` is not `NULL`.
#' @param point_size Size of spatial points in the plot.
#' @param point_alpha Initial opacity for spatial points.
#' @param show_zeros Show points whose `feature` value is zero or `NA`.
#'   Defaults to `TRUE`. Set to `FALSE`, or untick "Show zero points" in the
#'   app, to plot only feature-positive points -- on sparse Visium HD bins the
#'   zeros can bury the signal. This affects the DISPLAY only; clicks still
#'   snap against the full coordinate table, so a hidden point can still be
#'   chosen as a vertex.
#' @param launch_browser Passed to `shiny::runApp()`.
#' @param host Optional host passed to `shiny::runApp()`. This is useful when
#'   running from an interactive HPC session through an SSH tunnel.
#' @param port Optional port passed to `shiny::runApp()`. This is useful when
#'   running from an interactive HPC session through an SSH tunnel.
#' @param quiet Suppress warnings while the app is running.
#' @return Data frame with `spatial_id`, `x`, `y`, `num`, `polygon`, and `image` columns.
#' @export
selectSpatialPolygons <- function(object,
                                  feature = NULL,
                                  assay = "Spatial",
                                  layer = "counts",
                                  images = NULL,
                                  image_scale = NULL,
                                  swap_coords = NULL,
                                  save_path = NULL,
                                  plot_image = FALSE,
                                  downsample_n = NULL,
                                  downsample_seed = 1L,
                                  point_size = 0.8,
                                  point_alpha = 0.9,
                                  show_zeros = TRUE,
                                  launch_browser = TRUE,
                                  host = getOption("shiny.host", "127.0.0.1"),
  port = getOption("shiny.port", NULL),
  quiet = TRUE) {
  context <- .sirap_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  .require_namespace("SeuratObject")
  .require_namespace("shiny")
  if (!"package:shiny" %in% search()) {
    suppressPackageStartupMessages(base::attachNamespace("shiny"))
  }
  old_shiny_useragg <- getOption("shiny.useragg")
  options(shiny.useragg = FALSE)
  on.exit(options(shiny.useragg = old_shiny_useragg), add = TRUE)
  assay_choices <- tryCatch(SeuratObject::Assays(object), error = function(e) character())
  selected_assay <- assay %||% tryCatch(SeuratObject::DefaultAssay(object), error = function(e) NULL)
  if (!is.null(selected_assay)) {
    SeuratObject::DefaultAssay(object) <- selected_assay
  }
  if (length(assay_choices) == 0L) {
    assay_choices <- selected_assay %||% ""
  }
  layer_choices <- .sirap_assay_layers(object, selected_assay)

  coords <- .sirap_spatial_coordinate_table(object, images = images, image_scale = image_scale, swap_coords = swap_coords)
  if (nrow(coords) == 0L) {
    rlang::abort("No spatial coordinates were found in `object`.")
  }
  coords$spatial_id <- coords$cell
  coords$.sirap_plot_y <- -coords$y
  y_axis_offset <- min(coords$.sirap_plot_y, na.rm = TRUE)
  plot_coords <- .sirap_downsample_polygon_plot_coords(
    coords = coords,
    downsample_n = downsample_n,
    seed = downsample_seed
  )
  image_names <- unique(stats::na.omit(as.character(coords$image)))
  image_names <- setdiff(image_names, "metadata")
  plot_image_info <- NULL
  if (length(image_names) > 0L) {
    if (length(image_names) > 1L) {
      rlang::warn("Multiple spatial images are present; showing the first image in `selectSpatialPolygons()`.")
    }
    raster <- .sirap_image_raster(object[[image_names[[1]]]], alpha = 1)
    image_dim <- dim(raster)
    plot_image_info <- list(
      raster = raster,
      width = image_dim[[2]],
      height = image_dim[[1]]
    )
  }
  # One index for the whole session. Built here, not in the click observer.
  nn_index <- .sirap_nn_index(cbind(coords$x, coords$.sirap_plot_y))
  rlang::inform(paste0("Polygon selector nearest-neighbour backend: ", nn_index$backend))

  x_range <- diff(range(plot_coords$x, na.rm = TRUE))
  y_range <- diff(range(plot_coords$.sirap_plot_y, na.rm = TRUE))
  plot_height <- if (!is.finite(x_range) || !is.finite(y_range) || x_range <= 0 || y_range <= 0) {
    500
  } else {
    max(500, min(1200, round(900 * y_range / x_range)))
  }
  default_feature <- feature
  if (is.null(default_feature)) {
    features <- tryCatch(rownames(object), error = function(e) character())
    default_feature <- if (length(features) == 0L) NULL else features[[1]]
  }
  save_default <- save_path
  if (is.null(save_default)) {
    save_default <- ""
  }

  ui <- shiny::fluidPage(
    shiny::titlePanel("Select Spatial Polygon Vertices"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput("assay", "Assay", choices = assay_choices, selected = selected_assay),
        shiny::selectInput("layer", "Layer", choices = layer_choices, selected = layer),
        shiny::textInput("feature", "Feature", value = default_feature %||% ""),
        shiny::textInput("polygon_name", "Polygon name", value = "1"),
        shiny::numericInput("point_size", "Point size", value = point_size, min = 0.01, max = 10),
        shiny::sliderInput("point_alpha", "Point opacity", min = 0, max = 1, value = point_alpha, step = 0.05, ticks = FALSE),
        shiny::checkboxInput("plot_image", "Show image", value = isTRUE(plot_image)),
        shiny::checkboxInput("show_zeros", "Show zero points", value = isTRUE(show_zeros)),
        shiny::textOutput("feature_source"),
        shiny::textInput("save_path", "CSV path", value = save_default),
        shiny::actionButton("close_polygon", "Close polygon"),
        shiny::actionButton("new_polygon", "New polygon"),
        shiny::actionButton("remove_last_point", "Remove last point"),
        shiny::actionButton("reset_points", "Reset"),
        shiny::actionButton("done", "Done"),
        shiny::tableOutput("selected_points"),
        width = 4
      ),
      shiny::mainPanel(
        shiny::plotOutput(
          "spatial_plot",
          width = "100%",
          height = paste0(plot_height, "px"),
          click = shiny::clickOpts("spatial_plot_click")
        ),
        width = 8
      )
    )
  )

  server <- function(input, output, session) {
    old_warn <- getOption("warn")
    if (isTRUE(quiet)) {
      options(warn = -1)
      on.exit(options(warn = old_warn), add = TRUE)
    }

    selected_points <- shiny::reactiveVal(.sirap_empty_polygon_vertices())
    current_polygon_num <- shiny::reactiveVal(1L)
    current_polygon_name <- shiny::reactiveVal("1")

    shiny::observeEvent(input$assay, {
      shiny::updateSelectInput(
        session,
        "layer",
        choices = .sirap_assay_layers(object, input$assay),
        selected = input$layer %||% layer
      )
    }, ignoreInit = TRUE)

    feature_source <- shiny::reactiveVal("")

    plot_data <- shiny::reactive({
      dat <- plot_coords
      feature_name <- input$feature
      if (is.null(feature_name) || !nzchar(feature_name)) {
        feature_source("")
        return(dat)
      }
      expr <- .sirap_polygon_feature_values(
        object = object,
        feature = feature_name,
        assay = input$assay,
        layer = input$layer %||% layer
      )
      if (is.null(expr)) {
        feature_source(paste0("'", feature_name, "' not found in assay or metadata"))
        return(dat)
      }
      feature_source(paste0("'", feature_name, "' from ", attr(expr, "source")))
      dat <- merge(dat, expr, by.x = "cell", by.y = 0, all.x = TRUE)
      if (!isTRUE(input$show_zeros)) {
        # display filter only -- `coords`, which clicks snap against, is untouched
        keep <- !is.na(dat[[feature_name]]) & dat[[feature_name]] > 0
        if (any(keep)) {
          dat <- dat[keep, , drop = FALSE]
        }
      }
      dat[order(dat[[feature_name]], na.last = TRUE), , drop = FALSE]
    })

    output$feature_source <- shiny::renderText(feature_source())

    output$spatial_plot <- shiny::renderPlot({
      dat <- plot_data()
      feature_name <- input$feature
      has_feature <- !is.null(feature_name) && nzchar(feature_name) && feature_name %in% colnames(dat)
      aes_args <- if (has_feature) {
        ggplot2::aes(x = .data$x, y = .data$.sirap_plot_y, color = .data[[feature_name]], key = .data$cell)
      } else {
        ggplot2::aes(x = .data$x, y = .data$.sirap_plot_y, key = .data$cell)
      }

      p <- ggplot2::ggplot(dat, aes_args)
      if (isTRUE(input$plot_image) && !is.null(plot_image_info)) {
        p <- p +
          ggplot2::annotation_raster(
            plot_image_info$raster,
            xmin = 0,
            xmax = plot_image_info$width,
            ymin = -plot_image_info$height,
            ymax = 0,
            interpolate = FALSE
          )
      }
      p <- p +
        ggplot2::geom_point(size = input$point_size, alpha = input$point_alpha) +
        ggplot2::coord_fixed() +
        ggplot2::scale_y_continuous(name = "y", labels = function(x) round(x - y_axis_offset)) +
        ggplot2::labs(x = "x") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
      if (has_feature) {
        p <- p + ggplot2::scale_color_gradientn(colors = grDevices::hcl.colors(64, "viridis"))
      }

      selected <- selected_points()
      if (nrow(selected) > 0L) {
        selected_path <- selected[
          stats::ave(selected$num, selected$num, FUN = length) > 1L,
          ,
          drop = FALSE
        ]
        if (nrow(selected_path) > 0L) {
          p <- p +
            ggplot2::geom_path(
              data = selected_path,
              ggplot2::aes(x = .data$x, y = .data$.sirap_plot_y, group = .data$num),
              inherit.aes = FALSE,
              linewidth = 0.5,
              color = "black"
            )
        }
        p <- p +
          ggplot2::geom_point(
            data = selected,
            ggplot2::aes(x = .data$x, y = .data$.sirap_plot_y, group = .data$num),
            inherit.aes = FALSE,
            size = input$point_size * 2,
            shape = 21,
            fill = "white",
            color = "black",
            stroke = 0.8
          )
      }
      p
    })

    shiny::observeEvent(input$spatial_plot_click, {
      click_data <- input$spatial_plot_click
      if (is.null(click_data)) {
        return()
      }
      closest_idx <- nn_index$query(click_data$x, click_data$y)
      closest <- coords[closest_idx, , drop = FALSE]
      closest <- data.frame(
        spatial_id = closest$cell,
        x = closest$x,
        y = closest$y,
        .sirap_plot_y = closest$.sirap_plot_y,
        num = current_polygon_num(),
        polygon = current_polygon_name(),
        image = closest$image,
        stringsAsFactors = FALSE
      )
      selected_points(rbind(selected_points(), closest))
    })

    shiny::observeEvent(input$polygon_name, {
      polygon_name <- trimws(input$polygon_name %||% "")
      if (!nzchar(polygon_name)) {
        shiny::showNotification("Polygon name cannot be empty.", type = "error")
        return()
      }
      current <- selected_points()
      used_elsewhere <- nrow(current) > 0L &&
        any(current$num != current_polygon_num() & current$polygon == polygon_name)
      if (isTRUE(used_elsewhere)) {
        shiny::showNotification(
          "That polygon name is already used. Please choose a different name.",
          type = "error"
        )
        return()
      }
      current_polygon_name(polygon_name)
      if (nrow(current) > 0L) {
        current$polygon[current$num == current_polygon_num()] <- polygon_name
        selected_points(current)
      }
    }, ignoreInit = TRUE)

    output$selected_points <- shiny::renderTable({
      current <- selected_points()
      current <- current[current$num == current_polygon_num(), , drop = FALSE]
      current <- .sirap_public_polygon_vertices(current)
      current[, setdiff(colnames(current), "image"), drop = FALSE]
    })

    shiny::observeEvent(input$close_polygon, {
      selected_points(.sirap_close_current_polygon(selected_points(), current_polygon_num()))
    })

    shiny::observeEvent(input$new_polygon, {
      current <- selected_points()
      if (!any(current$num == current_polygon_num())) {
        return()
      }
      current <- .sirap_close_current_polygon(current, current_polygon_num())
      selected_points(current)
      next_num <- max(current$num, na.rm = TRUE) + 1L
      next_name <- as.character(next_num)
      current_polygon_num(next_num)
      current_polygon_name(next_name)
      shiny::updateTextInput(session, "polygon_name", value = next_name)
    })

    shiny::observeEvent(input$remove_last_point, {
      current <- selected_points()
      if (nrow(current) > 0L) {
        remove_idx <- utils::tail(which(current$num == current_polygon_num()), 1L)
        if (length(remove_idx) > 0L) {
          selected_points(current[-remove_idx, , drop = FALSE])
        }
      }
    })

    shiny::observeEvent(input$reset_points, {
      selected_points(.sirap_empty_polygon_vertices())
      current_polygon_num(1L)
      current_polygon_name("1")
      shiny::updateTextInput(session, "polygon_name", value = "1")
    })

    shiny::observeEvent(input$done, {
      out <- selected_points()
      for (polygon_num in unique(out$num)) {
        out <- .sirap_close_current_polygon(out, polygon_num)
      }
      out <- .sirap_public_polygon_vertices(out)
      path <- input$save_path
      if (!is.null(path) && nzchar(path)) {
        utils::write.csv(out, path, row.names = FALSE)
        shiny::showNotification(paste("Saved polygons to", path), type = "message")
      }
      shiny::stopApp(out)
    })
  }

  shiny::runApp(
    shiny::shinyApp(ui, server),
    launch.browser = launch_browser,
    host = host,
    port = port
  )
}

#' Format spatial polygon point sets as polygon vertices
#'
#' Converts a data frame containing all points inside one or more polygons into
#' the same vertex-table structure returned by `selectSpatialPolygons()`. By
#' default, interior point sets are converted to polygon outlines by buffering
#' and unioning spots with `sf`, then mapping the exterior ring back to the
#' nearest selected spatial_id coordinates. This preserves a valid boundary order
#' for concave polygons while returning Visium spatial_id coordinates as vertices.
#' If the input already contains boundary points, use `point_type = "boundary"`.
#'
#' @param object Seurat object with spatial coordinates, or a SIRAP object.
#' @param polygon_points Data frame containing all points inside each polygon,
#'   or boundary points when `point_type = "boundary"`.
#' @param spatial_id Column in `polygon_points` containing Seurat spatial IDs.
#'   Required when `x` and `y` are not present.
#' @param polygon Column in `polygon_points` identifying each polygon.
#' @param x X-coordinate column in `polygon_points`, if already present.
#' @param y Y-coordinate column in `polygon_points`, if already present.
#' @param image Image column in `polygon_points`, if already present.
#' @param point_type Are `polygon_points` all interior points or already
#'   boundary points?
#' @param buffer_radius Optional radius used to buffer points when
#'   converting interior points to a polygon outline. When `NULL`, the radius is
#'   estimated from the nearest-neighbor spacing in the full spatial coordinate
#'   grid.
#' @param order Optional column used to order points within each polygon
#'   before formatting. This is mainly useful with `point_type = "boundary"`.
#' @param images Optional Seurat image names to use for spatial coordinates.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate axes.
#' @param close_polygons Close each polygon by repeating the first vertex when
#'   needed.
#' @return Data frame with `spatial_id`, `x`, `y`, `polygon`, and `image` columns.
#'   Polygons with fewer than three points are excluded with a warning.
#' @export
formatSpatialPolygons <- function(object,
                                  polygon_points,
                                  spatial_id = "spatial_id",
                                  polygon = "polygon",
                                  x = "x",
                                  y = "y",
                                  image = "image",
                                  point_type = c("interior", "boundary"),
                                  buffer_radius = NULL,
                                  order = NULL,
                                  images = NULL,
                                  image_scale = NULL,
                                  swap_coords = NULL,
                                  close_polygons = TRUE) {
  context <- .sirap_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  .require_namespace("SeuratObject")
  point_type <- match.arg(point_type)
  polygon_points <- .as_data_frame(polygon_points)
  .check_required_cols(polygon_points, polygon, "polygon_points")
  if (nrow(polygon_points) == 0L) {
    rlang::abort("`polygon_points` must contain at least one row.")
  }
  if (any(is.na(polygon_points[[polygon]]))) {
    rlang::abort("`polygon_points` contains missing polygon IDs.")
  }
  if (!is.null(order)) {
    .check_required_cols(polygon_points, order, "polygon_points")
    polygon_points$.sirap_input_order <- seq_len(nrow(polygon_points))
    ordered_rows <- unlist(lapply(unique(polygon_points[[polygon]]), function(polygon_id) {
      idx <- which(polygon_points[[polygon]] == polygon_id)
      idx[order(polygon_points[[order]][idx], polygon_points$.sirap_input_order[idx])]
    }), use.names = FALSE)
    polygon_points <- polygon_points[ordered_rows, , drop = FALSE]
    polygon_points$.sirap_input_order <- NULL
  }

  has_spatial_id <- spatial_id %in% colnames(polygon_points)
  has_xy <- all(c(x, y) %in% colnames(polygon_points))
  if (!has_spatial_id && !has_xy) {
    rlang::abort("`polygon_points` must contain spatial_id values or x/y coordinate columns.")
  }

  coords <- .sirap_spatial_coordinate_table(
    object,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  coords$cell <- as.character(coords$cell)
  if (has_spatial_id) {
    polygon_points[[spatial_id]] <- as.character(polygon_points[[spatial_id]])
  }

  if (has_spatial_id && (!has_xy || !image %in% colnames(polygon_points))) {
    has_image <- image %in% colnames(polygon_points)
    if (anyDuplicated(coords$cell) > 0L && !has_image) {
      rlang::abort(
        "`object` has duplicated spatial IDs across images; provide `image` in `polygon_points`."
      )
    }
    if (has_image) {
      coord_key <- paste(coords$cell, coords$image, sep = "\r")
      point_key <- paste(polygon_points[[spatial_id]], polygon_points[[image]], sep = "\r")
      coord_lookup <- coords[match(point_key, coord_key), , drop = FALSE]
    } else {
      coord_lookup <- coords[match(polygon_points[[spatial_id]], coords$cell), , drop = FALSE]
    }
    missing_spatial_ids <- unique(polygon_points[[spatial_id]][is.na(coord_lookup$cell)])
    if (length(missing_spatial_ids) > 0L) {
      rlang::abort(paste0(
        "Some polygon point spatial IDs were not found in `object`: ",
        paste(utils::head(missing_spatial_ids, 10), collapse = ", ")
      ))
    }
    if (!has_xy) {
      polygon_points[[x]] <- coord_lookup$x
      polygon_points[[y]] <- coord_lookup$y
    }
    if (!image %in% colnames(polygon_points)) {
      polygon_points[[image]] <- coord_lookup$image
    }
  }

  if (!image %in% colnames(polygon_points)) {
    image_values <- unique(as.character(coords$image))
    if (length(image_values) != 1L) {
      rlang::abort("`polygon_points` must contain `image` when `object` has multiple spatial images.")
    }
    polygon_points[[image]] <- image_values[[1L]]
  }

  out <- data.frame(
    spatial_id = if (has_spatial_id) polygon_points[[spatial_id]] else NA_character_,
    x = as.numeric(polygon_points[[x]]),
    y = as.numeric(polygon_points[[y]]),
    polygon = polygon_points[[polygon]],
    image = as.character(polygon_points[[image]]),
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(out$x) | !is.finite(out$y))) {
    rlang::abort("`polygon_points` contains non-finite polygon coordinates.")
  }

  multi_image_polygons <- names(which(tapply(out$image, out$polygon, function(x) length(unique(x)) > 1L)))
  if (length(multi_image_polygons) > 0L) {
    rlang::abort(paste0(
      "Each polygon ID must belong to one image. These polygons span multiple images: ",
      paste(multi_image_polygons, collapse = ", ")
    ))
  }

  if (point_type == "interior") {
    .require_namespace("sf")
    if (is.null(buffer_radius)) {
      buffer_radius <- .sirap_default_polygon_buffer_radius(coords)
    }
    out <- .sirap_polygon_interior_to_sf_boundary(
      points = out,
      buffer_radius = buffer_radius
    )
  }

  if (isTRUE(close_polygons)) {
    out <- .sirap_close_formatted_polygons(out)
  }
  rownames(out) <- NULL
  out
}

#' Expand polygons outward from their centroids
#'
#' Grows each polygon by moving its vertices away from the polygon centroid.
#' Each vertex moves along the line from the centroid through the vertex, by a
#' fixed distance `delta`. The polygon keeps its shape and vertex count. The
#' polygons are expected to be convex.
#'
#' @param polygons Polygons to expand. Either a flat data frame of vertices
#'   (with `x`, `y`, and an optional polygon ID column), or a SIRAP object
#'   whose stored polygons are expanded.
#' @param delta Distance each vertex moves outward, in coordinate units. This
#'   is required. A negative value moves the vertices inward and shrinks the
#'   polygon.
#' @param x,y Vertex coordinate columns when `polygons` is a data frame.
#' @param polygon_id Column that separates polygons when `polygons` is a data
#'   frame. When the named column is absent, a `polygon` or `polygon_id` column
#'   is used if present; otherwise all rows form one polygon.
#' @return The expanded polygons in the same form as the input: a data frame
#'   for a data-frame input, or an updated SIRAP object for a SIRAP object
#'   input.
#' @export
expandPolygons <- function(polygons,
                           delta = NULL,
                           x = "x",
                           y = "y",
                           polygon_id = "polygon") {
  if (is.null(delta)) {
    rlang::abort("`delta` must be supplied (a single finite number).")
  }
  if (isSIRAP(polygons)) {
    object <- polygons
    stored <- .sirap_polygons(object)
    if (nrow(stored) == 0L) {
      rlang::abort("The SIRAP object does not contain stored polygons.")
    }
    flat <- .sirap_polygon_vertices_table(stored)
    # The stored vertices are moved coordinates, not spatial barcodes, so drop
    # the carried spatial IDs before re-storing the expanded polygons.
    flat$spatial_id <- NULL
    expanded <- .sirap_expand_polygon_df(flat, delta = delta, x = "x", y = "y", polygon_id = "polygon")
    return(addPolygons(
      object,
      polygons = expanded,
      include_boundary = isTRUE(object$metadata$polygons_include_boundary %||% TRUE)
    ))
  }
  polygons <- .as_data_frame(polygons)
  .sirap_expand_polygon_df(polygons, delta = delta, x = x, y = y, polygon_id = polygon_id)
}

# Expand every polygon in a flat vertex table. Only the `x`/`y` columns change;
# all other columns and the row order are kept. See `expandPolygons()`.
.sirap_expand_polygon_df <- function(polygon, delta, x = "x", y = "y", polygon_id = "polygon") {
  polygon <- .as_data_frame(polygon)
  .check_required_cols(polygon, c(x, y), "polygons")
  if (nrow(polygon) == 0L) {
    rlang::abort("`polygons` must contain at least one vertex.")
  }
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta)) {
    rlang::abort("`delta` must be a single finite number.")
  }
  # Resolve the polygon-id column. Fall back to the common names, or treat every
  # row as one polygon when no id column is present.
  if (!polygon_id %in% names(polygon)) {
    alt <- intersect(c("polygon", "polygon_id"), names(polygon))
    polygon_id <- if (length(alt) > 0L) alt[[1L]] else NA_character_
  }
  groups <- if (is.na(polygon_id)) {
    list(seq_len(nrow(polygon)))
  } else {
    ids <- as.character(polygon[[polygon_id]])
    unname(split(seq_len(nrow(polygon)), factor(ids, levels = unique(ids))))
  }
  out <- polygon
  out[[x]] <- as.numeric(out[[x]])
  out[[y]] <- as.numeric(out[[y]])
  for (idx in groups) {
    px <- as.numeric(polygon[[x]][idx])
    py <- as.numeric(polygon[[y]][idx])
    if (any(!is.finite(px) | !is.finite(py))) {
      rlang::abort("`polygons` contains non-finite coordinates.")
    }
    moved <- .sirap_expand_polygon_vertices(px, py, delta = delta)
    out[[x]][idx] <- moved$x
    out[[y]][idx] <- moved$y
  }
  out
}

# Move one polygon's vertices outward from its centroid by `delta`. A repeated
# closing vertex equals the first vertex, so it moves identically and the
# polygon stays closed. Returns a list with the moved `x`/`y` vectors.
.sirap_expand_polygon_vertices <- function(px, py, delta) {
  n_all <- length(px)
  closed <- n_all > 1L &&
    .sirap_same_xy(
      data.frame(x = px[[1L]], y = py[[1L]]),
      data.frame(x = px[[n_all]], y = py[[n_all]])
    )
  core <- if (closed) seq_len(n_all - 1L) else seq_len(n_all)
  if (length(core) < 3L) {
    rlang::abort("Each polygon must contain at least three vertices.")
  }
  centroid <- .sirap_polygon_centroid(px[core], py[core])
  vx <- px - centroid$x
  vy <- py - centroid$y
  dist <- sqrt(vx^2 + vy^2)
  # A vertex exactly at the centroid has no direction, so leave it in place.
  scale <- ifelse(dist > 0, delta / dist, 0)
  list(x = px + vx * scale, y = py + vy * scale)
}

# Area (geometric) centroid of a polygon, from its ordered vertices (without a
# repeated closing vertex). Falls back to the vertex mean for a zero-area
# (degenerate) polygon.
.sirap_polygon_centroid <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  n <- length(x)
  nxt <- c(seq_len(n)[-1L], 1L)
  cross <- x * y[nxt] - x[nxt] * y
  area2 <- sum(cross)
  if (abs(area2) < .Machine$double.eps^0.5) {
    return(list(x = mean(x), y = mean(y)))
  }
  list(
    x = sum((x + x[nxt]) * cross) / (3 * area2),
    y = sum((y + y[nxt]) * cross) / (3 * area2)
  )
}

#' Plot spatial polygons on a tissue image
#'
#' Plots one or more polygon vertex tables, such as the output of
#' `selectSpatialPolygons()`, over the Seurat tissue image.
#'
#' @param object Seurat object with spatial images, or a SIRAP object.
#' @param polygon Data frame of polygon vertices in order. Must contain `x` and
#'   `y` columns by default. If `object` is a SIRAP object and `polygon` is
#'   omitted, the object's stored polygons are used.
#' @param images Optional Seurat image names to show.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param polygon_val Optional polygon ID value(s) to select when `polygon_id`
#'   exists.
#' @param x Coordinate column in `polygon`.
#' @param y Coordinate column in `polygon`.
#' @param polygon_id Optional polygon grouping column in `polygon`.
#' @param polygon_alpha Fill opacity for polygons.
#' @param polygon_fill Fill color for polygons when no analysis is mapped.
#' @param fill_analysis Optional name of a stored SIRAP analysis containing
#'   one row per polygon.
#' @param fill_column Optional column used to color polygon interiors. When
#'   `fill_analysis` is `NULL`, this column is read from the polygon data. When
#'   `fill_analysis` is supplied, this column is read from that analysis, which
#'   must contain `polygon_id` or `polygon` as its polygon key.
#' @param analysis_fill Character vector of two colors used as the low and high
#'   ends of the numeric analysis gradient. Defaults to
#'   `c("#D3D3D3", "#FFFF00")`. Non-numeric analysis values use the Glasbey
#'   palette.
#' @param polygon_color Outline color for polygons.
#' @param linewidth Outline width for polygons.
#' @param image_alpha Tissue image opacity.
#' @param label Plot labels at polygon centers.
#' @param label_by Column in the polygon data to use for labels.
#' @param label_size Text size for polygon labels.
#' @param label_color Text color for polygon labels.
#' @param label_outside Place each label OUTSIDE its polygon instead of at the
#'   center. The label is pushed radially away from the centroid of all plotted
#'   polygons (so labels fan outward from a cluster of small polygons, e.g.
#'   crowded thymic eGCs), clearing the polygon by `label_offset`. A thin leader
#'   line links each label to its polygon unless `label_leader = FALSE`.
#' @param label_offset When `label_outside = TRUE`, the gap between the polygon
#'   edge and its label, as a fraction of the plotted coordinate span
#'   (default 0.08).
#' @param label_leader When `label_outside = TRUE`, draw a leader segment from
#'   each polygon to its displaced label.
#' @param label_leader_color Color of the leader segments.
#' @param zoom Zoom to the polygon extents instead of the full tissue.
#' @param combine Combine plots with patchwork when possible.
#' @return A ggplot/patchwork object or a named list of ggplots when
#'   `combine = FALSE`.
#' @export
plotSpatialPolygons <- function(object,
                                polygon = NULL,
                                images = NULL,
                                image_scale = NULL,
                                swap_coords = NULL,
                                polygon_val = NULL,
                                x = "x",
                                y = "y",
                                polygon_id = "polygon",
                                polygon_alpha = 0.2,
                                polygon_fill = "#0072B2",
                                fill_analysis = NULL,
                                fill_column = NULL,
                                analysis_fill = c("#D3D3D3", "#FFFF00"),
                                polygon_color = "#0072B2",
                                linewidth = 0.6,
                                image_alpha = 1,
                                label = TRUE,
                                label_by = "polygon",
                                label_size = 3,
                                label_color = "black",
                                label_outside = FALSE,
                                label_offset = 0.08,
                                label_leader = TRUE,
                                label_leader_color = "grey40",
                                zoom = FALSE,
                                combine = TRUE) {
  sirap_object <- if (isSIRAP(object)) object else NULL
  context <- .sirap_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  .require_namespace("SeuratObject")
  if (is.null(polygon)) {
    if (is.null(sirap_object)) {
      rlang::abort("`polygon` must be supplied when `object` is not a SIRAP object.")
    }
    polygon <- .sirap_polygon_vertices_table(.sirap_polygons(sirap_object))
    if (nrow(polygon) == 0L) {
      rlang::abort("The SIRAP object does not contain stored polygons.")
    }
  }
  polygon <- .as_data_frame(polygon)
  .check_required_cols(polygon, c(x, y), "polygon")
  if (!is.numeric(polygon_alpha) || length(polygon_alpha) != 1L || is.na(polygon_alpha) ||
      polygon_alpha < 0 || polygon_alpha > 1) {
    rlang::abort("`polygon_alpha` must be a number between 0 and 1.")
  }

  polygon_keep_cols <- label_by
  if (!is.null(fill_column) && is.null(fill_analysis)) {
    polygon_keep_cols <- unique(c(polygon_keep_cols, fill_column))
  }
  polygon_data <- .sirap_polygon_plot_data(
    polygon = polygon,
    polygon_id = polygon_val,
    polygon_col = polygon_id,
    x_col = x,
    y_col = y,
    keep_cols = polygon_keep_cols
  )
  if (!is.null(fill_analysis) || !is.null(fill_column)) {
    if (is.null(fill_column)) {
      rlang::abort("Supply `fill_column` when `fill_analysis` is supplied.")
    }
    if (is.null(fill_analysis)) {
      polygon_data <- .sirap_join_polygon_column_fill_data(
        polygon_data = polygon_data,
        fill_column = fill_column
      )
    } else {
      if (is.null(sirap_object)) {
        rlang::abort("`fill_analysis` requires a SIRAP object.")
      }
      polygon_data <- .sirap_join_polygon_fill_data(
        object = sirap_object,
        polygon_data = polygon_data,
        analysis_name = fill_analysis,
        fill_column = fill_column
      )
    }
  }
  if (!isTRUE(label) && !identical(label, FALSE)) {
    rlang::abort("`label` must be TRUE or FALSE.")
  }
  if (!isTRUE(label_outside) && !identical(label_outside, FALSE)) {
    rlang::abort("`label_outside` must be TRUE or FALSE.")
  }
  if (!isTRUE(label_leader) && !identical(label_leader, FALSE)) {
    rlang::abort("`label_leader` must be TRUE or FALSE.")
  }
  if (isTRUE(label)) {
    .check_required_cols(polygon_data, label_by, "polygon")
    label_idx <- !duplicated(polygon_data$polygon)
    label_lookup <- data.frame(
      polygon = polygon_data$polygon[label_idx],
      .sirap_label = as.character(polygon_data[[label_by]][label_idx]),
      stringsAsFactors = FALSE
    )
    label_data <- .sirap_polygon_centers(polygon_data, unique(polygon_data$polygon))
    label_data <- merge(label_data, label_lookup, by = "polygon", all.x = TRUE, sort = FALSE)
    if (isTRUE(label_outside)) {
      if (!is.numeric(label_offset) || length(label_offset) != 1L || is.na(label_offset) || label_offset < 0) {
        rlang::abort("`label_offset` must be a non-negative number.")
      }
      # Push each label outside its polygon, away from the centroid of all
      # polygons, then DODGE: search increasing distances and a few angle offsets
      # for the nearest spot that clears EVERY polygon (not just its own), so a
      # displaced label does not land on a neighbouring polygon.
      span <- max(
        diff(range(polygon_data$x, na.rm = TRUE)),
        diff(range(polygon_data$.sirap_y, na.rm = TRUE)),
        na.rm = TRUE
      )
      gx <- mean(label_data$x, na.rm = TRUE)
      gy <- mean(label_data$.sirap_y, na.rm = TRUE)
      poly_ids <- unique(polygon_data$polygon)
      poly_verts <- lapply(poly_ids, function(pg) polygon_data[polygon_data$polygon == pg, , drop = FALSE])
      # Clear the label FOOTPRINT, not just its anchor: test the anchor plus a
      # ring of points at `clearance` so a rendered digit cannot overlap a
      # neighbouring polygon even when its anchor sits just outside the edge.
      clearance <- 0.03 * span
      ring <- seq(0, 2 * pi, length.out = 13)[-13]
      in_any_poly <- function(px, py) {
        txs <- c(px, px + clearance * cos(ring))
        tys <- c(py, py + clearance * sin(ring))
        for (v in poly_verts) {
          if (nrow(v) < 4L) next
          hit <- tryCatch(
            .sirap_points_in_polygon(txs, tys, v$x, v$.sirap_y, include_boundary = TRUE),
            error = function(e) rep(FALSE, length(txs))
          )
          if (any(hit)) return(TRUE)
        }
        FALSE
      }
      rot <- function(v, deg) {
        th <- deg * pi / 180
        c(v[1] * cos(th) - v[2] * sin(th), v[1] * sin(th) + v[2] * cos(th))
      }
      angles <- c(0, 20, -20, 40, -40, 60, -60, 90, -90)
      dists <- (label_offset * span) * c(1, 1.6, 2.4, 3.4, 4.6, 6)
      lx <- numeric(nrow(label_data))
      ly <- numeric(nrow(label_data))
      ex <- numeric(nrow(label_data))
      ey <- numeric(nrow(label_data))
      for (i in seq_len(nrow(label_data))) {
        cx <- label_data$x[i]
        cy <- label_data$.sirap_y[i]
        v <- polygon_data[polygon_data$polygon == label_data$polygon[i], , drop = FALSE]
        radius <- max(sqrt((v$x - cx)^2 + (v$.sirap_y - cy)^2), na.rm = TRUE)
        dir0 <- c(cx - gx, cy - gy)
        n <- sqrt(sum(dir0^2))
        dir0 <- if (is.finite(n) && n > .Machine$double.eps) dir0 / n else c(0, 1)
        best <- NULL
        for (d in dists) {
          for (a in angles) {
            dd <- rot(dir0, a)
            px <- cx + dd[1] * (radius + d)
            py <- cy + dd[2] * (radius + d)
            if (!in_any_poly(px, py)) { best <- list(dd = dd, px = px, py = py); break }
          }
          if (!is.null(best)) break
        }
        if (is.null(best)) {  # nothing clear -> furthest radial position
          best <- list(dd = dir0,
                       px = cx + dir0[1] * (radius + max(dists)),
                       py = cy + dir0[2] * (radius + max(dists)))
        }
        ex[i] <- cx + best$dd[1] * radius
        ey[i] <- cy + best$dd[2] * radius
        lx[i] <- best$px
        ly[i] <- best$py
      }
      label_data$.sirap_lx <- lx
      label_data$.sirap_ly <- ly
      label_data$.sirap_ex <- ex
      label_data$.sirap_ey <- ey
    }
  }

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    # `plot_data` is the polygon geometry, so the spatial-barcode coordinates are
    # passed as `view_data`. This crops the plot to the sequenced spots (like
    # Seurat SpatialFeaturePlot) instead of the full tissue image.
    plots[[image]] <- .sirap_spatial_base(
      plot_data = polygon_data,
      image = spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = if (isTRUE(zoom)) NULL else spatial$coords
    ) +
      if (is.null(fill_column)) {
        ggplot2::geom_polygon(
          data = polygon_data,
          ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$polygon),
          inherit.aes = FALSE,
          fill = polygon_fill,
          color = polygon_color,
          alpha = polygon_alpha,
          linewidth = linewidth
        )
      } else {
        ggplot2::geom_polygon(
          data = polygon_data,
          ggplot2::aes(
            x = .data$x, y = .data$.sirap_y, group = .data$polygon,
            fill = .data$.sirap_fill_value
          ),
          inherit.aes = FALSE,
          color = polygon_color,
          alpha = polygon_alpha,
          linewidth = linewidth
        )
      }
    if (!is.null(fill_column)) {
      analysis_is_numeric <- is.numeric(polygon_data$.sirap_fill_value)
      if (analysis_is_numeric &&
          (!is.character(analysis_fill) || length(analysis_fill) != 2L ||
           anyNA(analysis_fill) || any(!nzchar(analysis_fill)))) {
        rlang::abort("`analysis_fill` must contain exactly two non-missing color values for numeric analysis data.")
      }
      fill_scale <- if (analysis_is_numeric) {
        ggplot2::scale_fill_gradient(
          low = analysis_fill[[1]],
          high = analysis_fill[[2]],
          na.value = "grey85"
        )
      } else {
        levels <- sort(unique(stats::na.omit(as.character(polygon_data$.sirap_fill_value))))
        fill_values <- .sirap_discrete_palette_values(levels)
        ggplot2::scale_fill_manual(values = fill_values, na.value = "grey85")
      }
      plots[[image]] <- plots[[image]] +
        fill_scale +
        ggplot2::labs(fill = fill_column)
    }
    if (isTRUE(label)) {
      if (isTRUE(label_outside)) {
        if (isTRUE(label_leader)) {
          plots[[image]] <- plots[[image]] +
            ggplot2::geom_segment(
              data = label_data,
              ggplot2::aes(x = .data$.sirap_ex, y = .data$.sirap_ey,
                           xend = .data$.sirap_lx, yend = .data$.sirap_ly),
              inherit.aes = FALSE,
              color = label_leader_color,
              linewidth = 0.25
            )
        }
        plots[[image]] <- plots[[image]] +
          ggplot2::geom_text(
            data = label_data,
            ggplot2::aes(x = .data$.sirap_lx, y = .data$.sirap_ly, label = .data$.sirap_label),
            inherit.aes = FALSE,
            size = label_size,
            color = label_color,
            fontface = "bold"
          )
      } else {
        plots[[image]] <- plots[[image]] +
          ggplot2::geom_text(
            data = label_data,
            ggplot2::aes(x = .data$x, y = .data$.sirap_y, label = .data$.sirap_label),
            inherit.aes = FALSE,
            size = label_size,
            color = label_color,
            fontface = "bold"
          )
      }
    }
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

#' Plot receptors inside spatial polygons
#'
#' Plots regional receptor rows from `spatialCloneDominance()` over the tissue
#' image, colored by a receptor-table column.
#'
#' @param dominance Object returned by `spatialCloneDominance()`.
#' @param object Seurat object with spatial images, or a SIRAP object.
#' @param color Column in `dominance$region_receptors` used to color receptor
#'   points. When omitted, the clone column recorded by
#'   `spatialCloneDominance()` is used.
#' @param images Optional Seurat image names to show.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param polygon_val Optional polygon ID value(s) to select from the dominance
#'   result.
#' @param point_size Size of receptor points.
#' @param point_alpha Opacity for receptor points.
#' @param polygon_alpha Fill opacity for polygons.
#' @param polygon_fill Fill color for polygons.
#' @param polygon_color Outline color for polygons.
#' @param image_alpha Tissue image opacity.
#' @param top_n Optional number of most frequent `color` levels to keep
#'   before grouping remaining levels as `other_label`.
#' @param other_label Label for grouped lower-frequency levels.
#' @param other_color Color used for `other_label`.
#' @param plot_other Plot receptors grouped into `other_label`.
#' @param zoom Zoom to the selected polygon extents. Defaults to `FALSE`, which
#'   shows the full tissue image.
#' @param combine Combine plots with patchwork when possible.
#' @return A ggplot/patchwork object or a named list of ggplots when
#'   `combine = FALSE`.
#' @export
plotSpatialRegionReceptors <- function(dominance,
                                       object,
                                       color = "clone_id",
                                       images = NULL,
                                       image_scale = NULL,
                                       swap_coords = NULL,
                                       polygon_val = NULL,
                                       point_size = 1,
                                       point_alpha = 0.5,
                                       polygon_alpha = 0.12,
                                       polygon_fill = "#FFFFFF",
                                       polygon_color = "black",
                                       image_alpha = 1,
                                       top_n = NULL,
                                       other_label = "Other",
                                       other_color = "#D3D3D3",
                                       plot_other = TRUE,
                                       zoom = FALSE,
                                       combine = TRUE) {
  color_missing <- missing(color)
  region <- .sirap_region_from_dominance(dominance, polygon_val = polygon_val)
  context <- .sirap_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (color_missing && !is.null(dominance$metadata$receptor_fields$clone)) {
    color <- dominance$metadata$receptor_fields$clone
  }
  .require_namespace("SeuratObject")
  if (nrow(region$receptors) == 0L) {
    rlang::abort("No receptor rows were found inside the selected polygon(s).")
  }
  .check_required_cols(region$receptors, color, "receptors")
  region$receptors$.sirap_color_by <- .sirap_lump_values(region$receptors[[color]], top_n = top_n, other_label = other_label)
  if (!isTRUE(plot_other)) {
    region$receptors <- region$receptors[as.character(region$receptors$.sirap_color_by) != other_label, , drop = FALSE]
  }
  if (nrow(region$receptors) == 0L) {
    rlang::abort("No receptor rows remain to plot after applying `top_n` and `plot_other`.")
  }
  region$receptors$.sirap_is_other <- as.character(region$receptors$.sirap_color_by) == other_label
  region$receptors <- region$receptors[order(!region$receptors$.sirap_is_other), , drop = FALSE]
  color_values <- .sirap_discrete_palette_values(
    levels(region$receptors$.sirap_color_by),
    other_label = other_label,
    other_color = other_color
  )

  plots <- list()
  full_image_data <- data.frame(x = numeric(), .sirap_y = numeric())
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    plot_receptors <- region$receptors[region$receptors$image == image | region$receptors$image == "metadata", , drop = FALSE]
    plots[[image]] <- .sirap_spatial_base(
      plot_data = if (isTRUE(zoom)) region$polygons else spatial$coords,
      image = spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = if (isTRUE(zoom)) NULL else full_image_data
    ) +
      ggplot2::geom_polygon(
        data = region$polygons,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$polygon),
        inherit.aes = FALSE,
        fill = polygon_fill,
        color = polygon_color,
        alpha = polygon_alpha,
        linewidth = 0.5
      ) +
      ggplot2::geom_point(
        data = plot_receptors,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, color = .data$.sirap_color_by),
        inherit.aes = FALSE,
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = color_values, drop = TRUE) +
      ggplot2::labs(color = color)
  }

  .sirap_return_spatial_plots(plots, combine = combine)
}

#' Plot 3D clone-density surfaces over a spatial image
#'
#' Uses the regional receptor rows returned by `spatialCloneDominance()` to draw
#' smoothed clone-abundance density surfaces above the original tissue image in
#' an interactive `rgl` scene.
#'
#' @param dominance Object returned by `spatialCloneDominance()`.
#' @param object Seurat object with spatial images, or a SIRAP object.
#' @param clone Column in `dominance$region_receptors` used to define clone
#'   groups. When omitted, the clone column recorded by
#'   `spatialCloneDominance()` is used.
#' @param weight Optional numeric column used as receptor weights. When
#'   `NULL`, each receptor row contributes one count.
#' @param images Optional Seurat image names to show.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param polygon_val Optional polygon ID value(s) to select from the dominance
#'   object.
#' @param top_n Optional number of most abundant clone groups to show
#'   individually. Remaining clone groups are combined as `other_label`.
#' @param other_label Label for grouped lower-abundance clone groups.
#' @param other_color Color used for `other_label`.
#' @param plot_other Plot grouped lower-abundance clone groups as
#'   `other_label`.
#' @param bandwidth Gaussian kernel bandwidth in spatial coordinate units. Use a
#'   single value for isotropic smoothing or two values for x/y bandwidths. When
#'   `NULL`, a bandwidth is estimated from the regional coordinates.
#' @param grid_n Number of grid points used for density surfaces. Use one value
#'   for both axes or two values for x/y grid sizes.
#' @param plot_type Plot smoothed density surfaces or binned 3D histogram bars.
#' @param hist_bins Number of x/y bins used when `plot_type = "histogram"`.
#' @param zoom Zoom the tissue image plane and 3D view to the selected region.
#'   Defaults to `FALSE`, which shows the full spatial image.
#' @param height_scale Maximum density height as a fraction of the plotted
#'   spatial span.
#' @param z_offset Vertical offset of the density surfaces above the tissue
#'   image, as a fraction of the plotted spatial span.
#' @param density_alpha Opacity for density surfaces.
#' @param density_floor Relative density threshold below which surface cells are
#'   hidden. This removes low-density Gaussian tails that otherwise look like a
#'   colored plane.
#' @param image_alpha Opacity for the tissue image plane.
#' @param show_points Show receptor locations above the tissue plane.
#' @param point_size Size of receptor-location points.
#' @param point_alpha Opacity for receptor-location points.
#' @param show_polygon Show polygon boundaries above the tissue plane.
#' @param polygon_color Color for polygon boundary lines.
#' @param polygon_lwd Width for polygon boundary lines.
#' @param label_density_peaks Label the highest point of the first few density
#'   surfaces.
#' @param label_top_n Number of density surfaces to label when
#'   `label_density_peaks = TRUE`.
#' @param label_cex Text scale for density peak labels.
#' @param show_axes Show the 3D bounding box and axis labels.
#' @param texture_dir Optional directory where the tissue texture PNG should be
#'   written. Use a persistent directory when saving an HTML widget.
#' @param as_widget Return an `rgl` htmlwidget. When `FALSE`, the scene is drawn
#'   and the current rgl device ID is returned.
#' @param return_data Return the computed density data instead of drawing an rgl
#'   scene. This does not require `rgl`.
#' @return An `rgl` widget/device ID, or a list of computed plot data when
#'   `return_data = TRUE`.
#' @export
plotSpatialCloneDensity3D <- function(dominance,
                                      object,
                                      clone = "clone_id",
                                      weight = NULL,
                                      images = NULL,
                                      image_scale = NULL,
                                      swap_coords = NULL,
                                      polygon_val = NULL,
                                      top_n = 10,
                                      other_label = "Other",
                                      other_color = "#D3D3D3",
                                      plot_other = TRUE,
                                      bandwidth = NULL,
                                      grid_n = 80,
                                      plot_type = c("density", "histogram"),
                                      hist_bins = 30,
                                      zoom = FALSE,
                                      height_scale = 0.22,
                                      z_offset = 0.04,
                                      density_alpha = 0.55,
                                      density_floor = 0.03,
                                      image_alpha = 1,
                                      show_points = FALSE,
                                      point_size = 5,
                                      point_alpha = 0.75,
                                      show_polygon = TRUE,
                                      polygon_color = "white",
                                      polygon_lwd = 3,
                                      label_density_peaks = FALSE,
                                      label_top_n = 5,
                                      label_cex = 0.8,
                                      show_axes = TRUE,
                                      texture_dir = NULL,
                                      as_widget = TRUE,
                                      return_data = FALSE) {
  clone_missing <- missing(clone)
  if (!inherits(dominance, "sirap_spatial_clone_dominance")) {
    rlang::abort("`dominance` must be an object returned by `spatialCloneDominance()`.")
  }
  if (is.null(dominance$region_receptors) || is.null(dominance$polygons)) {
    rlang::abort("`dominance` must contain `region_receptors` and `polygons` from `spatialCloneDominance()`.")
  }
  context <- .sirap_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (clone_missing && !is.null(dominance$metadata$receptor_fields$clone)) {
    clone <- dominance$metadata$receptor_fields$clone
  }
  plot_type <- match.arg(plot_type)
  .require_namespace("SeuratObject")

  region_receptors <- .as_data_frame(dominance$region_receptors)
  polygons <- .as_data_frame(dominance$polygons)
  if (!is.null(polygon_val)) {
    region_receptors <- region_receptors[as.character(region_receptors$polygon) %in% as.character(polygon_val), , drop = FALSE]
    polygons <- polygons[as.character(polygons$polygon) %in% as.character(polygon_val), , drop = FALSE]
  }
  if (nrow(region_receptors) == 0L) {
    rlang::abort("No regional receptor rows are available for the selected polygon(s).")
  }

  plot_data <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    image_receptors <- region_receptors[region_receptors$image == image | region_receptors$image == "metadata", , drop = FALSE]
    if (nrow(image_receptors) == 0L) {
      next
    }
    plot_data[[image]] <- .sirap_clone_density_3d_data(
      region_receptors = image_receptors,
      polygons = polygons,
      image = spatial$image,
      spatial_coords = spatial$coords,
      clone = clone,
      weight = weight,
      top_n = top_n,
      other_label = other_label,
      other_color = other_color,
      plot_other = plot_other,
      bandwidth = bandwidth,
      grid_n = grid_n,
      plot_type = plot_type,
      hist_bins = hist_bins,
      zoom = zoom,
      height_scale = height_scale,
      z_offset = z_offset,
      density_floor = density_floor
    )
  }
  plot_data <- Filter(Negate(is.null), plot_data)
  if (length(plot_data) == 0L) {
    rlang::abort("No regional receptor rows match the selected image(s).")
  }
  if (isTRUE(return_data)) {
    return(if (length(plot_data) == 1L) plot_data[[1L]] else plot_data)
  }

  if (isTRUE(as_widget)) {
    old_rgl_use_null <- getOption("rgl.useNULL")
    options(rgl.useNULL = TRUE)
    on.exit(options(rgl.useNULL = old_rgl_use_null), add = TRUE)
  }
  .require_namespace("rgl")
  plots <- lapply(plot_data, function(dat) {
    .sirap_draw_clone_density_3d(
      dat,
      density_alpha = density_alpha,
      image_alpha = image_alpha,
      show_points = show_points,
      point_size = point_size,
      point_alpha = point_alpha,
      show_polygon = show_polygon,
      polygon_color = polygon_color,
      polygon_lwd = polygon_lwd,
      label_density_peaks = label_density_peaks,
      label_top_n = label_top_n,
      label_cex = label_cex,
      show_axes = show_axes,
      texture_dir = texture_dir,
      as_widget = as_widget
    )
  })
  if (length(plots) == 1L) {
    plots[[1L]]
  } else {
    plots
  }
}

#' Plot spatial polygon pie summaries
#'
#' Summarizes regional receptor rows from `spatialCloneDominance()` for a
#' selected column. The default mode plots pie charts beside a numbered polygon
#' map. The `spatial_clipped` mode clips each pie chart to its polygon boundary
#' and requires the suggested `sf` package.
#'
#' @param dominance Object returned by `spatialCloneDominance()`.
#' @param object Seurat object with spatial images, or a SIRAP object.
#' @param fill Column in `dominance$region_receptors` summarized in each pie
#'   chart. When omitted, the clone column recorded by
#'   `spatialCloneDominance()` is used.
#' @param mode Plot style. `side_by_side` shows ordinary pies beside a polygon
#'   map. `spatial_clipped` clips pies to polygon boundaries on the tissue.
#' @param top_n Number of most frequent `fill` levels to identify within each
#'   polygon. The union of these polygon-level top levels is shown individually
#'   everywhere; remaining levels are grouped as `other_label`.
#' @param label_top_n Number of levels to label per polygon. Defaults to
#'   0 because this gets messy fast.
#' @param other_label Label for grouped lower-frequency levels.
#' @param other_color Color used for `other_label`.
#' @param plot_other Plot grouped lower-frequency BCRs as `other_label`.
#' @param images Optional Seurat image names to show.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param polygon_val Optional polygon ID value(s) to select from the dominance
#'   result.
#' @param polygon_alpha Fill opacity for polygons in the numbered map.
#' @param image_alpha Tissue image opacity.
#' @param pie_alpha Pie opacity.
#' @param label_size Text size for polygon and pie labels.
#' @param pies_per_page Maximum number of ordinary pie charts to show per page
#'   when `mode = "side_by_side"`. When more polygons are present, a named list
#'   of page plots is returned.
#' @param pie_order Order of the pie charts when `mode = "side_by_side"`.
#'   `"receptors"` (the default) shows the polygon with the most receptors
#'   first, then the next most, and so on. Each pie keeps its map label, so
#'   only the arrangement changes. `"polygon"` keeps the numbered map order
#'   (polygon 1, 2, 3, ...). Ignored when `mode = "spatial_clipped"`.
#' @param spatial_legend Show a shared fill legend for `spatial_clipped` pie
#'   plots. Defaults to `FALSE`.
#' @param zoom Zoom to the plotted polygon extents. Defaults to `FALSE`, which
#'   shows the full tissue image.
#' @return A ggplot, patchwork object, or named list of ggplots.
#' @export
plotSpatialPolygonPies <- function(dominance,
                                   object,
                                   fill = "clone_id",
                                   mode = c("side_by_side", "spatial_clipped", "dominance_only"),
                                   top_n = 10,
                                   label_top_n = 0,
                                   other_label = "Other",
                                   other_color = "#D3D3D3",
                                   plot_other = TRUE,
                                   images = NULL,
                                   image_scale = NULL,
                                   swap_coords = NULL,
                                   polygon_val = NULL,
                                   polygon_alpha = 0.2,
                                   image_alpha = 1,
                                   pie_alpha = 0.9,
                                   label_size = 3,
                                   pies_per_page = 6L,
                                   pie_order = c("receptors", "polygon"),
                                   spatial_legend = FALSE,
                                   zoom = FALSE) {
  fill_missing <- missing(fill)
  region <- .sirap_region_from_dominance(dominance, polygon_val = polygon_val)
  context <- .sirap_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (fill_missing && !is.null(dominance$metadata$receptor_fields$clone)) {
    fill <- dominance$metadata$receptor_fields$clone
  }
  .require_namespace("SeuratObject")
  mode <- match.arg(mode)
  pie_order <- match.arg(pie_order)
  if (is.null(label_top_n)) {
    label_top_n <- Inf
  }
  if (!is.numeric(label_top_n) || length(label_top_n) != 1L || is.na(label_top_n) || label_top_n < 0) {
    rlang::abort("`label_top_n` must be NULL or a non-negative number.")
  }
  if (!is.numeric(pies_per_page) || length(pies_per_page) != 1L || is.na(pies_per_page) || pies_per_page < 1) {
    rlang::abort("`pies_per_page` must be a positive number.")
  }
  pies_per_page <- as.integer(pies_per_page)
  if (nrow(region$receptors) == 0L) {
    rlang::abort("No receptor rows were found inside the selected polygon(s).")
  }
  .check_required_cols(region$receptors, fill, "receptors")
  pie_data <- .sirap_polygon_pie_summary(
    region$receptors,
    fill = fill,
    top_n = top_n,
    label_top_n = label_top_n,
    other_label = other_label,
    plot_other = plot_other
  )
  if (nrow(pie_data) == 0L) {
    rlang::abort("No receptor rows remain to plot after applying `top_n` and `plot_other`.")
  }
  fill_values <- .sirap_discrete_palette_values(
    unique(as.character(pie_data$level)),
    other_label = other_label,
    other_color = other_color
  )
  centers <- .sirap_polygon_centers(region$polygons, pie_data$polygon)

  if (mode == "dominance_only") {
    .sirap_dominance_only_polygon_pies(
      pie_data = pie_data, fill_values = fill_values,
      label_size = label_size, pie_order = pie_order
    )
  } else if (mode == "side_by_side") {
    .sirap_side_by_side_polygon_pies(
      object = object,
      images = images,
      image_scale = image_scale,
      swap_coords = swap_coords,
      polygon_plot_data = region$polygons,
      centers = centers,
      pie_data = pie_data,
      fill_values = fill_values,
      polygon_alpha = polygon_alpha,
      image_alpha = image_alpha,
      label_size = label_size,
      pies_per_page = pies_per_page,
      pie_order = pie_order,
      zoom = zoom
    )
  } else {
    .sirap_spatial_clipped_polygon_pies(
      object = object,
      images = images,
      image_scale = image_scale,
      swap_coords = swap_coords,
      polygon_plot_data = region$polygons,
      centers = centers,
      pie_data = pie_data,
      fill_values = fill_values,
      fill = fill,
      image_alpha = image_alpha,
      pie_alpha = pie_alpha,
      label_size = label_size,
      spatial_legend = spatial_legend,
      zoom = zoom
    )
  }
}

#' Subset a Seurat object to spatial IDs inside polygons
#'
#' Keeps only Seurat spatial IDs whose spatial coordinates fall inside one or more
#' polygon regions and records polygon membership in the returned object's
#' metadata.
#'
#' @param object Seurat object with spatial coordinates.
#' @param polygon Data frame of polygon vertices in order. Must contain `x` and
#'   `y` columns by default. If a `polygon` column is present, each polygon is
#'   evaluated independently.
#' @param images Optional Seurat image names to use for spatial coordinates.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param polygon_val Optional polygon ID value(s) to select when `polygon_id`
#'   exists.
#' @param x Coordinate column in `polygon`.
#' @param y Coordinate column in `polygon`.
#' @param polygon_id Optional polygon grouping column in `polygon`.
#' @param metadata Metadata column added to the returned Seurat object with the
#'   primary polygon assignment. Barcodes in overlapping polygons are assigned to
#'   the smallest polygon by area.
#' @param include_boundary Include Seurat spots on polygon edges.
#' @return A Seurat object subset to spatial IDs inside the selected polygon(s),
#'   with `metadata` appended to `@meta.data`.
#' @export
seuratSubsetPolygon <- function(object,
                                polygon,
                                images = NULL,
                                image_scale = NULL,
                                swap_coords = NULL,
                                polygon_val = NULL,
                                x = "x",
                                y = "y",
                                polygon_id = "polygon",
                                metadata = "polygon",
                                include_boundary = TRUE) {
  annotated <- assignSpatialPolygons(
    object = object,
    polygon = polygon,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords,
    polygon_val = polygon_val,
    x = x,
    y = y,
    polygon_id = polygon_id,
    metadata = metadata,
    include_boundary = include_boundary
  )
  keep_cells <- rownames(annotated@meta.data)[!is.na(annotated@meta.data[[metadata]])]
  keep_cells <- keep_cells[keep_cells %in% colnames(annotated)]
  if (length(keep_cells) == 0L) {
    rlang::abort("No Seurat spatial IDs were found inside the selected polygon(s).")
  }
  suppressWarnings(annotated[, keep_cells], classes = "validationWarning")
}

#' Assign polygon IDs to spatial metadata
#'
#' Adds polygon-membership columns to a Seurat object's metadata without
#' subsetting the object. Spatial IDs outside the selected polygon(s) receive
#' `NA` in the primary metadata column. Spatial IDs inside overlapping polygons
#' are assigned to the smallest polygon by area in the primary metadata column.
#' Additional logical columns named `polygon_<ID>` record membership in each
#' individual polygon. When `object` is a SIRAP object and `polygon` is omitted,
#' the stored polygons in the SIRAP object are used.
#'
#' @param object Seurat object with spatial coordinates, or a SIRAP object.
#' @param polygon Data frame of polygon vertices in order. Must contain `x` and
#'   `y` columns by default. If a `polygon` column is present, each polygon is
#'   evaluated independently. If `object` is a SIRAP object and `polygon` is
#'   omitted, the object's stored polygons are used.
#' @param images Optional Seurat image names to use for spatial coordinates.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param polygon_val Optional polygon ID value(s) to select when `polygon_id`
#'   exists.
#' @param x Coordinate column in `polygon`.
#' @param y Coordinate column in `polygon`.
#' @param polygon_id Optional polygon grouping column in `polygon`.
#' @param metadata Metadata column added to the returned Seurat object with the
#'   primary polygon assignment. Spatial IDs in overlapping polygons are assigned
#'   to the smallest polygon by area.
#' @param include_boundary Include Seurat spots on polygon edges.
#' @return A Seurat object with `metadata` added to `@meta.data`.
#' @export
assignSpatialPolygons <- function(object,
                                  polygon = NULL,
                                  images = NULL,
                                  image_scale = NULL,
                                  swap_coords = NULL,
                                  polygon_val = NULL,
                                  x = "x",
                                  y = "y",
                                  polygon_id = "polygon",
                                  metadata = "polygon_id",
                                  include_boundary = TRUE) {
  .require_namespace("SeuratObject")
  if (!.sirap_is_non_empty_string(metadata)) {
    rlang::abort("`metadata` must be a non-empty character string.")
  }
  if (isSIRAP(object) && is.null(polygon)) {
    polygon <- .sirap_polygon_vertices_table(.sirap_polygons(object))
    if (nrow(polygon) == 0L) {
      rlang::abort("The SIRAP object does not contain stored polygons.")
    }
  }
  context <- .sirap_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  membership <- .sirap_polygon_membership(
    object = object,
    polygon = polygon,
    images = images,
    image_scale = context$image_scale,
    swap_coords = context$swap_coords,
    polygon_val = polygon_val,
    x = x,
    y = y,
    polygon_id = polygon_id,
    include_boundary = include_boundary
  )

  polygon_lookup <- .sirap_polygon_lookup(membership)
  polygon_flags <- .sirap_polygon_flag_lookup(membership)
  object <- .sirap_update_seurat_image_slots(object)
  object@meta.data[[metadata]] <- unname(polygon_lookup[rownames(object@meta.data)])
  object@meta.data <- .sirap_append_polygon_flags(
    object@meta.data,
    rownames(object@meta.data),
    polygon_flags
  )
  object
}

#' Assign polygon IDs to a receptor table
#'
#' Adds polygon-membership columns to an AIRR-style receptor table without
#' subsetting the table. When `object` is supplied, receptor rows whose spatial
#' ID falls outside the selected polygon(s), or does not match the spatial
#' object, receive `NA`. When `object = NULL`, polygon membership is computed
#' directly from coordinate columns in `receptors`.
#' Receptors in spatial IDs inside overlapping polygons are assigned to the
#' smallest polygon by area in the primary polygon column. Additional logical
#' columns named `polygon_<ID>` record membership in each individual polygon.
#' When `receptors` is a SIRAP object and `object` and `polygon` are omitted,
#' the stored spatial object and polygons are used.
#'
#' @param receptors AIRR-style receptor data frame, or a SIRAP object when
#'   `object` is omitted.
#' @param object Seurat object with spatial coordinates, or a SIRAP object.
#' @param column Column added to the returned receptor table with polygon
#'   membership.
#' @param spatial_id Column in `receptors` matching Seurat spot/cell spatial IDs.
#'   If `NULL`, a spatial ID-like column is guessed or pulled from the SIRAP
#'   object.
#' @param polygon Data frame of polygon vertices in order. Must contain `x` and
#'   `y` columns by default. If a `polygon` column is present, each polygon is
#'   evaluated independently. If `object` is a SIRAP object and `polygon` is
#'   omitted, the object's stored polygons are used.
#' @param images Optional Seurat image names to use for spatial coordinates.
#' @param image_scale Image scale passed to
#'   `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Logical; if `TRUE`, swap the two Seurat coordinate
#'   axes before plotting or polygon membership checks. If `NULL`, the
#'   value must be supplied explicitly or by SIRAP object metadata.
#' @param polygon_val Optional polygon ID value(s) to select when `polygon_id`
#'   exists.
#' @param x Coordinate column in `polygon`.
#' @param y Coordinate column in `polygon`.
#' @param receptor_x,receptor_y Coordinate columns in `receptors` used when
#'   `object = NULL`.
#' @param polygon_id Optional polygon grouping column in `polygon`.
#' @param include_boundary Include Seurat spots on polygon edges.
#' @return A tibble with `column` added.
#' @export
assignReceptorPolygons <- function(receptors,
                                   object = NULL,
                                   polygon = NULL,
                                   column = "polygon_id",
                                   spatial_id = NULL,
                                   images = NULL,
                                   image_scale = NULL,
                                   swap_coords = NULL,
                                   polygon_val = NULL,
                                   x = "x",
                                   y = "y",
                                   receptor_x = "x",
                                   receptor_y = "y",
                                   polygon_id = "polygon",
                                   include_boundary = TRUE) {
  if (!.sirap_is_non_empty_string(column)) {
    rlang::abort("`column` must be a non-empty character string.")
  }
  if (isSIRAP(receptors) && is.null(object)) {
    object <- receptors
    receptors <- NULL
  }
  if (isSIRAP(object) && is.null(polygon)) {
    polygon <- .sirap_polygon_vertices_table(.sirap_polygons(object))
    if (nrow(polygon) == 0L) {
      rlang::abort("The SIRAP object does not contain stored polygons.")
    }
  }
  if (is.null(object)) {
    if (is.null(polygon)) {
      rlang::abort("`polygon` must be supplied when `object = NULL`.")
    }
    receptors <- .as_data_frame(receptors)
    .check_required_cols(receptors, c(receptor_x, receptor_y), "receptors")
    receptors$.sirap_receptor_row <- seq_len(nrow(receptors))
    coords <- data.frame(
      cell = as.character(receptors$.sirap_receptor_row),
      image = NA_character_,
      x = as.numeric(receptors[[receptor_x]]),
      y = as.numeric(receptors[[receptor_y]]),
      stringsAsFactors = FALSE
    )
    membership <- .sirap_polygon_membership_from_coords(
      coords = coords,
      polygon = polygon,
      polygon_val = polygon_val,
      x = x,
      y = y,
      polygon_id = polygon_id,
      include_boundary = include_boundary
    )
    polygon_lookup <- .sirap_polygon_lookup(membership)
    polygon_flags <- .sirap_polygon_flag_lookup(membership)
    out <- tibble::as_tibble(receptors)
    out[[column]] <- unname(polygon_lookup[as.character(out$.sirap_receptor_row)])
    out <- .sirap_append_polygon_flags(
      out,
      as.character(out$.sirap_receptor_row),
      polygon_flags
    )
    out$.sirap_receptor_row <- NULL
    return(out)
  }
  options <- .sirap_spatial_options(object, image_scale = image_scale, swap_coords = swap_coords)
  inputs <- .sirap_resolve_inputs(object, receptors, spatial_id)
  object <- inputs$object
  receptors <- .as_data_frame(inputs$receptors)
  spatial_id <- .guess_spatial_id_col(receptors, inputs$spatial_id_col)
  .check_required_cols(receptors, spatial_id, "receptors")

  membership <- .sirap_polygon_membership(
    object = object,
    polygon = polygon,
    images = images,
    image_scale = options$image_scale,
    swap_coords = options$swap_coords,
    polygon_val = polygon_val,
    x = x,
    y = y,
    polygon_id = polygon_id,
    include_boundary = include_boundary
  )
  polygon_lookup <- .sirap_polygon_lookup(membership)
  polygon_flags <- .sirap_polygon_flag_lookup(membership)
  out <- tibble::as_tibble(receptors)
  out[[column]] <- unname(polygon_lookup[as.character(out[[spatial_id]])])
  out <- .sirap_append_polygon_flags(
    out,
    as.character(out[[spatial_id]]),
    polygon_flags
  )
  out
}

#' Summarize receptor content inside spatial polygons
#'
#' Builds one receptor summary row per polygon from a SIRAP object or from an
#' already polygon-annotated receptor table. For SIRAP objects, the summary is
#' also stored as a named analysis and the updated object is returned. For a
#' receptor table input, the summary tibble is returned directly.
#'
#' @param object A SIRAP object or an AIRR-style receptor data frame.
#' @param polygons Optional polygon table used to define the output polygon
#'   order. When `NULL`, polygon IDs are taken from the SIRAP object or from
#'   `polygon_id` values already present in `receptors`.
#' @param spatial_id Column identifying spatial spots/cells. If `NULL`, the
#'   value is pulled from the SIRAP object or guessed from `receptors`.
#' @param polygon_id Primary polygon assignment column in `receptors`.
#' @param clone Clone column used for `n_clones` and `top_clone_fraction`.
#'   If this column is absent, those clone summaries are omitted.
#' @param weight Optional numeric receptor weight column. When supplied,
#'   `n_receptors`, `top_clone_fraction`, and numeric custom summaries are
#'   weighted by this column.
#' @param summary_columns Optional receptor columns to summarize.
#' @param summary_operations Operations to apply to `summary_columns`. Must be
#'   length 1 or the same length as `summary_columns`. Supported operations are
#'   `count`, `sum`, `mean`, `median`, `variance`, `var`, `sd`, `min`, `max`,
#'   `any`, `all`, and `cat`.
#' @param analysis_name Name used when storing the summary in a SIRAP object.
#' @param validate Whether to validate an updated SIRAP object.
#' @return An updated SIRAP object for SIRAP input, otherwise a summary tibble.
#' @export
polygonReceptorSummary <- function(object,
                                   polygons = NULL,
                                   spatial_id = NULL,
                                   polygon_id = "polygon",
                                   clone = "clone_id",
                                   weight = NULL,
                                   summary_columns = NULL,
                                   summary_operations = NULL,
                                   analysis_name = "polygon_receptor_summary",
                                   validate = TRUE) {
  is_sirap <- isSIRAP(object)
  x <- if (is_sirap) object else NULL
  polygon_id_missing <- missing(polygon_id)
  polygon_column <- if (is_sirap && isTRUE(polygon_id_missing)) "polygon_id" else polygon_id
  receptor_table <- if (is_sirap) getReceptors(x) else tibble::as_tibble(object)
  if (is_sirap) {
    fields <- .sirap_receptor_fields(x)
    spatial_id <- spatial_id %||% fields$spatial_id
    if (missing(clone) && !is.null(fields$clone)) {
      clone <- fields$clone
    }
    polygons <- polygons %||% getPolygons(x)
  } else {
    spatial_id <- tryCatch(
      .guess_spatial_id_col(receptor_table, spatial_id),
      error = function(e) spatial_id
    )
  }

  polygon_ids <- .sirap_polygon_ids(receptor_table, polygons = polygons, polygon_id = polygon_column)
  summary <- .sirap_polygon_receptor_summary_table(
    receptors = receptor_table,
    polygon_ids = polygon_ids,
    spatial_id = spatial_id,
    polygon_id = polygon_column,
    clone = clone,
    weight = weight,
    summary_columns = summary_columns,
    summary_operations = summary_operations
  )

  if (is_sirap) {
    return(setAnalysis(x, analysis_name, summary, validate = validate))
  }
  summary
}

#' Run a function on receptor rows from stored spatial polygons
#'
#' Splits the receptor table by polygon membership and applies a function to
#' each polygon's receptor rows, or applies the function once to all receptors
#' that fall in any polygon.
#'
#' @param receptors A SIRAP object or an AIRR-style receptor data frame.
#' @param .f Function to apply.
#' @param ... Additional arguments passed to `.f`.
#' @param polygons Optional polygon table used to define the output polygon
#'   order. When `NULL`, polygon IDs are taken from the SIRAP object or from
#'   `polygon_id` values already present in `receptors`.
#' @param polygon_id Primary polygon assignment column in `receptors`.
#' @param mode `"by_polygon"` runs `.f` once per polygon; `"all_polygons"` runs
#'   `.f` once on all receptors in any polygon.
#' @param pass_polygon_id Pass the current polygon ID to `.f` as `polygon_id`
#'   when `.f` can accept it. Defaults to `FALSE`.
#' @return A named list for `mode = "by_polygon"`, otherwise the function result.
#' @export
withPolygons <- function(receptors,
                         .f,
                         ...,
                         polygons = NULL,
                         polygon_id = "polygon",
                         mode = c("by_polygon", "all_polygons"),
                         pass_polygon_id = FALSE) {
  mode <- match.arg(mode)
  is_sirap <- isSIRAP(receptors)
  x <- if (is_sirap) receptors else NULL
  polygon_id_missing <- missing(polygon_id)
  polygon_column <- if (is_sirap && isTRUE(polygon_id_missing)) "polygon_id" else polygon_id
  receptor_table <- if (is_sirap) getReceptors(x) else tibble::as_tibble(receptors)
  if (is_sirap) {
    polygons <- polygons %||% getPolygons(x)
  }
  polygon_ids <- .sirap_polygon_ids(receptor_table, polygons = polygons, polygon_id = polygon_column)
  groups <- .sirap_polygon_receptor_groups(
    receptors = receptor_table,
    polygon_ids = polygon_ids,
    polygon_id = polygon_column,
    mode = mode
  )
  fun <- match.fun(.f)

  if (identical(mode, "all_polygons")) {
    return(.sirap_call_polygon_function(fun, groups[[1L]], polygon_id = NULL, pass_polygon_id = pass_polygon_id, ...))
  }
  stats::setNames(
    lapply(names(groups), function(id) .sirap_call_polygon_function(fun, groups[[id]], polygon_id = id, pass_polygon_id = pass_polygon_id, ...)),
    names(groups)
  )
}

#' Estimate receptor clonal abundance inside spatial polygons
#'
#' Runs [alakazam::estimateAbundance()] on receptor rows in each polygon, or on
#' all receptors that fall in any polygon. Polygons with fewer than `min_n`
#' receptor rows return `NULL`. SIRAP input stores the result as an analysis and
#' returns the updated object.
#'
#' @param receptors A SIRAP object or an AIRR-style receptor data frame.
#' @param polygons Optional polygon table used to define the output polygon
#'   order. When `NULL`, polygon IDs are taken from the SIRAP object or from
#'   `polygon_id` values already present in `receptors`.
#' @param polygon_id Primary polygon assignment column in `receptors`.
#' @param mode `"by_polygon"` runs the abundance calculation once per polygon;
#'   `"all_polygons"` runs it once on all receptors in any polygon.
#' @param analysis_name Name used when storing the abundance result in a SIRAP
#'   object.
#' @param clone alakazam: name of the data column containing clone identifiers.
#' @param copy alakazam: name of the data column containing copy numbers for
#'   each sequence. If copy=NULL (the default), then clone abundance is
#'   determined by the number of sequences. If a copy column is specified, then
#'   clone abundances is determined by the sum of copy numbers within each
#'   clonal group.
#' @param group alakazam: name of the data column containing group identifiers.
#'   If NULL then no grouping is performed and the group column of the output
#'   will contain the value NA for each row.
#' @param min_n alakazam: minimum number of observations to sample. A group with
#'   less observations than the minimum is excluded.
#' @param max_n alakazam: maximum number of observations to sample. If NULL then
#'   no maximum is set.
#' @param uniform alakazam: if TRUE then uniformly resample each group to the
#'   same number of observations. If FALSE then allow each group to be resampled
#'   to its original size or, if specified, max_size.
#' @param ci alakazam: confidence interval to calculate; the value must be
#'   between 0 and 1.
#' @param nboot alakazam: number of bootstrap realizations to generate.
#' @param cell_id alakazam: name of the data column containing cell identifiers.
#'   If cell_id=NULL then the function will assume bulk data.
#' @param progress alakazam: if TRUE show a progress bar.
#' @param validate Whether to validate an updated SIRAP object.
#' @return An updated SIRAP object for SIRAP input, otherwise abundance result(s).
#' @export
polygonAbundance <- function(receptors,
                             polygons = NULL,
                             polygon_id = "polygon",
                             mode = c("by_polygon", "all_polygons"),
                             analysis_name = "polygon_abundance",
                             clone = "clone_id",
                             copy = NULL,
                             group = NULL,
                             min_n = 5,
                             max_n = NULL,
                             uniform = TRUE,
                             ci = 0.95,
                             nboot = 2000,
                             cell_id = NULL,
                             progress = FALSE,
                             validate = TRUE) {
  .require_namespace("alakazam")
  clone_missing <- missing(clone)
  polygon_id_missing <- missing(polygon_id)
  if (isSIRAP(receptors) && isTRUE(clone_missing)) {
    clone <- .sirap_receptor_fields(receptors)$clone %||% clone
  }
  polygon_column <- if (isSIRAP(receptors) && isTRUE(polygon_id_missing)) "polygon_id" else polygon_id
  mode <- match.arg(mode)
  result <- .sirap_polygon_alakazam(
    receptors = receptors,
    polygons = polygons,
    polygon_id = polygon_column,
    mode = mode,
    threshold_n = min_n,
    fun = alakazam::estimateAbundance,
    clone = clone,
    copy = copy,
    group = group,
    min_n = min_n,
    max_n = max_n,
    uniform = uniform,
    ci = ci,
    nboot = nboot,
    cell_id = cell_id,
    progress = progress
  )
  if (isSIRAP(receptors)) {
    return(setAnalysis(receptors, analysis_name, result, validate = validate))
  }
  result
}

#' Estimate receptor clonal diversity inside spatial polygons
#'
#' Runs [alakazam::alphaDiversity()] on receptor rows in each polygon, or on all
#' receptors that fall in any polygon. Polygons with fewer than `min_n` receptor
#' rows return `NULL`. SIRAP input stores the result as an analysis and returns
#' the updated object.
#'
#' @param receptors A SIRAP object or an AIRR-style receptor data frame.
#' @param polygons Optional polygon table used to define the output polygon
#'   order. When `NULL`, polygon IDs are taken from the SIRAP object or from
#'   `polygon_id` values already present in `receptors`.
#' @param polygon_id Primary polygon assignment column in `receptors`.
#' @param mode `"by_polygon"` runs the diversity calculation once per polygon;
#'   `"all_polygons"` runs it once on all receptors in any polygon.
#' @param analysis_name Name used when storing the diversity result in a SIRAP
#'   object.
#' @param clone alakazam: name of the data column containing clone identifiers.
#' @param copy alakazam: name of the data column containing copy numbers for
#'   each sequence. If copy=NULL (the default), then clone abundance is
#'   determined by the number of sequences. If a copy column is specified, then
#'   clone abundances is determined by the sum of copy numbers within each
#'   clonal group.
#' @param group alakazam: name of the data column containing group identifiers.
#'   If NULL then no grouping is performed and the group column of the output
#'   will contain the value NA for each row.
#' @param min_n alakazam: minimum number of observations to sample. A group with
#'   less observations than the minimum is excluded.
#' @param min_q alakazam: minimum value of q.
#' @param max_q alakazam: maximum value of q.
#' @param step_q alakazam: value by which to increment q.
#' @param ci alakazam: confidence interval to calculate; the value must be
#'   between 0 and 1.
#' @param ... alakazam: additional arguments to pass to estimateAbundance.
#'   Additional arguments are ignored if a AbundanceCurve is provided as input.
#' @param validate Whether to validate an updated SIRAP object.
#' @return An updated SIRAP object for SIRAP input, otherwise diversity result(s).
#' @export
polygonDiversity <- function(receptors,
                             polygons = NULL,
                             polygon_id = "polygon",
                             mode = c("by_polygon", "all_polygons"),
                             analysis_name = "polygon_diversity",
                             clone = "clone_id",
                             copy = NULL,
                             group = NULL,
                             min_n = 30,
                             min_q = 0,
                             max_q = 4,
                             step_q = 0.1,
                             ci = 0.95,
                             ...,
                             validate = TRUE) {
  .require_namespace("alakazam")
  clone_missing <- missing(clone)
  polygon_id_missing <- missing(polygon_id)
  if (isSIRAP(receptors) && isTRUE(clone_missing)) {
    clone <- .sirap_receptor_fields(receptors)$clone %||% clone
  }
  polygon_column <- if (isSIRAP(receptors) && isTRUE(polygon_id_missing)) "polygon_id" else polygon_id
  mode <- match.arg(mode)
  result <- .sirap_polygon_alakazam(
    receptors = receptors,
    polygons = polygons,
    polygon_id = polygon_column,
    mode = mode,
    threshold_n = min_n,
    fun = alakazam::alphaDiversity,
    clone = clone,
    copy = copy,
    group = group,
    min_n = min_n,
    min_q = min_q,
    max_q = max_q,
    step_q = step_q,
    ci = ci,
    ...
  )
  if (isSIRAP(receptors)) {
    return(setAnalysis(receptors, analysis_name, result, validate = validate))
  }
  result
}

#' Summarize stored polygon diversity at one q value
#'
#' Extracts the diversity estimate for a specific diversity order from a stored
#' polygon diversity analysis, producing one table row per polygon and saving
#' that table as a new SIRAP analysis.
#'
#' @param object A SIRAP object.
#' @param q Diversity order to extract, such as `2` for Simpson diversity.
#' @param diversity_analysis Name of the stored analysis produced by
#'   [polygonDiversity()].
#' @param analysis_name Name used to store the summarized table.
#' @param tolerance Numeric tolerance used when matching `q`.
#' @param validate Whether to validate the updated SIRAP object.
#' @return An updated SIRAP object containing the diversity summary analysis.
#' @export
summarizePolygonDiversity <- function(object,
                                      q,
                                      diversity_analysis = "polygon_diversity",
                                      analysis_name = NULL,
                                      tolerance = sqrt(.Machine$double.eps),
                                      validate = TRUE) {
  validateSIRAP(object)
  if (!is.numeric(q) || length(q) != 1L || is.na(q)) {
    rlang::abort("`q` must be one numeric diversity order.")
  }
  if (is.null(analysis_name)) {
    analysis_name <- paste0("polygon_diversity_q", gsub("[^0-9A-Za-z]+", "_", as.character(q)))
  }
  diversity <- getAnalysis(object, diversity_analysis)
  if (!is.list(diversity)) {
    rlang::abort("`diversity_analysis` must be a polygon diversity analysis list.")
  }
  rows <- lapply(names(diversity), function(polygon) {
    curve <- diversity[[polygon]]
    if (is.null(curve)) {
      return(tibble::tibble(polygon_id = polygon, q = q, d = NA_real_, d_sd = NA_real_, d_lower = NA_real_, d_upper = NA_real_, e = NA_real_, e_lower = NA_real_, e_upper = NA_real_))
    }
    dat <- .sirap_diversity_curve_table(curve)
    if (!"q" %in% colnames(dat)) {
      rlang::abort("Diversity curve tables must contain a `q` column.")
    }
    idx <- which(abs(as.numeric(dat$q) - q) <= tolerance)
    if (length(idx) == 0L) {
      return(tibble::tibble(polygon_id = polygon, q = q, d = NA_real_, d_sd = NA_real_, d_lower = NA_real_, d_upper = NA_real_, e = NA_real_, e_lower = NA_real_, e_upper = NA_real_))
    }
    row <- dat[idx[[1L]], , drop = FALSE]
    row$polygon_id <- polygon
    row[, c("polygon_id", setdiff(colnames(row), "polygon_id")), drop = FALSE]
  })
  summary <- dplyr::bind_rows(rows)
  setAnalysis(object, analysis_name, summary, validate = validate)
}

# Return polygon IDs from stored polygons when available, otherwise from the
# receptor table annotation. This keeps empty polygons in summaries when their
# geometry is known.
.sirap_polygon_ids <- function(receptors, polygons = NULL, polygon_id = "polygon_id") {
  if (!is.null(polygons)) {
    polygons <- .sirap_prepare_polygons(polygons)
    return(as.character(polygons$polygon_id))
  }
  if (!polygon_id %in% colnames(receptors)) {
    rlang::abort(paste0(
      "`receptors` must contain `", polygon_id,
      "` when `polygons` are not supplied."
    ))
  }
  ids <- sort(unique(as.character(stats::na.omit(receptors[[polygon_id]]))))
  ids[nzchar(ids)]
}

# Split a receptor table by its primary polygon assignment. Per-polygon logical
# membership flags are intentionally ignored here; they are retained only so
# users can build a different primary polygon column before analysis.
.sirap_polygon_receptor_groups <- function(receptors,
                                           polygon_ids,
                                           polygon_id = "polygon_id",
                                           mode = c("by_polygon", "all_polygons")) {
  mode <- match.arg(mode)
  receptors <- tibble::as_tibble(receptors)
  if (!polygon_id %in% colnames(receptors)) {
    rlang::abort(paste0("`receptors` must contain `", polygon_id, "`."))
  }
  keep_list <- vector("list", length(polygon_ids))
  groups <- stats::setNames(vector("list", length(polygon_ids)), polygon_ids)
  for (i in seq_along(polygon_ids)) {
    keep <- as.character(receptors[[polygon_id]]) == as.character(polygon_ids[[i]])
    keep[is.na(keep)] <- FALSE
    keep_list[[i]] <- keep
    groups[[i]] <- receptors[keep, , drop = FALSE]
  }
  if (identical(mode, "by_polygon")) {
    return(groups)
  }
  if (length(groups) == 0L) {
    return(list(all_polygons = receptors[0, , drop = FALSE]))
  }
  keep_any <- Reduce(`|`, keep_list)
  list(all_polygons = receptors[keep_any, , drop = FALSE])
}

# Build one compact polygon summary table from split receptor rows.
.sirap_polygon_receptor_summary_table <- function(receptors,
                                                  polygon_ids,
                                                  spatial_id,
                                                  polygon_id,
                                                  clone,
                                                  weight,
                                                  summary_columns,
                                                  summary_operations) {
  if (length(polygon_ids) == 0L) {
    return(tibble::tibble(polygon_id = character(), n_receptors = integer()))
  }
  operations <- .sirap_summary_operations(summary_columns, summary_operations)
  missing_summary <- setdiff(summary_columns %||% character(), colnames(receptors))
  if (length(missing_summary) > 0L) {
    rlang::abort(paste0(
      "`summary_columns` were not found in `receptors`: ",
      paste(missing_summary, collapse = ", ")
    ))
  }
  groups <- .sirap_polygon_receptor_groups(
    receptors = receptors,
    polygon_ids = polygon_ids,
    polygon_id = polygon_id,
    mode = "by_polygon"
  )
  spatial_col <- if (!is.null(spatial_id) && spatial_id %in% colnames(receptors)) spatial_id else NULL
  clone_col <- if (!is.null(clone) && clone %in% colnames(receptors)) clone else NULL
  weight_col <- if (!is.null(weight) && weight %in% colnames(receptors)) weight else NULL
  if (!is.null(weight) && is.null(weight_col)) {
    rlang::abort(paste0("`weight` column was not found in `receptors`: ", weight))
  }
  if (!is.null(weight_col)) {
    rlang::warn("Weighted polygon summaries apply weights to numeric count/sum/mean/median/variance/sd summaries. Categorical/logical summaries such as `cat`, `any`, and `all` ignore weights and should be interpreted carefully.")
  }
  rows <- lapply(names(groups), function(id) {
    dat <- groups[[id]]
    weights <- .sirap_summary_weights(dat, weight_col)
    out <- list(
      polygon_id = id,
      n_receptors = if (is.null(weight_col)) nrow(dat) else sum(weights, na.rm = TRUE)
    )
    if (!is.null(spatial_col)) {
      vals <- as.character(dat[[spatial_col]])
      out$n_spatial_ids <- length(unique(vals[!is.na(vals) & nzchar(vals)]))
    }
    if (!is.null(clone_col)) {
      vals <- as.character(dat[[clone_col]])
      keep_vals <- !is.na(vals) & nzchar(vals)
      vals <- vals[keep_vals]
      out$n_clones <- length(unique(vals))
      out$top_clone_fraction <- if (length(vals) == 0L) {
        NA_real_
      } else if (is.null(weight_col)) {
        max(table(vals)) / length(vals)
      } else {
        clone_weights <- stats::aggregate(
          weights[keep_vals],
          by = list(clone_id = vals),
          FUN = sum,
          na.rm = TRUE
        )
        max(clone_weights$x) / sum(clone_weights$x)
      }
    }
    if (length(operations) > 0L) {
      for (i in seq_along(operations)) {
        column <- summary_columns[[i]]
        operation <- operations[[i]]
        out[[paste0(operation, "_", column)]] <- .sirap_apply_summary_operation(dat[[column]], operation, weights = weights)
      }
    }
    tibble::as_tibble(out)
  })
  dplyr::bind_rows(rows)
}

# Validate and recycle summary operations to match summary columns.
.sirap_summary_operations <- function(summary_columns, summary_operations) {
  if (is.null(summary_columns)) {
    if (!is.null(summary_operations)) {
      rlang::abort("`summary_columns` is required when `summary_operations` is supplied.")
    }
    return(character())
  }
  if (is.null(summary_operations)) {
    rlang::abort("`summary_operations` must be supplied when `summary_columns` is supplied.")
  }
  summary_columns <- as.character(summary_columns)
  summary_operations <- as.character(summary_operations)
  if (!length(summary_columns) || any(!nzchar(summary_columns))) {
    rlang::abort("`summary_columns` must contain non-empty column names.")
  }
  if (length(summary_operations) == 1L) {
    summary_operations <- rep(summary_operations, length(summary_columns))
  }
  if (length(summary_operations) != length(summary_columns)) {
    rlang::abort("`summary_operations` must have length 1 or match `summary_columns`.")
  }
  allowed <- c("count", "sum", "mean", "median", "variance", "var", "sd", "min", "max", "any", "all", "cat")
  bad <- setdiff(summary_operations, allowed)
  if (length(bad) > 0L) {
    rlang::abort(paste0(
      "Unsupported `summary_operations`: ",
      paste(bad, collapse = ", "),
      "."
    ))
  }
  summary_operations
}

# Apply one named summary operation to a vector from one polygon.
.sirap_apply_summary_operation <- function(x, operation, weights = NULL) {
  x_nonmissing <- x[!is.na(x)]
  w <- if (is.null(weights)) NULL else weights[!is.na(x)]
  if (identical(operation, "count")) {
    return(if (is.null(w)) length(x_nonmissing) else sum(w, na.rm = TRUE))
  }
  if (identical(operation, "cat")) {
    vals <- sort(unique(as.character(x_nonmissing)))
    vals <- vals[nzchar(vals)]
    return(if (length(vals) == 0L) NA_character_ else paste(vals, collapse = ","))
  }
  if (identical(operation, "any")) {
    return(if (length(x_nonmissing) == 0L) NA else any(as.logical(x_nonmissing)))
  }
  if (identical(operation, "all")) {
    return(if (length(x_nonmissing) == 0L) NA else all(as.logical(x_nonmissing)))
  }
  num <- suppressWarnings(as.numeric(x_nonmissing))
  if (!is.null(w)) {
    w <- w[!is.na(num)]
  }
  num <- num[!is.na(num)]
  if (length(num) == 0L) {
    return(NA_real_)
  }
  if (!is.null(w)) {
    w[is.na(w)] <- 0
    if (sum(w) <= 0) {
      return(NA_real_)
    }
    return(switch(
      operation,
      sum = sum(num * w),
      mean = stats::weighted.mean(num, w),
      median = .sirap_weighted_median(num, w),
      variance = .sirap_weighted_variance(num, w),
      var = .sirap_weighted_variance(num, w),
      sd = sqrt(.sirap_weighted_variance(num, w)),
      min = min(num),
      max = max(num)
    ))
  }
  switch(
    operation,
    sum = sum(num),
    mean = mean(num),
    median = stats::median(num),
    variance = stats::var(num),
    var = stats::var(num),
    sd = stats::sd(num),
    min = min(num),
    max = max(num)
  )
}

# Extract non-negative numeric weights for polygon summaries.
.sirap_summary_weights <- function(dat, weight_col) {
  if (is.null(weight_col)) {
    return(NULL)
  }
  weights <- suppressWarnings(as.numeric(dat[[weight_col]]))
  weights[is.na(weights) | weights < 0] <- 0
  weights
}

# Weighted median used by weighted polygon summaries.
.sirap_weighted_median <- function(x, w) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) >= sum(w) / 2)[[1L]]]
}

# Population-style weighted variance for weighted polygon summaries.
.sirap_weighted_variance <- function(x, w) {
  mu <- stats::weighted.mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

# Apply an Alakazam estimator to polygon-specific receptor tables, returning
# NULL for polygons that do not meet the requested receptor-count threshold.
.sirap_polygon_alakazam <- function(receptors,
                                    polygons,
                                    polygon_id,
                                    mode,
                                    threshold_n,
                                    fun,
                                    ...) {
  groups <- withPolygons(
    receptors = receptors,
    .f = identity,
    polygons = polygons,
    polygon_id = polygon_id,
    mode = mode
  )
  run_one <- function(dat) {
    if (nrow(dat) < threshold_n) {
      return(NULL)
    }
    args <- list(...)
    if ("cell_id" %in% names(args) && is.null(args$cell_id)) {
      args$cell_id <- .sirap_absent_column_name(dat, ".sirap_bulk_cell_id")
    }
    do.call(fun, c(list(data = dat), args))
  }
  if (identical(mode, "all_polygons")) {
    return(run_one(groups))
  }
  lapply(groups, run_one)
}

# Choose a column name guaranteed not to exist in a data frame.
.sirap_absent_column_name <- function(dat, prefix) {
  candidate <- prefix
  i <- 1L
  while (candidate %in% colnames(dat)) {
    candidate <- paste0(prefix, "_", i)
    i <- i + 1L
  }
  candidate
}

# Pass polygon_id only to functions that can accept it. This makes withPolygons
# useful with ordinary one-table functions and with polygon-aware callbacks.
.sirap_call_polygon_function <- function(fun, dat, polygon_id = NULL, pass_polygon_id = FALSE, ...) {
  if (!isTRUE(pass_polygon_id)) {
    return(fun(dat, ...))
  }
  fmls <- names(formals(fun))
  if ("polygon_id" %in% fmls || "..." %in% fmls) {
    return(fun(dat, polygon_id = polygon_id, ...))
  }
  fun(dat, ...)
}

# Return spatial coordinate rows for one stored polygon annotation.
.sirap_stored_polygon_coords <- function(coords, spatial, polygon_value, polygon_id = "polygon_id") {
  metadata <- spatial[[]]
  keep <- .sirap_polygon_membership_vector(metadata, rownames(metadata), polygon_value, polygon_id = polygon_id)
  coord_keep <- keep[as.character(coords$cell)]
  coord_keep[is.na(coord_keep)] <- FALSE
  coords[coord_keep, , drop = FALSE]
}

# Return receptor rows for one stored polygon annotation.
.sirap_receptors_in_polygon <- function(receptors, polygon_value, polygon_id = "polygon_id") {
  keep <- .sirap_polygon_membership_vector(receptors, seq_len(nrow(receptors)), polygon_value, polygon_id = polygon_id)
  receptors[keep, , drop = FALSE]
}

# Read a stored primary polygon membership vector. Per-polygon flags are not
# used by analyses; they are available for users to derive a different
# `polygon_id` column if they want overlapping/nested behavior.
.sirap_polygon_membership_vector <- function(dat, ids, polygon_value, polygon_id = "polygon_id") {
  if (!polygon_id %in% colnames(dat)) {
    rlang::abort("Stored polygon annotations were not found. Re-run `addPolygons()` on the SIRAP object.")
  }
  out <- as.character(dat[[polygon_id]]) == as.character(polygon_value)
  out[is.na(out)] <- FALSE
  stats::setNames(out, as.character(ids))
}

# Extract the tibble of diversity estimates from Alakazam DiversityCurve objects
# and simple data-frame/list fallbacks.
.sirap_diversity_curve_table <- function(curve) {
  if (is.data.frame(curve)) {
    return(tibble::as_tibble(curve))
  }
  if (methods::is(curve, "DiversityCurve")) {
    return(tibble::as_tibble(methods::slot(curve, "diversity")))
  }
  if (is.list(curve) && is.data.frame(curve$diversity)) {
    return(tibble::as_tibble(curve$diversity))
  }
  rlang::abort("Could not find a diversity table in the stored polygon diversity result.")
}

# Find which spatial IDs fall inside each polygon and return one row per
# spatial ID/polygon hit. The helper keeps membership calculation shared across
# Seurat annotation, receptor annotation, and subsetting.
.sirap_polygon_membership <- function(object,
                                      polygon,
                                      images,
                                      image_scale,
                                      swap_coords,
                                      polygon_val,
                                      x,
                                      y,
                                      polygon_id,
                                      include_boundary) {
  polygon <- .as_data_frame(polygon)
  .check_required_cols(polygon, c(x, y), "polygon")
  coords <- .sirap_spatial_coordinate_table(
    object,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  if (nrow(coords) == 0L) {
    rlang::abort("No finite spatial coordinates were found in `object`.")
  }
  .sirap_polygon_membership_from_coords(
    coords = coords,
    polygon = polygon,
    polygon_val = polygon_val,
    x = x,
    y = y,
    polygon_id = polygon_id,
    include_boundary = include_boundary
  )
}

.sirap_polygon_membership_from_coords <- function(coords,
                                                  polygon,
                                                  polygon_val,
                                                  x,
                                                  y,
                                                  polygon_id,
                                                  include_boundary) {
  polygon <- .as_data_frame(polygon)
  .check_required_cols(polygon, c(x, y), "polygon")
  .check_required_cols(coords, c("cell", "x", "y"), "coords")
  if (!"image" %in% colnames(coords)) {
    coords$image <- NA_character_
  }
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  if (nrow(coords) == 0L) {
    rlang::abort("No finite coordinates were found.")
  }
  polygon_specs <- .sirap_polygon_specs(polygon, polygon_id = polygon_val, polygon_col = polygon_id)
  polygon_values <- vapply(polygon_specs, function(spec) {
    as.character(if (isTRUE(spec$add_polygon_col)) spec$id else "polygon_1")
  }, character(1))
  polygon_columns <- .sirap_polygon_flag_columns(polygon_values)
  membership <- lapply(seq_along(polygon_specs), function(i) {
    spec <- polygon_specs[[i]]
    vertices <- .sirap_select_polygon_vertices(
      polygon = polygon,
      polygon_id = spec$id,
      polygon_col = polygon_id,
      x_col = x,
      y_col = y
    )
    inside <- .sirap_points_in_polygon(
      px = coords$x,
      py = coords$y,
      vx = vertices$x,
      vy = vertices$y,
      include_boundary = include_boundary
    )
    if (!any(inside)) {
      return(data.frame())
    }
    polygon_value <- polygon_values[[i]]
    data.frame(
      cell = as.character(coords$cell[inside]),
      polygon = rep(polygon_value, sum(inside)),
      polygon_area = rep(.sirap_polygon_area(vertices$x, vertices$y), sum(inside)),
      polygon_order = rep(i, sum(inside)),
      polygon_column = rep(polygon_columns[[i]], sum(inside)),
      image = as.character(coords$image[inside]),
      x = as.numeric(coords$x[inside]),
      y = as.numeric(coords$y[inside]),
      stringsAsFactors = FALSE
    )
  })
  out <- .sirap_bind_rows(membership)
  attr(out, "polygon_columns") <- polygon_columns
  out
}

# Collapse one-row-per-hit polygon membership into a named lookup vector from
# spatial ID to a primary polygon ID. Overlapping IDs resolve to the smallest
# polygon by area, then to polygon input order for exact ties.
.sirap_polygon_lookup <- function(membership) {
  if (nrow(membership) == 0L) {
    return(stats::setNames(character(), character()))
  }
  membership <- membership[order(membership$cell, membership$polygon_area, membership$polygon_order), , drop = FALSE]
  membership <- membership[!duplicated(membership$cell), , drop = FALSE]
  stats::setNames(as.character(membership$polygon), as.character(membership$cell))
}

# Return one logical flag column per polygon ID, keyed by spatial ID.
.sirap_polygon_flag_lookup <- function(membership) {
  all_polygon_columns <- attr(membership, "polygon_columns", exact = TRUE)
  if (nrow(membership) == 0L) {
    out <- data.frame(cell = character(), stringsAsFactors = FALSE)
    for (column in all_polygon_columns %||% character()) {
      out[[column]] <- logical()
    }
    return(out)
  }
  cells <- sort(unique(as.character(membership$cell)))
  polygon_columns <- all_polygon_columns %||% unique(as.character(membership$polygon_column)[order(membership$polygon_order)])
  out <- data.frame(cell = cells, stringsAsFactors = FALSE)
  for (column in polygon_columns) {
    out[[column]] <- cells %in% membership$cell[as.character(membership$polygon_column) == column]
  }
  out
}

.sirap_append_polygon_flags <- function(out, spatial_ids, polygon_flags) {
  flag_cols <- setdiff(colnames(polygon_flags), "cell")
  if (length(flag_cols) == 0L) {
    return(out)
  }
  row_idx <- match(as.character(spatial_ids), polygon_flags$cell)
  for (column in flag_cols) {
    out[[column]] <- !is.na(row_idx) & polygon_flags[[column]][row_idx]
  }
  out
}

.sirap_polygon_flag_columns <- function(polygon_id) {
  make.names(paste0("polygon_", as.character(polygon_id)), unique = TRUE)
}

.sirap_polygon_area <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) > 1L && .sirap_same_xy(data.frame(x = x[[1L]], y = y[[1L]]), data.frame(x = x[[length(x)]], y = y[[length(y)]]))) {
    x <- x[-length(x)]
    y <- y[-length(y)]
  }
  j <- c(seq_along(x)[-1L], 1L)
  abs(sum(x * y[j] - x[j] * y)) / 2
}

# Repair older Seurat image objects before subsetting or metadata updates so
# polygon helpers do not fail on images with legacy slot layouts.
.sirap_update_seurat_image_slots <- function(object) {
  update_slots <- utils::getFromNamespace("UpdateSlots", "SeuratObject")
  for (image_name in SeuratObject::Images(object)) {
    image <- object[[image_name]]
    has_missing_slots <- any(vapply(methods::slotNames(image), function(slot_name) {
      tryCatch(
        {
          methods::slot(image, slot_name)
          FALSE
        },
        error = function(...) TRUE
      )
    }, logical(1)))
    if (has_missing_slots) {
      methods::slot(object, "images")[[image_name]] <- suppressWarnings(update_slots(image))
    }
  }
  object
}

.sirap_select_polygon_vertices <- function(polygon,
                                           polygon_id,
                                           polygon_col,
                                           x_col,
                                           y_col) {
  if (!is.null(polygon_id) && polygon_col %in% colnames(polygon)) {
    polygon <- polygon[as.character(polygon[[polygon_col]]) == as.character(polygon_id), , drop = FALSE]
  } else if (is.null(polygon_id) && polygon_col %in% colnames(polygon)) {
    polygon_ids <- unique(stats::na.omit(polygon[[polygon_col]]))
    if (length(polygon_ids) > 1L) {
      rlang::abort("`polygon` contains multiple polygons; pass `polygon_id`.")
    }
  }
  if (nrow(polygon) < 3L) {
    rlang::abort("`polygon` must contain at least three rows.")
  }
  out <- data.frame(
    x = as.numeric(polygon[[x_col]]),
    y = as.numeric(polygon[[y_col]]),
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(out$x) | !is.finite(out$y))) {
    rlang::abort("`polygon` contains non-finite coordinates.")
  }
  closed <- .sirap_same_xy(out[1, , drop = FALSE], out[nrow(out), , drop = FALSE])
  if (!closed) {
    out <- rbind(out, out[1, , drop = FALSE])
  }
  out
}

.sirap_polygon_specs <- function(polygon, polygon_id, polygon_col) {
  has_polygon_col <- polygon_col %in% colnames(polygon)
  if (!has_polygon_col) {
    if (!is.null(polygon_id)) {
      rlang::abort("`polygon_id` can only be used when `polygon_col` exists in `polygon`.")
    }
    return(list(list(id = NULL, add_polygon_col = FALSE)))
  }

  ids <- if (is.null(polygon_id)) {
    unique(stats::na.omit(polygon[[polygon_col]]))
  } else {
    polygon_id
  }
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0L) {
    rlang::abort("No polygon IDs were found in `polygon`.")
  }
  missing_ids <- setdiff(as.character(ids), as.character(polygon[[polygon_col]]))
  if (length(missing_ids) > 0L) {
    rlang::abort(paste0(
      "`polygon_id` values were not found in `polygon`: ",
      paste(missing_ids, collapse = ", ")
    ))
  }
  lapply(ids, function(id) list(id = id, add_polygon_col = TRUE))
}

.sirap_rbind <- function(results, name) {
  tables <- lapply(results, `[[`, name)
  if (length(tables) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  out
}

.sirap_v_genes <- function(v_calls) {
  if (length(v_calls) == 0L) {
    return(character())
  }
  genes <- unlist(strsplit(as.character(v_calls), "[,;]"), use.names = FALSE)
  genes <- trimws(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  sort(unique(genes))
}

.sirap_v_gene_family <- function(v_genes) {
  sub("\\*.*$", "", v_genes)
}

.sirap_v_gene_superfamily <- function(v_gene_families) {
  sub("^([A-Za-z]+V[0-9]+).*$", "\\1", v_gene_families)
}

.sirap_points_in_polygon <- function(px,
                                     py,
                                     vx,
                                     vy,
                                     include_boundary = TRUE,
                                     tol = 1e-3) {
  # `tol` is the boundary band (in the plot's pixel coordinate units). The
  # previous default, sqrt(.Machine$double.eps) ~ 1.5e-8, was far too tight for
  # pixel coordinates in the hundreds-thousands: spots sitting exactly on a
  # polygon boundary/vertex could fall ~1e-6 outside when the polygon vertices
  # come from a rounded (saved) CSV, and get dropped. 1e-3 px is sub-pixel
  # (spot pitch is tens of px), so it absorbs coordinate round-trip error
  # without including any genuinely-outside spot.
  px <- as.numeric(px)
  py <- as.numeric(py)
  vx <- as.numeric(vx)
  vy <- as.numeric(vy)
  if (length(vx) != length(vy) || length(vx) < 4L) {
    rlang::abort("Polygon coordinates must contain at least four paired vertices.")
  }
  if (.sirap_same_xy(data.frame(x = vx[[1]], y = vy[[1]]), data.frame(x = vx[[length(vx)]], y = vy[[length(vy)]]))) {
    vx <- vx[-length(vx)]
    vy <- vy[-length(vy)]
  }

  inside <- rep(FALSE, length(px))
  j <- length(vx)
  for (i in seq_along(vx)) {
    crosses <- ((vy[[i]] > py) != (vy[[j]] > py)) &
      (px < (vx[[j]] - vx[[i]]) * (py - vy[[i]]) / (vy[[j]] - vy[[i]]) + vx[[i]])
    crosses[is.na(crosses)] <- FALSE
    inside <- xor(inside, crosses)
    j <- i
  }

  boundary <- rep(FALSE, length(px))
  j <- length(vx)
  for (i in seq_along(vx)) {
    boundary <- boundary | .sirap_points_on_segment(px, py, vx[[j]], vy[[j]], vx[[i]], vy[[i]], tol = tol)
    j <- i
  }
  if (isTRUE(include_boundary)) {
    inside <- inside | boundary
  } else {
    inside <- inside & !boundary
  }
  inside
}

.sirap_points_on_segment <- function(px, py, x1, y1, x2, y2, tol) {
  dx <- x2 - x1
  dy <- y2 - y1
  cross <- (px - x1) * dy - (py - y1) * dx
  scale <- pmax(abs(dx), abs(dy), 1)
  on_line <- abs(cross) <= tol * scale
  within_x <- px >= pmin(x1, x2) - tol & px <= pmax(x1, x2) + tol
  within_y <- py >= pmin(y1, y2) - tol & py <= pmax(y1, y2) + tol
  on_line & within_x & within_y
}

.sirap_same_xy <- function(a, b) {
  isTRUE(all.equal(as.numeric(a$x), as.numeric(b$x), tolerance = sqrt(.Machine$double.eps))) &&
    isTRUE(all.equal(as.numeric(a$y), as.numeric(b$y), tolerance = sqrt(.Machine$double.eps)))
}

.sirap_empty_polygon_vertices <- function() {
  data.frame(
    spatial_id = character(),
    x = numeric(),
    y = numeric(),
    .sirap_plot_y = numeric(),
    num = integer(),
    polygon = character(),
    image = character(),
    stringsAsFactors = FALSE
  )
}

.sirap_public_polygon_vertices <- function(points) {
  if (!"num" %in% colnames(points) && "polygon" %in% colnames(points)) {
    points$num <- match(as.character(points$polygon), unique(as.character(points$polygon)))
    points <- points[, c(setdiff(colnames(points), "polygon"), "polygon"), drop = FALSE]
  }
  points[, setdiff(colnames(points), ".sirap_plot_y"), drop = FALSE]
}

.sirap_assay_layers <- function(object, assay = NULL) {
  if (is.null(assay) || !nzchar(assay)) {
    return(c("counts", "data", "scale.data"))
  }
  layers <- tryCatch(
    SeuratObject::Layers(object[[assay]]),
    error = function(e) character()
  )
  layers <- unique(as.character(layers))
  layers <- layers[nzchar(layers)]
  if (length(layers) == 0L) {
    layers <- c("counts", "data", "scale.data")
  }
  layers
}

.sirap_polygon_plot_data <- function(polygon,
                                     polygon_id,
                                     polygon_col,
                                     x_col,
                                     y_col,
                                     keep_cols = character()) {
  specs <- .sirap_polygon_specs(polygon, polygon_id = polygon_id, polygon_col = polygon_col)
  keep_cols <- unique(keep_cols[keep_cols %in% colnames(polygon)])
  out <- lapply(specs, function(spec) {
    polygon_rows <- if (isTRUE(spec$add_polygon_col)) {
      polygon[as.character(polygon[[polygon_col]]) == as.character(spec$id), , drop = FALSE]
    } else {
      polygon
    }
    vertices <- .sirap_select_polygon_vertices(
      polygon = polygon,
      polygon_id = spec$id,
      polygon_col = polygon_col,
      x_col = x_col,
      y_col = y_col
    )
    vertices$polygon <- if (isTRUE(spec$add_polygon_col)) {
      as.character(spec$id)
    } else {
      "polygon_1"
    }
    for (keep_col in setdiff(keep_cols, c(x_col, y_col, "polygon"))) {
      values <- unique(polygon_rows[[keep_col]])
      values <- values[!is.na(values)]
      vertices[[keep_col]] <- if (length(values) == 0L) NA_character_ else values[[1L]]
    }
    vertices
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out$.sirap_y <- -out$y
  out
}

.sirap_polygon_vertices_table <- function(polygons) {
  .sirap_validate_polygons(polygons)
  if (nrow(polygons) == 0L) {
    return(data.frame(x = numeric(), y = numeric(), polygon = character(), stringsAsFactors = FALSE))
  }
  out <- lapply(seq_len(nrow(polygons)), function(i) {
    vertices <- as.data.frame(polygons$vertices[[i]], stringsAsFactors = FALSE)
    vertices$polygon <- as.character(polygons$polygon_id[[i]])
    extra <- setdiff(names(polygons), c("polygon_id", "vertices"))
    for (column in extra) {
      if (!column %in% names(vertices)) vertices[[column]] <- polygons[[column]][[i]]
    }
    vertices
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

.sirap_join_polygon_column_fill_data <- function(polygon_data, fill_column) {
  if (!fill_column %in% names(polygon_data)) {
    rlang::abort(paste0(
      "`fill_column` was not found in polygons. `plotSpatialPolygons()` looks ",
      "for `fill_column` in the SIRAP object's polygons slot, or in the supplied ",
      "`polygon` data frame. To color polygons by a stored analysis instead, ",
      "specify `fill_analysis`; the selected analysis must contain a `polygon_id` ",
      "or `polygon` column that matches the plotted polygons and a `",
      fill_column, "` column."
    ))
  }
  polygon_data$.sirap_fill_value <- polygon_data[[fill_column]]
  polygon_data
}

.sirap_join_polygon_fill_data <- function(object,
                                          polygon_data,
                                          analysis_name,
                                          fill_column) {
  analysis <- getAnalysis(object, analysis_name)
  if (is.list(analysis) && !is.data.frame(analysis) &&
      is.data.frame(analysis$dominance)) {
    analysis <- analysis$dominance
  }
  if (!is.data.frame(analysis)) {
    rlang::abort("The selected polygon analysis must be a data frame or contain a `dominance` data frame.")
  }
  key <- if ("polygon_id" %in% names(analysis)) "polygon_id" else "polygon"
  if (!all(c(key, fill_column) %in% names(analysis))) {
    rlang::abort(paste0(
      "The selected polygon analysis must contain `", key, "` and `",
      fill_column, "`."
    ))
  }
  values <- data.frame(
    polygon = as.character(analysis[[key]]),
    .sirap_fill_value = analysis[[fill_column]],
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(values$polygon)) {
    rlang::abort("The selected polygon analysis must contain one row per polygon.")
  }
  polygon_data$polygon <- as.character(polygon_data$polygon)
  matched <- match(polygon_data$polygon, values$polygon)
  polygon_data$.sirap_fill_value <- values$.sirap_fill_value[matched]
  polygon_data
}

.sirap_region_from_dominance <- function(dominance, polygon_val = NULL) {
  if (!inherits(dominance, "sirap_spatial_clone_dominance")) {
    rlang::abort("`dominance` must be an object returned by `spatialCloneDominance()`.")
  }
  if (is.null(dominance$region_receptors) || is.null(dominance$polygons)) {
    rlang::abort("`dominance` must contain `region_receptors` and `polygons` from `spatialCloneDominance()`.")
  }
  region_receptors <- .as_data_frame(dominance$region_receptors)
  polygon_plot_data <- .as_data_frame(dominance$polygons)
  if (!is.null(polygon_val)) {
    region_receptors <- region_receptors[
      as.character(region_receptors$polygon) %in% as.character(polygon_val),
      ,
      drop = FALSE
    ]
    polygon_plot_data <- polygon_plot_data[
      as.character(polygon_plot_data$polygon) %in% as.character(polygon_val),
      ,
      drop = FALSE
    ]
  }
  list(
    receptors = region_receptors,
    polygons = polygon_plot_data
  )
}

.sirap_drop_plot_coordinate_cols <- function(dat) {
  dat[, setdiff(colnames(dat), c("x", "y", ".sirap_y", "image")), drop = FALSE]
}

.sirap_add_v_gene_columns <- function(dat, v_call_col = "v_call") {
  if (nrow(dat) == 0L) {
    if (!"v_gene" %in% colnames(dat)) {
      dat$v_gene <- character()
    }
    if (!"v_gene_family" %in% colnames(dat)) {
      dat$v_gene_family <- character()
    }
    if (!"v_gene_superfamily" %in% colnames(dat)) {
      dat$v_gene_superfamily <- character()
    }
    return(dat)
  }
  if (!"v_gene" %in% colnames(dat)) {
    dat$v_gene <- if (v_call_col %in% colnames(dat)) {
      vapply(as.character(dat[[v_call_col]]), function(x) paste(.sirap_v_genes(x), collapse = ","), character(1))
    } else {
      rep(NA_character_, nrow(dat))
    }
    dat$v_gene[dat$v_gene == ""] <- NA_character_
  }
  if (!"v_gene_family" %in% colnames(dat)) {
    dat$v_gene_family <- if (v_call_col %in% colnames(dat)) {
      vapply(as.character(dat[[v_call_col]]), function(v_calls) {
        families <- unique(stats::na.omit(.sirap_v_gene_family(.sirap_v_genes(v_calls))))
        families <- sort(families)
        if (length(families) == 0L) NA_character_ else paste(families, collapse = ",")
      }, character(1))
    } else {
      .sirap_v_gene_family(dat$v_gene)
    }
  }
  if (!"v_gene_superfamily" %in% colnames(dat)) {
    dat$v_gene_superfamily <- if (v_call_col %in% colnames(dat)) {
      vapply(as.character(dat[[v_call_col]]), function(v_calls) {
        superfamilies <- .sirap_v_genes(v_calls)
        superfamilies <- .sirap_v_gene_family(superfamilies)
        superfamilies <- unique(stats::na.omit(.sirap_v_gene_superfamily(superfamilies)))
        superfamilies <- sort(superfamilies)
        if (length(superfamilies) == 0L) NA_character_ else paste(superfamilies, collapse = ",")
      }, character(1))
    } else {
      .sirap_v_gene_superfamily(dat$v_gene_family)
    }
  }
  dat
}

.sirap_bind_rows <- function(dots) {
  dots <- Filter(function(x) !is.null(x) && nrow(x) > 0L, dots)
  if (length(dots) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, dots)
  rownames(out) <- NULL
  out
}

.sirap_lump_values <- function(x, top_n = NULL, other_label = "Other") {
  values <- as.character(x)
  values[is.na(values) | !nzchar(values)] <- NA_character_
  if (is.null(top_n)) {
    return(factor(values, levels = sort(unique(stats::na.omit(values)))))
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
    rlang::abort("`top_n` must be NULL or a positive number.")
  }
  counts <- sort(table(values, useNA = "no"), decreasing = TRUE)
  keep <- names(counts)[seq_len(min(as.integer(top_n), length(counts)))]
  values[!is.na(values) & !values %in% keep] <- other_label
  factor(values, levels = c(keep, if (any(values == other_label, na.rm = TRUE)) other_label))
}

.sirap_discrete_palette_values <- function(levels, other_label = "Other", other_color = "#D3D3D3") {
  levels <- as.character(levels)
  levels <- levels[!is.na(levels) & nzchar(levels)]
  has_other <- other_label %in% levels
  main_levels <- setdiff(levels, other_label)
  .require_namespace("pals")
  values <- if (length(main_levels) == 0L) {
    character()
  } else if (length(main_levels) <= 32L) {
    pals::glasbey(length(main_levels))
  } else {
    grDevices::colorRampPalette(pals::glasbey(32L))(length(main_levels))
  }
  names(values) <- main_levels
  if (has_other) {
    values <- c(values, stats::setNames(other_color, other_label))
  }
  values[levels]
}

.sirap_clone_density_3d_data <- function(region_receptors,
                                         polygons,
                                         image,
                                         spatial_coords,
                                         clone,
                                         weight,
                                         top_n,
                                         other_label,
                                         other_color,
                                         plot_other,
                                         bandwidth,
                                         grid_n,
                                         plot_type,
                                         hist_bins,
                                         zoom,
                                         height_scale,
                                         z_offset,
                                         density_floor) {
  .check_required_cols(region_receptors, c(clone, "x", "y", ".sirap_y", "polygon"), "dominance$region_receptors")
  .check_required_cols(polygons, c("x", "y", ".sirap_y", "polygon"), "dominance$polygons")
  if (!is.null(weight)) {
    .check_required_cols(region_receptors, weight, "dominance$region_receptors")
  }
  if (!is.numeric(height_scale) || length(height_scale) != 1L || !is.finite(height_scale) || height_scale <= 0) {
    rlang::abort("`height_scale` must be a positive number.")
  }
  if (!is.numeric(z_offset) || length(z_offset) != 1L || !is.finite(z_offset) || z_offset < 0) {
    rlang::abort("`z_offset` must be a non-negative number.")
  }
  if (!is.numeric(density_floor) || length(density_floor) != 1L || !is.finite(density_floor) || density_floor < 0 || density_floor >= 1) {
    rlang::abort("`density_floor` must be a number greater than or equal to 0 and less than 1.")
  }

  region_receptors <- region_receptors[
    is.finite(region_receptors$x) & is.finite(region_receptors$y) & is.finite(region_receptors$.sirap_y),
    ,
    drop = FALSE
  ]
  if (nrow(region_receptors) == 0L) {
    rlang::abort("No finite regional receptor coordinates are available.")
  }
  region_receptors$.sirap_density_group <- .sirap_lump_values(region_receptors[[clone]], top_n = top_n, other_label = other_label)
  if (!isTRUE(plot_other)) {
    region_receptors <- region_receptors[as.character(region_receptors$.sirap_density_group) != other_label, , drop = FALSE]
  }
  if (nrow(region_receptors) == 0L) {
    rlang::abort("No receptor rows remain to plot after applying `top_n` and `plot_other`.")
  }

  weights <- if (is.null(weight)) {
    rep(1, nrow(region_receptors))
  } else {
    as.numeric(region_receptors[[weight]])
  }
  weights[!is.finite(weights) | weights < 0] <- 0
  if (!any(weights > 0)) {
    rlang::abort("All receptor weights are zero or non-finite.")
  }
  region_receptors$.sirap_density_weight <- weights

  image_dim <- if (is.null(image)) {
    .sirap_image_dim(image, spatial_coords = spatial_coords)
  } else {
    .sirap_image_dim(image)
  }
  region_xlim <- .sirap_density_plot_range(c(region_receptors$x, polygons$x), c(0, image_dim$width))
  region_ylim <- .sirap_density_plot_range(c(region_receptors$.sirap_y, polygons$.sirap_y), c(-image_dim$height, 0))
  object_xlim <- .sirap_density_plot_range(spatial_coords$x, c(0, image_dim$width))
  object_ylim <- .sirap_density_plot_range(spatial_coords$.sirap_y, c(-image_dim$height, 0))
  image_xlim <- if (isTRUE(zoom)) region_xlim else object_xlim
  image_ylim <- if (isTRUE(zoom)) region_ylim else object_ylim
  image_extent <- if (is.null(image)) {
    list(rows = 1:2, cols = 1:2, xlim = image_xlim, ylim = image_ylim)
  } else {
    .sirap_image_texture_extent(image_dim, image_xlim, image_ylim)
  }
  span <- max(diff(region_xlim), diff(region_ylim))
  z_base <- z_offset * span
  max_height <- height_scale * span
  levels <- levels(region_receptors$.sirap_density_group)
  levels <- levels[levels %in% as.character(region_receptors$.sirap_density_group)]
  color_values <- .sirap_discrete_palette_values(levels, other_label = other_label, other_color = other_color)

  if (identical(plot_type, "histogram")) {
    surfaces <- list()
    bandwidth <- NULL
    histograms <- .sirap_clone_histogram_3d_bars(
      region_bcrs = region_receptors,
      levels = levels,
      color_values = color_values,
      xlim = region_xlim,
      ylim = region_ylim,
      hist_bins = hist_bins,
      z_base = z_base,
      max_height = max_height
    )
  } else {
    histograms <- list()
    grid <- .sirap_density_grid(region_xlim, region_ylim, grid_n)
    bandwidth <- .sirap_density_bandwidth(region_receptors$x, region_receptors$.sirap_y, bandwidth)
    raw_surfaces <- lapply(levels, function(level) {
      dat <- region_receptors[as.character(region_receptors$.sirap_density_group) == level, , drop = FALSE]
      .sirap_weighted_gaussian_grid(
        x = dat$x,
        y = dat$.sirap_y,
        weights = dat$.sirap_density_weight,
        x_grid = grid$x,
        y_grid = grid$y,
        bandwidth = bandwidth
      )
    })
    names(raw_surfaces) <- levels
    max_density <- max(vapply(raw_surfaces, function(z) max(z, na.rm = TRUE), numeric(1)), na.rm = TRUE)
    surfaces <- lapply(levels, function(level) {
      z <- raw_surfaces[[level]]
      surface_max <- max(z, na.rm = TRUE)
      if (is.finite(surface_max) && surface_max > 0 && density_floor > 0) {
        z[z < surface_max * density_floor] <- NA_real_
      }
      if (is.finite(max_density) && max_density > 0) {
        z <- z / max_density * max_height + z_base
      } else {
        z[] <- z_base
      }
      list(
        level = level,
        x = grid$x,
        y = grid$y,
        z = z,
        color = unname(color_values[level])
      )
    })
    names(surfaces) <- levels
  }

  list(
    image = image,
    image_extent = image_extent,
    region_xlim = region_xlim,
    region_ylim = region_ylim,
    plot_type = plot_type,
    surfaces = surfaces,
    histograms = histograms,
    receptors = region_receptors,
    polygons = polygons,
    colors = color_values,
    z_base = z_base,
    z_top = z_base + max_height,
    bandwidth = bandwidth,
    clone = clone
  )
}

.sirap_draw_clone_density_3d <- function(dat,
                                         density_alpha,
                                         image_alpha,
                                         show_points,
                                         point_size,
                                         point_alpha,
                                         show_polygon,
                                         polygon_color,
                                         polygon_lwd,
                                         label_density_peaks,
                                         label_top_n,
                                         label_cex,
                                         show_axes,
                                         texture_dir,
                                         as_widget) {
  if (isTRUE(as_widget)) {
    old_rgl_use_null <- getOption("rgl.useNULL")
    options(rgl.useNULL = TRUE)
    on.exit(options(rgl.useNULL = old_rgl_use_null), add = TRUE)
  }
  rgl::open3d(useNULL = isTRUE(as_widget))
  rgl::clear3d()

  texture <- .sirap_image_texture_file(dat$image, dat$image_extent, alpha = image_alpha, texture_dir = texture_dir)
  x_plane <- dat$image_extent$xlim
  y_plane <- dat$image_extent$ylim
  z_plane <- matrix(0, nrow = 2, ncol = 2)
  texture_s <- matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2)
  texture_t <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
  rgl::surface3d(
    x = x_plane,
    y = y_plane,
    z = z_plane,
    texture_s = texture_s,
    texture_t = texture_t,
    texture = texture,
    color = "white",
    alpha = image_alpha,
    lit = FALSE,
    back = "filled"
  )

  for (surface in dat$surfaces) {
    rgl::surface3d(
      x = surface$x,
      y = surface$y,
      z = surface$z,
      color = surface$color,
      alpha = density_alpha,
      lit = FALSE,
      back = "filled"
    )
  }
  if (identical(dat$plot_type, "histogram")) {
    for (histogram in dat$histograms) {
      .sirap_draw_histogram_bars_3d(histogram, alpha = density_alpha)
    }
  }

  if (isTRUE(show_polygon) && nrow(dat$polygons) > 0L) {
    for (polygon_id in unique(dat$polygons$polygon)) {
      polygon <- dat$polygons[dat$polygons$polygon == polygon_id, , drop = FALSE]
      rgl::lines3d(
        polygon$x,
        polygon$.sirap_y,
        rep(dat$z_base * 0.35, nrow(polygon)),
        color = polygon_color,
        lwd = polygon_lwd
      )
    }
  }

  if (isTRUE(label_density_peaks)) {
    labeled_items <- if (identical(dat$plot_type, "histogram")) {
      utils::head(dat$histograms, label_top_n)
    } else {
      utils::head(dat$surfaces, label_top_n)
    }
    for (item in labeled_items) {
      peak <- if (identical(dat$plot_type, "histogram")) {
        .sirap_histogram_bar_peak(item)
      } else {
        .sirap_density_surface_peak(item)
      }
      if (!is.null(peak)) {
        rgl::text3d(
          peak$x,
          peak$y,
          peak$z + (dat$z_top - dat$z_base) * 0.05,
          texts = item$level,
          color = item$color,
          cex = label_cex,
          font = 2
        )
      }
    }
  }

  if (isTRUE(show_points)) {
    point_colors <- unname(dat$colors[as.character(dat$receptors$.sirap_density_group)])
    point_colors[is.na(point_colors)] <- "#222222"
    rgl::points3d(
      dat$receptors$x,
      dat$receptors$.sirap_y,
      rep(dat$z_base * 0.65, nrow(dat$receptors)),
      color = grDevices::adjustcolor(point_colors, alpha.f = point_alpha),
      size = point_size
    )
  }

  if (isTRUE(show_axes)) {
    xlim <- range(c(dat$image_extent$xlim, dat$region_xlim), na.rm = TRUE)
    ylim <- range(c(dat$image_extent$ylim, dat$region_ylim), na.rm = TRUE)
    rgl::bbox3d(
      xat = pretty(xlim, n = 4),
      yat = pretty(ylim, n = 4),
      zat = pretty(c(0, dat$z_top), n = 3),
      color = "#666666",
      emission = "#333333"
    )
    rgl::decorate3d(xlab = "x", ylab = "y", zlab = "clone density")
  }
  rgl::aspect3d(1, 1, 0.55)
  rgl::view3d(theta = 40, phi = 25, zoom = 0.72)

  if (isTRUE(as_widget)) {
    rgl::rglwidget()
  } else {
    rgl::cur3d()
  }
}

.sirap_image_dim <- function(image, spatial_coords = NULL) {
  if (is.null(image)) {
    if (!is.null(spatial_coords) && nrow(spatial_coords) > 0L) {
      xlim <- range(spatial_coords$x, na.rm = TRUE)
      ylim <- range(spatial_coords$.sirap_y, na.rm = TRUE)
      width <- ceiling(max(abs(xlim), diff(xlim), 1, na.rm = TRUE))
      height <- ceiling(max(abs(ylim), diff(ylim), 1, na.rm = TRUE))
      return(list(height = height, width = width))
    }
    return(list(height = 1, width = 1))
  }
  img <- SeuratObject::GetImage(image, mode = "raw")
  dims <- dim(img)
  list(height = dims[[1]], width = dims[[2]])
}

.sirap_density_surface_peak <- function(surface) {
  if (is.null(surface$z) || !any(is.finite(surface$z))) {
    return(NULL)
  }
  idx <- arrayInd(which.max(surface$z), dim(surface$z))
  list(
    x = surface$x[[idx[[1]]]],
    y = surface$y[[idx[[2]]]],
    z = surface$z[idx[[1]], idx[[2]]]
  )
}

.sirap_histogram_bar_peak <- function(histogram) {
  bars <- histogram$bars
  if (is.null(bars) || nrow(bars) == 0L || !any(is.finite(bars$z1))) {
    return(NULL)
  }
  i <- which.max(bars$z1)
  list(
    x = (bars$x0[[i]] + bars$x1[[i]]) / 2,
    y = (bars$y0[[i]] + bars$y1[[i]]) / 2,
    z = bars$z1[[i]]
  )
}

.sirap_density_plot_range <- function(x, fallback) {
  out <- range(x, na.rm = TRUE)
  if (!all(is.finite(out))) {
    return(fallback)
  }
  if (out[[1]] == out[[2]]) {
    out <- out + c(-1, 1)
  }
  pad <- max(diff(out) * 0.08, 1)
  out <- out + c(-pad, pad)
  c(max(min(fallback), out[[1]]), min(max(fallback), out[[2]]))
}

.sirap_density_grid <- function(xlim, ylim, grid_n) {
  if (!is.numeric(grid_n) || length(grid_n) > 2L || any(!is.finite(grid_n)) || any(grid_n < 3)) {
    rlang::abort("`grid_n` must be one or two numbers greater than or equal to 3.")
  }
  grid_n <- as.integer(grid_n)
  if (length(grid_n) == 1L) {
    grid_n <- rep(grid_n, 2L)
  }
  list(
    x = seq(xlim[[1]], xlim[[2]], length.out = grid_n[[1]]),
    y = seq(ylim[[1]], ylim[[2]], length.out = grid_n[[2]])
  )
}

.sirap_density_bandwidth <- function(x, y, bandwidth) {
  if (!is.null(bandwidth)) {
    if (!is.numeric(bandwidth) || length(bandwidth) > 2L || any(!is.finite(bandwidth)) || any(bandwidth <= 0)) {
      rlang::abort("`bandwidth` must be NULL or one or two positive numbers.")
    }
    bandwidth <- as.numeric(bandwidth)
    if (length(bandwidth) == 1L) {
      bandwidth <- rep(bandwidth, 2L)
    }
    return(bandwidth)
  }
  estimate <- function(v) {
    v <- v[is.finite(v)]
    if (length(unique(v)) < 2L) {
      return(1)
    }
    bw <- stats::bw.nrd0(v)
    if (!is.finite(bw) || bw <= 0) {
      bw <- diff(range(v)) / 8
    }
    max(bw, 1)
  }
  c(estimate(x), estimate(y))
}

.sirap_weighted_gaussian_grid <- function(x, y, weights, x_grid, y_grid, bandwidth) {
  z <- matrix(0, nrow = length(x_grid), ncol = length(y_grid))
  if (length(x) == 0L) {
    return(z)
  }
  keep <- is.finite(x) & is.finite(y) & is.finite(weights) & weights > 0
  x <- x[keep]
  y <- y[keep]
  weights <- weights[keep]
  if (length(x) == 0L) {
    return(z)
  }
  for (i in seq_along(x)) {
    dx <- stats::dnorm((x_grid - x[[i]]) / bandwidth[[1]])
    dy <- stats::dnorm((y_grid - y[[i]]) / bandwidth[[2]])
    z <- z + weights[[i]] * tcrossprod(dx, dy)
  }
  z / sum(weights)
}

.sirap_clone_histogram_3d_bars <- function(region_bcrs,
                                           levels,
                                           color_values,
                                           xlim,
                                           ylim,
                                           hist_bins,
                                           z_base,
                                           max_height) {
  grid <- .sirap_density_grid(xlim, ylim, hist_bins + 1L)
  x_breaks <- grid$x
  y_breaks <- grid$y
  x_bins <- length(x_breaks) - 1L
  y_bins <- length(y_breaks) - 1L
  raw_counts <- lapply(levels, function(level) {
    dat <- region_bcrs[as.character(region_bcrs$.sirap_density_group) == level, , drop = FALSE]
    counts <- matrix(0, nrow = x_bins, ncol = y_bins)
    if (nrow(dat) == 0L) {
      return(counts)
    }
    xi <- findInterval(dat$x, x_breaks, all.inside = TRUE)
    yi <- findInterval(dat$.sirap_y, y_breaks, all.inside = TRUE)
    xi <- pmin(xi, x_bins)
    yi <- pmin(yi, y_bins)
    for (i in seq_len(nrow(dat))) {
      counts[xi[[i]], yi[[i]]] <- counts[xi[[i]], yi[[i]]] + dat$.sirap_density_weight[[i]]
    }
    counts
  })
  names(raw_counts) <- levels
  max_count <- max(vapply(raw_counts, function(x) max(x, na.rm = TRUE), numeric(1)), na.rm = TRUE)
  lapply(levels, function(level) {
    counts <- raw_counts[[level]]
    idx <- which(counts > 0, arr.ind = TRUE)
    if (nrow(idx) == 0L || !is.finite(max_count) || max_count <= 0) {
      bars <- data.frame(x0 = numeric(), x1 = numeric(), y0 = numeric(), y1 = numeric(), z0 = numeric(), z1 = numeric())
    } else {
      bars <- data.frame(
        x0 = x_breaks[idx[, 1]],
        x1 = x_breaks[idx[, 1] + 1L],
        y0 = y_breaks[idx[, 2]],
        y1 = y_breaks[idx[, 2] + 1L],
        z0 = z_base,
        z1 = z_base + counts[idx] / max_count * max_height,
        stringsAsFactors = FALSE
      )
    }
    list(
      level = level,
      bars = bars,
      color = unname(color_values[level])
    )
  }) |> stats::setNames(levels)
}

.sirap_draw_histogram_bars_3d <- function(histogram, alpha) {
  bars <- histogram$bars
  if (is.null(bars) || nrow(bars) == 0L) {
    return(invisible(NULL))
  }
  color <- grDevices::adjustcolor(histogram$color, alpha.f = alpha)
  for (i in seq_len(nrow(bars))) {
    vertices <- .sirap_histogram_bar_vertices(bars[i, , drop = FALSE])
    rgl::quads3d(
      vertices$x,
      vertices$y,
      vertices$z,
      color = color,
      lit = FALSE
    )
  }
  invisible(NULL)
}

.sirap_histogram_bar_vertices <- function(bar) {
  x0 <- bar$x0[[1]]
  x1 <- bar$x1[[1]]
  y0 <- bar$y0[[1]]
  y1 <- bar$y1[[1]]
  z0 <- bar$z0[[1]]
  z1 <- bar$z1[[1]]
  points <- rbind(
    c(x0, y0, z1), c(x1, y0, z1), c(x1, y1, z1), c(x0, y1, z1),
    c(x0, y0, z0), c(x0, y1, z0), c(x1, y1, z0), c(x1, y0, z0),
    c(x0, y0, z0), c(x1, y0, z0), c(x1, y0, z1), c(x0, y0, z1),
    c(x1, y0, z0), c(x1, y1, z0), c(x1, y1, z1), c(x1, y0, z1),
    c(x1, y1, z0), c(x0, y1, z0), c(x0, y1, z1), c(x1, y1, z1),
    c(x0, y1, z0), c(x0, y0, z0), c(x0, y0, z1), c(x0, y1, z1)
  )
  data.frame(x = points[, 1], y = points[, 2], z = points[, 3])
}

.sirap_image_texture_extent <- function(image_dim, xlim, ylim) {
  xlim <- range(xlim, na.rm = TRUE)
  ylim <- range(ylim, na.rm = TRUE)
  xlim <- c(max(0, xlim[[1]]), min(image_dim$width, xlim[[2]]))
  ylim <- c(max(-image_dim$height, ylim[[1]]), min(0, ylim[[2]]))
  if (diff(xlim) <= 0 || diff(ylim) <= 0) {
    xlim <- c(0, image_dim$width)
    ylim <- c(-image_dim$height, 0)
  }
  col_start <- max(1L, floor(xlim[[1]]) + 1L)
  col_end <- min(image_dim$width, ceiling(xlim[[2]]))
  row_start <- max(1L, floor(-ylim[[2]]) + 1L)
  row_end <- min(image_dim$height, ceiling(-ylim[[1]]))
  list(
    cols = seq.int(col_start, col_end),
    rows = seq.int(row_start, row_end),
    xlim = c(col_start - 1, col_end),
    ylim = c(-row_end, -(row_start - 1))
  )
}

.sirap_image_texture_file <- function(image, extent, alpha = 1, texture_dir = NULL) {
  if (is.null(image)) {
    crop <- array(1, dim = c(2, 2, 4))
    crop[, , 4] <- alpha
  } else {
    img <- SeuratObject::GetImage(image, mode = "raw")
    crop <- img[extent$rows, extent$cols, , drop = FALSE]
  }
  if (alpha < 1) {
    rgba <- array(NA_real_, dim = c(dim(crop)[1:2], 4))
    rgba[, , 1:3] <- crop[, , 1:3]
    rgba[, , 4] <- alpha
    crop <- rgba
  }
  if (!is.null(texture_dir)) {
    dir.create(texture_dir, recursive = TRUE, showWarnings = FALSE)
    alpha_label <- gsub("[^0-9A-Za-z]+", "_", format(alpha, scientific = FALSE, trim = TRUE))
    texture <- file.path(
      texture_dir,
      paste0(
        "sirap_spatial_clone_density_texture",
        "_r", min(extent$rows), "-", max(extent$rows),
        "_c", min(extent$cols), "-", max(extent$cols),
        "_a", alpha_label,
        ".png"
      )
    )
    if (file.exists(texture)) {
      return(texture)
    }
  } else {
    texture <- tempfile(
      pattern = "sirap_spatial_clone_density_texture_",
      tmpdir = tempdir(),
      fileext = ".png"
    )
  }
  grDevices::png(texture, width = dim(crop)[[2]], height = dim(crop)[[1]], units = "px", bg = "transparent")
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.raster(grDevices::as.raster(crop), interpolate = FALSE)
  texture
}

.sirap_polygon_pie_summary <- function(region_bcrs,
                                       fill,
                                       top_n,
                                       label_top_n,
                                       other_label,
                                       plot_other = TRUE) {
  if (nrow(region_bcrs) == 0L) {
    return(data.frame(
      polygon = character(),
      polygon_label = character(),
      level = character(),
      n = integer(),
      proportion = numeric(),
      ymin = numeric(),
      ymax = numeric(),
      label = character(),
      stringsAsFactors = FALSE
    ))
  }
  region_bcrs$.sirap_fill_by <- as.character(region_bcrs[[fill]])
  region_bcrs$.sirap_fill_by[is.na(region_bcrs$.sirap_fill_by) | !nzchar(region_bcrs$.sirap_fill_by)] <- NA_character_
  if (is.null(top_n)) {
    region_bcrs$.sirap_level <- region_bcrs$.sirap_fill_by
  } else {
    if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
      rlang::abort("`top_n` must be NULL or a positive number.")
    }
    top_levels <- unlist(lapply(split(region_bcrs, region_bcrs$polygon), function(dat) {
      counts <- sort(table(dat$.sirap_fill_by, useNA = "no"), decreasing = TRUE)
      names(counts)[seq_len(min(as.integer(top_n), length(counts)))]
    }), use.names = FALSE)
    top_levels <- unique(top_levels[!is.na(top_levels) & nzchar(top_levels)])
    region_bcrs$.sirap_level <- region_bcrs$.sirap_fill_by
    region_bcrs$.sirap_level[
      !is.na(region_bcrs$.sirap_level) &
        nzchar(region_bcrs$.sirap_level) &
        !region_bcrs$.sirap_level %in% top_levels
    ] <- other_label
  }
  pieces <- lapply(split(region_bcrs, region_bcrs$polygon), function(dat) {
    if (!isTRUE(plot_other)) {
      dat <- dat[dat$.sirap_level != other_label, , drop = FALSE]
    }
    if (nrow(dat) == 0L) {
      return(data.frame(level = character(), n = integer(), polygon = character()))
    }
    counts <- as.data.frame(table(dat$.sirap_level, useNA = "no"), stringsAsFactors = FALSE)
    names(counts) <- c("level", "n")
    counts <- counts[counts$n > 0L, , drop = FALSE]
    other <- counts[counts$level == other_label, , drop = FALSE]
    counts <- counts[counts$level != other_label, , drop = FALSE]
    counts <- counts[order(-counts$n, counts$level), , drop = FALSE]
    if (nrow(other) > 0L) {
      counts <- rbind(counts, other)
    }
    counts$.sirap_order <- seq_len(nrow(counts))
    counts$polygon <- unique(dat$polygon)[[1]]
    counts
  })
  out <- .sirap_bind_rows(pieces)
  if (nrow(out) == 0L) {
    return(out)
  }
  # Label pies by the actual polygon id (not a sequential 1..N index) so pie
  # titles match the stored polygon_id and the map labels.
  out$polygon_label <- as.character(out$polygon)
  out <- out[order(out$polygon_label, out$.sirap_order), , drop = FALSE]
  out <- do.call(rbind, lapply(split(out, out$polygon), function(dat) {
    total <- sum(dat$n)
    dat$proportion <- dat$n / total
    dat$ymax <- cumsum(dat$proportion)
    dat$ymin <- c(0, utils::head(dat$ymax, -1))
    dat$label_rank <- seq_len(nrow(dat))
    dat$label <- ifelse(dat$label_rank <= label_top_n, dat$level, "")
    dat
  }))
  out$.sirap_order <- NULL
  rownames(out) <- NULL
  out
}

# Decide the order of the side-by-side pie charts. `"receptors"` orders the
# polygons by their total receptor count (most first); ties break by the polygon
# map label. `"polygon"` keeps the numbered map order. Returns the polygon
# labels in the chosen order.
.sirap_pie_polygon_order <- function(pie_data, pie_order = c("receptors", "polygon")) {
  pie_order <- match.arg(pie_order)
  if (identical(pie_order, "polygon")) {
    return(unique(pie_data$polygon_label))
  }
  totals <- tapply(pie_data$n, pie_data$polygon_label, sum)
  names(totals)[order(-totals, suppressWarnings(as.integer(names(totals))))]
}

.sirap_polygon_centers <- function(polygon_plot_data, polygons) {
  centers <- lapply(split(polygon_plot_data, polygon_plot_data$polygon), function(dat) {
    dat <- dat[!duplicated(paste(dat$x, dat$y)), , drop = FALSE]
    data.frame(
      polygon = unique(dat$polygon)[[1]],
      x = mean(range(dat$x, na.rm = TRUE)),
      y = mean(range(dat$y, na.rm = TRUE)),
      .sirap_y = -mean(range(dat$y, na.rm = TRUE)),
      radius = max(diff(range(dat$x, na.rm = TRUE)), diff(range(dat$y, na.rm = TRUE))) / 2,
      stringsAsFactors = FALSE
    )
  })
  centers <- .sirap_bind_rows(centers)
  centers <- centers[centers$polygon %in% polygons, , drop = FALSE]
  # Label polygons by their actual id, not a sequential 1..N index, so pie titles
  # and map labels match the stored polygon_id.
  centers$polygon_label <- as.character(centers$polygon)
  centers
}

.sirap_dominance_only_polygon_pies <- function(pie_data, fill_values, label_size, pie_order) {
  ordered_labels <- .sirap_pie_polygon_order(pie_data, pie_order)
  pie_groups <- split(pie_data, factor(pie_data$polygon_label, levels = ordered_labels))
  pie_plots <- lapply(pie_groups, function(dat) {
    local_levels <- unique(as.character(dat$level))
    ggplot2::ggplot(dat, ggplot2::aes(fill = .data$level)) +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = 0.5, xmax = 1.5, ymin = .data$ymin, ymax = .data$ymax),
        color = "white", linewidth = 0.2
      ) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::geom_text(
        data = dat[dat$label != "", , drop = FALSE],
        ggplot2::aes(y = (.data$ymin + .data$ymax) / 2, label = .data$label),
        x = 1.25, inherit.aes = TRUE, size = label_size
      ) +
      ggplot2::scale_fill_manual(values = fill_values[local_levels], limits = local_levels, drop = TRUE) +
      ggplot2::theme_void(base_size = 8) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 8, hjust = 0.5),
        plot.margin = ggplot2::margin(2, 2, 2, 2)
      ) +
      ggplot2::labs(title = paste("Polygon", unique(dat$polygon_label)))
  })
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    return(pie_plots)
  }
  patchwork::wrap_plots(pie_plots, ncol = max(1L, ceiling(sqrt(length(pie_plots)))))
}

.sirap_side_by_side_polygon_pies <- function(object,
                                             images,
                                             image_scale,
                                             swap_coords,
                                             polygon_plot_data,
                                             centers,
                                             pie_data,
                                             fill_values,
                                             polygon_alpha,
                                             image_alpha,
                                             label_size,
                                             pies_per_page,
                                             pie_order = c("receptors", "polygon"),
                                             zoom) {
  pie_order <- match.arg(pie_order)
  map_polygon_data <- polygon_plot_data
  polygon_labels <- stats::setNames(centers$polygon_label, centers$polygon)
  map_polygon_data$polygon_label <- unname(polygon_labels[as.character(map_polygon_data$polygon)])
  map_plots <- plotSpatialPolygons(
    object = object,
    polygon = map_polygon_data,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords,
    polygon_id = "polygon",
    polygon_alpha = polygon_alpha,
    image_alpha = image_alpha,
    label = TRUE,
    label_by = "polygon_label",
    label_size = label_size,
    zoom = zoom,
    combine = FALSE
  )
  map_plot <- if (is.list(map_plots)) map_plots[[1L]] else map_plots

  ordered_labels <- .sirap_pie_polygon_order(pie_data, pie_order)
  pie_groups <- split(
    pie_data,
    factor(pie_data$polygon_label, levels = ordered_labels)
  )
  pie_plots <- lapply(pie_groups, function(dat) {
    local_levels <- unique(as.character(dat$level))
    legend_ncol <- max(1L, min(4L, ceiling(length(local_levels) / 4)))
    ggplot2::ggplot(dat, ggplot2::aes(fill = .data$level)) +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = 0.5, xmax = 1.5, ymin = .data$ymin, ymax = .data$ymax),
        color = "white",
        linewidth = 0.2
      ) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::geom_text(
        data = dat[dat$label != "", , drop = FALSE],
        ggplot2::aes(y = (.data$ymin + .data$ymax) / 2, label = .data$label),
        x = 1.25,
        inherit.aes = TRUE,
        size = label_size
      ) +
      ggplot2::scale_fill_manual(values = fill_values[local_levels], limits = local_levels, drop = TRUE) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, ncol = legend_ncol, byrow = TRUE)) +
      ggplot2::theme_void(base_size = 8) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 5),
        legend.key.size = grid::unit(0.12, "in"),
        legend.spacing.x = grid::unit(0.03, "in"),
        legend.spacing.y = grid::unit(0.01, "in"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        plot.title = ggplot2::element_text(size = 8, hjust = 0.5),
        plot.margin = ggplot2::margin(2, 2, 2, 2)
      ) +
      ggplot2::labs(title = paste("Polygon", unique(dat$polygon_label)))
  })

  page_indices <- split(seq_along(pie_plots), ceiling(seq_along(pie_plots) / pies_per_page))
  page_names <- paste0("page_", seq_along(page_indices))
  names(page_indices) <- page_names

  if (requireNamespace("patchwork", quietly = TRUE)) {
    pages <- lapply(seq_along(page_indices), function(page_i) {
      idx <- page_indices[[page_i]]
      page_pies <- pie_plots[idx]
      pie_panel <- patchwork::wrap_plots(page_pies, ncol = max(1L, ceiling(sqrt(length(page_pies)))))
      page <- map_plot + pie_panel + patchwork::plot_layout(widths = c(1, 1))
      if (length(page_indices) > 1L) {
        page <- page + patchwork::plot_annotation(title = paste0("Page ", page_i, " of ", length(page_indices)))
      }
      page
    })
    names(pages) <- page_names
    if (length(pages) == 1L) pages[[1L]] else pages
  } else {
    pages <- lapply(page_indices, function(idx) {
      list(pies = pie_plots[idx], polygons = map_plot)
    })
    if (length(pages) == 1L) pages[[1L]] else pages
  }
}

.sirap_spatial_clipped_polygon_pies <- function(object,
                                                images,
                                                image_scale,
                                                swap_coords,
                                                polygon_plot_data,
                                                centers,
                                                pie_data,
                                                fill_values,
                                                fill,
                                                image_alpha,
                                                pie_alpha,
                                                label_size,
                                                spatial_legend,
                                                zoom) {
  .require_namespace("sf")
  clipped <- .sirap_clipped_pie_polygons(polygon_plot_data, centers, pie_data)
  labels <- .sirap_spatial_pie_labels(centers, pie_data)
  legend_ncol <- max(1L, min(4L, ceiling(length(fill_values) / 8)))
  plotted_polygons <- unique(as.character(pie_data$polygon))
  zoom_data <- polygon_plot_data[as.character(polygon_plot_data$polygon) %in% plotted_polygons, , drop = FALSE]
  full_image_data <- data.frame(x = numeric(), .sirap_y = numeric())

  plots <- list()
  for (image in .sirap_spatial_images(object, images = images)) {
    spatial <- .sirap_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    plots[[image]] <- .sirap_spatial_base(
      plot_data = if (isTRUE(zoom)) zoom_data else full_image_data,
      image = spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = NULL
    ) +
      ggplot2::geom_polygon(
        data = clipped,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$piece, fill = .data$level),
        inherit.aes = FALSE,
        color = "white",
        linewidth = 0.15,
        alpha = pie_alpha
      ) +
      ggplot2::geom_path(
        data = polygon_plot_data,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, group = .data$polygon),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.6
      ) +
      ggplot2::geom_text(
        data = labels,
        ggplot2::aes(x = .data$x, y = .data$.sirap_y, label = .data$label),
        inherit.aes = FALSE,
        size = label_size
      ) +
      ggplot2::scale_fill_manual(values = fill_values, limits = unique(as.character(pie_data$level)), drop = TRUE) +
      ggplot2::labs(fill = fill)
    if (isTRUE(spatial_legend)) {
      plots[[image]] <- plots[[image]] +
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE)) +
        ggplot2::theme(
          legend.position = "right",
          legend.text = ggplot2::element_text(size = 6),
          legend.title = ggplot2::element_text(size = 7),
          legend.key.size = grid::unit(0.14, "in"),
          legend.spacing.x = grid::unit(0.04, "in"),
          legend.spacing.y = grid::unit(0.02, "in")
        )
    } else {
      plots[[image]] <- plots[[image]] + ggplot2::guides(fill = "none")
    }
  }

  .sirap_return_spatial_plots(plots, combine = TRUE)
}

.sirap_clipped_pie_polygons <- function(polygon_plot_data, centers, pie_data) {
  pieces <- list()
  idx <- 1L
  for (polygon_id in unique(pie_data$polygon)) {
    polygon_vertices <- polygon_plot_data[polygon_plot_data$polygon == polygon_id, , drop = FALSE]
    center <- centers[centers$polygon == polygon_id, , drop = FALSE]
    if (nrow(polygon_vertices) == 0L || nrow(center) == 0L) {
      next
    }
    boundary <- .sirap_repair_polygon_sfc(
      sf::st_sfc(sf::st_polygon(list(as.matrix(polygon_vertices[, c("x", "y")]))))
    )
    if (is.null(boundary)) {
      rlang::warn(paste0("Polygon ", polygon_id, " is invalid and could not be repaired for clipped pie plotting."))
      next
    }
    this_pie <- pie_data[pie_data$polygon == polygon_id, , drop = FALSE]
    this_pie <- this_pie[order(this_pie$ymin), , drop = FALSE]
    for (i in seq_len(nrow(this_pie))) {
      sector <- .sirap_pie_sector(center$x[[1]], center$y[[1]], center$radius[[1]] * 1.25, this_pie$ymin[[i]], this_pie$ymax[[i]])
      sector_sf <- .sirap_repair_polygon_sfc(
        sf::st_sfc(sf::st_polygon(list(as.matrix(sector[, c("x", "y")]))))
      )
      if (is.null(sector_sf)) {
        rlang::warn(paste0(
          "Could not create a valid pie slice for polygon ",
          polygon_id,
          " and level ",
          this_pie$level[[i]],
          "."
        ))
        next
      }
      intersection <- tryCatch(
        suppressWarnings(.sirap_polygon_intersection(sector_sf, boundary)),
        error = function(e) {
          rlang::warn(paste0(
            "Could not clip pie slice for polygon ",
            polygon_id,
            " and level ",
            this_pie$level[[i]],
            ": ",
            conditionMessage(e)
          ))
          NULL
        }
      )
      if (is.null(intersection)) {
        next
      }
      intersection <- .sirap_repair_polygon_sfc(intersection)
      if (is.null(intersection) || length(intersection) == 0L || all(sf::st_is_empty(intersection))) {
        next
      }
      coords <- as.data.frame(sf::st_coordinates(intersection))
      geometry_levels <- grep("^L[0-9]+$", colnames(coords), value = TRUE)
      if (length(geometry_levels) == 0L) {
        coordinate_piece <- rep("1", nrow(coords))
      } else {
        coordinate_piece <- as.character(interaction(coords[, geometry_levels, drop = FALSE], drop = TRUE, lex.order = TRUE))
      }
      coords$piece <- paste(polygon_id, this_pie$level[[i]], coordinate_piece, idx, sep = "__")
      coords$polygon <- polygon_id
      coords$level <- this_pie$level[[i]]
      coords$x <- coords$X
      coords$y <- coords$Y
      coords$.sirap_y <- -coords$y
      pieces[[idx]] <- coords[, c("polygon", "level", "piece", "x", "y", ".sirap_y"), drop = FALSE]
      idx <- idx + 1L
    }
  }
  .sirap_bind_rows(pieces)
}

.sirap_repair_polygon_sfc <- function(geometry) {
  repaired <- tryCatch({
    if ("st_make_valid" %in% getNamespaceExports("sf")) {
      geometry <- sf::st_make_valid(geometry)
    } else if (!all(sf::st_is_valid(geometry))) {
      geometry <- sf::st_buffer(geometry, 0)
    }
    geometry <- suppressWarnings(sf::st_collection_extract(geometry, "POLYGON"))
    geometry <- suppressWarnings(sf::st_cast(geometry, "POLYGON"))
    if (length(geometry) == 0L || all(sf::st_is_empty(geometry))) {
      return(NULL)
    }
    if (!all(sf::st_is_valid(geometry))) {
      geometry <- suppressWarnings(sf::st_buffer(geometry, 0))
      geometry <- suppressWarnings(sf::st_collection_extract(geometry, "POLYGON"))
      geometry <- suppressWarnings(sf::st_cast(geometry, "POLYGON"))
    }
    if (length(geometry) == 0L || all(sf::st_is_empty(geometry)) || !all(sf::st_is_valid(geometry))) {
      return(NULL)
    }
    geometry
  }, error = function(e) NULL)
  repaired
}

.sirap_polygon_intersection <- function(x, y) {
  tryCatch(
    suppressWarnings(sf::st_intersection(x, y)),
    error = function(e) {
      x <- .sirap_repair_polygon_sfc(suppressWarnings(sf::st_buffer(x, 0)))
      y <- .sirap_repair_polygon_sfc(suppressWarnings(sf::st_buffer(y, 0)))
      if (is.null(x) || is.null(y)) {
        stop(e)
      }
      suppressWarnings(sf::st_intersection(x, y))
    }
  )
}

.sirap_pie_sector <- function(center_x, center_y, radius, ymin, ymax, n = 80L) {
  start <- (ymin * 2 * pi) - (pi / 2)
  end <- (ymax * 2 * pi) - (pi / 2)
  angles <- seq(start, end, length.out = max(3L, ceiling(n * (ymax - ymin))))
  data.frame(
    x = c(center_x, center_x + radius * cos(angles), center_x),
    y = c(center_y, center_y + radius * sin(angles), center_y),
    stringsAsFactors = FALSE
  )
}

.sirap_spatial_pie_labels <- function(centers, pie_data) {
  labels <- pie_data[pie_data$label != "", , drop = FALSE]
  if (nrow(labels) == 0L) {
    return(data.frame(x = numeric(), y = numeric(), .sirap_y = numeric(), label = character()))
  }
  labels <- merge(labels, centers[, c("polygon", "x", "y", "radius"), drop = FALSE], by = "polygon", sort = FALSE)
  mid <- ((labels$ymin + labels$ymax) / 2) * 2 * pi - (pi / 2)
  labels$x <- labels$x + labels$radius * 0.45 * cos(mid)
  labels$y <- labels$y + labels$radius * 0.45 * sin(mid)
  labels$.sirap_y <- -labels$y
  labels[, c("x", "y", ".sirap_y", "label"), drop = FALSE]
}

.sirap_downsample_polygon_plot_coords <- function(coords, downsample_n = NULL, seed = 1L) {
  if (is.null(downsample_n)) {
    return(coords)
  }
  if (!is.numeric(downsample_n) || length(downsample_n) != 1L || is.na(downsample_n) || downsample_n < 1) {
    rlang::abort("`downsample_n` must be NULL or a positive number.")
  }
  downsample_n <- as.integer(downsample_n)
  if (nrow(coords) <= downsample_n) {
    return(coords)
  }
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  coords[sample.int(nrow(coords), downsample_n), , drop = FALSE]
}

.sirap_close_current_polygon <- function(points, polygon_num) {
  id_col <- if ("num" %in% colnames(points)) "num" else "polygon"
  current <- points[points[[id_col]] == polygon_num, , drop = FALSE]
  if (nrow(current) == 0L) {
    return(points)
  }
  if (nrow(current) < 3L) {
    rlang::warn("A polygon needs at least three distinct vertices before it can be closed.")
    return(points)
  }
  if (.sirap_same_xy(current[1, , drop = FALSE], current[nrow(current), , drop = FALSE])) {
    return(points)
  }
  rbind(points, current[1, , drop = FALSE])
}

.sirap_default_polygon_buffer_radius <- function(coords) {
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  coords <- unique(coords[, c("x", "y"), drop = FALSE])
  if (nrow(coords) < 2L) {
    rlang::abort("At least two spatial coordinates are required to estimate `buffer_radius`.")
  }
  nn <- RANN::nn2(
    data = as.matrix(coords[, c("x", "y"), drop = FALSE]),
    query = as.matrix(coords[, c("x", "y"), drop = FALSE]),
    k = 2L,
    treetype = "kd",
    searchtype = "standard"
  )
  distances <- nn$nn.dists[, 2L]
  distances <- distances[is.finite(distances) & distances > 0]
  if (length(distances) == 0L) {
    rlang::abort("Could not estimate `buffer_radius` from the spatial coordinates.")
  }
  stats::median(distances) * 0.55
}

.sirap_polygon_interior_to_sf_boundary <- function(points, buffer_radius) {
  if (!is.numeric(buffer_radius) || length(buffer_radius) != 1L ||
      is.na(buffer_radius) || buffer_radius <= 0) {
    rlang::abort("`buffer_radius` must be a positive number.")
  }
  boundary_tables <- lapply(unique(points$polygon), function(polygon_id) {
    current <- points[points$polygon == polygon_id, , drop = FALSE]
    if (nrow(current) < 3L) {
      rlang::warn(paste0("Polygon ", polygon_id, " has fewer than three points and will be excluded."))
      return(NULL)
    }

    point_sf <- sf::st_as_sf(
      current,
      coords = c("x", "y"),
      crs = NA,
      remove = FALSE
    )
    buffered <- suppressWarnings(sf::st_buffer(point_sf, dist = buffer_radius))
    unioned <- suppressWarnings(sf::st_union(sf::st_geometry(buffered)))
    polygons <- suppressWarnings(sf::st_collection_extract(unioned, "POLYGON"))
    polygons <- suppressWarnings(sf::st_cast(polygons, "POLYGON"))

    if (length(polygons) == 0L || all(sf::st_is_empty(polygons))) {
      rlang::warn(paste0(
        "No polygon boundary could be created for polygon ",
        polygon_id,
        "; it will be excluded."
      ))
      return(NULL)
    }

    if (length(polygons) > 1L) {
      polygon_areas <- as.numeric(sf::st_area(polygons))
      polygons <- polygons[which.max(polygon_areas)]
      rlang::warn(paste0(
        "Polygon ",
        polygon_id,
        " produced multiple disconnected outlines; keeping the largest outline."
      ))
    }

    boundary <- as.data.frame(sf::st_coordinates(polygons))
    if ("L1" %in% colnames(boundary)) {
      boundary <- boundary[boundary$L1 == 1L, , drop = FALSE]
    }
    boundary <- data.frame(
      x = as.numeric(boundary$X),
      y = as.numeric(boundary$Y),
      stringsAsFactors = FALSE
    )
    if (nrow(boundary) < 3L) {
      rlang::warn(paste0("Polygon ", polygon_id, " has fewer than three boundary vertices and will be excluded."))
      return(NULL)
    }

    nearest <- RANN::nn2(
      data = as.matrix(current[, c("x", "y"), drop = FALSE]),
      query = as.matrix(boundary[, c("x", "y"), drop = FALSE]),
      k = 1L,
      treetype = "kd",
      searchtype = "standard"
    )
    nearest_idx <- nearest$nn.idx[, 1L]
    nearest_idx <- nearest_idx[!is.na(nearest_idx)]
    if (length(nearest_idx) == 0L) {
      rlang::warn(paste0("No boundary spatial IDs could be identified for polygon ", polygon_id, "."))
      return(NULL)
    }
    nearest_idx <- nearest_idx[c(TRUE, diff(nearest_idx) != 0)]
    nearest_idx <- nearest_idx[!duplicated(nearest_idx)]
    if (length(nearest_idx) < 3L) {
      rlang::warn(paste0("Fewer than three boundary spatial IDs were found for polygon ", polygon_id, "."))
      return(NULL)
    }
    boundary_points <- current[nearest_idx, , drop = FALSE]
    boundary_points <- .sirap_clean_boundary_barcode_vertices(boundary_points)
    if (nrow(boundary_points) < 3L) {
      rlang::warn(paste0("Fewer than three non-degenerate boundary spatial IDs were found for polygon ", polygon_id, "."))
      return(NULL)
    }
    if (!.sirap_boundary_barcode_vertices_valid(boundary_points)) {
      boundary_points <- .sirap_repair_boundary_barcode_vertices(
        boundary_points = boundary_points,
        candidate_points = current
      )
    }
    if (!.sirap_boundary_barcode_vertices_valid(boundary_points)) {
      rlang::warn(paste0(
        "Polygon ",
        polygon_id,
        " could not be represented as a valid spatial_id-coordinate boundary and will be excluded."
      ))
      return(NULL)
    }

    data.frame(
      spatial_id = boundary_points$spatial_id,
      x = boundary_points$x,
      y = boundary_points$y,
      polygon = boundary_points$polygon,
      image = boundary_points$image,
      stringsAsFactors = FALSE
    )
  })
  boundary_tables <- Filter(Negate(is.null), boundary_tables)
  if (length(boundary_tables) == 0L) {
    rlang::abort("No polygons with at least three points remain.")
  }
  out <- do.call(rbind, boundary_tables)
  rownames(out) <- NULL
  out
}

.sirap_clean_boundary_barcode_vertices <- function(points, tol = 1e-8) {
  points <- points[is.finite(points$x) & is.finite(points$y), , drop = FALSE]
  if (nrow(points) < 3L) {
    return(points)
  }

  changed <- TRUE
  while (isTRUE(changed) && nrow(points) >= 3L) {
    old_n <- nrow(points)
    if (nrow(points) > 1L && .sirap_same_xy(points[1, , drop = FALSE], points[nrow(points), , drop = FALSE])) {
      points <- points[-nrow(points), , drop = FALSE]
    }
    if (nrow(points) > 1L) {
      duplicate_previous <- c(
        FALSE,
        abs(diff(points$x)) <= tol & abs(diff(points$y)) <= tol
      )
      points <- points[!duplicate_previous, , drop = FALSE]
    }
    if (nrow(points) >= 4L) {
      remove <- logical(nrow(points))
      for (i in seq_len(nrow(points))) {
        prev_i <- if (i == 1L) nrow(points) else i - 1L
        next_i <- if (i == nrow(points)) 1L else i + 1L
        remove[[i]] <- .sirap_degenerate_boundary_vertex(
          prev = points[prev_i, , drop = FALSE],
          current = points[i, , drop = FALSE],
          next_point = points[next_i, , drop = FALSE],
          tol = tol
        )
      }
      points <- points[!remove, , drop = FALSE]
    }
    changed <- nrow(points) != old_n
  }
  rownames(points) <- NULL
  points
}

.sirap_degenerate_boundary_vertex <- function(prev, current, next_point, tol = 1e-8) {
  if (.sirap_same_xy(prev, next_point)) {
    return(TRUE)
  }
  v1 <- c(prev$x[[1]] - current$x[[1]], prev$y[[1]] - current$y[[1]])
  v2 <- c(next_point$x[[1]] - current$x[[1]], next_point$y[[1]] - current$y[[1]])
  cross <- v1[[1]] * v2[[2]] - v1[[2]] * v2[[1]]
  dot <- sum(v1 * v2)
  isTRUE(abs(cross) <= tol && dot > tol)
}

.sirap_closed_boundary_matrix <- function(points) {
  xy <- as.matrix(points[, c("x", "y"), drop = FALSE])
  if (nrow(xy) == 0L) {
    return(xy)
  }
  if (!all(xy[1L, ] == xy[nrow(xy), ])) {
    xy <- rbind(xy, xy[1L, ])
  }
  xy
}

.sirap_boundary_barcode_vertices_valid <- function(points) {
  points <- .sirap_clean_boundary_barcode_vertices(points)
  if (nrow(unique(points[, c("x", "y"), drop = FALSE])) < 3L) {
    return(FALSE)
  }
  xy <- .sirap_closed_boundary_matrix(points)
  polygon <- tryCatch(
    sf::st_sfc(sf::st_polygon(list(xy))),
    error = function(e) NULL
  )
  if (is.null(polygon) || length(polygon) == 0L || any(sf::st_is_empty(polygon))) {
    return(FALSE)
  }
  isTRUE(all(sf::st_is_valid(polygon))) && isTRUE(as.numeric(sf::st_area(polygon)) > 0)
}

.sirap_repair_boundary_barcode_vertices <- function(boundary_points, candidate_points) {
  if (nrow(boundary_points) < 3L || nrow(candidate_points) < 3L) {
    return(boundary_points)
  }
  polygon <- tryCatch(
    sf::st_sfc(sf::st_polygon(list(.sirap_closed_boundary_matrix(boundary_points)))),
    error = function(e) NULL
  )
  if (is.null(polygon)) {
    return(boundary_points)
  }
  repaired <- if ("st_make_valid" %in% getNamespaceExports("sf")) {
    sf::st_make_valid(polygon)
  } else {
    sf::st_buffer(polygon, 0)
  }
  repaired <- suppressWarnings(sf::st_collection_extract(repaired, "POLYGON"))
  repaired <- suppressWarnings(sf::st_cast(repaired, "POLYGON"))
  if (length(repaired) == 0L || all(sf::st_is_empty(repaired))) {
    return(boundary_points)
  }
  if (length(repaired) > 1L) {
    repaired <- repaired[which.max(as.numeric(sf::st_area(repaired)))]
  }
  boundary <- as.data.frame(sf::st_coordinates(repaired))
  if ("L1" %in% colnames(boundary)) {
    boundary <- boundary[boundary$L1 == 1L, , drop = FALSE]
  }
  boundary <- data.frame(
    x = as.numeric(boundary$X),
    y = as.numeric(boundary$Y),
    stringsAsFactors = FALSE
  )
  nearest <- RANN::nn2(
    data = as.matrix(candidate_points[, c("x", "y"), drop = FALSE]),
    query = as.matrix(boundary[, c("x", "y"), drop = FALSE]),
    k = 1L,
    treetype = "kd",
    searchtype = "standard"
  )
  nearest_idx <- nearest$nn.idx[, 1L]
  nearest_idx <- nearest_idx[!is.na(nearest_idx)]
  nearest_idx <- nearest_idx[c(TRUE, diff(nearest_idx) != 0)]
  nearest_idx <- nearest_idx[!duplicated(nearest_idx)]
  .sirap_clean_boundary_barcode_vertices(candidate_points[nearest_idx, , drop = FALSE])
}

.sirap_close_formatted_polygons <- function(points) {
  closed <- lapply(unique(points$polygon), function(polygon_id) {
    current <- points[points$polygon == polygon_id, , drop = FALSE]
    if (nrow(current) < 3L) {
      rlang::warn(paste0("Polygon ", polygon_id, " has fewer than three points and will be excluded."))
      return(NULL)
    }
    if (.sirap_same_xy(current[1, , drop = FALSE], current[nrow(current), , drop = FALSE])) {
      return(current)
    }
    rbind(current, current[1, , drop = FALSE])
  })
  closed <- Filter(Negate(is.null), closed)
  if (length(closed) == 0L) {
    rlang::abort("No polygons with at least three points remain.")
  }
  out <- do.call(rbind, closed)
  rownames(out) <- NULL
  out
}
