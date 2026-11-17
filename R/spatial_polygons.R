#' Calculate clonal dominance inside polygons
#'
#' `spatialCloneDominance` calculates the share of each clone among the
#' receptors inside spatial polygons.
#'
#' @details
#' The function finds the spatial IDs whose coordinates lie inside each polygon.
#' It keeps the receptor rows with those spatial IDs. Then it calculates the
#' fraction of these receptors that belong to each clone.
#'
#' If you give a SIRUP object and omit `clone`, `locus`, or `count`, the
#' function uses the field names stored in the object. If you omit `polygon`, it
#' uses the stored polygons.
#'
#' @param receptors AIRR-style receptor table with clone IDs and a spatial ID
#'   column. You can give a SIRUP object here if you omit `object`.
#' @param object Seurat object or SIRUP object. If it is a SIRUP object, you can
#'   omit `receptors`.
#' @param polygon Data frame of polygon vertices in order. It needs the columns
#'   set by `x` and `y`. If the column set by `polygon_id` exists, the function
#'   treats each polygon on its own. If `object` is a SIRUP object and `polygon`
#'   is `NULL`, the function uses the stored polygons.
#' @param spatial_id Name of the `receptors` column that holds the spatial IDs.
#'   If `NULL`, the function uses the stored field of a SIRUP object, or looks
#'   for a standard spatial ID column.
#' @param clone Name of the `receptors` column that holds the clone IDs (default
#'   `"clone_id"`).
#' @param locus_val Locus values to keep, for example `"IGH"`. If `NULL`, the
#'   function keeps all loci.
#' @param locus Name of the `receptors` column that holds the locus (default
#'   `"locus"`).
#' @param count Name of the numeric `receptors` column that holds the counts. If
#'   `NULL`, each receptor row counts as one.
#' @param v_call Name of the `receptors` column that holds the V-gene calls
#'   (default `"v_call"`). If the column exists, the function also returns the
#'   V genes of each clone at three levels: with allele, without allele, and
#'   V-gene superfamily.
#' @param exclude_singletons Exclude clones that have only one receptor row
#'   after the locus filter (default `FALSE`)?
#' @param images Names of the Seurat images to use.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param polygon_val IDs of the polygons to use. If `NULL`, the function uses
#'   all polygons.
#' @param x Name of the x-coordinate column in `polygon` (default `"x"`).
#' @param y Name of the y-coordinate column in `polygon` (default `"y"`).
#' @param polygon_id Name of the column in `polygon` that groups the vertices
#'   by polygon (default `"polygon"`).
#' @param include_boundary Count spots on a polygon edge as inside the polygon
#'   (default `TRUE`)? The function ignores this argument if
#'   `object` is a SIRUP object with stored polygon annotations and `polygon` is
#'   `NULL`.
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#'
#' @return An object of class `sirup_spatial_clone_dominance`. It is a list with
#'   these items: `dominance` (the dominance table), `region_spatial_ids`,
#'   `region_coordinates`, `region_receptors`, `polygons` (the vertices prepared
#'   for plotting), and `metadata` (the receptor field names that the function
#'   used).
#' @seealso [plotSpatialRegionReceptors()], [plotSpatialPolygonPies()],
#'   [plotSpatialCloneDensity3D()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Use the polygons stored in the SIRUP object.
#' dominance <- spatialCloneDominance(x)
#' dominance$dominance[, c("polygon_id", "clone_id", "n_receptors", "proportion")]
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
                                  include_boundary = TRUE,
                                  slot = NULL) {
  clone_missing <- missing(clone)
  locus_missing <- missing(locus)
  count_missing <- missing(count)
  if (isSIRUP(receptors) && is.null(object)) {
    object <- receptors
    receptors <- NULL
  }
  options <- .sirup_spatial_options(object, image_scale = image_scale, swap_coords = swap_coords)
  image_scale <- options$image_scale
  swap_coords <- options$swap_coords
  is_sirup <- isSIRUP(object)
  sirup_container <- if (is_sirup) object else NULL
  use_stored_polygons <- is_sirup && is.null(polygon)
  if (is.null(polygon)) {
    if (!is_sirup) {
      rlang::abort("`polygon` must be supplied when `object` is not a SIRUP object.")
    }
    polygon <- .sirup_polygon_vertices_table(.sirup_polygons(object))
    if (nrow(polygon) == 0L) {
      rlang::abort("The SIRUP object does not contain stored polygons.")
    }
  }
  inputs <- .sirup_resolve_inputs(object, receptors, spatial_id, slot)
  object <- inputs$object
  receptors <- inputs$receptors
  spatial_id <- inputs$spatial_id_col
  sample_col <- .sirup_receptor_sample_col(receptors %||% data.frame(), inputs$receptor_fields, object)
  if (is.null(receptors)) {
    rlang::abort("`receptors` is required when `object` is not a SIRUP object.")
  }
  if (is_sirup && clone_missing && !is.null(inputs$receptor_fields$clone)) {
    clone <- inputs$receptor_fields$clone
  }
  if (is_sirup && locus_missing && !is.null(inputs$receptor_fields$locus)) {
    locus <- inputs$receptor_fields$locus
  }
  if (is_sirup && count_missing && !is.null(inputs$receptor_fields$count)) {
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

  coords <- .sirup_spatial_coordinate_table(object, images = images, image_scale = image_scale, swap_coords = swap_coords)
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  coords <- .sirup_add_coord_spatial_keys(coords)
  receptors <- .sirup_add_receptor_spatial_keys(receptors, spatial_id, sample_col)
  polygon_specs <- .sirup_polygon_specs(polygon, polygon_id = polygon_val, polygon_col = polygon_id)
  polygon_plot_data <- .sirup_polygon_plot_data(
    polygon = polygon,
    polygon_id = polygon_val,
    polygon_col = polygon_id,
    x_col = x,
    y_col = y,
    keep_cols = .sirup_polygon_sample_col(polygon)
  )

  results <- lapply(polygon_specs, function(spec) {
    polygon_value <- if (isTRUE(spec$add_polygon_col)) spec$id else "polygon_1"
    polygon_display <- .sirup_polygon_display_id(polygon_plot_data, polygon_value)
    polygon_plot_rows <- polygon_plot_data[
      as.character(polygon_plot_data$polygon) == as.character(polygon_value),
      ,
      drop = FALSE
    ]
    polygon_sample <- .sirup_polygon_sample_values(polygon_plot_rows)
    polygon_sample <- if (length(polygon_sample) == 0L) NA_character_ else polygon_sample[[1L]]
    if (isTRUE(use_stored_polygons)) {
      region_coords <- .sirup_stored_polygon_coords(
        coords = coords,
        spatial = getSpatial(sirup_container),
        polygon_value = polygon_value
      )
    } else {
      polygon_rows <- .sirup_polygon_spec_rows(polygon, spec, polygon_col = polygon_id)
      vertices <- .sirup_select_polygon_vertices(
        polygon = polygon,
        polygon_id = spec$id,
        polygon_col = polygon_id,
        x_col = x,
        y_col = y
      )
      inside <- .sirup_points_in_polygon(
        px = coords$x,
        py = coords$y,
        vx = vertices$x,
        vy = vertices$y,
        include_boundary = include_boundary
      )
      polygon_samples <- .sirup_polygon_sample_values(polygon_rows)
      if (length(polygon_samples) > 0L && "sample" %in% colnames(coords)) {
        inside <- inside & !is.na(coords$sample) & as.character(coords$sample) %in% polygon_samples
      }
      region_coords <- coords[inside, , drop = FALSE]
    }
    region_coords <- .sirup_add_coord_spatial_keys(region_coords)
    region_key_cols <- .sirup_spatial_key_cols(region_coords)
    region_spatial_keys <- unique(region_coords[, region_key_cols, drop = FALSE])
    region_spatial_ids <- unique(as.character(region_coords$cell))

    regional_receptors <- if (isTRUE(use_stored_polygons)) {
      .sirup_receptors_in_polygon(receptors, polygon_value, polygon_id = "polygon_id")
    } else {
      receptor_keys <- receptors[, .sirup_spatial_key_cols(receptors, region_spatial_keys), drop = FALSE]
      receptors[.sirup_key_vector(receptor_keys) %in% .sirup_key_vector(region_spatial_keys), , drop = FALSE]
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
          genes <- .sirup_v_genes(calls)
          families <- unique(stats::na.omit(.sirup_v_gene_family(genes)))
          families <- sort(families)
          superfamilies <- unique(stats::na.omit(.sirup_v_gene_superfamily(.sirup_v_gene_family(genes))))
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
    out$polygon <- rep(as.character(polygon_display), nrow(out))
    out$polygon_id <- rep(as.character(polygon_value), nrow(out))
    out$polygon_label <- rep(as.character(polygon_display), nrow(out))
    if (!is.na(polygon_sample)) {
      out$sample <- rep(as.character(polygon_sample), nrow(out))
    }

    region_spatial_id_table <- unique(region_coords[, c(intersect("sample", colnames(region_coords)), "cell"), drop = FALSE])
    names(region_spatial_id_table)[names(region_spatial_id_table) == "cell"] <- "spatial_id"
    region_spatial_id_table$polygon <- rep(as.character(polygon_display), nrow(region_spatial_id_table))
    region_spatial_id_table$polygon_id <- rep(as.character(polygon_value), nrow(region_spatial_id_table))
    region_coords$polygon <- rep(as.character(polygon_display), nrow(region_coords))
    region_coords$polygon_id <- rep(as.character(polygon_value), nrow(region_coords))
    region_coords$polygon_label <- rep(as.character(polygon_display), nrow(region_coords))
    if (nrow(regional_receptors) == 0L || nrow(region_coords) == 0L) {
      region_receptor_table <- regional_receptors[0, , drop = FALSE]
      region_receptor_table$.sirup_spatial_id <- character()
      region_receptor_table$x <- numeric()
      region_receptor_table$y <- numeric()
      region_receptor_table$.sirup_y <- numeric()
      region_receptor_table$image <- character()
      region_receptor_table$polygon <- character()
      region_receptor_table$polygon_id <- character()
      region_receptor_table$polygon_label <- character()
      region_receptor_table <- .sirup_add_v_gene_columns(region_receptor_table, v_call_col = v_call)
    } else {
      region_coord_table <- region_coords[, unique(c(region_key_cols, "cell", "x", "y", "image")), drop = FALSE]
      region_coord_table$.sirup_y <- -region_coord_table$y
      region_receptor_table <- .sirup_drop_plot_coordinate_cols(regional_receptors)
      region_receptor_table <- .sirup_add_receptor_spatial_keys(region_receptor_table, spatial_id, sample_col)
      join_by <- .sirup_spatial_key_cols(region_receptor_table, region_coord_table)
      region_receptor_table <- merge(
        region_receptor_table,
        region_coord_table,
        by = join_by,
        all.x = FALSE,
        all.y = FALSE,
        sort = FALSE
      )
      region_receptor_table$polygon <- as.character(polygon_display)
      region_receptor_table$polygon_id <- as.character(polygon_value)
      region_receptor_table$polygon_label <- as.character(polygon_display)
      region_receptor_table <- .sirup_add_v_gene_columns(region_receptor_table, v_call_col = v_call)
    }
    list(
      dominance = out,
      region_spatial_ids = region_spatial_id_table,
      region_coordinates = region_coords,
      region_receptors = region_receptor_table
    )
  })

  out <- .sirup_rbind(results, "dominance")

  out <- list(
    dominance = out,
    region_spatial_ids = .sirup_rbind(results, "region_spatial_ids"),
    region_coordinates = .sirup_rbind(results, "region_coordinates"),
    region_receptors = .sirup_rbind(results, "region_receptors"),
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
  class(out) <- c("sirup_spatial_clone_dominance", class(out))
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
.sirup_polygon_feature_values <- function(object, feature, assay = NULL, layer = "counts") {
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
    .sirup_fetch_data(fetch_object, vars = feature, layer = layer),
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
.sirup_nn_index <- function(xy) {
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

#' Select polygon vertices in an interactive app
#'
#' `selectSpatialPolygons` opens a Shiny app. You click polygon vertices on a
#' plot of the spatial coordinates.
#'
#' @details
#' Each click snaps to the nearest spatial ID. Click "Close polygon" or "New
#' polygon" before you finish. This repeats the first vertex at the end of each
#' polygon.
#'
#' The vertices can describe overlapping or nested polygons. The function that
#' receives the vertices decides the membership: [assignSpatialPolygons()],
#' [assignReceptorPolygons()], [createSIRUP()], or [addPolygons()]. The primary
#' `polygon_id` gets the ID of the smallest matching polygon by area. The
#' logical `polygon_<ID>` columns show membership in each polygon.
#'
#' @param object Seurat object or SIRUP object.
#' @param feature Feature to show as a color scale under the vertices. You can
#'   change it in the app. The function looks in the assay first. If the feature
#'   is not there, it looks in the object metadata. This lets you show module
#'   scores (from `AddModuleScore()`, `UCell`, or any numeric column) without
#'   adding them to an assay.
#' @param assay Assay to set before the function gets the feature values
#'   (default `"Spatial"`).
#' @param layer Expression layer to use for `feature` (default `"counts"`).
#' @param images Names of the Seurat images to use.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param save_path Path of a CSV file to save the vertices to.
#' @param plot_image Show the Seurat tissue image under the points (default
#'   `FALSE`)?
#' @param downsample_n Number of points to show in the app. If `NULL`, the app
#'   shows all points. A click still snaps to the nearest spatial ID in the full
#'   coordinate table.
#' @param downsample_seed Random seed for the downsampling (default 1).
#' @param point_size Size of the points (default 0.8).
#' @param point_alpha Starting opacity of the points (default 0.9).
#' @param show_zeros Show the points where `feature` is zero or `NA` (default
#'   `TRUE`)? On sparse Visium HD data, the zeros can hide the signal. Use
#'   `FALSE`, or clear "Show zero points" in the app, to show only the positive
#'   points. This changes the display only. A click can still select a hidden
#'   point as a vertex.
#' @param launch_browser Passed to `shiny::runApp()` (default `TRUE`).
#' @param host Host for `shiny::runApp()`. Set it when you run the function in
#'   an interactive HPC session through an SSH tunnel.
#' @param port Port for `shiny::runApp()`. Set it when you run the function in
#'   an interactive HPC session through an SSH tunnel.
#' @param quiet Hide warnings while the app runs (default `TRUE`)?
#'
#' @return A data frame with the columns `spatial_id`, `x`, `y`, `num`,
#'   `polygon`, and `image`.
#' @seealso [formatSpatialPolygons()], [addPolygons()], [plotSpatialPolygons()]
#'
#' @examplesIf interactive()
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Click the vertices in the app. Use "Close polygon" when you finish.
#' polygons <- selectSpatialPolygons(
#'   getSpatial(x),
#'   image_scale = "hires",
#'   swap_coords = FALSE
#' )
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
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  selected_images <- .sirup_spatial_images(object, images = images)
  if (length(selected_images) > 1L) {
    rlang::warn("Multiple spatial images are present; showing the first image in `selectSpatialPolygons()`.")
    selected_images <- selected_images[[1L]]
  }
  feature_object <- object
  if (.sirup_is_spatial_list(object)) {
    parsed_image <- .sirup_parse_spatial_image_label(selected_images[[1L]])
    if (is.null(parsed_image) || !parsed_image$sample %in% names(object)) {
      rlang::abort("Multi-spatial image labels must use '<sample>::<image>'.")
    }
    feature_object <- object[[parsed_image$sample]]
  }
  if (!"package:shiny" %in% search()) {
    suppressPackageStartupMessages(base::attachNamespace("shiny"))
  }
  old_shiny_useragg <- getOption("shiny.useragg")
  options(shiny.useragg = FALSE)
  on.exit(options(shiny.useragg = old_shiny_useragg), add = TRUE)
  assay_choices <- tryCatch(SeuratObject::Assays(feature_object), error = function(e) character())
  selected_assay <- assay %||% tryCatch(SeuratObject::DefaultAssay(feature_object), error = function(e) NULL)
  if (!is.null(selected_assay)) {
    SeuratObject::DefaultAssay(feature_object) <- selected_assay
  }
  if (length(assay_choices) == 0L) {
    assay_choices <- selected_assay %||% ""
  }
  layer_choices <- .sirup_assay_layers(feature_object, selected_assay)

  coords <- .sirup_spatial_coordinate_table(object, images = selected_images, image_scale = image_scale, swap_coords = swap_coords)
  if (nrow(coords) == 0L) {
    rlang::abort("No spatial coordinates were found in `object`.")
  }
  coords$spatial_id <- coords$cell
  coords$.sirup_plot_y <- -coords$y
  y_axis_offset <- min(coords$.sirup_plot_y, na.rm = TRUE)
  plot_coords <- .sirup_downsample_polygon_plot_coords(
    coords = coords,
    downsample_n = downsample_n,
    seed = downsample_seed
  )
  image_names <- unique(stats::na.omit(as.character(coords$image)))
  image_names <- setdiff(image_names, "metadata")
  plot_image_info <- NULL
  if (length(image_names) > 0L) {
    spatial_image <- .sirup_spatial_data(
      object,
      image = image_names[[1L]],
      image_scale = image_scale,
      swap_coords = swap_coords
    )
    raster <- .sirup_image_raster(spatial_image$image, alpha = 1)
    image_dim <- dim(raster)
    plot_image_info <- list(
      raster = raster,
      width = image_dim[[2]],
      height = image_dim[[1]]
    )
  }
  # One index for the whole session. Built here, not in the click observer.
  nn_index <- .sirup_nn_index(cbind(coords$x, coords$.sirup_plot_y))
  rlang::inform(paste0("Polygon selector nearest-neighbour backend: ", nn_index$backend))

  x_range <- diff(range(plot_coords$x, na.rm = TRUE))
  y_range <- diff(range(plot_coords$.sirup_plot_y, na.rm = TRUE))
  plot_height <- if (!is.finite(x_range) || !is.finite(y_range) || x_range <= 0 || y_range <= 0) {
    500
  } else {
    max(500, min(1200, round(900 * y_range / x_range)))
  }
  default_feature <- feature
  if (is.null(default_feature)) {
    features <- tryCatch(rownames(feature_object), error = function(e) character())
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

    selected_points <- shiny::reactiveVal(.sirup_empty_polygon_vertices())
    current_polygon_num <- shiny::reactiveVal(1L)
    current_polygon_name <- shiny::reactiveVal("1")

    shiny::observeEvent(input$assay, {
      shiny::updateSelectInput(
        session,
        "layer",
        choices = .sirup_assay_layers(feature_object, input$assay),
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
      expr <- .sirup_polygon_feature_values(
        object = feature_object,
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
        ggplot2::aes(x = .data$x, y = .data$.sirup_plot_y, color = .data[[feature_name]], key = .data$cell)
      } else {
        ggplot2::aes(x = .data$x, y = .data$.sirup_plot_y, key = .data$cell)
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
              ggplot2::aes(x = .data$x, y = .data$.sirup_plot_y, group = .data$num),
              inherit.aes = FALSE,
              linewidth = 0.5,
              color = "black"
            )
        }
        p <- p +
          ggplot2::geom_point(
            data = selected,
            ggplot2::aes(x = .data$x, y = .data$.sirup_plot_y, group = .data$num),
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
        .sirup_plot_y = closest$.sirup_plot_y,
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
      current <- .sirup_public_polygon_vertices(current)
      current[, setdiff(colnames(current), "image"), drop = FALSE]
    })

    shiny::observeEvent(input$close_polygon, {
      selected_points(.sirup_close_current_polygon(selected_points(), current_polygon_num()))
    })

    shiny::observeEvent(input$new_polygon, {
      current <- selected_points()
      if (!any(current$num == current_polygon_num())) {
        return()
      }
      current <- .sirup_close_current_polygon(current, current_polygon_num())
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
      selected_points(.sirup_empty_polygon_vertices())
      current_polygon_num(1L)
      current_polygon_name("1")
      shiny::updateTextInput(session, "polygon_name", value = "1")
    })

    shiny::observeEvent(input$done, {
      out <- selected_points()
      for (polygon_num in unique(out$num)) {
        out <- .sirup_close_current_polygon(out, polygon_num)
      }
      out <- .sirup_public_polygon_vertices(out)
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

#' Format point sets as polygon vertices
#'
#' `formatSpatialPolygons` converts a table of points into the vertex table that
#' [selectSpatialPolygons()] returns.
#'
#' @details
#' By default, `polygon_points` holds all points inside each polygon. The
#' function merges the points with `sf`. Then it maps the outer ring
#' back to the nearest selected spatial ID coordinates. The boundary keeps a
#' valid order for concave polygons, and the vertices are spatial ID
#' coordinates. If the input already holds boundary points, use
#' `point_type = "boundary"`.
#'
#' @param object Seurat object or SIRUP object.
#' @param polygon_points Data frame of points. With `point_type = "interior"`,
#'   it holds all points inside each polygon. With `point_type = "boundary"`, it
#'   holds the boundary points.
#' @param spatial_id Name of the column in `polygon_points` that holds the
#'   spatial IDs (default `"spatial_id"`). You need it if `polygon_points` has
#'   no `x` and `y` columns.
#' @param polygon Name of the column in `polygon_points` that identifies the
#'   polygons (default `"polygon"`).
#' @param x Name of the x-coordinate column in `polygon_points`, if it exists
#'   (default `"x"`).
#' @param y Name of the y-coordinate column in `polygon_points`, if it exists
#'   (default `"y"`).
#' @param image Name of the image column in `polygon_points`, if it exists
#'   (default `"image"`).
#' @param point_type Type of the points in `polygon_points`: `interior` (all
#'   points inside) or `boundary`.
#' @param buffer_radius Radius for the buffer around the points. If `NULL`, the
#'   function estimates it from the spacing between nearest neighbors in the
#'   full coordinate grid.
#' @param order Name of a column that orders the points in each polygon. It is
#'   most useful with `point_type = "boundary"`.
#' @param images Names of the Seurat images to use.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object.
#' @param close_polygons Close each polygon by repeating the first vertex, if
#'   needed (default `TRUE`)?
#'
#' @return A data frame with the columns `spatial_id`, `x`, `y`, `polygon`, and
#'   `image`. The function drops polygons with fewer than three points and gives
#'   a warning.
#' @seealso [selectSpatialPolygons()], [addPolygons()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Four spots that lie inside one polygon.
#' points <- data.frame(
#'   spatial_id = c("spot04", "spot05", "spot09", "spot10"),
#'   polygon = "block"
#' )
#' formatSpatialPolygons(x, points, image_scale = "hires", swap_coords = FALSE)
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
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
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
    polygon_points$.sirup_input_order <- seq_len(nrow(polygon_points))
    ordered_rows <- unlist(lapply(unique(polygon_points[[polygon]]), function(polygon_id) {
      idx <- which(polygon_points[[polygon]] == polygon_id)
      idx[order(polygon_points[[order]][idx], polygon_points$.sirup_input_order[idx])]
    }), use.names = FALSE)
    polygon_points <- polygon_points[ordered_rows, , drop = FALSE]
    polygon_points$.sirup_input_order <- NULL
  }

  has_spatial_id <- spatial_id %in% colnames(polygon_points)
  has_xy <- all(c(x, y) %in% colnames(polygon_points))
  if (!has_spatial_id && !has_xy) {
    rlang::abort("`polygon_points` must contain spatial_id values or x/y coordinate columns.")
  }

  coords <- .sirup_spatial_coordinate_table(
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
    if (is.null(buffer_radius)) {
      buffer_radius <- .sirup_default_polygon_buffer_radius(coords)
    }
    out <- .sirup_polygon_interior_to_sf_boundary(
      points = out,
      buffer_radius = buffer_radius
    )
  }

  if (isTRUE(close_polygons)) {
    out <- .sirup_close_formatted_polygons(out)
  }
  rownames(out) <- NULL
  out
}

#' Expand polygons outward from their centroids
#'
#' `expandPolygons` makes each polygon larger by moving its vertices away from
#' the polygon centroid.
#'
#' @details
#' Each vertex moves along the line from the centroid through the vertex. It
#' moves by the distance `delta`. The polygon keeps its shape and its number of
#' vertices. The polygons must be convex.
#'
#' @param polygons Polygons to expand. Give a flat data frame of vertices (with
#'   `x`, `y`, and, if you want, a polygon ID column). You can also pass a SIRUP
#'   object. The function then expands the stored polygons.
#' @param delta Distance that each vertex moves outward, in coordinate units.
#'   Required. A negative value moves the vertices inward, so the
#'   polygon becomes smaller.
#' @param x Name of the x-coordinate column, if `polygons` is a data frame
#'   (default `"x"`).
#' @param y Name of the y-coordinate column, if `polygons` is a data frame
#'   (default `"y"`).
#' @param polygon_id Name of the column that separates the polygons, if
#'   `polygons` is a data frame (default `"polygon"`). If the column does not
#'   exist, the function uses a `polygon` or `polygon_id` column. If neither
#'   exists, all rows form one polygon.
#'
#' @return The expanded polygons in the form of the input: a data frame for a
#'   data frame, or the updated SIRUP object for a SIRUP object.
#' @seealso [addPolygons()], [plotSpatialPolygons()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Make the stored polygons larger by 0.25 coordinate units.
#' x2 <- expandPolygons(x, delta = 0.25)
#'
#' # Use a data frame of vertices.
#' flat <- data.frame(
#'   polygon = "block",
#'   x = c(3.5, 5.5, 5.5, 3.5),
#'   y = c(0.5, 0.5, 2.5, 2.5)
#' )
#' expandPolygons(flat, delta = 0.25)
#' @export
expandPolygons <- function(polygons,
                           delta = NULL,
                           x = "x",
                           y = "y",
                           polygon_id = "polygon") {
  if (is.null(delta)) {
    rlang::abort("`delta` must be supplied (a single finite number).")
  }
  if (isSIRUP(polygons)) {
    object <- polygons
    stored <- .sirup_polygons(object)
    if (nrow(stored) == 0L) {
      rlang::abort("The SIRUP object does not contain stored polygons.")
    }
    flat <- .sirup_polygon_vertices_table(stored)
    # The stored vertices are moved coordinates, not spatial barcodes, so drop
    # the carried spatial IDs before re-storing the expanded polygons.
    flat$spatial_id <- NULL
    expanded <- .sirup_expand_polygon_df(flat, delta = delta, x = "x", y = "y", polygon_id = "polygon")
    return(addPolygons(
      object,
      polygons = expanded,
      include_boundary = isTRUE(object$metadata$polygons_include_boundary %||% TRUE)
    ))
  }
  polygons <- .as_data_frame(polygons)
  .sirup_expand_polygon_df(polygons, delta = delta, x = x, y = y, polygon_id = polygon_id)
}

# Expand every polygon in a flat vertex table. Only the `x`/`y` columns change;
# all other columns and the row order are kept. See `expandPolygons()`.
.sirup_expand_polygon_df <- function(polygon, delta, x = "x", y = "y", polygon_id = "polygon") {
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
    moved <- .sirup_expand_polygon_vertices(px, py, delta = delta)
    out[[x]][idx] <- moved$x
    out[[y]][idx] <- moved$y
  }
  out
}

# Move one polygon's vertices outward from its centroid by `delta`. A repeated
# closing vertex equals the first vertex, so it moves identically and the
# polygon stays closed. Returns a list with the moved `x`/`y` vectors.
.sirup_expand_polygon_vertices <- function(px, py, delta) {
  n_all <- length(px)
  closed <- n_all > 1L &&
    .sirup_same_xy(
      data.frame(x = px[[1L]], y = py[[1L]]),
      data.frame(x = px[[n_all]], y = py[[n_all]])
    )
  core <- if (closed) seq_len(n_all - 1L) else seq_len(n_all)
  if (length(core) < 3L) {
    rlang::abort("Each polygon must contain at least three vertices.")
  }
  centroid <- .sirup_polygon_centroid(px[core], py[core])
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
.sirup_polygon_centroid <- function(x, y) {
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

#' Plot polygons on a tissue image
#'
#' @details
#' `plotSpatialPolygons` plots polygons over the Seurat tissue image. If there
#' are multiple samples or images, the function creates separate plots for each.
#'
#' @param object Seurat object or SIRUP object.
#' @param polygon Data frame of polygon vertices in order. It needs the columns
#'   set by `x` and `y`. If the column set by `polygon_id` exists, the function
#'   treats each polygon on its own. If `object` is a SIRUP object and `polygon`
#'   is `NULL`, the function uses the stored polygons.
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param polygon_val IDs of the polygons to use. If `NULL`, the function uses
#'   all polygons.
#' @param x Name of the x-coordinate column in `polygon` (default `"x"`).
#' @param y Name of the y-coordinate column in `polygon` (default `"y"`).
#' @param polygon_id Name of the column in `polygon` that groups the vertices
#'   by polygon (default `"polygon"`).
#' @param polygon_alpha Opacity of the polygon fill (default 0.2).
#' @param polygon_fill Fill color of the polygons, if no analysis is mapped
#'   (default `"#0072B2"`).
#' @param fill_analysis Name of a stored SIRUP analysis that has one row for
#'   each polygon. The function colors the polygons with it.
#' @param fill_column Name of the column that sets the fill color. If
#'   `fill_analysis` is `NULL`, the function reads the column from the polygon
#'   data. Otherwise, it reads the column from the analysis. The analysis must
#'   have a `polygon_id` or `polygon` column as the polygon key.
#' @param analysis_fill Two colors for the low and high end of the gradient for
#'   numeric analysis values (default `c("#D3D3D3", "#FFFF00")`). The function
#'   uses the Glasbey palette for other values.
#' @param polygon_color Color of the polygon outlines (default `"#0072B2"`).
#' @param linewidth Width of the polygon outlines (default 0.6).
#' @param image_alpha Opacity of the tissue image.
#' @param label Write a label at the center of each polygon (default `TRUE`)?
#' @param label_by Name of the column in the polygon data that holds the labels
#'   (default `"polygon"`).
#' @param label_size Text size of the labels (default 3).
#' @param label_color Text color of the labels (default `"black"`).
#' @param label_outside Put each label outside its polygon (default `FALSE`)?
#'   The function pushes each label away from the centroid of all plotted
#'   polygons, so that labels spread outward from a group of small polygons. The
#'   label is `label_offset` away from the polygon.
#' @param label_offset Gap between a polygon edge and its label, if
#'   `label_outside` is `TRUE` (default 0.08). The value is a fraction of the
#'   plotted coordinate range.
#' @param label_leader Draw a line from each polygon to its label, if
#'   `label_outside` is `TRUE` (default `TRUE`)?
#' @param label_leader_color Color of the leader lines (default `"grey40"`).
#' @param zoom Zoom to the extent of the polygons, not the full tissue (default
#'   `FALSE`)?
#' @param combine Combine the plots with patchwork, if possible (default
#'   `TRUE`)?
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#'
#' @return A ggplot or patchwork object. If `combine = FALSE`, a named list of
#'   ggplot objects.
#' @seealso [selectSpatialPolygons()], [spatialCloneDominance()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' plotSpatialPolygons(x)
#' plotSpatialPolygons(x, label_outside = TRUE)
#' @export
plotSpatialPolygons <- function(object,
                                polygon = NULL,
                                images = NULL,
                                sample_id = NULL,
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
                                combine = TRUE,
                                combine_sample_plots = FALSE) {
  sirup_object <- if (isSIRUP(object)) object else NULL
  if (!is.null(sample_id)) {
    sample_id <- .sirup_plot_sample_id(object, sample_id)
    images <- unique(c(images %||% character(), as.character(sample_id)))
  }
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (is.null(polygon)) {
    if (is.null(sirup_object)) {
      rlang::abort("`polygon` must be supplied when `object` is not a SIRUP object.")
    }
    polygon <- .sirup_polygon_vertices_table(.sirup_polygons(sirup_object))
    if (nrow(polygon) == 0L) {
      rlang::abort("The SIRUP object does not contain stored polygons.")
    }
  }
  polygon <- .as_data_frame(polygon)
  .check_required_cols(polygon, c(x, y), "polygon")
  if (!is.numeric(polygon_alpha) || length(polygon_alpha) != 1L || is.na(polygon_alpha) ||
      polygon_alpha < 0 || polygon_alpha > 1) {
    rlang::abort("`polygon_alpha` must be a number between 0 and 1.")
  }

  polygon_keep_cols <- unique(c(label_by, ".sirup_polygon_label", .sirup_polygon_sample_col(polygon)))
  if (!is.null(fill_column) && is.null(fill_analysis)) {
    polygon_keep_cols <- unique(c(polygon_keep_cols, fill_column))
  }
  polygon_data <- .sirup_polygon_plot_data(
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
      polygon_data <- .sirup_join_polygon_column_fill_data(
        polygon_data = polygon_data,
        fill_column = fill_column
      )
    } else {
      if (is.null(sirup_object)) {
        rlang::abort("`fill_analysis` requires a SIRUP object.")
      }
      polygon_data <- .sirup_join_polygon_fill_data(
        object = sirup_object,
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
    label_source <- label_by
    if (identical(label_by, "polygon") && ".sirup_polygon_label" %in% colnames(polygon_data)) {
      label_source <- ".sirup_polygon_label"
    }
    .check_required_cols(polygon_data, label_source, "polygon")
    label_idx <- !duplicated(polygon_data$polygon)
    label_lookup <- data.frame(
      polygon = polygon_data$polygon[label_idx],
      .sirup_label = as.character(polygon_data[[label_source]][label_idx]),
      stringsAsFactors = FALSE
    )
    label_data <- .sirup_polygon_centers(polygon_data, unique(polygon_data$polygon))
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
        diff(range(polygon_data$.sirup_y, na.rm = TRUE)),
        na.rm = TRUE
      )
      gx <- mean(label_data$x, na.rm = TRUE)
      gy <- mean(label_data$.sirup_y, na.rm = TRUE)
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
            .sirup_points_in_polygon(txs, tys, v$x, v$.sirup_y, include_boundary = TRUE),
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
        cy <- label_data$.sirup_y[i]
        v <- polygon_data[polygon_data$polygon == label_data$polygon[i], , drop = FALSE]
        radius <- max(sqrt((v$x - cx)^2 + (v$.sirup_y - cy)^2), na.rm = TRUE)
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
      label_data$.sirup_lx <- lx
      label_data$.sirup_ly <- ly
      label_data$.sirup_ex <- ex
      label_data$.sirup_ey <- ey
    }
  }

  plots <- list()
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    image_polygon_data <- .sirup_filter_polygon_plot_sample(polygon_data, image)
    image_label_data <- if (isTRUE(label)) {
      label_data[as.character(label_data$polygon) %in% unique(as.character(image_polygon_data$polygon)), , drop = FALSE]
    } else {
      NULL
    }
    # `plot_data` is the polygon geometry, so the spatial-barcode coordinates are
    # passed as `view_data`. This crops the plot to the sequenced spots (like
    # Seurat SpatialFeaturePlot) instead of the full tissue image.
    plots[[image]] <- .sirup_spatial_base(
      plot_data = image_polygon_data,
      image = spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = if (isTRUE(zoom)) NULL else spatial$coords
    ) +
      if (is.null(fill_column)) {
        ggplot2::geom_polygon(
          data = image_polygon_data,
          ggplot2::aes(x = .data$x, y = .data$.sirup_y, group = .data$polygon),
          inherit.aes = FALSE,
          fill = polygon_fill,
          color = polygon_color,
          alpha = polygon_alpha,
          linewidth = linewidth
        )
      } else {
        ggplot2::geom_polygon(
          data = image_polygon_data,
          ggplot2::aes(
            x = .data$x, y = .data$.sirup_y, group = .data$polygon,
            fill = .data$.sirup_fill_value
          ),
          inherit.aes = FALSE,
          color = polygon_color,
          alpha = polygon_alpha,
          linewidth = linewidth
        )
      }
    if (!is.null(fill_column)) {
      analysis_is_numeric <- is.numeric(polygon_data$.sirup_fill_value)
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
        levels <- sort(unique(stats::na.omit(as.character(polygon_data$.sirup_fill_value))))
        fill_values <- .sirup_discrete_palette_values(levels)
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
              data = image_label_data,
              ggplot2::aes(x = .data$.sirup_ex, y = .data$.sirup_ey,
                           xend = .data$.sirup_lx, yend = .data$.sirup_ly),
              inherit.aes = FALSE,
              color = label_leader_color,
              linewidth = 0.25
            )
        }
        plots[[image]] <- plots[[image]] +
          ggplot2::geom_text(
            data = image_label_data,
            ggplot2::aes(x = .data$.sirup_lx, y = .data$.sirup_ly, label = .data$.sirup_label),
            inherit.aes = FALSE,
            size = label_size,
            color = label_color,
            fontface = "bold"
          )
      } else {
        plots[[image]] <- plots[[image]] +
          ggplot2::geom_text(
            data = image_label_data,
            ggplot2::aes(x = .data$x, y = .data$.sirup_y, label = .data$.sirup_label),
            inherit.aes = FALSE,
            size = label_size,
            color = label_color,
            fontface = "bold"
          )
      }
    }
  }

  .sirup_return_spatial_plots(
    plots,
    combine = combine,
    combine_sample_plots = combine_sample_plots
  )
}

#' Plot receptors inside polygons using the results of [spatialCloneDominance()]
#'
#' @details
#' `plotSpatialRegionReceptors` plots the receptors inside the polygons on the
#' tissue image. The result of [spatialCloneDominance()] supplies the receptors.
#' If there are multiple samples or images, the function creates separate plots for each.
#'
#' @param dominance Result of [spatialCloneDominance()].
#' @param object Seurat object or SIRUP object.
#' @param color Name of the column in `dominance$region_receptors` that sets the
#'   color of the points. If you omit it, the function uses the clone column
#'   that [spatialCloneDominance()] recorded (default `"clone_id"`).
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param polygon_val IDs of the polygons to select from `dominance`.
#' @param point_size Size of the receptor points (default 1).
#' @param point_alpha Opacity of the receptor points (default 0.5).
#' @param polygon_alpha Opacity of the polygon fill (default 0.12).
#' @param polygon_fill Fill color of the polygons (default `"#FFFFFF"`).
#' @param polygon_color Color of the polygon outlines (default `"black"`).
#' @param image_alpha Opacity of the tissue image.
#' @param top_n Number of the most frequent `color` levels to show. The function
#'   groups the other levels under `other_label`. If `NULL`, it shows all
#'   levels.
#' @param other_label Label for the grouped levels (default `"Other"`).
#' @param other_color Color of the grouped levels (default `"#D3D3D3"`).
#' @param plot_other Plot the grouped levels (default `TRUE`)?
#' @param zoom Zoom to the extent of the selected polygons (default `FALSE`)? If
#'   `FALSE`, the plot shows the full tissue image.
#' @param combine Combine the plots with patchwork, if possible (default
#'   `TRUE`)?
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#'
#' @return A ggplot or patchwork object. If `combine = FALSE`, a named list of
#'   ggplot objects.
#' @seealso [spatialCloneDominance()], [plotSpatialPolygonPies()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' dominance <- spatialCloneDominance(x)
#' plotSpatialRegionReceptors(dominance, x)
#' @export
plotSpatialRegionReceptors <- function(dominance,
                                       object,
                                       color = "clone_id",
                                       images = NULL,
                                       sample_id = NULL,
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
                                       combine = TRUE,
                                       combine_sample_plots = FALSE) {
  color_missing <- missing(color)
  region <- .sirup_region_from_dominance(dominance, polygon_val = polygon_val)
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (color_missing && !is.null(dominance$metadata$receptor_fields$clone)) {
    color <- dominance$metadata$receptor_fields$clone
  }
  if (nrow(region$receptors) == 0L) {
    rlang::abort("No receptor rows were found inside the selected polygon(s).")
  }
  .check_required_cols(region$receptors, color, "receptors")
  region$receptors$.sirup_color_by <- .sirup_lump_values(region$receptors[[color]], top_n = top_n, other_label = other_label)
  if (!isTRUE(plot_other)) {
    region$receptors <- region$receptors[as.character(region$receptors$.sirup_color_by) != other_label, , drop = FALSE]
  }
  if (nrow(region$receptors) == 0L) {
    rlang::abort("No receptor rows remain to plot after applying `top_n` and `plot_other`.")
  }
  region$receptors$.sirup_is_other <- as.character(region$receptors$.sirup_color_by) == other_label
  region$receptors <- region$receptors[order(!region$receptors$.sirup_is_other), , drop = FALSE]
  color_values <- .sirup_discrete_palette_values(
    levels(region$receptors$.sirup_color_by),
    other_label = other_label,
    other_color = other_color
  )

  plots <- list()
  full_image_data <- data.frame(x = numeric(), .sirup_y = numeric())
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    plot_receptors <- region$receptors[region$receptors$image == image | region$receptors$image == "metadata", , drop = FALSE]
    plot_polygons <- .sirup_filter_polygon_plot_sample(region$polygons, image)
    plots[[image]] <- .sirup_spatial_base(
      plot_data = if (isTRUE(zoom) && nrow(plot_polygons) > 0L) plot_polygons else spatial$coords,
      image = spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = if (isTRUE(zoom)) NULL else full_image_data
    ) +
      ggplot2::geom_polygon(
        data = plot_polygons,
        ggplot2::aes(x = .data$x, y = .data$.sirup_y, group = .data$polygon),
        inherit.aes = FALSE,
        fill = polygon_fill,
        color = polygon_color,
        alpha = polygon_alpha,
        linewidth = 0.5
      ) +
      ggplot2::geom_point(
        data = plot_receptors,
        ggplot2::aes(x = .data$x, y = .data$.sirup_y, color = .data$.sirup_color_by),
        inherit.aes = FALSE,
        size = point_size,
        alpha = point_alpha,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = color_values, drop = TRUE) +
      ggplot2::labs(color = color)
  }

  .sirup_return_spatial_plots(plots, combine = combine, combine_sample_plots = combine_sample_plots)
}

#' Plot 3D clone-density surfaces on a tissue image
#'
#' `plotSpatialCloneDensity3D` draws a smoothed density surface for each clone
#' above the tissue image in an interactive rgl scene.
#'
#' @details
#' The function uses the receptor rows in the result of
#' [spatialCloneDominance()]. It needs the `rgl` package, except if you set
#' `return_data = TRUE`.
#'
#' @param dominance Result of [spatialCloneDominance()].
#' @param object Seurat object or SIRUP object.
#' @param clone Name of the column in `dominance$region_receptors` that defines
#'   the clone groups. If you omit it, the function uses the clone column that
#'   [spatialCloneDominance()] recorded (default `"clone_id"`).
#' @param weight Name of a numeric column that gives the weight of each
#'   receptor. If `NULL`, each receptor row counts as one.
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param polygon_val IDs of the polygons to select from `dominance`.
#' @param top_n Number of the most abundant clone groups to show on their own
#'   (default 10). The function groups the other clone groups under
#'   `other_label`.
#' @param other_label Label for the grouped clone groups (default `"Other"`).
#' @param other_color Color of the grouped clone groups (default `"#D3D3D3"`).
#' @param plot_other Plot the grouped clone groups (default `TRUE`)?
#' @param bandwidth Bandwidth of the Gaussian kernel, in coordinate units. Give
#'   one value for equal smoothing in x and y, or two values (x, y). If `NULL`,
#'   the function estimates it from the coordinates in the regions.
#' @param grid_n Number of grid points for the density surfaces (default 80).
#'   Give one value for both axes, or two values (x, y).
#' @param plot_type Type of plot: `density` (smoothed surfaces) or `histogram`
#'   (binned 3D bars).
#' @param hist_bins Number of x and y bins, if `plot_type = "histogram"`
#'   (default 30).
#' @param zoom Zoom the image plane and the 3D view to the selected region
#'   (default `FALSE`)? If `FALSE`, the plot shows the full image.
#' @param height_scale Maximum height of the density, as a fraction of the
#'   plotted range (default 0.22).
#' @param z_offset Vertical gap between the density surfaces and the tissue
#'   image, as a fraction of the plotted range (default 0.04).
#' @param density_alpha Opacity of the density surfaces (default 0.55).
#' @param density_floor Relative density below which the function hides surface
#'   cells (default 0.03). This removes the low Gaussian tails, which can look
#'   like a colored plane.
#' @param image_alpha Opacity of the tissue image plane.
#' @param show_points Show the receptor locations above the tissue plane
#'   (default `FALSE`)?
#' @param point_size Size of the receptor points (default 5).
#' @param point_alpha Opacity of the receptor points (default 0.75).
#' @param show_polygon Show the polygon boundaries above the tissue plane
#'   (default `TRUE`)?
#' @param polygon_color Color of the polygon boundary lines (default `"white"`).
#' @param polygon_lwd Width of the polygon boundary lines (default 3).
#' @param label_density_peaks Label the highest point of the first density
#'   surfaces (default `FALSE`)?
#' @param label_top_n Number of surfaces to label, if `label_density_peaks` is
#'   `TRUE` (default 5).
#' @param label_cex Text scale of the peak labels (default 0.8).
#' @param show_axes Show the 3D bounding box and the axis labels (default
#'   `TRUE`)?
#' @param texture_dir Directory to write the tissue texture PNG file to. Use a
#'   permanent directory if you save an HTML widget.
#' @param as_widget Return an rgl htmlwidget (default `TRUE`)? If `FALSE`, the
#'   function draws the scene and returns the ID of the current rgl device.
#' @param return_data Return the density data and do not draw the scene (default
#'   `FALSE`)? This does not need `rgl`.
#'
#' @return An rgl widget or an rgl device ID. If `return_data = TRUE`, a list of
#'   the plot data.
#' @seealso [spatialCloneDominance()], [plotSpatialRegionReceptors()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' dominance <- spatialCloneDominance(x)
#'
#' # Return the density data. This needs no rgl window.
#' density <- plotSpatialCloneDensity3D(dominance, x, return_data = TRUE)
#' names(density)
#' @export
plotSpatialCloneDensity3D <- function(dominance,
                                      object,
                                      clone = "clone_id",
                                      weight = NULL,
                                      images = NULL,
                                      sample_id = NULL,
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
  if (!inherits(dominance, "sirup_spatial_clone_dominance")) {
    rlang::abort("`dominance` must be an object returned by `spatialCloneDominance()`.")
  }
  if (is.null(dominance$region_receptors) || is.null(dominance$polygons)) {
    rlang::abort("`dominance` must contain `region_receptors` and `polygons` from `spatialCloneDominance()`.")
  }
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (clone_missing && !is.null(dominance$metadata$receptor_fields$clone)) {
    clone <- dominance$metadata$receptor_fields$clone
  }
  plot_type <- match.arg(plot_type)

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
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    image_receptors <- region_receptors[region_receptors$image == image | region_receptors$image == "metadata", , drop = FALSE]
    if (nrow(image_receptors) == 0L) {
      next
    }
    image_polygons <- .sirup_filter_polygon_plot_sample(polygons, image)
    plot_data[[image]] <- .sirup_clone_density_3d_data(
      region_receptors = image_receptors,
      polygons = image_polygons,
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
  if (!rlang::is_installed("rgl")) {
    rlang::abort("The optional `rgl` package is required for 3D clone-density plotting.")
  }
  plots <- lapply(plot_data, function(dat) {
    .sirup_draw_clone_density_3d(
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

#' Plot clonal dominance results within polygons as pies or bars
#'
#' `plotSpatialPolygonPies` summarizes one receptor column for each polygon. It
#' plots the result as pie charts or stacked bars.
#'
#' @details
#' The function uses the receptor rows in the result of
#' [spatialCloneDominance()]. The `mode` argument sets the plot style:
#' \itemize{
#'   \item `side_by_side`: pie charts next to a numbered map of the polygons.
#'   \item `spatial_clipped`: pie charts inside the polygon boundaries on
#'     the tissue. It needs the `sf` package.
#'   \item `dominance_only`: the pie charts only, with no map.
#'   \item `stacked_bar`: one stacked bar for each polygon in one panel, with no
#'     map. Each bar stacks the `fill` levels in dominance order.
#' }
#'
#' If there are multiple samples or images, the function creates separate plots for each.
#'
#' @param dominance Result of [spatialCloneDominance()].
#' @param object Seurat object or SIRUP object.
#' @param fill Name of the column in `dominance$region_receptors` to summarize.
#'   If you omit it, the function uses the clone column that
#'   [spatialCloneDominance()] recorded (default `"clone_id"`).
#' @param mode Plot style: `side_by_side`, `spatial_clipped`, `dominance_only`,
#'   or `stacked_bar`. See Details.
#' @param top_n Number of the most frequent `fill` levels to find in each
#'   polygon (default 10). The function shows the union of these levels on their
#'   own in all polygons. It groups the other levels under `other_label`.
#' @param label_top_n Number of levels to label in each polygon (default 0).
#' @param other_label Label for the grouped levels (default `"Other"`).
#' @param other_color Color of the grouped levels (default `"#D3D3D3"`).
#' @param plot_other Plot the grouped levels (default `TRUE`)?
#' @param images Names of the Seurat images to use.
#' @param sample_id Names of the samples to plot. Use it for multi-sample SIRUP
#'   objects. If `NULL`, the function plots all samples.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param polygon_val IDs of the polygons to select from `dominance`.
#' @param polygon_alpha Opacity of the polygon fill in the numbered map (default
#'   0.2).
#' @param image_alpha Opacity of the tissue image.
#' @param pie_alpha Opacity of the pies (default 0.9).
#' @param label_size Text size of the polygon labels and pie labels (default 3).
#' @param pies_per_page Largest number of pie charts on one page, if
#'   `mode = "side_by_side"` (default 6). If there are more polygons, the
#'   function returns a named list with one plot for each page.
#' @param pie_order Order of the pie charts, if `mode = "side_by_side"`:
#'   `"receptors"` (default) puts the polygon with the most receptors first,
#'   then the next, and so on. `"polygon"` keeps the order of the numbered map.
#'   Each pie keeps its map label. The function ignores this argument if `mode =
#'   "spatial_clipped"`.
#' @param spatial_legend Show one shared fill legend for the `spatial_clipped`
#'   plots (default `FALSE`)?
#' @param zoom Zoom to the extent of the plotted polygons (default `FALSE`)? If
#'   `FALSE`, the plot shows the full tissue image.
#' @param combine_sample_plots Combine the plots of all samples into one
#'   patchwork (default `FALSE`)? The function uses it only for multi-sample
#'   objects. If `FALSE`, it returns a named list with one plot for each sample.
#'
#' @return A ggplot object, a patchwork object, or a named list of ggplot
#'   objects.
#' @seealso [spatialCloneDominance()], [plotSpatialRegionReceptors()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' dominance <- spatialCloneDominance(x)
#' plotSpatialPolygonPies(dominance, x)
#' plotSpatialPolygonPies(dominance, x, mode = "stacked_bar")
#' @export
plotSpatialPolygonPies <- function(dominance,
                                   object,
                                   fill = "clone_id",
                                   mode = c("side_by_side", "spatial_clipped", "dominance_only", "stacked_bar"),
                                   top_n = 10,
                                   label_top_n = 0,
                                   other_label = "Other",
                                   other_color = "#D3D3D3",
                                   plot_other = TRUE,
                                   images = NULL,
                                   sample_id = NULL,
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
                                   zoom = FALSE,
                                   combine_sample_plots = FALSE) {
  fill_missing <- missing(fill)
  region <- .sirup_region_from_dominance(dominance, polygon_val = polygon_val)
  sample_id <- .sirup_plot_sample_id(object, sample_id)
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  image_scale <- context$image_scale
  swap_coords <- context$swap_coords
  if (fill_missing && !is.null(dominance$metadata$receptor_fields$clone)) {
    fill <- dominance$metadata$receptor_fields$clone
  }
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
  pie_data <- .sirup_polygon_pie_summary(
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
  fill_values <- .sirup_discrete_palette_values(
    unique(as.character(pie_data$level)),
    other_label = other_label,
    other_color = other_color
  )
  centers <- .sirup_polygon_centers(region$polygons, pie_data$polygon)

  if (mode == "dominance_only") {
    .sirup_dominance_only_polygon_pies(
      pie_data = pie_data, fill_values = fill_values,
      label_size = label_size, pie_order = pie_order
    )
  } else if (mode == "stacked_bar") {
    .sirup_stacked_bar_polygon_dominance(
      pie_data = pie_data, fill_values = fill_values,
      label_size = label_size, pie_order = pie_order
    )
  } else if (mode == "side_by_side") {
    .sirup_side_by_side_polygon_pies(
      object = object,
      images = images,
      sample_id = sample_id,
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
      zoom = zoom,
      combine_sample_plots = combine_sample_plots
    )
  } else {
    .sirup_spatial_clipped_polygon_pies(
      object = object,
      images = images,
      sample_id = sample_id,
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
      zoom = zoom,
      combine_sample_plots = combine_sample_plots
    )
  }
}

#' Subset a Seurat object to the spots inside polygons
#'
#' `seuratSubsetPolygon` keeps only the spatial IDs inside the polygons. It also
#' records the polygon of each spot in the metadata.
#'
#' @param object Seurat object, or a named list of Seurat objects.
#' @param polygon Data frame of polygon vertices in order. It needs the columns
#'   set by `x` and `y`. If the column set by `polygon_id` exists, the function
#'   treats each polygon on its own.
#' @param images Names of the Seurat images to use.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#' @param swap_coords Swap the two Seurat coordinate axes? You must give a
#'   value.
#' @param polygon_val IDs of the polygons to use. If `NULL`, the function uses
#'   all polygons.
#' @param x Name of the x-coordinate column in `polygon` (default `"x"`).
#' @param y Name of the y-coordinate column in `polygon` (default `"y"`).
#' @param polygon_id Name of the column in `polygon` that groups the vertices
#'   by polygon (default `"polygon"`).
#' @param metadata Name of the metadata column that holds the primary polygon of
#'   each spot (default `"polygon"`). A spot in overlapping polygons gets the
#'   smallest polygon by area.
#' @param include_boundary Count spots on a polygon edge as inside the polygon
#'   (default `TRUE`)?
#'
#' @return A Seurat object with only the spatial IDs inside the selected
#'   polygons. For a list of Seurat objects, a list of the subset objects. The
#'   function drops an object with no spot inside a polygon. It adds the
#'   `metadata` column to `@meta.data`.
#' @seealso [assignSpatialPolygons()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' flat <- data.frame(
#'   polygon = "block",
#'   x = c(3.5, 5.5, 5.5, 3.5),
#'   y = c(0.5, 0.5, 2.5, 2.5)
#' )
#'
#' # Keep the four spots inside the polygon.
#' subset <- seuratSubsetPolygon(
#'   getSpatial(x),
#'   polygon = flat,
#'   image_scale = "hires",
#'   swap_coords = FALSE
#' )
#' ncol(subset)
#' @export
seuratSubsetPolygon <- function(object,
                                polygon = NULL,
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
  subset_one <- function(spatial) {
    keep_cells <- rownames(spatial@meta.data)[!is.na(spatial@meta.data[[metadata]])]
    keep_cells <- keep_cells[keep_cells %in% colnames(spatial)]
    if (length(keep_cells) == 0L) {
      return(NULL)
    }
    suppressWarnings(spatial[, keep_cells], classes = "validationWarning")
  }
  if (.sirup_is_spatial_list(annotated)) {
    out <- lapply(annotated, subset_one)
    out <- out[!vapply(out, is.null, logical(1))]
    if (length(out) == 0L) {
      rlang::abort("No Seurat spatial IDs were found inside the selected polygon(s).")
    }
    return(out)
  }
  subset <- subset_one(annotated)
  if (is.null(subset)) {
    rlang::abort("No Seurat spatial IDs were found inside the selected polygon(s).")
  }
  subset
}

#' Add polygon IDs to the spatial metadata
#'
#' `assignSpatialPolygons` adds polygon columns to the Seurat metadata. It keeps
#' all spots.
#'
#' @details
#' A spot outside the selected polygons gets `NA` in the primary column. A spot
#' in overlapping polygons gets the smallest polygon by area. The function also
#' adds one logical column `polygon_<ID>` for each polygon. These columns show
#' membership in every polygon.
#'
#' @param object Seurat object or SIRUP object.
#' @param polygon Data frame of polygon vertices in order. It needs the columns
#'   set by `x` and `y`. If the column set by `polygon_id` exists, the function
#'   treats each polygon on its own. If `object` is a SIRUP object and `polygon`
#'   is `NULL`, the function uses the stored polygons.
#' @param images Names of the Seurat images to use.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param polygon_val IDs of the polygons to use. If `NULL`, the function uses
#'   all polygons.
#' @param x Name of the x-coordinate column in `polygon` (default `"x"`).
#' @param y Name of the y-coordinate column in `polygon` (default `"y"`).
#' @param polygon_id Name of the column in `polygon` that groups the vertices
#'   by polygon (default `"polygon"`).
#' @param metadata Name of the metadata column that holds the primary polygon of
#'   each spot (default `"polygon_id"`).
#' @param include_boundary Count spots on a polygon edge as inside the polygon
#'   (default `TRUE`)?
#'
#' @return A Seurat object with the `metadata` column added to `@meta.data`.
#' @seealso [seuratSubsetPolygon()], [assignReceptorPolygons()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' flat <- data.frame(
#'   polygon = "block",
#'   x = c(3.5, 5.5, 5.5, 3.5),
#'   y = c(0.5, 0.5, 2.5, 2.5)
#' )
#' spatial <- assignSpatialPolygons(
#'   getSpatial(x),
#'   polygon = flat,
#'   image_scale = "hires",
#'   swap_coords = FALSE
#' )
#' table(spatial$polygon_id, useNA = "ifany")
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
  if (!.sirup_is_non_empty_string(metadata)) {
    rlang::abort("`metadata` must be a non-empty character string.")
  }
  if (isSIRUP(object) && is.null(polygon)) {
    polygon <- .sirup_polygon_vertices_table(.sirup_polygons(object))
    if (nrow(polygon) == 0L) {
      rlang::abort("The SIRUP object does not contain stored polygons.")
    }
  }
  context <- .sirup_spatial_context(object, image_scale = image_scale, swap_coords = swap_coords)
  object <- context$object
  if (.sirup_is_spatial_list(object)) {
    out <- lapply(names(object), function(sample) {
      sample_polygon <- .sirup_filter_polygon_sample(polygon, sample)
      if (!.sirup_polygon_has_vertices(sample_polygon, x = x, y = y)) {
        object[[sample]]
      } else {
        assignSpatialPolygons(
          object = object[[sample]],
          polygon = sample_polygon,
          images = images,
          image_scale = context$image_scale,
          swap_coords = context$swap_coords,
          polygon_val = polygon_val,
          x = x,
          y = y,
          polygon_id = polygon_id,
          metadata = metadata,
          include_boundary = include_boundary
        )
      }
    })
    names(out) <- names(object)
    return(out)
  }
  membership <- .sirup_polygon_membership(
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

  polygon_lookup <- .sirup_polygon_lookup(membership)
  polygon_flags <- .sirup_polygon_flag_lookup(membership)
  object <- .sirup_update_seurat_image_slots(object)
  object@meta.data[[metadata]] <- unname(polygon_lookup[rownames(object@meta.data)])
  object@meta.data <- .sirup_append_polygon_flags(
    object@meta.data,
    rownames(object@meta.data),
    polygon_flags
  )
  object
}

#' Add polygon IDs to a receptor table
#'
#' `assignReceptorPolygons` adds polygon columns to a receptor table. It keeps
#' all rows.
#'
#' @details
#' If you give `object`, a receptor row gets `NA` if its spatial ID is outside
#' the selected polygons or does not match the spatial object. If `object` is
#' `NULL`, the function finds the polygons from the coordinate columns in
#' `receptors`.
#'
#' A row in overlapping polygons gets the smallest polygon by area in the
#' primary column. The function also adds one logical column `polygon_<ID>` for
#' each polygon. These columns show membership in every polygon.
#'
#' If `receptors` is a SIRUP object and you omit `object` and `polygon`, the
#' function uses the stored spatial object and polygons.
#'
#' @param receptors AIRR-style receptor table. You can give a SIRUP object here
#'   if you omit `object`.
#' @param object Seurat object or SIRUP object.
#' @param column Name of the column to add to the receptor table (default
#'   `"polygon_id"`).
#' @param spatial_id Name of the `receptors` column that holds the spatial IDs.
#'   If `NULL`, the function uses the stored field of a SIRUP object, or looks
#'   for a standard spatial ID column.
#' @param polygon Data frame of polygon vertices in order. It needs the columns
#'   set by `x` and `y`. If the column set by `polygon_id` exists, the function
#'   treats each polygon on its own. If `object` is a SIRUP object and `polygon`
#'   is `NULL`, the function uses the stored polygons.
#' @param images Names of the Seurat images to use.
#' @param image_scale Image scale for `SeuratObject::GetTissueCoordinates()`.
#'   If `NULL`, the function uses the value stored in the SIRUP object.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in the SIRUP object. Give a value if
#'   `object` is a Seurat object.
#' @param polygon_val IDs of the polygons to use. If `NULL`, the function uses
#'   all polygons.
#' @param x Name of the x-coordinate column in `polygon` (default `"x"`).
#' @param y Name of the y-coordinate column in `polygon` (default `"y"`).
#' @param receptor_x Name of the x-coordinate column in `receptors`, for
#'   `object = NULL` (default `"x"`).
#' @param receptor_y Name of the y-coordinate column in `receptors`, for
#'   `object = NULL` (default `"y"`).
#' @param polygon_id Name of the column in `polygon` that groups the vertices
#'   by polygon (default `"polygon"`).
#' @param include_boundary Count spots on a polygon edge as inside the polygon
#'   (default `TRUE`)?
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#'
#' @return A tibble with the `column` added.
#' @seealso [assignSpatialPolygons()], [withPolygons()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Use the spatial object and the polygons stored in the SIRUP object.
#' receptors <- assignReceptorPolygons(x)
#' table(receptors$polygon_id, useNA = "ifany")
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
                                   include_boundary = TRUE,
                                   slot = NULL) {
  if (!.sirup_is_non_empty_string(column)) {
    rlang::abort("`column` must be a non-empty character string.")
  }
  if (isSIRUP(receptors) && is.null(object)) {
    object <- receptors
    receptors <- NULL
  }
  if (isSIRUP(object) && is.null(polygon)) {
    polygon <- .sirup_polygon_vertices_table(.sirup_polygons(object))
    if (nrow(polygon) == 0L) {
      rlang::abort("The SIRUP object does not contain stored polygons.")
    }
  }
  if (is.null(object)) {
    if (is.null(polygon)) {
      rlang::abort("`polygon` must be supplied when `object = NULL`.")
    }
    receptors <- .as_data_frame(receptors)
    .check_required_cols(receptors, c(receptor_x, receptor_y), "receptors")
    receptors$.sirup_receptor_row <- seq_len(nrow(receptors))
    coords <- data.frame(
      cell = as.character(receptors$.sirup_receptor_row),
      image = NA_character_,
      x = as.numeric(receptors[[receptor_x]]),
      y = as.numeric(receptors[[receptor_y]]),
      stringsAsFactors = FALSE
    )
    membership <- .sirup_polygon_membership_from_coords(
      coords = coords,
      polygon = polygon,
      polygon_val = polygon_val,
      x = x,
      y = y,
      polygon_id = polygon_id,
      include_boundary = include_boundary
    )
    polygon_lookup <- .sirup_polygon_lookup(membership)
    polygon_flags <- .sirup_polygon_flag_lookup(membership)
    out <- tibble::as_tibble(receptors)
    out[[column]] <- unname(polygon_lookup[as.character(out$.sirup_receptor_row)])
    out <- .sirup_append_polygon_flags(
      out,
      as.character(out$.sirup_receptor_row),
      polygon_flags
    )
    out$.sirup_receptor_row <- NULL
    return(out)
  }
  options <- .sirup_spatial_options(object, image_scale = image_scale, swap_coords = swap_coords)
  inputs <- .sirup_resolve_inputs(object, receptors, spatial_id, slot)
  object <- inputs$object
  receptors <- .as_data_frame(inputs$receptors)
  spatial_id <- .guess_spatial_id_col(receptors, inputs$spatial_id_col)
  .check_required_cols(receptors, spatial_id, "receptors")

  if (.sirup_is_spatial_list(object)) {
    sample_col <- inputs$receptor_fields$sample %||% "sample"
    .check_required_cols(receptors, sample_col, "receptors")
    out <- tibble::as_tibble(receptors)
    out[[column]] <- NA_character_
    for (sample in names(object)) {
      keep <- !is.na(out[[sample_col]]) & as.character(out[[sample_col]]) == sample
      if (!any(keep)) next
      sample_polygon <- .sirup_filter_polygon_sample(polygon, sample)
      if (!.sirup_polygon_has_vertices(sample_polygon, x = x, y = y)) next
      assigned <- assignReceptorPolygons(
        receptors = out[keep, , drop = FALSE],
        object = object[[sample]],
        polygon = sample_polygon,
        column = column,
        spatial_id = spatial_id,
        images = images,
        image_scale = options$image_scale,
        swap_coords = options$swap_coords,
        polygon_val = polygon_val,
        x = x,
        y = y,
        receptor_x = receptor_x,
        receptor_y = receptor_y,
        polygon_id = polygon_id,
        include_boundary = include_boundary
      )
      missing_cols <- setdiff(colnames(assigned), colnames(out))
      for (missing_col in missing_cols) {
        out[[missing_col]] <- if (is.logical(assigned[[missing_col]])) FALSE else NA
      }
      extra_cols <- setdiff(colnames(out), colnames(assigned))
      for (extra_col in extra_cols) {
        assigned[[extra_col]] <- if (is.logical(out[[extra_col]])) FALSE else NA
      }
      assigned <- assigned[, colnames(out), drop = FALSE]
      out[keep, ] <- assigned
    }
    return(out)
  }

  membership <- .sirup_polygon_membership(
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
  polygon_lookup <- .sirup_polygon_lookup(membership)
  polygon_flags <- .sirup_polygon_flag_lookup(membership)
  out <- tibble::as_tibble(receptors)
  out[[column]] <- unname(polygon_lookup[as.character(out[[spatial_id]])])
  out <- .sirup_append_polygon_flags(
    out,
    as.character(out[[spatial_id]]),
    polygon_flags
  )
  out
}

.sirup_filter_polygon_sample <- function(polygon, sample_id) {
  polygon <- .as_data_frame(polygon)
  sample_col <- .sirup_polygon_sample_col(polygon)
  if (is.null(sample_col)) {
    return(polygon)
  }
  keep <- is.na(polygon[[sample_col]]) | as.character(polygon[[sample_col]]) == sample_id
  polygon[keep, , drop = FALSE]
}

.sirup_polygon_sample_col <- function(polygon) {
  sample_col <- intersect(c("sample", "spatial_sample"), colnames(polygon))
  if (length(sample_col) == 0L) NULL else sample_col[[1L]]
}

.sirup_polygon_sample_values <- function(polygon) {
  sample_col <- .sirup_polygon_sample_col(polygon)
  if (is.null(sample_col)) {
    return(character())
  }
  samples <- unique(as.character(stats::na.omit(polygon[[sample_col]])))
  samples[nzchar(samples)]
}

.sirup_polygon_display_id <- function(polygon_data, polygon_value) {
  if (".sirup_polygon_label" %in% colnames(polygon_data)) {
    labels <- unique(stats::na.omit(as.character(polygon_data$.sirup_polygon_label[
      as.character(polygon_data$polygon) == as.character(polygon_value)
    ])))
    labels <- labels[nzchar(labels)]
    if (length(labels) > 0L) {
      return(labels[[1L]])
    }
  }
  parsed <- .sirup_parse_spatial_image_label(as.character(polygon_value))
  if (!is.null(parsed)) {
    return(parsed$image)
  }
  as.character(polygon_value)
}

.sirup_polygon_spec_rows <- function(polygon, spec, polygon_col) {
  if (isTRUE(spec$add_polygon_col)) {
    polygon[as.character(polygon[[polygon_col]]) == as.character(spec$id), , drop = FALSE]
  } else {
    polygon
  }
}

.sirup_filter_polygon_plot_sample <- function(polygons, image) {
  sample_col <- .sirup_polygon_sample_col(polygons)
  parsed <- .sirup_parse_spatial_image_label(image)
  if (is.null(sample_col)) {
    return(polygons)
  }
  if (is.null(parsed)) {
    sample <- as.character(image)
    sample_values <- .sirup_polygon_sample_values(polygons)
    if (!sample %in% sample_values) {
      return(polygons)
    }
  } else {
    sample <- parsed$sample
  }
  keep <- is.na(polygons[[sample_col]]) | as.character(polygons[[sample_col]]) == sample
  polygons[keep, , drop = FALSE]
}

.sirup_polygon_has_vertices <- function(polygon, x = "x", y = "y") {
  is.data.frame(polygon) && all(c(x, y) %in% colnames(polygon)) && nrow(polygon) >= 3L
}

#' Summarize columns in receptors data inside polygons
#'
#' `polygonReceptorSummary` makes one summary row for each polygon.
#'
#' @details
#' The input can be a SIRUP object or a receptor table that already has polygon
#' columns. For a SIRUP object, the function stores the summary as an analysis
#' and returns the updated object. For a receptor table, it returns the summary
#' tibble.
#'
#' @param object SIRUP object, or an AIRR-style receptor table.
#' @param polygons Polygon table that sets the order of the output polygons. If
#'   `NULL`, the function takes the polygon IDs from the SIRUP object, or from
#'   the `polygon_id` values in `object`.
#' @param spatial_id Name of the column that holds the spatial IDs. If `NULL`,
#'   the function uses the stored field of a SIRUP object, or looks for a
#'   standard spatial ID column in the receptor table.
#' @param polygon_id Name of the column that holds the primary polygon of each
#'   receptor (default `"polygon"`).
#' @param clone Name of the column that holds the clone IDs (default
#'   `"clone_id"`). The function uses it for `n_clones` and
#'   `top_clone_fraction`. If the column does not exist, the function leaves out
#'   these two summaries.
#' @param weight Name of a numeric column that holds the receptor weights. If
#'   you give it, `n_receptors`, `top_clone_fraction`, and the numeric custom
#'   summaries use the weights.
#' @param summary_columns Names of the receptor columns to summarize.
#' @param summary_operations Operations for `summary_columns`. Give one
#'   operation, or one for each column. The operations are `count`, `sum`,
#'   `mean`, `median`, `variance`, `var`, `sd`, `min`, `max`, `any`, `all`, and
#'   `cat`.
#' @param include_total_spots Add `n_total_spatial_ids` (default `FALSE`)? This
#'   column holds the number of all spots in each polygon, including spots with
#'   no receptor rows. It works only for SIRUP objects.
#' @param analysis_name Name for the stored summary in a SIRUP object (default
#'   `"polygon_receptor_summary"`). If the object has several receptor slots, the default also ends
#'   with `_` and the slot name.
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#' @param validate Validate the updated SIRUP object (default `TRUE`)?
#'
#' @return For a SIRUP object, the updated SIRUP object. Otherwise, a summary
#'   tibble.
#' @seealso [assignReceptorPolygons()], [polygonDiversity()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- polygonReceptorSummary(x)
#' getAnalysis(x, "polygon_receptor_summary")
#' @export
polygonReceptorSummary <- function(object,
                                   polygons = NULL,
                                   spatial_id = NULL,
                                   polygon_id = "polygon",
                                   clone = "clone_id",
                                   weight = NULL,
                                   summary_columns = NULL,
                                   summary_operations = NULL,
                                   include_total_spots = FALSE,
                                   analysis_name = NULL,
                                   slot = NULL,
                                   validate = TRUE) {
  is_sirup <- isSIRUP(object)
  x <- if (is_sirup) object else NULL
  polygon_id_missing <- missing(polygon_id)
  polygon_column <- if (is_sirup && isTRUE(polygon_id_missing)) "polygon_id" else polygon_id
  receptor_table <- if (is_sirup) getReceptors(x, slot = slot) else tibble::as_tibble(object)
  if (is_sirup) {
    fields <- .sirup_receptor_fields(x)
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

  polygon_ids <- .sirup_polygon_ids(receptor_table, polygons = polygons, polygon_id = polygon_column)
  summary <- .sirup_polygon_receptor_summary_table(
    receptors = receptor_table,
    polygon_ids = polygon_ids,
    spatial_id = spatial_id,
    polygon_id = polygon_column,
    clone = clone,
    weight = weight,
    summary_columns = summary_columns,
    summary_operations = summary_operations
  )
  if (is_sirup && isTRUE(include_total_spots)) {
    spot_counts <- .sirup_polygon_spot_counts(x, polygon_ids = polygon_ids)
    summary <- merge(summary, spot_counts, by = "polygon_id", all.x = TRUE, sort = FALSE)
    summary$n_total_spatial_ids[is.na(summary$n_total_spatial_ids)] <- 0L
  } else if (!is_sirup && isTRUE(include_total_spots)) {
    rlang::warn("`include_total_spots = TRUE` requires a SIRUP object; ignoring it for receptor-table input.")
  }

  if (is_sirup) {
    analysis_name <- analysis_name %||% .sirup_slot_analysis_name(x, "polygon_receptor_summary", slot)
    return(setAnalysis(x, analysis_name, summary, validate = validate))
  }
  summary
}

.sirup_polygon_spot_counts <- function(x, polygon_ids) {
  spatial <- getSpatial(x)
  rows <- lapply(.sirup_spatial_names(spatial), function(sample_name) {
    spatial_object <- if (.sirup_is_spatial_list(spatial)) spatial[[sample_name]] else spatial
    metadata <- spatial_object[[]]
    if (!"polygon_id" %in% colnames(metadata)) {
      return(data.frame(polygon_id = character(), n_total_spatial_ids = integer()))
    }
    dat <- data.frame(
      polygon_id = as.character(metadata$polygon_id),
      spatial_id = rownames(metadata),
      sample = sample_name,
      stringsAsFactors = FALSE
    )
    dat <- dat[!is.na(dat$polygon_id) & nzchar(dat$polygon_id), , drop = FALSE]
    if (nrow(dat) == 0L) {
      return(data.frame(polygon_id = character(), n_total_spatial_ids = integer()))
    }
    agg <- stats::aggregate(
      paste(dat$sample, dat$spatial_id, sep = "\r"),
      by = list(polygon_id = dat$polygon_id),
      FUN = function(v) length(unique(v))
    )
    stats::setNames(agg, c("polygon_id", "n_total_spatial_ids"))
  })
  counts <- .sirup_bind_rows(rows)
  out <- data.frame(polygon_id = as.character(polygon_ids), stringsAsFactors = FALSE)
  out <- merge(out, counts, by = "polygon_id", all.x = TRUE, sort = FALSE)
  out$n_total_spatial_ids[is.na(out$n_total_spatial_ids)] <- 0L
  out$n_total_spatial_ids <- as.integer(out$n_total_spatial_ids)
  out
}

#' Calculate polygon centroids
#'
#' `polygonCentroids` returns the centroid of each polygon.
#'
#' @param object SIRUP object with stored polygons, or a table of polygon
#'   vertices.
#' @param polygons Table of polygon vertices. If you give it, you can omit
#'   `object`, or use `object` only for context.
#' @param polygon_id Name of the polygon ID column, for flat polygon input
#'   (default `"polygon"`).
#' @param x Name of the x-coordinate column, for flat polygon input (default
#'   `"x"`).
#' @param y Name of the y-coordinate column, for flat polygon input (default
#'   `"y"`).
#'
#' @return A tibble with one row for each polygon. The columns are `sample`,
#'   `polygon_id`, `polygon_label`, `x`, `y`, `.sirup_y`, and `n_vertices`.
#' @seealso [polygonPairwiseDistance()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' polygonCentroids(x)
#' @export
polygonCentroids <- function(object = NULL,
                             polygons = NULL,
                             polygon_id = "polygon",
                             x = "x",
                             y = "y") {
  if (isSIRUP(object)) {
    polygons <- polygons %||% .sirup_polygon_vertices_table(getPolygons(object))
    if (identical(polygon_id, "polygon") && "polygon" %in% colnames(polygons)) {
      polygon_id <- "polygon"
    }
  } else {
    polygons <- polygons %||% object
  }
  if (is.null(polygons)) {
    rlang::abort("Supply a SIRUP object with stored polygons or a polygon vertex table.")
  }
  polygons <- .as_data_frame(polygons)
  .check_required_cols(polygons, c(polygon_id, x, y), "polygons")
  sample_col <- .sirup_polygon_sample_col(polygons)
  label_col <- if (".sirup_polygon_label" %in% colnames(polygons)) ".sirup_polygon_label" else polygon_id
  groups <- split(polygons, as.character(polygons[[polygon_id]]), drop = TRUE)
  if (length(groups) == 0L) {
    return(tibble::tibble(sample = character(), polygon_id = character(), polygon_label = character(), x = numeric(), y = numeric(), .sirup_y = numeric(), n_vertices = integer()))
  }
  rows <- lapply(names(groups), function(id) {
    dat <- groups[[id]]
    centroid <- .sirup_polygon_centroid_xy(dat[[x]], dat[[y]])
    sample_value <- if (!is.null(sample_col)) {
      values <- unique(stats::na.omit(as.character(dat[[sample_col]])))
      if (length(values) > 0L) values[[1L]] else NA_character_
    } else {
      parsed <- .sirup_parse_spatial_image_label(id)
      if (is.null(parsed)) NA_character_ else parsed$sample
    }
    data.frame(
      sample = sample_value,
      polygon_id = id,
      polygon_label = as.character(dat[[label_col]][[1L]]),
      x = centroid[[1L]],
      y = centroid[[2L]],
      .sirup_y = -centroid[[2L]],
      n_vertices = nrow(dat),
      stringsAsFactors = FALSE
    )
  })
  tibble::as_tibble(dplyr::bind_rows(rows))
}

#' Calculate distances between polygon centroids
#'
#' `polygonPairwiseDistance` calculates the distance between the centroids of
#' each pair of polygons in a sample.
#'
#' @param object SIRUP object with stored polygons.
#' @param polygons Table of polygon vertices.
#' @param normalize_by_diagonal Divide each distance by the diagonal of the
#'   coordinate range of its sample (default `TRUE`)?
#' @param images Names of the images or spatial samples to use for the
#'   coordinate range.
#' @param image_scale Image scale for the Seurat image coordinates.
#' @param swap_coords Swap the x and y coordinates before the function
#'   calculates the coordinate range?
#'
#' @return A tibble with one row for each pair of polygons in a sample. The
#'   columns are `sample`, `polygon_id_1`, `polygon_id_2`, `polygon_label_1`,
#'   `polygon_label_2`, `distance`, `normalized_distance`, and
#'   `section_diagonal`.
#' @seealso [polygonCentroids()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' polygonPairwiseDistance(x)
#' @export
polygonPairwiseDistance <- function(object,
                                    polygons = NULL,
                                    normalize_by_diagonal = TRUE,
                                    images = NULL,
                                    image_scale = NULL,
                                    swap_coords = NULL) {
  centroids <- polygonCentroids(object = object, polygons = polygons)
  centroids$sample[is.na(centroids$sample)] <- "spatial"
  extents <- if (isTRUE(normalize_by_diagonal)) {
    spatialCoordinateExtents(
      object,
      by = "sample",
      images = images,
      image_scale = image_scale,
      swap_coords = swap_coords
    )
  } else {
    NULL
  }
  rows <- lapply(split(centroids, centroids$sample, drop = TRUE), function(dat) {
    if (nrow(dat) < 2L) {
      return(NULL)
    }
    pairs <- utils::combn(seq_len(nrow(dat)), 2L)
    out <- lapply(seq_len(ncol(pairs)), function(i) {
      a <- dat[pairs[1L, i], , drop = FALSE]
      b <- dat[pairs[2L, i], , drop = FALSE]
      distance <- sqrt((a$x - b$x)^2 + (a$y - b$y)^2)
      diagonal <- NA_real_
      if (!is.null(extents)) {
        match_idx <- match(a$sample, extents$sample)
        if (!is.na(match_idx)) {
          diagonal <- extents$diagonal[[match_idx]]
        }
      }
      data.frame(
        sample = a$sample,
        polygon_id_1 = a$polygon_id,
        polygon_id_2 = b$polygon_id,
        polygon_label_1 = a$polygon_label,
        polygon_label_2 = b$polygon_label,
        distance = distance,
        normalized_distance = if (is.finite(diagonal) && diagonal > 0) distance / diagonal else NA_real_,
        section_diagonal = diagonal,
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(out)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(tibble::tibble(sample = character(), polygon_id_1 = character(), polygon_id_2 = character(), polygon_label_1 = character(), polygon_label_2 = character(), distance = numeric(), normalized_distance = numeric(), section_diagonal = numeric()))
  }
  tibble::as_tibble(dplyr::bind_rows(rows))
}

.sirup_polygon_centroid_xy <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) >= 3L) {
    xy <- cbind(x, y)
    if (!all(xy[1L, ] == xy[nrow(xy), ])) {
      xy <- rbind(xy, xy[1L, ])
    }
    centroid <- tryCatch({
      geom <- .sirup_repair_polygon_sfc(sf::st_sfc(sf::st_polygon(list(xy))))
      if (is.null(geom)) {
        NULL
      } else {
        as.numeric(sf::st_coordinates(sf::st_centroid(geom))[1L, 1:2])
      }
    }, error = function(e) NULL)
    if (!is.null(centroid) && all(is.finite(centroid))) {
      return(centroid)
    }
  }
  c(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE))
}

#' Apply a function to the receptors of each polygon
#'
#' `withPolygons` splits the receptor table by polygon and applies a function to
#' the rows of each polygon. Or it applies the function once to all receptors in
#' any polygon.
#'
#' @param receptors SIRUP object, or an AIRR-style receptor table.
#' @param .f Function to apply.
#' @param ... Other arguments for `.f`.
#' @param polygons Polygon table that sets the order of the output polygons. If
#'   `NULL`, the function takes the polygon IDs from the SIRUP object, or from
#'   the `polygon_id` values in `receptors`.
#' @param polygon_id Name of the column that holds the primary polygon of each
#'   receptor (default `"polygon"`).
#' @param mode `"by_polygon"` runs `.f` once for each polygon. `"all_polygons"`
#'   runs `.f` once on all receptors in any polygon.
#' @param pass_polygon_id Pass the current polygon ID to `.f` as the argument
#'   `polygon_id`, if `.f` accepts it (default `FALSE`)?
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#'
#' @return For `mode = "by_polygon"`, a named list with one result for each
#'   polygon. Otherwise, the result of `.f`.
#' @seealso [assignReceptorPolygons()], [polygonAbundance()],
#'   [polygonDiversity()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Count the receptors in each polygon.
#' withPolygons(x, nrow)
#'
#' # Count the receptors in all polygons together.
#' withPolygons(x, nrow, mode = "all_polygons")
#' @export
withPolygons <- function(receptors,
                         .f,
                         ...,
                         polygons = NULL,
                         polygon_id = "polygon",
                         mode = c("by_polygon", "all_polygons"),
                         pass_polygon_id = FALSE,
                         slot = NULL) {
  mode <- match.arg(mode)
  is_sirup <- isSIRUP(receptors)
  x <- if (is_sirup) receptors else NULL
  polygon_id_missing <- missing(polygon_id)
  polygon_column <- if (is_sirup && isTRUE(polygon_id_missing)) "polygon_id" else polygon_id
  receptor_table <- if (is_sirup) getReceptors(x, slot = slot) else tibble::as_tibble(receptors)
  if (is_sirup) {
    polygons <- polygons %||% getPolygons(x)
  }
  polygon_ids <- .sirup_polygon_ids(receptor_table, polygons = polygons, polygon_id = polygon_column)
  groups <- .sirup_polygon_receptor_groups(
    receptors = receptor_table,
    polygon_ids = polygon_ids,
    polygon_id = polygon_column,
    mode = mode
  )
  fun <- match.fun(.f)

  if (identical(mode, "all_polygons")) {
    return(.sirup_call_polygon_function(fun, groups[[1L]], polygon_id = NULL, pass_polygon_id = pass_polygon_id, ...))
  }
  stats::setNames(
    lapply(names(groups), function(id) .sirup_call_polygon_function(fun, groups[[id]], polygon_id = id, pass_polygon_id = pass_polygon_id, ...)),
    names(groups)
  )
}

#' Estimate clonal abundance inside polygons
#'
#' `polygonAbundance` runs [alakazam::estimateAbundance()] once on the receptors
#' of all polygons, with the polygon as the group. alakazam then resamples every
#' polygon to the same number of receptors, so the curves of different polygons
#' are comparable. With `common_size = FALSE`, the function analyses each polygon
#' alone instead.
#'
#' @details
#' The function puts the receptors of all polygons in one table and passes the
#' polygon column to [alakazam::estimateAbundance()] as `group`. With
#' `uniform = TRUE` (the default), alakazam resamples each polygon to the size
#' of the smallest polygon that has at least `min_n` receptors. Use `max_n` to
#' set a smaller size. The function leaves out polygons with fewer than `min_n`
#' receptors. If no polygon has enough receptors, the function gives an error.
#'
#' The result is one `AbundanceCurve` object. Its `abundance` table has one
#' column of polygon IDs. [alakazam::plotAbundanceCurve()] plots all polygons
#' together. For a SIRUP object, the function stores the result as an analysis
#' and returns the updated object.
#'
#' With `common_size = FALSE`, the function runs
#' [alakazam::estimateAbundance()] once for each polygon. Each polygon is then
#' resampled to its own size, so polygons with more receptors can look more
#' diverse. A polygon with fewer than `min_n` receptors gives `NULL`. The result
#' is a named list with one `AbundanceCurve` for each polygon.
#'
#' The arguments from `clone` to `progress` go to
#' [alakazam::estimateAbundance()].
#'
#' @param receptors SIRUP object, or an AIRR-style receptor table.
#' @param polygons Polygon table that sets the polygons to use. If `NULL`, the
#'   function takes the polygon IDs from the SIRUP object, or from the
#'   `polygon_id` values in `receptors`.
#' @param polygon_id Name of the column that holds the primary polygon of each
#'   receptor (default `"polygon"`). For a SIRUP object, the default is
#'   `"polygon_id"`.
#' @param analysis_name Name for the stored result in a SIRUP object (default
#'   `"polygon_abundance"`). If the object has several receptor slots, the default also ends
#'   with `_` and the slot name.
#' @param clone Name of the column that holds the clone IDs (default
#'   `"clone_id"`).
#' @param copy Name of the column that holds the copy number of each sequence
#'   (default `NULL`). If `NULL`, the function measures clone abundance as the
#'   number of sequences. Otherwise, it uses the sum of the copy numbers in each
#'   clone.
#' @param common_size Resample all polygons to the same size (default `TRUE`)?
#'   If `FALSE`, the function analyses each polygon alone, as in earlier
#'   versions of SIRUP, and returns a list with one result for each polygon.
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#' @param min_n Smallest number of receptors that a polygon needs (default 10).
#'   The function drops a polygon with fewer receptors. With
#'   `common_size = TRUE`, the smallest polygon that stays sets the size of the
#'   resampling, so a larger `min_n` gives less noisy curves and keeps fewer
#'   polygons.
#' @param max_n Largest number of receptors to sample from each polygon
#'   (default `NULL`). If `NULL`, the size is the size of the smallest polygon
#'   that stays.
#' @param uniform Resample each polygon to the same number of receptors
#'   (default `TRUE`)? Passed to [alakazam::estimateAbundance()].
#' @param ci Confidence interval, between 0 and 1 (default 0.95).
#' @param nboot Number of bootstrap samples (default 2000).
#' @param cell_id Name of the column that holds the cell IDs (default `NULL`).
#'   If `NULL`, the function assumes bulk data.
#' @param progress Show a progress bar (default `FALSE`)?
#' @param validate Validate the updated SIRUP object (default `TRUE`)?
#'
#' @return For a SIRUP object, the updated SIRUP object. Otherwise, the
#'   abundance results: one `AbundanceCurve` object, or a list with one object
#'   for each polygon if `common_size = FALSE`.
#' @seealso [polygonDiversity()], [withPolygons()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # The example polygons have few receptors, so lower min_n.
#' x <- polygonAbundance(x, min_n = 2, nboot = 50)
#' listAnalyses(x)
#' @export
polygonAbundance <- function(receptors,
                             polygons = NULL,
                             polygon_id = "polygon",
                             analysis_name = NULL,
                             clone = "clone_id",
                             copy = NULL,
                             common_size = TRUE,
                             slot = NULL,
                             min_n = 10,
                             max_n = NULL,
                             uniform = TRUE,
                             ci = 0.95,
                             nboot = 2000,
                             cell_id = NULL,
                             progress = FALSE,
                             validate = TRUE) {
  clone_missing <- missing(clone)
  polygon_id_missing <- missing(polygon_id)
  if (isSIRUP(receptors) && isTRUE(clone_missing)) {
    clone <- .sirup_receptor_fields(receptors)$clone %||% clone
  }
  polygon_column <- if (isSIRUP(receptors) && isTRUE(polygon_id_missing)) "polygon_id" else polygon_id
  result <- .sirup_polygon_alakazam(
    receptors = receptors,
    polygons = polygons,
    polygon_id = polygon_column,
    common_size = common_size,
    threshold_n = min_n,
    fun = alakazam::estimateAbundance,
    slot = slot,
    clone = clone,
    copy = copy,
    min_n = min_n,
    max_n = max_n,
    uniform = uniform,
    ci = ci,
    nboot = nboot,
    cell_id = cell_id,
    progress = progress
  )
  if (isSIRUP(receptors)) {
    analysis_name <- analysis_name %||% .sirup_slot_analysis_name(receptors, "polygon_abundance", slot)
    return(setAnalysis(receptors, analysis_name, result, validate = validate))
  }
  result
}

#' Estimate clonal diversity inside polygons
#'
#' `polygonDiversity` runs [alakazam::alphaDiversity()] once on the receptors of
#' all polygons, with the polygon as the group. alakazam then resamples every
#' polygon to the same number of receptors, so the diversity of different
#' polygons is comparable. With `common_size = FALSE`, the function analyses each
#' polygon alone instead.
#'
#' @details
#' The function puts the receptors of all polygons in one table and passes the
#' polygon column to [alakazam::alphaDiversity()] as `group`. alakazam resamples
#' each polygon to the size of the smallest polygon that has at least `min_n`
#' receptors. Use `max_n` (through `...`) to set a smaller size. The function
#' leaves out polygons with fewer than `min_n` receptors. If no polygon has
#' enough receptors, the function gives an error.
#'
#' The result is one `DiversityCurve` object. Its `diversity` table has one
#' column of polygon IDs. [alakazam::plotDiversityCurve()] plots all polygons
#' together. For a SIRUP object, the function stores the result as an analysis
#' and returns the updated object.
#'
#' With `common_size = FALSE`, the function runs [alakazam::alphaDiversity()]
#' once for each polygon. Each polygon is then resampled to its own size, so
#' polygons with more receptors can look more diverse. A polygon with fewer than
#' `min_n` receptors gives `NULL`. The result is a named list with one
#' `DiversityCurve` for each polygon.
#'
#' The arguments from `clone` to `...` go to [alakazam::alphaDiversity()].
#'
#' @param receptors SIRUP object, or an AIRR-style receptor table.
#' @param polygons Polygon table that sets the polygons to use. If `NULL`, the
#'   function takes the polygon IDs from the SIRUP object, or from the
#'   `polygon_id` values in `receptors`.
#' @param polygon_id Name of the column that holds the primary polygon of each
#'   receptor (default `"polygon"`). For a SIRUP object, the default is
#'   `"polygon_id"`.
#' @param analysis_name Name for the stored result in a SIRUP object (default
#'   `"polygon_diversity"`). If the object has several receptor slots, the default also ends
#'   with `_` and the slot name.
#' @param clone Name of the column that holds the clone IDs (default
#'   `"clone_id"`).
#' @param copy Name of the column that holds the copy number of each sequence
#'   (default `NULL`). If `NULL`, the function measures clone abundance as the
#'   number of sequences. Otherwise, it uses the sum of the copy numbers in each
#'   clone.
#' @param common_size Resample all polygons to the same size (default `TRUE`)?
#'   If `FALSE`, the function analyses each polygon alone, as in earlier
#'   versions of SIRUP, and returns a list with one result for each polygon.
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. The function ignores `slot` if you give a plain receptor table.
#' @param min_n Smallest number of receptors that a polygon needs (default 10).
#'   The function drops a polygon with fewer receptors. With
#'   `common_size = TRUE`, the smallest polygon that stays sets the size of the
#'   resampling, so a larger `min_n` gives less noisy estimates and keeps fewer
#'   polygons.
#' @param min_q Smallest diversity order q (default 0).
#' @param max_q Largest diversity order q (default 4).
#' @param step_q Step between the values of q (default 0.1).
#' @param ci Confidence interval, between 0 and 1 (default 0.95).
#' @param ... Other arguments for [alakazam::alphaDiversity()], such as
#'   `nboot`, `max_n` and `uniform`. You cannot give `group`, because the
#'   function uses the polygon as the group.
#' @param validate Validate the updated SIRUP object (default `TRUE`)?
#'
#' @return For a SIRUP object, the updated SIRUP object. Otherwise, the
#'   diversity results: one `DiversityCurve` object, or a list with one object
#'   for each polygon if `common_size = FALSE`.
#' @seealso [summarizePolygonDiversity()], [polygonAbundance()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- polygonDiversity(x, min_n = 2, max_q = 2, step_q = 1)
#' listAnalyses(x)
#' @export
polygonDiversity <- function(receptors,
                             polygons = NULL,
                             polygon_id = "polygon",
                             analysis_name = NULL,
                             clone = "clone_id",
                             copy = NULL,
                             common_size = TRUE,
                             slot = NULL,
                             min_n = 10,
                             min_q = 0,
                             max_q = 4,
                             step_q = 0.1,
                             ci = 0.95,
                             ...,
                             validate = TRUE) {
  if ("group" %in% names(list(...))) {
    rlang::abort("`group` is not an argument of `polygonDiversity()`. The function uses the polygon as the group.")
  }
  clone_missing <- missing(clone)
  polygon_id_missing <- missing(polygon_id)
  if (isSIRUP(receptors) && isTRUE(clone_missing)) {
    clone <- .sirup_receptor_fields(receptors)$clone %||% clone
  }
  polygon_column <- if (isSIRUP(receptors) && isTRUE(polygon_id_missing)) "polygon_id" else polygon_id
  result <- .sirup_polygon_alakazam(
    receptors = receptors,
    polygons = polygons,
    polygon_id = polygon_column,
    common_size = common_size,
    threshold_n = min_n,
    fun = alakazam::alphaDiversity,
    slot = slot,
    clone = clone,
    copy = copy,
    min_n = min_n,
    min_q = min_q,
    max_q = max_q,
    step_q = step_q,
    ci = ci,
    ...
  )
  if (isSIRUP(receptors)) {
    analysis_name <- analysis_name %||% .sirup_slot_analysis_name(receptors, "polygon_diversity", slot)
    return(setAnalysis(receptors, analysis_name, result, validate = validate))
  }
  result
}

#' Summarize polygon diversity at one value of q
#'
#' `summarizePolygonDiversity` takes the diversity estimate for one diversity
#' order q from a stored analysis. It makes one table row for each polygon and
#' stores the table as a new analysis.
#'
#' @details
#' Polygons that [polygonDiversity()] left out, because they had fewer than
#' `min_n` receptors, get a row with `NA` values. [plotSpatialPolygons()] draws
#' these polygons in grey. The function reads both kinds of result: one
#' `DiversityCurve` for all polygons, and the list that
#' `polygonDiversity(common_size = FALSE)` makes.
#'
#' @param object SIRUP object.
#' @param q Diversity order to take, for example `2` for Simpson diversity.
#' @param diversity_analysis Name of the stored analysis from
#'   [polygonDiversity()]. If `NULL`, the function uses `"polygon_diversity"`,
#'   or `"polygon_diversity_<slot>"` if the object has several receptor slots.
#' @param analysis_name Name for the stored table. If `NULL`, the function uses
#'   `polygon_diversity_q` followed by the value of `q`, and the slot name if
#'   the object has several receptor slots.
#' @param tolerance Tolerance for the match of `q` (default
#'   `sqrt(.Machine$double.eps)`).
#' @param slot Name of the receptor slot whose diversity you summarize. The
#'   function uses it only to find the default analysis names. Leave it `NULL`
#'   if the object has one receptor table. You must give it if the object has
#'   several slots and you use the default names.
#' @param validate Validate the updated SIRUP object (default `TRUE`)?
#'
#' @return The updated SIRUP object, with the diversity summary as an analysis.
#' @seealso [polygonDiversity()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- polygonDiversity(x, min_n = 2, max_q = 2, step_q = 1)
#' x <- summarizePolygonDiversity(x, q = 1)
#' listAnalyses(x)
#' @export
summarizePolygonDiversity <- function(object,
                                      q,
                                      diversity_analysis = NULL,
                                      analysis_name = NULL,
                                      tolerance = sqrt(.Machine$double.eps),
                                      slot = NULL,
                                      validate = TRUE) {
  validateSIRUP(object)
  diversity_analysis <- diversity_analysis %||% .sirup_slot_analysis_name(object, "polygon_diversity", slot)
  if (!is.numeric(q) || length(q) != 1L || is.na(q)) {
    rlang::abort("`q` must be one numeric diversity order.")
  }
  if (is.null(analysis_name)) {
    analysis_name <- .sirup_slot_analysis_name(
      object,
      paste0("polygon_diversity_q", gsub("[^0-9A-Za-z]+", "_", as.character(q))),
      slot
    )
  }
  diversity <- getAnalysis(object, diversity_analysis)
  if (is.list(diversity) && !methods::is(diversity, "DiversityCurve")) {
    return(.sirup_summarize_diversity_list(object, diversity, q, analysis_name, tolerance, validate))
  }
  if (!methods::is(diversity, "DiversityCurve")) {
    rlang::abort("`diversity_analysis` must be a result of `polygonDiversity()`.")
  }
  dat <- .sirup_diversity_curve_table(diversity)
  group_column <- methods::slot(diversity, "group_by")
  if (!all(c("q", group_column) %in% colnames(dat))) {
    rlang::abort("The diversity table must contain a `q` column and a column of polygon IDs.")
  }
  dat <- dat[abs(as.numeric(dat$q) - q) <= tolerance, , drop = FALSE]
  dat$polygon_id <- as.character(dat[[group_column]])
  if (!identical(group_column, "polygon_id")) {
    dat[[group_column]] <- NULL
  }
  dat <- dat[!duplicated(dat$polygon_id), , drop = FALSE]
  # Keep polygons that polygonDiversity() left out, with NA values.
  all_ids <- as.character(getPolygons(object)$polygon_id)
  missing_ids <- setdiff(all_ids, dat$polygon_id)
  if (length(missing_ids) > 0L) {
    empty <- tibble::tibble(polygon_id = missing_ids, q = q)
    dat <- dplyr::bind_rows(dat, empty)
  }
  dat <- dat[match(unique(c(all_ids, dat$polygon_id)), dat$polygon_id), , drop = FALSE]
  dat <- dat[!is.na(dat$polygon_id), c("polygon_id", setdiff(colnames(dat), "polygon_id")), drop = FALSE]
  setAnalysis(object, analysis_name, tibble::as_tibble(dat), validate = validate)
}

# Summarize a per-polygon list of diversity curves (polygonDiversity(common_size = FALSE)).
.sirup_summarize_diversity_list <- function(object, diversity, q, analysis_name, tolerance, validate) {
  rows <- lapply(names(diversity), function(polygon) {
    curve <- diversity[[polygon]]
    if (is.null(curve)) {
      return(tibble::tibble(polygon_id = polygon, q = q, d = NA_real_, d_sd = NA_real_, d_lower = NA_real_, d_upper = NA_real_, e = NA_real_, e_lower = NA_real_, e_upper = NA_real_))
    }
    dat <- .sirup_diversity_curve_table(curve)
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
  setAnalysis(object, analysis_name, dplyr::bind_rows(rows), validate = validate)
}

# Return polygon IDs from stored polygons when available, otherwise from the
# receptor table annotation. This keeps empty polygons in summaries when their
# geometry is known.
.sirup_polygon_ids <- function(receptors, polygons = NULL, polygon_id = "polygon_id") {
  if (!is.null(polygons)) {
    polygons <- .sirup_prepare_polygons(polygons)
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
.sirup_polygon_receptor_groups <- function(receptors,
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
.sirup_polygon_receptor_summary_table <- function(receptors,
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
  operations <- .sirup_summary_operations(summary_columns, summary_operations)
  missing_summary <- setdiff(summary_columns %||% character(), colnames(receptors))
  if (length(missing_summary) > 0L) {
    rlang::abort(paste0(
      "`summary_columns` were not found in `receptors`: ",
      paste(missing_summary, collapse = ", ")
    ))
  }
  groups <- .sirup_polygon_receptor_groups(
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
    weights <- .sirup_summary_weights(dat, weight_col)
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
        out[[paste0(operation, "_", column)]] <- .sirup_apply_summary_operation(dat[[column]], operation, weights = weights)
      }
    }
    tibble::as_tibble(out)
  })
  dplyr::bind_rows(rows)
}

# Validate and recycle summary operations to match summary columns.
.sirup_summary_operations <- function(summary_columns, summary_operations) {
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
.sirup_apply_summary_operation <- function(x, operation, weights = NULL) {
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
      median = .sirup_weighted_median(num, w),
      variance = .sirup_weighted_variance(num, w),
      var = .sirup_weighted_variance(num, w),
      sd = sqrt(.sirup_weighted_variance(num, w)),
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
.sirup_summary_weights <- function(dat, weight_col) {
  if (is.null(weight_col)) {
    return(NULL)
  }
  weights <- suppressWarnings(as.numeric(dat[[weight_col]]))
  weights[is.na(weights) | weights < 0] <- 0
  weights
}

# Weighted median used by weighted polygon summaries.
.sirup_weighted_median <- function(x, w) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) >= sum(w) / 2)[[1L]]]
}

# Population-style weighted variance for weighted polygon summaries.
.sirup_weighted_variance <- function(x, w) {
  mu <- stats::weighted.mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

# Run an alakazam function on the receptors of polygons.
# common_size = TRUE: one call on all polygons, with the polygon column as the group.
#   alakazam resamples every group to the same size (the smallest group with at least
#   `threshold_n` rows, or `max_n`), so the results of different polygons are comparable.
# common_size = FALSE: one call for each polygon, as in earlier versions. A polygon with fewer
#   than `threshold_n` rows gives NULL, and the result is a named list.
.sirup_polygon_alakazam <- function(receptors,
                                    polygons,
                                    polygon_id,
                                    common_size = TRUE,
                                    threshold_n,
                                    fun,
                                    slot = NULL,
                                    ...) {
  is_sirup <- isSIRUP(receptors)
  receptor_table <- if (is_sirup) getReceptors(receptors, slot = slot) else tibble::as_tibble(receptors)
  if (is_sirup) {
    polygons <- polygons %||% getPolygons(receptors)
  }
  polygon_ids <- .sirup_polygon_ids(receptor_table, polygons = polygons, polygon_id = polygon_id)
  args <- list(...)
  if (!isTRUE(common_size)) {
    groups <- .sirup_polygon_receptor_groups(
      receptors = receptor_table,
      polygon_ids = polygon_ids,
      polygon_id = polygon_id,
      mode = "by_polygon"
    )
    run_one <- function(dat) {
      if (nrow(dat) < threshold_n) {
        return(NULL)
      }
      one_args <- args
      if ("cell_id" %in% names(one_args) && is.null(one_args$cell_id)) {
        one_args$cell_id <- .sirup_absent_column_name(dat, ".sirup_bulk_cell_id")
      }
      do.call(fun, c(list(data = dat), one_args))
    }
    return(lapply(groups, run_one))
  }
  dat <- .sirup_polygon_receptor_groups(
    receptors = receptor_table,
    polygon_ids = polygon_ids,
    polygon_id = polygon_id,
    mode = "all_polygons"
  )[[1L]]
  dat[[polygon_id]] <- as.character(dat[[polygon_id]])
  sizes <- table(dat[[polygon_id]])
  if (!any(sizes >= threshold_n)) {
    rlang::abort(paste0(
      "No polygon has at least ", threshold_n, " receptors. Lower `min_n`."
    ))
  }
  if ("cell_id" %in% names(args) && is.null(args$cell_id)) {
    args$cell_id <- .sirup_absent_column_name(dat, ".sirup_bulk_cell_id")
  }
  do.call(fun, c(list(data = dat, group = polygon_id), args))
}

# Choose a column name guaranteed not to exist in a data frame.
.sirup_absent_column_name <- function(dat, prefix) {
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
.sirup_call_polygon_function <- function(fun, dat, polygon_id = NULL, pass_polygon_id = FALSE, ...) {
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
.sirup_stored_polygon_coords <- function(coords, spatial, polygon_value, polygon_id = "polygon_id") {
  if (.sirup_is_spatial_list(spatial)) {
    out <- lapply(names(spatial), function(sample) {
      sample_coords <- coords
      if ("sample" %in% colnames(sample_coords)) {
        sample_coords <- sample_coords[as.character(sample_coords$sample) == sample, , drop = FALSE]
      }
      .sirup_stored_polygon_coords(
        coords = sample_coords,
        spatial = spatial[[sample]],
        polygon_value = polygon_value,
        polygon_id = polygon_id
      )
    })
    return(.sirup_bind_rows(out))
  }
  metadata <- spatial[[]]
  keep <- .sirup_polygon_membership_vector(metadata, rownames(metadata), polygon_value, polygon_id = polygon_id)
  coord_keep <- keep[as.character(coords$cell)]
  coord_keep[is.na(coord_keep)] <- FALSE
  coords[coord_keep, , drop = FALSE]
}

# Return receptor rows for one stored polygon annotation.
.sirup_receptors_in_polygon <- function(receptors, polygon_value, polygon_id = "polygon_id") {
  keep <- .sirup_polygon_membership_vector(receptors, seq_len(nrow(receptors)), polygon_value, polygon_id = polygon_id)
  receptors[keep, , drop = FALSE]
}

# Read a stored primary polygon membership vector. Per-polygon flags are not
# used by analyses; they are available for users to derive a different
# `polygon_id` column if they want overlapping/nested behavior.
.sirup_polygon_membership_vector <- function(dat, ids, polygon_value, polygon_id = "polygon_id") {
  if (!polygon_id %in% colnames(dat)) {
    rlang::abort("Stored polygon annotations were not found. Re-run `addPolygons()` on the SIRUP object.")
  }
  out <- as.character(dat[[polygon_id]]) == as.character(polygon_value)
  out[is.na(out)] <- FALSE
  stats::setNames(out, as.character(ids))
}

# Extract the tibble of diversity estimates from Alakazam DiversityCurve objects
# and simple data-frame/list fallbacks.
.sirup_diversity_curve_table <- function(curve) {
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
.sirup_polygon_membership <- function(object,
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
  coords <- .sirup_spatial_coordinate_table(
    object,
    images = images,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  coords <- coords[is.finite(coords$x) & is.finite(coords$y), , drop = FALSE]
  if (nrow(coords) == 0L) {
    rlang::abort("No finite spatial coordinates were found in `object`.")
  }
  .sirup_polygon_membership_from_coords(
    coords = coords,
    polygon = polygon,
    polygon_val = polygon_val,
    x = x,
    y = y,
    polygon_id = polygon_id,
    include_boundary = include_boundary
  )
}

.sirup_polygon_membership_from_coords <- function(coords,
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
  polygon_specs <- .sirup_polygon_specs(polygon, polygon_id = polygon_val, polygon_col = polygon_id)
  polygon_values <- vapply(polygon_specs, function(spec) {
    as.character(if (isTRUE(spec$add_polygon_col)) spec$id else "polygon_1")
  }, character(1))
  polygon_columns <- .sirup_polygon_flag_columns(polygon_values)
  membership <- lapply(seq_along(polygon_specs), function(i) {
    spec <- polygon_specs[[i]]
    vertices <- .sirup_select_polygon_vertices(
      polygon = polygon,
      polygon_id = spec$id,
      polygon_col = polygon_id,
      x_col = x,
      y_col = y
    )
    inside <- .sirup_points_in_polygon(
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
      sample = if ("sample" %in% colnames(coords)) as.character(coords$sample[inside]) else NA_character_,
      polygon = rep(polygon_value, sum(inside)),
      polygon_area = rep(.sirup_polygon_area(vertices$x, vertices$y), sum(inside)),
      polygon_order = rep(i, sum(inside)),
      polygon_column = rep(polygon_columns[[i]], sum(inside)),
      image = as.character(coords$image[inside]),
      x = as.numeric(coords$x[inside]),
      y = as.numeric(coords$y[inside]),
      stringsAsFactors = FALSE
    )
  })
  out <- .sirup_bind_rows(membership)
  attr(out, "polygon_columns") <- polygon_columns
  out
}

# Collapse one-row-per-hit polygon membership into a named lookup vector from
# spatial ID to a primary polygon ID. Overlapping IDs resolve to the smallest
# polygon by area, then to polygon input order for exact ties.
.sirup_polygon_lookup <- function(membership) {
  if (nrow(membership) == 0L) {
    return(stats::setNames(character(), character()))
  }
  membership <- membership[order(membership$cell, membership$polygon_area, membership$polygon_order), , drop = FALSE]
  membership <- membership[!duplicated(membership$cell), , drop = FALSE]
  stats::setNames(as.character(membership$polygon), as.character(membership$cell))
}

# Return one logical flag column per polygon ID, keyed by spatial ID.
.sirup_polygon_flag_lookup <- function(membership) {
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

.sirup_append_polygon_flags <- function(out, spatial_ids, polygon_flags) {
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

.sirup_polygon_flag_columns <- function(polygon_id) {
  make.names(paste0("polygon_", as.character(polygon_id)), unique = TRUE)
}

.sirup_polygon_area <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) > 1L && .sirup_same_xy(data.frame(x = x[[1L]], y = y[[1L]]), data.frame(x = x[[length(x)]], y = y[[length(y)]]))) {
    x <- x[-length(x)]
    y <- y[-length(y)]
  }
  j <- c(seq_along(x)[-1L], 1L)
  abs(sum(x * y[j] - x[j] * y)) / 2
}

# Repair older Seurat image objects before subsetting or metadata updates so
# polygon helpers do not fail on images with legacy slot layouts.
.sirup_update_seurat_image_slots <- function(object) {
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

.sirup_select_polygon_vertices <- function(polygon,
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
  closed <- .sirup_same_xy(out[1, , drop = FALSE], out[nrow(out), , drop = FALSE])
  if (!closed) {
    out <- rbind(out, out[1, , drop = FALSE])
  }
  out
}

.sirup_polygon_specs <- function(polygon, polygon_id, polygon_col) {
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

.sirup_rbind <- function(results, name) {
  tables <- lapply(results, `[[`, name)
  if (length(tables) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  out
}

.sirup_v_genes <- function(v_calls) {
  if (length(v_calls) == 0L) {
    return(character())
  }
  genes <- unlist(strsplit(as.character(v_calls), "[,;]"), use.names = FALSE)
  genes <- trimws(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  sort(unique(genes))
}

.sirup_v_gene_family <- function(v_genes) {
  sub("\\*.*$", "", v_genes)
}

.sirup_v_gene_superfamily <- function(v_gene_families) {
  sub("^([A-Za-z]+V[0-9]+).*$", "\\1", v_gene_families)
}

.sirup_points_in_polygon <- function(px,
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
  if (.sirup_same_xy(data.frame(x = vx[[1]], y = vy[[1]]), data.frame(x = vx[[length(vx)]], y = vy[[length(vy)]]))) {
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
    boundary <- boundary | .sirup_points_on_segment(px, py, vx[[j]], vy[[j]], vx[[i]], vy[[i]], tol = tol)
    j <- i
  }
  if (isTRUE(include_boundary)) {
    inside <- inside | boundary
  } else {
    inside <- inside & !boundary
  }
  inside
}

.sirup_points_on_segment <- function(px, py, x1, y1, x2, y2, tol) {
  dx <- x2 - x1
  dy <- y2 - y1
  cross <- (px - x1) * dy - (py - y1) * dx
  scale <- pmax(abs(dx), abs(dy), 1)
  on_line <- abs(cross) <= tol * scale
  within_x <- px >= pmin(x1, x2) - tol & px <= pmax(x1, x2) + tol
  within_y <- py >= pmin(y1, y2) - tol & py <= pmax(y1, y2) + tol
  on_line & within_x & within_y
}

.sirup_same_xy <- function(a, b) {
  isTRUE(all.equal(as.numeric(a$x), as.numeric(b$x), tolerance = sqrt(.Machine$double.eps))) &&
    isTRUE(all.equal(as.numeric(a$y), as.numeric(b$y), tolerance = sqrt(.Machine$double.eps)))
}

.sirup_empty_polygon_vertices <- function() {
  data.frame(
    spatial_id = character(),
    x = numeric(),
    y = numeric(),
    .sirup_plot_y = numeric(),
    num = integer(),
    polygon = character(),
    image = character(),
    stringsAsFactors = FALSE
  )
}

.sirup_public_polygon_vertices <- function(points) {
  if (!"num" %in% colnames(points) && "polygon" %in% colnames(points)) {
    points$num <- match(as.character(points$polygon), unique(as.character(points$polygon)))
    points <- points[, c(setdiff(colnames(points), "polygon"), "polygon"), drop = FALSE]
  }
  points[, setdiff(colnames(points), ".sirup_plot_y"), drop = FALSE]
}

.sirup_assay_layers <- function(object, assay = NULL) {
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

.sirup_polygon_plot_data <- function(polygon,
                                     polygon_id,
                                     polygon_col,
                                     x_col,
                                     y_col,
                                     keep_cols = character()) {
  specs <- .sirup_polygon_specs(polygon, polygon_id = polygon_id, polygon_col = polygon_col)
  keep_cols <- unique(keep_cols[keep_cols %in% colnames(polygon)])
  out <- lapply(specs, function(spec) {
    polygon_rows <- if (isTRUE(spec$add_polygon_col)) {
      polygon[as.character(polygon[[polygon_col]]) == as.character(spec$id), , drop = FALSE]
    } else {
      polygon
    }
    vertices <- .sirup_select_polygon_vertices(
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
  out$.sirup_y <- -out$y
  out
}

.sirup_polygon_vertices_table <- function(polygons) {
  .sirup_validate_polygons(polygons)
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

.sirup_join_polygon_column_fill_data <- function(polygon_data, fill_column) {
  if (!fill_column %in% names(polygon_data)) {
    rlang::abort(paste0(
      "`fill_column` was not found in polygons. `plotSpatialPolygons()` looks ",
      "for `fill_column` in the SIRUP object's polygons slot, or in the supplied ",
      "`polygon` data frame. To color polygons by a stored analysis instead, ",
      "specify `fill_analysis`; the selected analysis must contain a `polygon_id` ",
      "or `polygon` column that matches the plotted polygons and a `",
      fill_column, "` column."
    ))
  }
  polygon_data$.sirup_fill_value <- polygon_data[[fill_column]]
  polygon_data
}

.sirup_join_polygon_fill_data <- function(object,
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
    .sirup_fill_value = analysis[[fill_column]],
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(values$polygon)) {
    rlang::abort("The selected polygon analysis must contain one row per polygon.")
  }
  polygon_data$polygon <- as.character(polygon_data$polygon)
  matched <- match(polygon_data$polygon, values$polygon)
  polygon_data$.sirup_fill_value <- values$.sirup_fill_value[matched]
  polygon_data
}

.sirup_region_from_dominance <- function(dominance, polygon_val = NULL) {
  if (!inherits(dominance, "sirup_spatial_clone_dominance")) {
    rlang::abort("`dominance` must be an object returned by `spatialCloneDominance()`.")
  }
  if (is.null(dominance$region_receptors) || is.null(dominance$polygons)) {
    rlang::abort("`dominance` must contain `region_receptors` and `polygons` from `spatialCloneDominance()`.")
  }
  region_receptors <- .as_data_frame(dominance$region_receptors)
  polygon_plot_data <- .as_data_frame(dominance$polygons)
  if (!is.null(polygon_val)) {
    # Match polygon_val against the display `polygon` label OR the section-unique
    # `polygon_id` (composite "<sample>::<id>"). Multi-sample objects share clean
    # polygon labels across sections, so selecting by polygon_id is what keeps a
    # per-section selection unambiguous; single-sample callers still match `polygon`.
    pv <- as.character(polygon_val)
    sel_rows <- function(df) {
      keep <- as.character(df$polygon) %in% pv
      if ("polygon_id" %in% names(df)) keep <- keep | (as.character(df$polygon_id) %in% pv)
      df[keep, , drop = FALSE]
    }
    region_receptors <- sel_rows(region_receptors)
    polygon_plot_data <- sel_rows(polygon_plot_data)
  }
  list(
    receptors = region_receptors,
    polygons = polygon_plot_data
  )
}

.sirup_drop_plot_coordinate_cols <- function(dat) {
  dat[, setdiff(colnames(dat), c("x", "y", ".sirup_y", "image")), drop = FALSE]
}

.sirup_add_v_gene_columns <- function(dat, v_call_col = "v_call") {
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
      vapply(as.character(dat[[v_call_col]]), function(x) paste(.sirup_v_genes(x), collapse = ","), character(1))
    } else {
      rep(NA_character_, nrow(dat))
    }
    dat$v_gene[dat$v_gene == ""] <- NA_character_
  }
  if (!"v_gene_family" %in% colnames(dat)) {
    dat$v_gene_family <- if (v_call_col %in% colnames(dat)) {
      vapply(as.character(dat[[v_call_col]]), function(v_calls) {
        families <- unique(stats::na.omit(.sirup_v_gene_family(.sirup_v_genes(v_calls))))
        families <- sort(families)
        if (length(families) == 0L) NA_character_ else paste(families, collapse = ",")
      }, character(1))
    } else {
      .sirup_v_gene_family(dat$v_gene)
    }
  }
  if (!"v_gene_superfamily" %in% colnames(dat)) {
    dat$v_gene_superfamily <- if (v_call_col %in% colnames(dat)) {
      vapply(as.character(dat[[v_call_col]]), function(v_calls) {
        superfamilies <- .sirup_v_genes(v_calls)
        superfamilies <- .sirup_v_gene_family(superfamilies)
        superfamilies <- unique(stats::na.omit(.sirup_v_gene_superfamily(superfamilies)))
        superfamilies <- sort(superfamilies)
        if (length(superfamilies) == 0L) NA_character_ else paste(superfamilies, collapse = ",")
      }, character(1))
    } else {
      .sirup_v_gene_superfamily(dat$v_gene_family)
    }
  }
  dat
}

.sirup_bind_rows <- function(dots) {
  dots <- Filter(function(x) !is.null(x) && nrow(x) > 0L, dots)
  if (length(dots) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, dots)
  rownames(out) <- NULL
  out
}

.sirup_lump_values <- function(x, top_n = NULL, other_label = "Other") {
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

.sirup_discrete_palette_values <- function(levels, other_label = "Other", other_color = "#D3D3D3") {
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

.sirup_clone_density_3d_data <- function(region_receptors,
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
  .check_required_cols(region_receptors, c(clone, "x", "y", ".sirup_y", "polygon"), "dominance$region_receptors")
  .check_required_cols(polygons, c("x", "y", ".sirup_y", "polygon"), "dominance$polygons")
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
    is.finite(region_receptors$x) & is.finite(region_receptors$y) & is.finite(region_receptors$.sirup_y),
    ,
    drop = FALSE
  ]
  if (nrow(region_receptors) == 0L) {
    rlang::abort("No finite regional receptor coordinates are available.")
  }
  region_receptors$.sirup_density_group <- .sirup_lump_values(region_receptors[[clone]], top_n = top_n, other_label = other_label)
  if (!isTRUE(plot_other)) {
    region_receptors <- region_receptors[as.character(region_receptors$.sirup_density_group) != other_label, , drop = FALSE]
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
  region_receptors$.sirup_density_weight <- weights

  image_dim <- if (is.null(image)) {
    .sirup_image_dim(image, spatial_coords = spatial_coords)
  } else {
    .sirup_image_dim(image)
  }
  region_xlim <- .sirup_density_plot_range(c(region_receptors$x, polygons$x), c(0, image_dim$width))
  region_ylim <- .sirup_density_plot_range(c(region_receptors$.sirup_y, polygons$.sirup_y), c(-image_dim$height, 0))
  object_xlim <- .sirup_density_plot_range(spatial_coords$x, c(0, image_dim$width))
  object_ylim <- .sirup_density_plot_range(spatial_coords$.sirup_y, c(-image_dim$height, 0))
  image_xlim <- if (isTRUE(zoom)) region_xlim else object_xlim
  image_ylim <- if (isTRUE(zoom)) region_ylim else object_ylim
  image_extent <- if (is.null(image)) {
    list(rows = 1:2, cols = 1:2, xlim = image_xlim, ylim = image_ylim)
  } else {
    .sirup_image_texture_extent(image_dim, image_xlim, image_ylim)
  }
  span <- max(diff(region_xlim), diff(region_ylim))
  z_base <- z_offset * span
  max_height <- height_scale * span
  levels <- levels(region_receptors$.sirup_density_group)
  levels <- levels[levels %in% as.character(region_receptors$.sirup_density_group)]
  color_values <- .sirup_discrete_palette_values(levels, other_label = other_label, other_color = other_color)

  if (identical(plot_type, "histogram")) {
    surfaces <- list()
    bandwidth <- NULL
    histograms <- .sirup_clone_histogram_3d_bars(
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
    grid <- .sirup_density_grid(region_xlim, region_ylim, grid_n)
    bandwidth <- .sirup_density_bandwidth(region_receptors$x, region_receptors$.sirup_y, bandwidth)
    raw_surfaces <- lapply(levels, function(level) {
      dat <- region_receptors[as.character(region_receptors$.sirup_density_group) == level, , drop = FALSE]
      .sirup_weighted_gaussian_grid(
        x = dat$x,
        y = dat$.sirup_y,
        weights = dat$.sirup_density_weight,
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

.sirup_draw_clone_density_3d <- function(dat,
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
  .sirup_rgl("open3d")(useNULL = isTRUE(as_widget))
  .sirup_rgl("clear3d")()

  texture <- .sirup_image_texture_file(dat$image, dat$image_extent, alpha = image_alpha, texture_dir = texture_dir)
  x_plane <- dat$image_extent$xlim
  y_plane <- dat$image_extent$ylim
  z_plane <- matrix(0, nrow = 2, ncol = 2)
  texture_s <- matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2)
  texture_t <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
  .sirup_rgl("surface3d")(
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
    .sirup_rgl("surface3d")(
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
      .sirup_draw_histogram_bars_3d(histogram, alpha = density_alpha)
    }
  }

  if (isTRUE(show_polygon) && nrow(dat$polygons) > 0L) {
    for (polygon_id in unique(dat$polygons$polygon)) {
      polygon <- dat$polygons[dat$polygons$polygon == polygon_id, , drop = FALSE]
      .sirup_rgl("lines3d")(
        polygon$x,
        polygon$.sirup_y,
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
        .sirup_histogram_bar_peak(item)
      } else {
        .sirup_density_surface_peak(item)
      }
      if (!is.null(peak)) {
        .sirup_rgl("text3d")(
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
    point_colors <- unname(dat$colors[as.character(dat$receptors$.sirup_density_group)])
    point_colors[is.na(point_colors)] <- "#222222"
    .sirup_rgl("points3d")(
      dat$receptors$x,
      dat$receptors$.sirup_y,
      rep(dat$z_base * 0.65, nrow(dat$receptors)),
      color = grDevices::adjustcolor(point_colors, alpha.f = point_alpha),
      size = point_size
    )
  }

  if (isTRUE(show_axes)) {
    xlim <- range(c(dat$image_extent$xlim, dat$region_xlim), na.rm = TRUE)
    ylim <- range(c(dat$image_extent$ylim, dat$region_ylim), na.rm = TRUE)
    .sirup_rgl("bbox3d")(
      xat = pretty(xlim, n = 4),
      yat = pretty(ylim, n = 4),
      zat = pretty(c(0, dat$z_top), n = 3),
      color = "#666666",
      emission = "#333333"
    )
    .sirup_rgl("decorate3d")(xlab = "x", ylab = "y", zlab = "clone density")
  }
  .sirup_rgl("aspect3d")(1, 1, 0.55)
  .sirup_rgl("view3d")(theta = 40, phi = 25, zoom = 0.72)

  if (isTRUE(as_widget)) {
    .sirup_rgl("rglwidget")()
  } else {
    .sirup_rgl("cur3d")()
  }
}

.sirup_rgl <- function(name) {
  getExportedValue("rgl", name)
}

.sirup_image_dim <- function(image, spatial_coords = NULL) {
  if (is.null(image)) {
    if (!is.null(spatial_coords) && nrow(spatial_coords) > 0L) {
      xlim <- range(spatial_coords$x, na.rm = TRUE)
      ylim <- range(spatial_coords$.sirup_y, na.rm = TRUE)
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

.sirup_density_surface_peak <- function(surface) {
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

.sirup_histogram_bar_peak <- function(histogram) {
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

.sirup_density_plot_range <- function(x, fallback) {
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

.sirup_density_grid <- function(xlim, ylim, grid_n) {
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

.sirup_density_bandwidth <- function(x, y, bandwidth) {
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

.sirup_weighted_gaussian_grid <- function(x, y, weights, x_grid, y_grid, bandwidth) {
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

.sirup_clone_histogram_3d_bars <- function(region_bcrs,
                                           levels,
                                           color_values,
                                           xlim,
                                           ylim,
                                           hist_bins,
                                           z_base,
                                           max_height) {
  grid <- .sirup_density_grid(xlim, ylim, hist_bins + 1L)
  x_breaks <- grid$x
  y_breaks <- grid$y
  x_bins <- length(x_breaks) - 1L
  y_bins <- length(y_breaks) - 1L
  raw_counts <- lapply(levels, function(level) {
    dat <- region_bcrs[as.character(region_bcrs$.sirup_density_group) == level, , drop = FALSE]
    counts <- matrix(0, nrow = x_bins, ncol = y_bins)
    if (nrow(dat) == 0L) {
      return(counts)
    }
    xi <- findInterval(dat$x, x_breaks, all.inside = TRUE)
    yi <- findInterval(dat$.sirup_y, y_breaks, all.inside = TRUE)
    xi <- pmin(xi, x_bins)
    yi <- pmin(yi, y_bins)
    for (i in seq_len(nrow(dat))) {
      counts[xi[[i]], yi[[i]]] <- counts[xi[[i]], yi[[i]]] + dat$.sirup_density_weight[[i]]
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

.sirup_draw_histogram_bars_3d <- function(histogram, alpha) {
  bars <- histogram$bars
  if (is.null(bars) || nrow(bars) == 0L) {
    return(invisible(NULL))
  }
  color <- grDevices::adjustcolor(histogram$color, alpha.f = alpha)
  for (i in seq_len(nrow(bars))) {
    vertices <- .sirup_histogram_bar_vertices(bars[i, , drop = FALSE])
    .sirup_rgl("quads3d")(
      vertices$x,
      vertices$y,
      vertices$z,
      color = color,
      lit = FALSE
    )
  }
  invisible(NULL)
}

.sirup_histogram_bar_vertices <- function(bar) {
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

.sirup_image_texture_extent <- function(image_dim, xlim, ylim) {
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

.sirup_image_texture_file <- function(image, extent, alpha = 1, texture_dir = NULL) {
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
        "sirup_spatial_clone_density_texture",
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
      pattern = "sirup_spatial_clone_density_texture_",
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

.sirup_polygon_pie_summary <- function(region_bcrs,
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
  region_bcrs$.sirup_fill_by <- as.character(region_bcrs[[fill]])
  region_bcrs$.sirup_fill_by[is.na(region_bcrs$.sirup_fill_by) | !nzchar(region_bcrs$.sirup_fill_by)] <- NA_character_
  if (is.null(top_n)) {
    region_bcrs$.sirup_level <- region_bcrs$.sirup_fill_by
  } else {
    if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
      rlang::abort("`top_n` must be NULL or a positive number.")
    }
    top_levels <- unlist(lapply(split(region_bcrs, region_bcrs$polygon), function(dat) {
      counts <- sort(table(dat$.sirup_fill_by, useNA = "no"), decreasing = TRUE)
      names(counts)[seq_len(min(as.integer(top_n), length(counts)))]
    }), use.names = FALSE)
    top_levels <- unique(top_levels[!is.na(top_levels) & nzchar(top_levels)])
    region_bcrs$.sirup_level <- region_bcrs$.sirup_fill_by
    region_bcrs$.sirup_level[
      !is.na(region_bcrs$.sirup_level) &
        nzchar(region_bcrs$.sirup_level) &
        !region_bcrs$.sirup_level %in% top_levels
    ] <- other_label
  }
  pieces <- lapply(split(region_bcrs, region_bcrs$polygon), function(dat) {
    if (!isTRUE(plot_other)) {
      dat <- dat[dat$.sirup_level != other_label, , drop = FALSE]
    }
    if (nrow(dat) == 0L) {
      return(data.frame(level = character(), n = integer(), polygon = character()))
    }
    counts <- as.data.frame(table(dat$.sirup_level, useNA = "no"), stringsAsFactors = FALSE)
    names(counts) <- c("level", "n")
    counts <- counts[counts$n > 0L, , drop = FALSE]
    other <- counts[counts$level == other_label, , drop = FALSE]
    counts <- counts[counts$level != other_label, , drop = FALSE]
    counts <- counts[order(-counts$n, counts$level), , drop = FALSE]
    if (nrow(other) > 0L) {
      counts <- rbind(counts, other)
    }
    counts$.sirup_order <- seq_len(nrow(counts))
    polygon_value <- unique(dat$polygon)[[1]]
    counts$polygon <- polygon_value
    if ("polygon_label" %in% colnames(dat)) {
      labels <- unique(stats::na.omit(as.character(dat$polygon_label)))
      labels <- labels[nzchar(labels)]
      counts$polygon_label <- if (length(labels) > 0L) labels[[1L]] else as.character(polygon_value)
    } else {
      counts$polygon_label <- as.character(polygon_value)
    }
    counts
  })
  out <- .sirup_bind_rows(pieces)
  if (nrow(out) == 0L) {
    return(out)
  }
  out <- out[order(out$polygon_label, out$.sirup_order), , drop = FALSE]
  out <- do.call(rbind, lapply(split(out, out$polygon), function(dat) {
    total <- sum(dat$n)
    dat$proportion <- dat$n / total
    dat$ymax <- cumsum(dat$proportion)
    dat$ymin <- c(0, utils::head(dat$ymax, -1))
    dat$label_rank <- seq_len(nrow(dat))
    dat$label <- ifelse(dat$label_rank <= label_top_n, dat$level, "")
    dat
  }))
  out$.sirup_order <- NULL
  rownames(out) <- NULL
  out
}

# Decide the order of the polygon pie charts. `"receptors"` orders the polygons
# by their total receptor count (most first); ties break by the displayed map
# label. `"polygon"` keeps the numbered map order. Returns stable polygon IDs,
# not display labels, so equal labels from different samples remain distinct.
.sirup_pie_polygon_order <- function(pie_data, pie_order = c("receptors", "polygon")) {
  pie_order <- match.arg(pie_order)
  if (identical(pie_order, "polygon")) {
    return(unique(as.character(pie_data$polygon)))
  }
  totals <- tapply(pie_data$n, pie_data$polygon, sum)
  label_lookup <- stats::setNames(as.character(pie_data$polygon_label), as.character(pie_data$polygon))
  labels <- unname(label_lookup[names(totals)])
  names(totals)[order(-totals, suppressWarnings(as.integer(labels)), labels, names(totals))]
}

.sirup_polygon_centers <- function(polygon_plot_data, polygons) {
  centers <- lapply(split(polygon_plot_data, polygon_plot_data$polygon), function(dat) {
    dat <- dat[!duplicated(paste(dat$x, dat$y)), , drop = FALSE]
    data.frame(
      polygon = unique(dat$polygon)[[1]],
      x = mean(range(dat$x, na.rm = TRUE)),
      y = mean(range(dat$y, na.rm = TRUE)),
      .sirup_y = -mean(range(dat$y, na.rm = TRUE)),
      radius = max(diff(range(dat$x, na.rm = TRUE)), diff(range(dat$y, na.rm = TRUE))) / 2,
      stringsAsFactors = FALSE
    )
  })
  centers <- .sirup_bind_rows(centers)
  centers <- centers[centers$polygon %in% polygons, , drop = FALSE]
  centers$polygon_label <- vapply(centers$polygon, function(polygon_value) {
    .sirup_polygon_display_id(polygon_plot_data, polygon_value)
  }, character(1))
  centers
}

.sirup_dominance_only_polygon_pies <- function(pie_data, fill_values, label_size, pie_order) {
  ordered_polygons <- .sirup_pie_polygon_order(pie_data, pie_order)
  pie_groups <- split(pie_data, factor(as.character(pie_data$polygon), levels = ordered_polygons))
  pie_plots <- lapply(pie_groups, function(dat) {
    local_levels <- unique(as.character(dat$level))
    legend_ncol <- max(1L, min(4L, ceiling(length(local_levels) / 4)))
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
  patchwork::wrap_plots(pie_plots, ncol = max(1L, ceiling(sqrt(length(pie_plots)))))
}

# One combined stacked-bar panel: a bar per polygon (ordered by `pie_order`),
# stacked by `fill` level in dominance order. Uses the cumulative ymin/ymax from
# .sirup_polygon_pie_summary so each bar keeps its own per-polygon ordering, and a
# numeric x position so geom_rect can place one bar per polygon.
.sirup_stacked_bar_polygon_dominance <- function(pie_data, fill_values, label_size, pie_order) {
  ordered_polygons <- .sirup_pie_polygon_order(pie_data, pie_order)
  label_lookup <- stats::setNames(as.character(pie_data$polygon_label), as.character(pie_data$polygon))
  ordered_labels <- unname(label_lookup[ordered_polygons])
  pie_data$.sirup_polygon <- factor(as.character(pie_data$polygon), levels = ordered_polygons)
  pie_data$.sirup_xpos <- as.integer(pie_data$.sirup_polygon)
  local_levels <- unique(as.character(pie_data$level))
  legend_ncol <- max(1L, min(4L, ceiling(length(local_levels) / 4)))
  ggplot2::ggplot(pie_data) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$.sirup_xpos - 0.45, xmax = .data$.sirup_xpos + 0.45,
        ymin = .data$ymin, ymax = .data$ymax, fill = .data$level
      ),
      color = "white", linewidth = 0.15
    ) +
    ggplot2::geom_text(
      data = pie_data[pie_data$label != "", , drop = FALSE],
      ggplot2::aes(x = .data$.sirup_xpos, y = (.data$ymin + .data$ymax) / 2, label = .data$label),
      inherit.aes = FALSE, size = label_size
    ) +
    ggplot2::scale_fill_manual(values = fill_values[local_levels], limits = local_levels, drop = TRUE) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, ncol = legend_ncol, byrow = TRUE)) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(ordered_labels), labels = ordered_labels,
      expand = ggplot2::expansion(add = 0.6)
    ) +
    ggplot2::scale_y_continuous(
      labels = function(z) paste0(z * 100, "%"),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(x = NULL, y = "% of BCRs", fill = NULL) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 6),
      legend.key.size = grid::unit(0.14, "in"),
      legend.spacing.x = grid::unit(0.04, "in"),
      legend.spacing.y = grid::unit(0.01, "in"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      panel.grid = ggplot2::element_blank()
    )
}

.sirup_side_by_side_polygon_pies <- function(object,
                                             images,
                                             sample_id,
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
                                             zoom,
                                             combine_sample_plots = FALSE) {
  pie_order <- match.arg(pie_order)
  map_polygon_data <- polygon_plot_data
  polygon_labels <- stats::setNames(centers$polygon_label, centers$polygon)
  map_polygon_data$polygon_label <- unname(polygon_labels[as.character(map_polygon_data$polygon)])
  plot_object <- object
  plot_sample_id <- sample_id
  if (!is.null(sample_id) && .sirup_is_spatial_list(object)) {
    missing_samples <- setdiff(as.character(sample_id), names(object))
    if (length(missing_samples) > 0L) {
      rlang::abort(paste0(
        "`sample_id` not found in spatial samples: ",
        paste(missing_samples, collapse = ", ")
      ))
    }
    plot_object <- object[as.character(sample_id)]
    plot_sample_id <- NULL
  }
  map_plots <- plotSpatialPolygons(
    object = plot_object,
    polygon = map_polygon_data,
    images = images,
    sample_id = plot_sample_id,
    image_scale = image_scale,
    swap_coords = swap_coords,
    polygon_id = "polygon",
    polygon_alpha = polygon_alpha,
    image_alpha = image_alpha,
    label = TRUE,
    label_by = "polygon_label",
    label_size = label_size,
    zoom = zoom,
    combine = TRUE,
    combine_sample_plots = combine_sample_plots
  )

  ordered_polygons <- .sirup_pie_polygon_order(pie_data, pie_order)

  build_pie_plots <- function(dat) {
    local_order <- ordered_polygons[ordered_polygons %in% unique(as.character(dat$polygon))]
    pie_groups <- split(dat, factor(as.character(dat$polygon), levels = local_order))
    lapply(pie_groups, function(dat) {
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
  }

  build_pages <- function(map_plot, dat) {
    pie_plots <- build_pie_plots(dat)
    page_indices <- split(seq_along(pie_plots), ceiling(seq_along(pie_plots) / pies_per_page))
    page_names <- paste0("page_", seq_along(page_indices))
    names(page_indices) <- page_names

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
  }

  if (is.list(map_plots) && !inherits(map_plots, "ggplot")) {
    out <- list()
    for (image in names(map_plots)) {
      image_polygons <- .sirup_filter_polygon_plot_sample(map_polygon_data, image)
      image_ids <- unique(as.character(image_polygons$polygon))
      image_pies <- pie_data[as.character(pie_data$polygon) %in% image_ids, , drop = FALSE]
      if (nrow(image_pies) == 0L) {
        next
      }
      out[[image]] <- build_pages(map_plots[[image]], image_pies)
    }
    return(out)
  }

  build_pages(map_plots, pie_data)
}

.sirup_spatial_clipped_polygon_pies <- function(object,
                                                images,
                                                sample_id,
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
                                                zoom,
                                                combine_sample_plots = FALSE) {
  clipped <- .sirup_clipped_pie_polygons(polygon_plot_data, centers, pie_data)
  labels <- .sirup_spatial_pie_labels(centers, pie_data)
  legend_ncol <- max(1L, min(4L, ceiling(length(fill_values) / 8)))
  plotted_polygons <- unique(as.character(pie_data$polygon))
  zoom_data <- polygon_plot_data[as.character(polygon_plot_data$polygon) %in% plotted_polygons, , drop = FALSE]
  full_image_data <- data.frame(x = numeric(), .sirup_y = numeric())

  plots <- list()
  for (image in .sirup_spatial_images(object, images = images, sample_id = sample_id)) {
    spatial <- .sirup_spatial_data(object, image = image, image_scale = image_scale, swap_coords = swap_coords)
    image_polygons <- .sirup_filter_polygon_plot_sample(polygon_plot_data, image)
    image_polygon_ids <- unique(as.character(image_polygons$polygon))
    image_clipped <- clipped[as.character(clipped$polygon) %in% image_polygon_ids, , drop = FALSE]
    image_labels <- labels[as.character(labels$polygon) %in% image_polygon_ids, , drop = FALSE]
    image_zoom_data <- zoom_data[as.character(zoom_data$polygon) %in% image_polygon_ids, , drop = FALSE]
    plots[[image]] <- .sirup_spatial_base(
      plot_data = if (isTRUE(zoom)) image_zoom_data else full_image_data,
      image = spatial$image,
      image_alpha = image_alpha,
      zoom = zoom,
      view_data = NULL
    ) +
      ggplot2::geom_polygon(
        data = image_clipped,
        ggplot2::aes(x = .data$x, y = .data$.sirup_y, group = .data$piece, fill = .data$level),
        inherit.aes = FALSE,
        color = "white",
        linewidth = 0.15,
        alpha = pie_alpha
      ) +
      ggplot2::geom_path(
        data = image_polygons,
        ggplot2::aes(x = .data$x, y = .data$.sirup_y, group = .data$polygon),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.6
      ) +
      ggplot2::geom_text(
        data = image_labels,
        ggplot2::aes(x = .data$x, y = .data$.sirup_y, label = .data$label),
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

  .sirup_return_spatial_plots(plots, combine = TRUE, combine_sample_plots = combine_sample_plots)
}

.sirup_clipped_pie_polygons <- function(polygon_plot_data, centers, pie_data) {
  pieces <- list()
  idx <- 1L
  for (polygon_id in unique(pie_data$polygon)) {
    polygon_vertices <- polygon_plot_data[polygon_plot_data$polygon == polygon_id, , drop = FALSE]
    center <- centers[centers$polygon == polygon_id, , drop = FALSE]
    if (nrow(polygon_vertices) == 0L || nrow(center) == 0L) {
      next
    }
    boundary <- .sirup_repair_polygon_sfc(
      sf::st_sfc(sf::st_polygon(list(as.matrix(polygon_vertices[, c("x", "y")]))))
    )
    if (is.null(boundary)) {
      rlang::warn(paste0("Polygon ", polygon_id, " is invalid and could not be repaired for clipped pie plotting."))
      next
    }
    this_pie <- pie_data[pie_data$polygon == polygon_id, , drop = FALSE]
    this_pie <- this_pie[order(this_pie$ymin), , drop = FALSE]
    for (i in seq_len(nrow(this_pie))) {
      sector <- .sirup_pie_sector(center$x[[1]], center$y[[1]], center$radius[[1]] * 1.25, this_pie$ymin[[i]], this_pie$ymax[[i]])
      sector_sf <- .sirup_repair_polygon_sfc(
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
        suppressWarnings(.sirup_polygon_intersection(sector_sf, boundary)),
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
      intersection <- .sirup_repair_polygon_sfc(intersection)
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
      coords$.sirup_y <- -coords$y
      pieces[[idx]] <- coords[, c("polygon", "level", "piece", "x", "y", ".sirup_y"), drop = FALSE]
      idx <- idx + 1L
    }
  }
  .sirup_bind_rows(pieces)
}

.sirup_repair_polygon_sfc <- function(geometry) {
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

.sirup_polygon_intersection <- function(x, y) {
  tryCatch(
    suppressWarnings(sf::st_intersection(x, y)),
    error = function(e) {
      x <- .sirup_repair_polygon_sfc(suppressWarnings(sf::st_buffer(x, 0)))
      y <- .sirup_repair_polygon_sfc(suppressWarnings(sf::st_buffer(y, 0)))
      if (is.null(x) || is.null(y)) {
        stop(e)
      }
      suppressWarnings(sf::st_intersection(x, y))
    }
  )
}

.sirup_pie_sector <- function(center_x, center_y, radius, ymin, ymax, n = 80L) {
  start <- (ymin * 2 * pi) - (pi / 2)
  end <- (ymax * 2 * pi) - (pi / 2)
  angles <- seq(start, end, length.out = max(3L, ceiling(n * (ymax - ymin))))
  data.frame(
    x = c(center_x, center_x + radius * cos(angles), center_x),
    y = c(center_y, center_y + radius * sin(angles), center_y),
    stringsAsFactors = FALSE
  )
}

.sirup_spatial_pie_labels <- function(centers, pie_data) {
  labels <- pie_data[pie_data$label != "", , drop = FALSE]
  if (nrow(labels) == 0L) {
    return(data.frame(x = numeric(), y = numeric(), .sirup_y = numeric(), label = character()))
  }
  labels <- merge(labels, centers[, c("polygon", "x", "y", "radius"), drop = FALSE], by = "polygon", sort = FALSE)
  mid <- ((labels$ymin + labels$ymax) / 2) * 2 * pi - (pi / 2)
  labels$x <- labels$x + labels$radius * 0.45 * cos(mid)
  labels$y <- labels$y + labels$radius * 0.45 * sin(mid)
  labels$.sirup_y <- -labels$y
  labels[, c("x", "y", ".sirup_y", "label"), drop = FALSE]
}

.sirup_downsample_polygon_plot_coords <- function(coords, downsample_n = NULL, seed = 1L) {
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

.sirup_close_current_polygon <- function(points, polygon_num) {
  id_col <- if ("num" %in% colnames(points)) "num" else "polygon"
  current <- points[points[[id_col]] == polygon_num, , drop = FALSE]
  if (nrow(current) == 0L) {
    return(points)
  }
  if (nrow(current) < 3L) {
    rlang::warn("A polygon needs at least three distinct vertices before it can be closed.")
    return(points)
  }
  if (.sirup_same_xy(current[1, , drop = FALSE], current[nrow(current), , drop = FALSE])) {
    return(points)
  }
  rbind(points, current[1, , drop = FALSE])
}

.sirup_default_polygon_buffer_radius <- function(coords) {
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

.sirup_polygon_interior_to_sf_boundary <- function(points, buffer_radius) {
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
    boundary_points <- .sirup_clean_boundary_barcode_vertices(boundary_points)
    if (nrow(boundary_points) < 3L) {
      rlang::warn(paste0("Fewer than three non-degenerate boundary spatial IDs were found for polygon ", polygon_id, "."))
      return(NULL)
    }
    if (!.sirup_boundary_barcode_vertices_valid(boundary_points)) {
      boundary_points <- .sirup_repair_boundary_barcode_vertices(
        boundary_points = boundary_points,
        candidate_points = current
      )
    }
    if (!.sirup_boundary_barcode_vertices_valid(boundary_points)) {
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

.sirup_clean_boundary_barcode_vertices <- function(points, tol = 1e-8) {
  points <- points[is.finite(points$x) & is.finite(points$y), , drop = FALSE]
  if (nrow(points) < 3L) {
    return(points)
  }

  changed <- TRUE
  while (isTRUE(changed) && nrow(points) >= 3L) {
    old_n <- nrow(points)
    if (nrow(points) > 1L && .sirup_same_xy(points[1, , drop = FALSE], points[nrow(points), , drop = FALSE])) {
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
        remove[[i]] <- .sirup_degenerate_boundary_vertex(
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

.sirup_degenerate_boundary_vertex <- function(prev, current, next_point, tol = 1e-8) {
  if (.sirup_same_xy(prev, next_point)) {
    return(TRUE)
  }
  v1 <- c(prev$x[[1]] - current$x[[1]], prev$y[[1]] - current$y[[1]])
  v2 <- c(next_point$x[[1]] - current$x[[1]], next_point$y[[1]] - current$y[[1]])
  cross <- v1[[1]] * v2[[2]] - v1[[2]] * v2[[1]]
  dot <- sum(v1 * v2)
  isTRUE(abs(cross) <= tol && dot > tol)
}

.sirup_closed_boundary_matrix <- function(points) {
  xy <- as.matrix(points[, c("x", "y"), drop = FALSE])
  if (nrow(xy) == 0L) {
    return(xy)
  }
  if (!all(xy[1L, ] == xy[nrow(xy), ])) {
    xy <- rbind(xy, xy[1L, ])
  }
  xy
}

.sirup_boundary_barcode_vertices_valid <- function(points) {
  points <- .sirup_clean_boundary_barcode_vertices(points)
  if (nrow(unique(points[, c("x", "y"), drop = FALSE])) < 3L) {
    return(FALSE)
  }
  xy <- .sirup_closed_boundary_matrix(points)
  polygon <- tryCatch(
    sf::st_sfc(sf::st_polygon(list(xy))),
    error = function(e) NULL
  )
  if (is.null(polygon) || length(polygon) == 0L || any(sf::st_is_empty(polygon))) {
    return(FALSE)
  }
  isTRUE(all(sf::st_is_valid(polygon))) && isTRUE(as.numeric(sf::st_area(polygon)) > 0)
}

.sirup_repair_boundary_barcode_vertices <- function(boundary_points, candidate_points) {
  if (nrow(boundary_points) < 3L || nrow(candidate_points) < 3L) {
    return(boundary_points)
  }
  polygon <- tryCatch(
    sf::st_sfc(sf::st_polygon(list(.sirup_closed_boundary_matrix(boundary_points)))),
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
  .sirup_clean_boundary_barcode_vertices(candidate_points[nearest_idx, , drop = FALSE])
}

.sirup_close_formatted_polygons <- function(points) {
  closed <- lapply(unique(points$polygon), function(polygon_id) {
    current <- points[points$polygon == polygon_id, , drop = FALSE]
    if (nrow(current) < 3L) {
      rlang::warn(paste0("Polygon ", polygon_id, " has fewer than three points and will be excluded."))
      return(NULL)
    }
    if (.sirup_same_xy(current[1, , drop = FALSE], current[nrow(current), , drop = FALSE])) {
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
