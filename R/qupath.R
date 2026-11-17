#' Read QuPath annotations as SIRUP polygons
#'
#' `readQuPathPolygons` reads a GeoJSON file exported from QuPath. It returns a
#' polygon table. Pass the table to [addPolygons()].
#'
#' @details
#' In QuPath, use `File > Object data... > Export as GeoJSON`. Select the
#' annotations first, or export all objects. A script can also call
#' `exportObjectsToGeoJson(getAnnotationObjects(), path, "FEATURE_COLLECTION")`.
#'
#' QuPath writes coordinates in pixels of the full-resolution image. The origin
#' is the top-left corner, and y increases downward. SIRUP polygons must use the
#' same coordinates as `SeuratObject::GetTissueCoordinates()` at the
#' `image_scale` of the object. If the object uses the `hires` or `lowres`
#' image, set `scale` to the matching scale factor from
#' `scalefactors_json.json`. Use `scale = 1` for full-resolution coordinates.
#'
#' The function accepts a FeatureCollection, a single Feature, and the bare list
#' of Features that older QuPath versions write. It keeps features with
#' `objectType` equal to `"annotation"`, or with no `objectType`. It reads only `Polygon` and `MultiPolygon` shapes.
#' A `MultiPolygon` becomes one SIRUP polygon for each part. The IDs of the
#' parts end in `_p1`, `_p2`, and so on.
#'
#' A SIRUP polygon has one ring and cannot have holes. If an annotation has
#' holes, the function drops them and gives a warning. Spots in a hole then
#' count as inside the polygon.
#'
#' Polygon IDs come from the annotation name. If an annotation has no name, the
#' ID comes from the class name and a number, for example `ENDO_01`. An
#' annotation with neither gets the label `annotation`.
#'
#' @param path Path to the GeoJSON file.
#' @param sample_id Sample name to add to all polygons. A multi-sample SIRUP
#'   object needs it. Use `NULL` for a single-sample object.
#' @param class Class names to keep, for example `"GC"`. If `NULL`, the function
#'   keeps all classes.
#' @param keep_unclassified Keep annotations that have no class (default
#'   `TRUE`)? The function ignores this if you set `class`.
#' @param scale Number that multiplies all x and y values. See Details.
#'
#' @return A tibble with the columns `polygon_id`, `class`, `vertices`, and
#'   `sample` (only if you set `sample_id`). `vertices` is a list-column of
#'   tables with the columns `x` and `y`.
#' @seealso [addPolygons()], [plotSpatialPolygons()]
#'
#' @examples
#' geojson <- tempfile(fileext = ".geojson")
#' writeLines(paste0(
#'   '{"type":"FeatureCollection","features":[{"type":"Feature",',
#'   '"geometry":{"type":"Polygon","coordinates":',
#'   '[[[2.5,2.5],[3.5,2.5],[3.5,3.5],[2.5,3.5],[2.5,2.5]]]},',
#'   '"properties":{"objectType":"annotation",',
#'   '"classification":{"name":"GC","color":[255,0,0]}}}]}'
#' ), geojson)
#' polygons <- readQuPathPolygons(geojson)
#' polygons
#'
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- addPolygons(x, polygons)
#' getPolygons(x)$polygon_id
#' @export
readQuPathPolygons <- function(path,
                               sample_id = NULL,
                               class = NULL,
                               keep_unclassified = TRUE,
                               scale = 1) {
  if (!.sirup_is_non_empty_string(path)) {
    rlang::abort("`path` must be one non-empty character path.")
  }
  if (!file.exists(path)) {
    rlang::abort(paste0("`path` does not exist: ", path))
  }
  if (!is.null(sample_id) && !.sirup_is_non_empty_string(sample_id)) {
    rlang::abort("`sample_id` must be NULL or one non-empty character string.")
  }
  if (!is.null(class) && (!is.character(class) || anyNA(class))) {
    rlang::abort("`class` must be NULL or a character vector.")
  }
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    rlang::abort("`scale` must be one positive number.")
  }

  geojson <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  # Older QuPath versions write a bare list of features. Newer versions write a
  # FeatureCollection. A single Feature is also possible.
  features <- if (is.null(names(geojson))) {
    geojson
  } else if (identical(geojson[["type"]], "FeatureCollection")) {
    geojson[["features"]]
  } else if (identical(geojson[["type"]], "Feature")) {
    list(geojson)
  } else {
    rlang::abort("`path` must hold a GeoJSON FeatureCollection, a Feature, or a list of Features.")
  }

  # Return a property as one string. Return NA if it is missing or malformed.
  # `[[` matches names exactly. `$` would match partial names.
  property_string <- function(value) {
    if (is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)) value else NA_character_
  }

  # Keep annotations with a polygon shape.
  keep <- vapply(features, function(feature) {
    object_type <- feature[["properties"]][["objectType"]]
    shape <- property_string(feature[["geometry"]][["type"]])
    (is.null(object_type) || identical(object_type, "annotation")) &&
      shape %in% c("Polygon", "MultiPolygon")
  }, logical(1))
  features <- features[keep]

  feature_class <- vapply(features, function(feature) {
    property_string(feature[["properties"]][["classification"]][["name"]])
  }, character(1))
  feature_name <- vapply(features, function(feature) {
    property_string(feature[["properties"]][["name"]])
  }, character(1))

  keep <- if (!is.null(class)) {
    feature_class %in% class
  } else if (isTRUE(keep_unclassified)) {
    rep(TRUE, length(features))
  } else {
    !is.na(feature_class)
  }
  features <- features[keep]
  feature_class <- feature_class[keep]
  feature_name <- feature_name[keep]
  if (length(features) == 0L) {
    rlang::abort(paste0(
      "No QuPath annotations matched. The function reads only annotations with a ",
      "Polygon or MultiPolygon shape. Check `class`."
    ))
  }

  # Build one ID for each annotation. A name that is unique is used as it is.
  label <- ifelse(!is.na(feature_name), feature_name,
                  ifelse(!is.na(feature_class), feature_class, "annotation"))
  label_total <- as.integer(table(label)[label])
  label_number <- stats::ave(seq_along(label), label, FUN = seq_along)
  width <- max(2L, nchar(as.character(max(label_total))))
  annotation_id <- ifelse(
    !is.na(feature_name) & label_total == 1L,
    label,
    sprintf(paste0("%s_%0", width, "d"), label, label_number)
  )

  ring_to_vertices <- function(ring) {
    tibble::tibble(
      x = vapply(ring, function(point) as.numeric(point[[1L]]), numeric(1)) * scale,
      y = vapply(ring, function(point) as.numeric(point[[2L]]), numeric(1)) * scale
    )
  }

  rows <- list()
  holes <- character()
  for (i in seq_along(features)) {
    geometry <- features[[i]][["geometry"]]
    parts <- if (identical(geometry[["type"]], "Polygon")) {
      list(geometry[["coordinates"]])
    } else {
      geometry[["coordinates"]]
    }
    for (j in seq_along(parts)) {
      rings <- parts[[j]]
      if (length(rings) > 1L) {
        holes <- c(holes, annotation_id[[i]])
      }
      polygon_id <- if (length(parts) > 1L) {
        paste0(annotation_id[[i]], "_p", j)
      } else {
        annotation_id[[i]]
      }
      rows[[length(rows) + 1L]] <- tibble::tibble(
        polygon_id = polygon_id,
        class = feature_class[[i]],
        vertices = list(ring_to_vertices(rings[[1L]]))
      )
    }
  }
  out <- dplyr::bind_rows(rows)
  out$polygon_id <- make.unique(out$polygon_id, sep = "_")

  if (length(holes) > 0L) {
    holes <- unique(holes)
    rlang::warn(paste0(
      "SIRUP polygons cannot have holes. The function dropped the holes of ",
      length(holes), " annotation(s): ",
      paste(utils::head(holes, 10L), collapse = ", "),
      if (length(holes) > 10L) ", ..." else "",
      ". Spots in a hole count as inside the polygon."
    ))
  }
  if (!is.null(sample_id)) {
    out$sample <- sample_id
    out <- out[, c("polygon_id", "sample", "class", "vertices")]
  }
  out
}

#' Add QuPath annotations to a SIRUP object as polygons
#'
#' `addQuPathPolygons` reads a GeoJSON file that QuPath exported. It converts
#' the annotations to the coordinates of a spatial sample in a SIRUP object. It
#' adds a buffer to each polygon, so that every spot whose circle touches the
#' polygon becomes a member of it. Then it stores the polygons in the object.
#'
#' @details
#' QuPath writes coordinates in pixels of the full-resolution image. The
#' function finds the factor that converts them to the coordinates of the
#' object. It compares `SeuratObject::GetTissueCoordinates()` at the
#' `image_scale` of the object with the unscaled coordinates. This works for
#' `"hires"` and `"lowres"`. The image that you annotated in QuPath must be the
#' full-resolution image that the spots were registered to.
#'
#' `image_scale` and `swap_coords` come from the SIRUP object. If the object
#' has no `image_scale`, you must give it. The function swaps the polygon axes
#' in the same way as the spots when `swap_coords` is `TRUE`.
#'
#' A spot is a circle. By default, the function grows each polygon by the spot
#' radius, which it reads from the scale factors of the Seurat image. If the
#' sample has no Seurat image, the function uses the QuPath coordinates as they
#' are and sets the default buffer to 0. Then
#' [assignSpatialPolygons()] counts a spot as inside when its centre is inside
#' the grown polygon. That is the same as counting every spot whose circle
#' touches the original polygon. Set `buffer = 0` to count only spots with the
#' centre inside the QuPath polygon. For Visium HD, check that the stored spot
#' size suits your bins, or give `buffer`.
#'
#' A SIRUP polygon cannot have holes. The function drops holes and gives a
#' warning. See [readQuPathPolygons()] for the rules about shapes, classes, and
#' IDs.
#'
#' @param x A SIRUP object.
#' @param path Path to the QuPath GeoJSON file.
#' @param sample_id Name of the spatial sample that the annotations belong to.
#'   A SIRUP object with one spatial sample does not need it. A multi-sample
#'   object does.
#' @param class Class names to keep. If `NULL`, the function keeps all classes.
#' @param keep_unclassified Keep annotations that have no class (default
#'   `TRUE`)? The function ignores this if you set `class`.
#' @param image_scale Image scale of the polygon coordinates, for example
#'   `"hires"`. If `NULL`, the function uses the value stored in `x`. If `x` has
#'   no value, you must give one, and the function stores it in `x`. A value
#'   that differs from the stored one gives an error.
#' @param swap_coords Swap the two Seurat coordinate axes? If `NULL`, the
#'   function uses the value stored in `x`. The rules are the same as for
#'   `image_scale`.
#' @param image Name of the Seurat image to use. If `NULL`, the function uses
#'   the only image of the sample.
#' @param buffer Distance to grow each polygon, in the coordinates of the
#'   object. If `NULL` (default), the function uses the spot radius. Use `0` for
#'   no buffer.
#' @param replace Replace the stored polygons of this sample (default `FALSE`)?
#'   If `FALSE`, the function adds the new polygons to the stored ones, and the
#'   polygon IDs must not repeat.
#' @param include_boundary Count spots on a polygon edge as members of the
#'   polygon?
#' @param validate Validate the updated object?
#'
#' @return The updated SIRUP object. Spatial metadata and receptor records get
#'   the new `polygon_id` values.
#' @seealso [readQuPathPolygons()], [addPolygons()], [plotSpatialPolygons()]
#'
#' @examples
#' # Write a small GeoJSON file with one annotation. QuPath makes files like it.
#' geojson <- tempfile(fileext = ".geojson")
#' writeLines(paste0(
#'   '{"type":"FeatureCollection","features":[{"type":"Feature",',
#'   '"geometry":{"type":"Polygon","coordinates":',
#'   '[[[3.5,0.5],[5.5,0.5],[5.5,2.5],[3.5,2.5],[3.5,0.5]]]},',
#'   '"properties":{"objectType":"annotation","name":"GC1"}}]}'
#' ), geojson)
#'
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- addQuPathPolygons(x, geojson, replace = TRUE)
#' getPolygons(x)
#' table(getSpatial(x)$polygon_id, useNA = "ifany")
#' @export
addQuPathPolygons <- function(x,
                              path,
                              sample_id = NULL,
                              class = NULL,
                              keep_unclassified = TRUE,
                              image_scale = NULL,
                              swap_coords = NULL,
                              image = NULL,
                              buffer = NULL,
                              replace = FALSE,
                              include_boundary = TRUE,
                              validate = TRUE) {
  validateSIRUP(x)
  multi <- .sirup_is_spatial_list(x$spatial)
  if (multi) {
    .sirup_validate_sample_name(x, sample_id)
    sp <- x$spatial[[sample_id]]
  } else {
    if (!is.null(sample_id)) {
      rlang::abort("`sample_id` can only be used when this SIRUP object contains a spatial list.")
    }
    sp <- x$spatial
  }
  options <- .sirup_require_spatial_options(
    image_scale %||% x$metadata$image_scale,
    swap_coords %||% x$metadata$swap_coords,
    sample_id = sample_id
  )
  # Spot assignment reads these two options from the object. Store a value that
  # the user gave if the object has none. Refuse a value that conflicts.
  for (option in c("image_scale", "swap_coords")) {
    given <- if (option == "image_scale") image_scale else swap_coords
    stored_value <- x$metadata[[option]]
    if (is.null(stored_value)) {
      x$metadata[[option]] <- given
    } else if (!is.null(given)) {
      resolved <- .sirup_sample_option(
        if (option == "image_scale") .sirup_resolve_image_scale(stored_value) else .sirup_resolve_swap_coords(stored_value),
        sample_id = sample_id, argument = option
      )
      if (!identical(unname(as.vector(given)), unname(as.vector(resolved)))) {
        rlang::abort(paste0(
          "`", option, "` differs from the value stored in `x`. Update the stored ",
          "value first, or leave `", option, "` as `NULL`."
        ))
      }
    }
  }
  if (!is.null(buffer) && (!is.numeric(buffer) || length(buffer) != 1L || !is.finite(buffer) || buffer < 0)) {
    rlang::abort("`buffer` must be `NULL` or one number that is 0 or more.")
  }
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    rlang::abort("`replace` must be `TRUE` or `FALSE`.")
  }

  image_names <- SeuratObject::Images(sp)
  if (length(image_names) == 0L) {
    # No Seurat image: the coordinates are in the metadata. Use them as they are.
    px_factor <- 1
    if (is.null(buffer)) buffer <- 0
  } else {
    if (is.null(image)) {
      if (length(image_names) != 1L) {
        rlang::abort("The sample has more than one Seurat image. Give `image`.")
      }
      image <- image_names[[1L]]
    } else if (!image %in% image_names) {
      rlang::abort(paste0("Unknown Seurat image: ", image))
    }
    tissue <- sp[[image]]

    # Factor from full-resolution pixels (QuPath) to the coordinates of the object.
    pick_xy <- function(coords) {
      coords <- as.data.frame(coords, stringsAsFactors = FALSE)
      if (!"cell" %in% colnames(coords)) coords$cell <- rownames(coords)
      x_col <- intersect(c("x", "imagecol", "pxl_col_in_fullres", "col"), colnames(coords))[[1L]]
      list(cell = as.character(coords$cell), x = as.numeric(coords[[x_col]]))
    }
    scaled <- pick_xy(SeuratObject::GetTissueCoordinates(tissue, scale = options$image_scale))
    unscaled <- pick_xy(SeuratObject::GetTissueCoordinates(tissue, scale = NULL))
    ratios <- scaled$x / unscaled$x[match(scaled$cell, unscaled$cell)]
    ratios <- ratios[is.finite(ratios) & ratios > 0]
    if (length(ratios) == 0L) {
      rlang::abort(paste0(
        "Could not find the factor between full-resolution pixels and image scale \"",
        options$image_scale, "\"."
      ))
    }
    px_factor <- stats::median(ratios)

    # Buffer distance: the spot radius, in the coordinates of the object.
    if (is.null(buffer)) {
      spot <- tryCatch(unname(SeuratObject::ScaleFactors(tissue)[["spot"]]), error = function(e) NULL)
      if (is.null(spot)) spot <- tryCatch(tissue@scale.factors$spot, error = function(e) NULL)
      if (is.null(spot) || !is.numeric(spot) || !is.finite(spot)) {
        rlang::abort("Could not read the spot size from the Seurat image. Give `buffer`.")
      }
      buffer <- spot / 2 * px_factor
    }
  }

  qupath <- readQuPathPolygons(path, class = class, keep_unclassified = keep_unclassified, scale = px_factor)
  if (isTRUE(options$swap_coords)) {
    qupath$vertices <- lapply(qupath$vertices, function(v) tibble::tibble(x = v$y, y = v$x))
  }

  # Grow each polygon by the buffer. Keep the outer ring of the result.
  if (buffer > 0) {
    qupath$vertices <- lapply(seq_len(nrow(qupath)), function(i) {
      v <- qupath$vertices[[i]]
      ring <- sf::st_make_valid(sf::st_polygon(list(cbind(v$x, v$y))))
      grown <- sf::st_buffer(ring, dist = buffer, nQuadSegs = 8L)
      polygons <- sf::st_cast(sf::st_sfc(grown), "POLYGON")
      if (length(polygons) > 1L) {
        rlang::warn(paste0(
          "Polygon ", qupath$polygon_id[[i]], " has more than one part after the buffer. ",
          "The function kept the largest part."
        ))
        polygons <- polygons[which.max(as.numeric(sf::st_area(polygons)))]
      }
      xy <- sf::st_coordinates(polygons[[1L]])[, 1:2, drop = FALSE]
      tibble::tibble(x = xy[, 1L], y = xy[, 2L])
    })
  }

  # Build the new rows in the same layout that addPolygons() makes.
  flat <- do.call(rbind, lapply(seq_len(nrow(qupath)), function(i) {
    data.frame(
      polygon = qupath$polygon_id[[i]],
      x = qupath$vertices[[i]]$x,
      y = qupath$vertices[[i]]$y,
      stringsAsFactors = FALSE
    )
  }))
  if (multi) flat$sample <- sample_id
  new_rows <- .sirup_prepare_polygons(flat)
  new_rows$class <- qupath$class[match(vapply(new_rows$vertices, function(v) v$.sirup_polygon_label[[1L]], ""), qupath$polygon_id)]

  # Combine with the stored polygons.
  stored <- .sirup_polygons(x)
  if (nrow(stored) > 0L && isTRUE(replace)) {
    if (multi) {
      other <- setdiff(names(x$spatial), sample_id)
      stored <- .sirup_subset_polygons_component(stored, sample_id = other)
    } else {
      stored <- stored[0, , drop = FALSE]
    }
  }
  if (nrow(stored) > 0L) {
    clash <- new_rows$polygon_id %in% stored$polygon_id
    if (multi && any(clash)) {
      new_rows$polygon_id[clash] <- paste(sample_id, new_rows$polygon_id[clash], sep = "::")
      clash <- new_rows$polygon_id %in% stored$polygon_id
    }
    if (any(clash)) {
      rlang::abort(paste0(
        "These polygon IDs are already stored: ",
        paste(new_rows$polygon_id[clash], collapse = ", "),
        ". Use `replace = TRUE` or give the QuPath annotations different names."
      ))
    }
  }
  combined <- dplyr::bind_rows(stored, new_rows)
  # Give every polygon the same vertex columns. Stored polygons can have fewer.
  vertex_columns <- unique(unlist(lapply(combined$vertices, names)))
  combined$vertices <- lapply(seq_len(nrow(combined)), function(i) {
    vertices <- tibble::as_tibble(combined$vertices[[i]])
    for (column in setdiff(vertex_columns, names(vertices))) {
      vertices[[column]] <- if (column == ".sirup_polygon_label") {
        as.character(combined$polygon_id[[i]])
      } else if (column %in% c("sample", "spatial_sample") && column %in% names(combined)) {
        as.character(combined[[column]][[i]])
      } else {
        NA_character_
      }
    }
    vertices[, vertex_columns, drop = FALSE]
  })
  addPolygons(x, combined, include_boundary = include_boundary, validate = validate)
}
