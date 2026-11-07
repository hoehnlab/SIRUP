#' Create a SIRUP object
#'
#' `createSIRUP` joins a Seurat spatial object and an AIRR receptor table into
#' one SIRUP object.
#'
#' @details
#' `spatial` can be one Seurat object or a named list of Seurat objects. If you
#' give a list, the `sample` column of `receptors` must name the corresponding
#' list entry for each row. The function keeps rows with no matching entry as
#' nonspatial records. For example, these rows can hold bulk or single-cell
#' data.
#'
#' The function stores key receptor column names in `metadata$receptor_fields`
#' and converts `receptors` to a [tibble::tibble()]. The optional `analyses`
#' list keeps derived results.
#'
#' The examples in this package use a small synthetic SIRUP object. It has 25
#' spots, 10 genes, 5 clones, and 2 polygons. The spots and gene counts are
#' simulated. The receptor rows are real BCR sequences from a human thymus
#' section, assigned to synthetic spots. Load the object with
#' `readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))`.
#' The script `data-raw/make_synthetic_example.R` builds it.
#'
#' @param spatial Seurat object with spatial expression data, or a named list of
#'   Seurat objects (one per spatial sample).
#' @param receptors AIRR-style receptor table (a data frame or similar).
#' @param spatial_id Name of the `receptors` column that holds the spatial cell
#'   or spot IDs. If `NULL`, the function uses the first column found from this
#'   list: `spatial_id`, `spatial_id_8um`, `spatial_id_2um`, `spatial_id_16um`,
#'   `barcode_8um`, `barcode_2um`, `cell`, `cell_barcode`. It stores the name in
#'   the metadata and gives an error if none exists.
#' @param image Name of the `receptors` column that identifies the Seurat image
#'   (default `NULL`).
#' @param sample Name of the `receptors` column that holds the sample name. If
#'   `spatial` is a named list, this argument is required. If a Seurat object
#'   has no such metadata column, the function fills it with the list name.
#' @param clone Name of the `receptors` column that holds the clone ID (default
#'   `NULL`). The function gives an error if the column does not exist.
#' @param locus Name of the `receptors` column that holds the receptor locus
#'   (default `"locus"`). If you use the default and the column does not exist,
#'   the function stores `NULL`. If you give another name that does not exist,
#'   the function gives an error.
#' @param count Name of the `receptors` column that holds the observation count
#'   or weight (default `NULL`).
#' @param image_scale Default Seurat image scale (default `NULL`). Plotting and
#'   spatial statistics functions use it when their own `image_scale` is `NULL`.
#'   For multi-sample objects, give a named vector or list with one value for
#'   each sample, or one unnamed value that applies to all samples.
#' @param swap_coords Swap the two Seurat coordinate axes by default (default
#'   `NULL`)? Plotting and spatial statistics functions use this value when their
#'   own `swap_coords` is `NULL`. For multi-sample objects, give a named vector
#'   or list with one value for each sample, or one unnamed value that applies to
#'   all samples.
#' @param metadata Named list of dataset-level metadata.
#' @param polygons Polygon table (default `NULL`). The table has one row per
#'   polygon, with a `polygon_id` column and a list column `vertices`. Each
#'   `vertices` table needs columns `x` and `y`. It can also have `spatial_id`
#'   and `image`. The function adds a `polygon_id` column and one logical
#'   `polygon_<ID>` column per polygon to both the Seurat metadata and the
#'   receptor table.
#' @param polygons_path Path to a flat CSV file of polygon vertices (default
#'   `NULL`). The file needs `x` and `y` columns and a `polygon_id` or `polygon`
#'   column. The function annotates these polygons in the same way as
#'   `polygons`.
#' @param include_boundary Count spots on a polygon edge as members of the
#'   polygon (default `TRUE`)?
#' @param analyses Named list of derived analysis results.
#' @param validate Validate the new object (default `TRUE`)?
#'
#' @return A SIRUP object (S3 class `sirup`).
#' @seealso [validateSIRUP()], [mergeSIRUP()], [saveSIRUP()]
#'
#' @examples
#' # Use the small synthetic object that comes with the package.
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x
#'
#' # Build a SIRUP object from a Seurat object and a receptor table.
#' sirup <- createSIRUP(
#'   spatial = getSpatial(x),
#'   receptors = getReceptors(x),
#'   spatial_id = "spatial_id",
#'   clone = "clone_id",
#'   image_scale = "hires",
#'   swap_coords = FALSE
#' )
#' sirup
#' @export
createSIRUP <- function(spatial,
                        receptors,
                        spatial_id = NULL,
                        image = NULL,
                        sample = NULL,
                        clone = NULL,
                        locus = "locus",
                        count = NULL,
                        image_scale = NULL,
                        swap_coords = NULL,
                        metadata = list(),
                        polygons = NULL,
                        polygons_path = NULL,
                        include_boundary = TRUE,
                        analyses = list(),
                        validate = TRUE) {
  clone_missing <- missing(clone)
  locus_missing <- missing(locus)
  receptors <- tibble::as_tibble(receptors)
  spatial_id <- .guess_spatial_id_col(receptors, spatial_id)

  clone <- .sirup_resolve_optional_col(
    receptors,
    clone,
    "clone",
    allow_missing = clone_missing
  )
  locus <- .sirup_resolve_optional_col(
    receptors,
    locus,
    "locus",
    allow_missing = locus_missing
  )
  image <- .sirup_resolve_optional_col(receptors, image, "image")
  sample <- .sirup_resolve_optional_col(receptors, sample, "sample")
  if (.sirup_is_spatial_list(spatial) && is.null(sample)) {
    rlang::abort("`sample` is required when `spatial` is a named list of Seurat objects.")
  }
  count <- .sirup_resolve_optional_col(receptors, count, "count")
  image_scale <- .sirup_resolve_image_scale(image_scale)
  swap_coords <- .sirup_resolve_swap_coords(swap_coords)
  spatial <- .sirup_populate_spatial_sample_metadata(spatial, sample)

  if (!is.list(metadata)) {
    rlang::abort("`metadata` must be a named list.")
  }
  if (length(metadata) > 0L && (is.null(names(metadata)) || any(!nzchar(names(metadata))))) {
    rlang::abort("`metadata` must be a named list.")
  }
  reserved_metadata <- c("receptor_fields", "image_scale", "swap_coords", "sample_names")
  reserved_present <- intersect(reserved_metadata, names(metadata))
  if (length(reserved_present) > 0L) {
    rlang::abort(paste0(
      "`metadata$", reserved_present[[1L]],
      "` is reserved; use the corresponding SIRUP object arguments or spatial sample names."
    ))
  }
  if (!is.null(polygons) && !is.null(polygons_path)) {
    rlang::abort("Supply only one of `polygons` and `polygons_path`.")
  }
  if (!is.null(polygons_path)) {
    polygons <- .sirup_read_polygons_paths(polygons_path)
    metadata$polygons_path <- as.character(unlist(polygons_path, recursive = FALSE, use.names = TRUE))
  }
  .sirup_validate_multisample_polygon_input(spatial, polygons)
  polygons <- .sirup_prepare_polygons(polygons)
  .sirup_validate_analyses(analyses)
  metadata$receptor_fields <- list(
    spatial_id = spatial_id,
    image = image,
    sample = sample,
    clone = clone,
    locus = locus,
    count = count
  )
  metadata$sample_names <- .sirup_spatial_names(spatial)
  metadata$image_scale <- image_scale
  metadata$swap_coords <- swap_coords
  metadata$polygons_include_boundary <- isTRUE(include_boundary)

  out <- structure(
    list(
      spatial = spatial,
      receptors = receptors,
      metadata = metadata,
      polygons = polygons,
      analyses = analyses
    ),
    class = c("sirup", "list")
  )

  out <- .sirup_annotate_stored_polygons(out)

  if (isTRUE(validate)) {
    validateSIRUP(out)
  }
  out
}

#' Validate a SIRUP object
#'
#' `validateSIRUP` checks that the parts of a SIRUP object agree with each other.
#'
#' @details
#' The function gives an error if one of these checks fails:
#' \itemize{
#'   \item The object does not have class `sirup`.
#'   \item No receptor spatial ID matches a cell or spot name in the Seurat
#'     object.
#'   \item A column named in `metadata$receptor_fields` does not exist in the
#'     receptor table.
#'   \item The `sample` field is set, but the column is missing from the receptor
#'     table or from the Seurat metadata.
#' }
#'
#' @param x Object to validate.
#' @param verbose Print a short summary of the result (default `FALSE`)?
#'
#' @return `x`, invisibly, if the object is valid. If not, the function gives an
#'   error.
#' @seealso [createSIRUP()], [isSIRUP()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' validateSIRUP(x)
#' validateSIRUP(x, verbose = TRUE)
#' @export
validateSIRUP <- function(x, verbose = FALSE) {
  if (!isSIRUP(x)) {
    rlang::abort("`x` must inherit from class 'sirup'.")
  }
  legacy_components <- c("spatial", "receptors", "metadata")
  legacy_with_polygons <- c(legacy_components, "polygons")
  legacy_with_analyses <- c(legacy_components, "analyses")
  required_components <- c(legacy_components, "polygons", "analyses")
  if (!identical(names(x), required_components) &&
      !identical(names(x), legacy_with_polygons) &&
      !identical(names(x), legacy_with_analyses) &&
      !identical(names(x), legacy_components)) {
    rlang::abort(paste0(
      "A SIRUP object must contain, in order: ",
      paste(required_components, collapse = ", "),
      "."
    ))
  }
  .sirup_validate_spatial_component(x$spatial)
  if (.sirup_is_spatial_list(x$spatial) && is.null(.sirup_receptor_fields(x)$sample)) {
    rlang::abort("SIRUP objects with a spatial list require a non-NULL sample receptor field.")
  }
  if (!tibble::is_tibble(x$receptors)) {
    rlang::abort("`x$receptors` must be a tibble.")
  }
  if (!is.list(x$metadata)) {
    rlang::abort("`x$metadata` must be a list.")
  }
  if (!is.null(x$polygons)) {
    .sirup_validate_polygons(x$polygons)
  }
  if (length(x$metadata) > 0L && (is.null(names(x$metadata)) || any(!nzchar(names(x$metadata))))) {
    rlang::abort("`x$metadata` must be a named list.")
  }
  .sirup_validate_analyses(.sirup_analyses(x))
  .sirup_validate_spatial_metadata(x$metadata)
  .sirup_validate_metadata_sample_names(x)
  receptor_fields <- .sirup_receptor_fields(x)
  field_names <- c("spatial_id", "image", "sample", "clone", "locus", "count")
  if (!identical(names(receptor_fields), field_names)) {
    rlang::abort(paste0(
      "SIRUP receptor fields must contain, in order: ",
      paste(field_names, collapse = ", "),
      "."
    ))
  }

  spatial_id_col <- receptor_fields$spatial_id
  if (!is.character(spatial_id_col) || length(spatial_id_col) != 1L || is.na(spatial_id_col) || !nzchar(spatial_id_col)) {
    rlang::abort("The SIRUP spatial_id field must name one receptor column.")
  }
  .check_required_cols(x$receptors, spatial_id_col, "x$receptors")
  for (field in setdiff(field_names, "spatial_id")) {
    column <- receptor_fields[[field]]
    if (!is.null(column)) {
      if (!is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column)) {
        rlang::abort(paste0("The SIRUP ", field, " field must be `NULL` or one column name."))
      }
      .check_required_cols(x$receptors, column, "x$receptors")
    }
  }

  if (!is.null(receptor_fields$sample)) {
    sample_field <- receptor_fields$sample
    if (inherits(x$spatial, "Seurat")) {
      spatial_metadata <- x$spatial[[]]
      if (sample_field %in% colnames(spatial_metadata)) {
        .check_required_cols(spatial_metadata, sample_field, "x$spatial metadata")
      }
    } else if (.sirup_is_spatial_list(x$spatial)) {
      for (sample_name in names(x$spatial)) {
        spatial_metadata <- x$spatial[[sample_name]][[]]
        if (!sample_field %in% colnames(spatial_metadata)) {
          next
        }
        values <- unique(stats::na.omit(as.character(spatial_metadata[[sample_field]])))
        values <- values[nzchar(values)]
        if (length(values) > 0L && any(values != sample_name)) {
          rlang::abort(paste0(
            "Spatial metadata sample field `", sample_field, "` in spatial sample `",
            sample_name, "` contains value(s) that do not match the spatial list name: ",
            paste(values[values != sample_name], collapse = ", "), "."
          ))
        }
      }
    }
  }

  diagnostics <- .sirup_join_diagnostics(x)
  if (diagnostics$n_receptor_ids == 0L) {
    rlang::abort("`x$receptors` contains no non-missing spatial IDs.")
  }
  if (diagnostics$n_matched_ids == 0L) {
    rlang::abort("No receptor IDs match the spatial Seurat object.")
  }
  if (isTRUE(verbose)) {
    print(x)
  }
  invisible(x)
}

#' Test whether an object is a SIRUP object
#'
#' `isSIRUP` checks if an object has class `sirup`. It does not validate the
#' object. Use [validateSIRUP()] for that.
#'
#' @param x Object to test.
#'
#' @return `TRUE` if `x` has class `sirup`. Otherwise, `FALSE`.
#' @seealso [validateSIRUP()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' isSIRUP(x)
#' isSIRUP(getReceptors(x))
#' @export
isSIRUP <- function(x) {
  inherits(x, "sirup")
}

#' Merge SIRUP objects into one multi-sample object
#'
#' `mergeSIRUP` combines several SIRUP objects into one multi-sample SIRUP
#' object.
#'
#' @details
#' Each input object becomes one spatial sample, named by `samples`. The
#' function does these steps:
#' \itemize{
#'   \item It sets `sample_field` in each receptor table, then binds the rows.
#'   \item It stores the Seurat objects in a named list.
#'   \item It adds the sample name to each stored polygon (if any) and keeps
#'     the polygons.
#'   \item It stores each analysis as a list with one entry per sample, if
#'     `keep_analyses` is `TRUE`.
#' }
#'
#' @param ... Two or more SIRUP objects to merge.
#' @param samples Character vector of sample names. Give one unique, non-empty
#'   name for each SIRUP object.
#' @param sample_field Name of the receptor-table column that holds the sample
#'   name (default `"sample_id"`). The function creates the column or checks it.
#' @param keep_analyses Keep the stored analyses (default `TRUE`)? Each analysis
#'   becomes a named list with one entry for each name in `samples`.
#' @param validate Validate the merged object (default `TRUE`)?
#'
#' @return A multi-sample SIRUP object.
#' @seealso [createSIRUP()], [subsetSIRUP()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' merged <- mergeSIRUP(x, x, samples = c("a", "b"))
#' spatialSamples(merged)
#' @export
mergeSIRUP <- function(...,
                       samples,
                       sample_field = "sample_id",
                       keep_analyses = TRUE,
                       validate = TRUE) {
  objects <- list(...)
  if (length(objects) == 1L && is.list(objects[[1L]]) && !isSIRUP(objects[[1L]])) {
    objects <- objects[[1L]]
  }
  if (length(objects) < 2L) {
    rlang::abort("`mergeSIRUP()` requires at least two SIRUP objects.")
  }
  if (!is.character(samples) || length(samples) != length(objects) ||
      any(is.na(samples) | !nzchar(samples)) || anyDuplicated(samples)) {
    rlang::abort("`samples` must contain one unique, non-empty sample name per SIRUP object.")
  }
  if (!.sirup_is_non_empty_string(sample_field)) {
    rlang::abort("`sample_field` must be one non-empty receptor column name.")
  }
  objects <- lapply(objects, function(object) {
    validateSIRUP(object)
    object
  })
  fields <- .sirup_check_merge_receptor_fields(objects)

  spatial <- list()
  receptor_tables <- vector("list", length(objects))
  polygon_tables <- vector("list", length(objects))
  analyses <- list()
  for (i in seq_along(objects)) {
    sample_name <- samples[[i]]
    object <- objects[[i]]
    object_spatial <- getSpatial(object)
    if (.sirup_is_spatial_list(object_spatial)) {
      if (length(object_spatial) != 1L) {
        rlang::abort("Each input SIRUP object must contain exactly one spatial sample. Use `filterSample()` before merging multi-sample objects.")
      }
      object_spatial <- object_spatial[[1L]]
    }
    spatial[[sample_name]] <- object_spatial
    receptor_tables[[i]] <- .sirup_merge_receptors_for_sample(
      object,
      sample = sample_name,
      sample_field = sample_field
    )
    polygon_tables[[i]] <- .sirup_polygons_for_merge(object, sample = sample_name)
    if (isTRUE(keep_analyses)) {
      analyses <- .sirup_merge_analyses(analyses, .sirup_analyses(object), sample = sample_name)
    }
  }

  receptors <- Reduce(.sirup_bind_receptor_tables, receptor_tables)
  metadata <- .sirup_merge_metadata(objects, samples = samples)
  metadata$receptor_fields <- fields
  metadata$receptor_fields$sample <- sample_field
  metadata$sample_names <- samples
  metadata$image_scale <- .sirup_merge_sample_options(objects, samples, "image_scale")
  metadata$swap_coords <- .sirup_merge_sample_options(objects, samples, "swap_coords")
  metadata$polygons_include_boundary <- all(vapply(objects, function(object) {
    isTRUE(object$metadata$polygons_include_boundary %||% TRUE)
  }, logical(1)))
  polygons <- .sirup_bind_polygons_for_merge(polygon_tables)

  out <- structure(
    list(
      spatial = spatial,
      receptors = tibble::as_tibble(receptors),
      metadata = metadata,
      polygons = polygons,
      analyses = analyses
    ),
    class = c("sirup", "list")
  )
  out$spatial <- .sirup_populate_spatial_sample_metadata(out$spatial, sample_field)
  out <- .sirup_annotate_stored_polygons(out)
  if (isTRUE(validate)) {
    validateSIRUP(out)
  }
  out
}

#' Extract data from a SIRUP object
#'
#' `getSpatial` returns the Seurat object(s). `spatialSamples` returns the
#' spatial sample names. `getReceptors` returns the receptor table.
#'
#' @param x A SIRUP object.
#' @param sample Name of one spatial sample. Use it only if the object holds a
#'   named list of Seurat objects.
#'
#' @return `getSpatial` returns the Seurat object. For a multi-sample object, it
#'   returns the named list of Seurat objects, or one Seurat object if you give
#'   `sample`. `spatialSamples` returns a character vector of sample names.
#'   `getReceptors` returns the receptor tibble. If you give `sample`, it
#'   returns only the rows for that sample.
#' @seealso [getPolygons()], [getAnalysis()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' getSpatial(x)
#' head(getReceptors(x))
#' spatialSamples(x)
#' @name sirup-accessors
NULL

# The roxygen NULL above intentionally creates one shared help topic for the
# accessor functions below.

#' @rdname sirup-accessors
#' @export
getSpatial <- function(x, sample = NULL) {
  validateSIRUP(x)
  if (is.null(sample)) {
    return(x$spatial)
  }
  if (!.sirup_is_spatial_list(x$spatial)) {
    rlang::abort("`sample` can only be used when this SIRUP object contains a spatial list.")
  }
  if (!.sirup_is_non_empty_string(sample)) {
    rlang::abort("`sample` must be one non-empty spatial sample name.")
  }
  if (!sample %in% names(x$spatial)) {
    rlang::abort(paste0("Unknown spatial sample: ", sample))
  }
  x$spatial[[sample]]
}

#' @rdname sirup-accessors
#' @export
spatialSamples <- function(x) {
  validateSIRUP(x)
  .sirup_spatial_names(x$spatial)
}

#' @rdname sirup-accessors
#' @export
getReceptors <- function(x, sample = NULL) {
  validateSIRUP(x)
  if (is.null(sample)) {
    return(x$receptors)
  }
  .sirup_receptors_for_sample(x, sample)
}

#' Subset a SIRUP object
#'
#' `subsetSIRUP` keeps selected samples, spatial IDs, and receptor rows of a
#' SIRUP object. `filterSample` keeps only the samples you name.
#'
#' @details
#' The function subsets the Seurat data, the receptor rows, the stored polygons,
#' and the analyses. It applies the stored polygons again to the subset, so the
#' polygon columns stay consistent.
#'
#' @param x A SIRUP object.
#' @param sample Name or names of the spatial samples to keep. Use it only if
#'   the object holds a named list of Seurat objects.
#' @param spatial_id Spatial IDs (barcodes) to keep. In a multi-sample object,
#'   the function matches the IDs in each selected sample.
#' @param receptor_rows Logical vector with length `nrow(getReceptors(x))`, or a
#'   function. The function takes the receptor tibble and returns such a vector.
#'   The function keeps the rows where the vector is `TRUE`.
#' @param keep_nonspatial_receptors Keep receptor rows that have no spatial ID
#'   after the other filters (default `FALSE`)? The function uses this argument
#'   only if you give `spatial_id`.
#' @param keep_analyses Keep the stored analyses (default `FALSE`)? Most
#'   analyses depend on the original spatial and receptor data.
#' @param drop_empty_spatial Drop samples that have no cells left after the
#'   `spatial_id` filter (default `TRUE`)? The function uses this argument only
#'   for multi-sample objects.
#' @param validate Validate the subset object (default `TRUE`)?
#'
#' @return A SIRUP object that contains the subset.
#' @seealso [mergeSIRUP()], [createSIRUP()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Keep three spots.
#' subsetSIRUP(x, spatial_id = c("spot01", "spot02", "spot03"))
#'
#' # Keep the receptor rows of clone 2.
#' subsetSIRUP(x, receptor_rows = function(receptors) receptors$clone_id == "2")
#'
#' # Keep one sample of a multi-sample object.
#' merged <- mergeSIRUP(x, x, samples = c("a", "b"))
#' filterSample(merged, "a")
#' @export
subsetSIRUP <- function(x,
                        sample = NULL,
                        spatial_id = NULL,
                        receptor_rows = NULL,
                        keep_nonspatial_receptors = FALSE,
                        keep_analyses = FALSE,
                        drop_empty_spatial = TRUE,
                        validate = TRUE) {
  validateSIRUP(x)
  receptor_fields <- .sirup_receptor_fields(x)
  spatial_id_col <- receptor_fields$spatial_id
  sample_col <- receptor_fields$sample

  if (!is.null(sample)) {
    sample <- unique(as.character(sample))
    sample <- sample[!is.na(sample) & nzchar(sample)]
    if (length(sample) == 0L) {
      rlang::abort("`sample` must contain at least one non-empty sample name.")
    }
    if (!.sirup_is_spatial_list(x$spatial)) {
      rlang::abort("`sample` can only be used when this SIRUP object contains a spatial list.")
    }
    missing_samples <- setdiff(sample, names(x$spatial))
    if (length(missing_samples) > 0L) {
      rlang::abort(paste0("Unknown spatial sample(s): ", paste(missing_samples, collapse = ", ")))
    }
  }
  if (!is.null(spatial_id)) {
    spatial_id <- unique(as.character(spatial_id))
    spatial_id <- spatial_id[!is.na(spatial_id) & nzchar(spatial_id)]
    if (length(spatial_id) == 0L) {
      rlang::abort("`spatial_id` must contain at least one non-missing spatial ID.")
    }
  }

  out <- x
  out$spatial <- .sirup_subset_spatial_component(
    x$spatial,
    sample = sample,
    spatial_id = spatial_id,
    drop_empty_spatial = drop_empty_spatial
  )
  out$metadata$sample_names <- .sirup_spatial_names(out$spatial)

  keep <- rep(TRUE, nrow(x$receptors))
  if (!is.null(sample)) {
    if (is.null(sample_col) || !sample_col %in% colnames(x$receptors)) {
      rlang::abort("Cannot subset receptors by `sample`; this SIRUP object has no receptor sample field.")
    }
    keep <- keep & !is.na(x$receptors[[sample_col]]) & as.character(x$receptors[[sample_col]]) %in% sample
  }
  if (!is.null(spatial_id)) {
    receptor_ids <- as.character(x$receptors[[spatial_id_col]])
    if (.sirup_is_spatial_list(out$spatial) && !is.null(sample_col) && sample_col %in% colnames(x$receptors)) {
      allowed_spatial <- .sirup_spatial_cells(out$spatial)
      allowed_keys <- paste(allowed_spatial$.sirup_sample, allowed_spatial$spatial_id, sep = "\r")
      receptor_keys <- paste(as.character(x$receptors[[sample_col]]), receptor_ids, sep = "\r")
      spatial_keep <- !is.na(receptor_ids) & receptor_keys %in% allowed_keys
    } else {
      spatial_keep <- !is.na(receptor_ids) & receptor_ids %in% spatial_id
    }
    if (isTRUE(keep_nonspatial_receptors)) {
      spatial_keep <- spatial_keep | is.na(receptor_ids) | !nzchar(receptor_ids)
    }
    keep <- keep & spatial_keep
  }
  if (!is.null(receptor_rows)) {
    row_keep <- if (is.function(receptor_rows)) {
      receptor_rows(x$receptors)
    } else {
      receptor_rows
    }
    if (!is.logical(row_keep) || length(row_keep) != nrow(x$receptors) || anyNA(row_keep)) {
      rlang::abort("`receptor_rows` must be a non-missing logical vector with one value per receptor row, or a function returning one.")
    }
    keep <- keep & row_keep
  }
  out$receptors <- tibble::as_tibble(x$receptors[keep, , drop = FALSE])
  out$polygons <- .sirup_subset_polygons_component(
    .sirup_polygons(x),
    sample = sample,
    relabel_single_sample = !is.null(sample) && length(sample) == 1L
  )
  if (!isTRUE(keep_analyses)) {
    out$analyses <- list()
  }
  out <- .sirup_annotate_stored_polygons(out)
  if (isTRUE(validate)) {
    validateSIRUP(out)
  }
  out
}

#' @rdname subsetSIRUP
#' @export
filterSample <- function(x,
                         sample,
                         keep_analyses = FALSE,
                         validate = TRUE) {
  out <- subsetSIRUP(
    x,
    sample = sample,
    keep_analyses = keep_analyses,
    validate = FALSE
  )
  if (.sirup_is_spatial_list(out$spatial) && length(out$spatial) == 1L) {
    out$spatial <- out$spatial[[1L]]
    out$metadata$sample_names <- .sirup_spatial_names(out$spatial)
  }
  if (isTRUE(validate)) {
    validateSIRUP(out)
  }
  out
}

#' Add or retrieve polygon geometry from a SIRUP object
#'
#' `addPolygons` stores a polygon table in a SIRUP object. It replaces any stored
#' polygons. `getPolygons` returns the stored table.
#'
#' @details
#' The input to `addPolygons` can be a nested tibble with `polygon_id` and
#' `vertices` columns. It can also be a flat vertex table, such as the CSV file
#' from [selectSpatialPolygons()].
#'
#' The function adds a `polygon_id` column to both the Seurat metadata and the
#' receptor table. If a barcode lies in overlapping or nested polygons, it gets
#' the ID of the smallest polygon by area. The function also adds one logical
#' `polygon_<ID>` column for each polygon. These columns show membership in
#' every polygon as booleans.
#'
#' @param x A SIRUP object.
#' @param polygons Polygon table.
#' @param polygons_path Path to a flat CSV file of polygon vertices. For a
#'   multi-sample object, give a named vector or list of paths, with one path
#'   for each sample.
#' @param sample Sample name or names to add to the polygons from
#'   `polygons_path`. Use it only if the paths have no names.
#' @param include_boundary Count spots on a polygon edge as members of the
#'   polygon?
#' @param validate Validate the updated object?
#' @return `getPolygons` returns a polygon tibble. `addPolygons` returns the
#'   updated SIRUP object.
#' @seealso [selectSpatialPolygons()], [createSIRUP()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' getPolygons(x)
#'
#' # Replace the stored polygons with one polygon in the middle of the data.
#' center <- tibble::tibble(
#'   polygon_id = "center",
#'   vertices = list(data.frame(x = c(2.5, 3.5, 3.5, 2.5), y = c(2.5, 2.5, 3.5, 3.5)))
#' )
#' x2 <- addPolygons(x, center)
#' getPolygons(x2)$polygon_id
#' @name sirup-polygons
NULL

#' @rdname sirup-polygons
#' @export
getPolygons <- function(x, sample = NULL) {
  validateSIRUP(x)
  polygons <- .sirup_polygons(x)
  if (is.null(sample) || nrow(polygons) == 0L) {
    return(polygons)
  }
  if (!.sirup_is_spatial_list(x$spatial)) {
    rlang::abort("`sample` can only be used when this SIRUP object contains a spatial list.")
  }
  sample <- unique(as.character(sample))
  missing_samples <- setdiff(sample, names(x$spatial))
  if (length(missing_samples) > 0L) {
    rlang::abort(paste0("Unknown spatial sample(s): ", paste(missing_samples, collapse = ", ")))
  }
  .sirup_subset_polygons_component(
    polygons,
    sample = sample,
    relabel_single_sample = length(sample) == 1L
  )
}

#' @rdname sirup-polygons
#' @export
addPolygons <- function(x,
                        polygons = NULL,
                        polygons_path = NULL,
                        sample = NULL,
                        include_boundary = TRUE,
                        validate = TRUE) {
  validateSIRUP(x)
  if (!is.null(polygons) && !is.null(polygons_path)) {
    rlang::abort("Supply only one of `polygons` and `polygons_path`.")
  }
  if (!is.null(polygons_path)) {
    polygons <- .sirup_read_polygons_paths(polygons_path, sample = sample)
    x$metadata$polygons_path <- as.character(unlist(polygons_path, recursive = FALSE, use.names = TRUE))
  }
  .sirup_validate_multisample_polygons(x, polygons)
  x$polygons <- .sirup_prepare_polygons(polygons)
  x$metadata$polygons_include_boundary <- isTRUE(include_boundary)
  x <- .sirup_annotate_stored_polygons(x)
  if (isTRUE(validate)) validateSIRUP(x)
  x
}

#' Store and retrieve analyses in a SIRUP object
#'
#' `getAnalysis` returns one stored result. `setAnalysis` adds or replaces a
#' result. `listAnalyses` returns the names of the stored results.
#' `removeAnalysis` deletes a result.
#'
#' @details
#' The `analyses` list of a SIRUP object holds any R object. Use it to keep
#' derived results with the data that produced them. Examples are Dowser clones
#' and trees, spatial statistics, and polygon summaries.
#'
#' @param x A SIRUP object.
#' @param name Name of the analysis. Give one non-empty character string.
#' @param value R object to store. It can be `NULL`.
#' @param sample Name of a spatial sample. Use it to store or get a result for
#'   one sample.
#' @param validate Validate the updated SIRUP object?
#'
#' @return `getAnalysis` returns the stored object. `listAnalyses` returns a
#'   character vector. `setAnalysis` and `removeAnalysis` return the updated
#'   SIRUP object.
#' @seealso [updateAnalysis()], [mergeSIRUP()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- setAnalysis(x, "n_clones", 5)
#' listAnalyses(x)
#' getAnalysis(x, "n_clones")
#' x <- removeAnalysis(x, "n_clones")
#' listAnalyses(x)
#' @name sirup-analyses
NULL

#' @rdname sirup-analyses
#' @export
getAnalysis <- function(x, name, sample = NULL) {
  validateSIRUP(x)
  name <- .sirup_analysis_name(name)
  analyses <- .sirup_analyses(x)
  if (!name %in% names(analyses)) {
    rlang::abort(paste0("Analysis '", name, "' is not stored in this SIRUP object."))
  }
  value <- analyses[[name]]
  if (!is.null(sample)) {
    .sirup_validate_sample_name(x, sample)
    if (!is.list(value) || is.null(names(value)) || !sample %in% names(value)) {
      rlang::abort(paste0("Analysis '", name, "' does not contain sample '", sample, "'."))
    }
    value <- value[[sample]]
  }
  value
}

#' @rdname sirup-analyses
#' @export
setAnalysis <- function(x, name, value, sample = NULL, validate = TRUE) {
  validateSIRUP(x)
  name <- .sirup_analysis_name(name)
  analyses <- .sirup_analyses(x)
  if (is.null(sample)) {
    analyses[name] <- list(value)
  } else {
    .sirup_validate_sample_name(x, sample)
    current <- analyses[[name]] %||% list()
    if (!is.list(current) || (length(current) > 0L && is.null(names(current)))) {
      rlang::abort(paste0("Analysis '", name, "' is not a named per-sample analysis list."))
    }
    current[[sample]] <- value
    analyses[name] <- list(current)
  }
  x$analyses <- analyses
  if (isTRUE(validate)) {
    validateSIRUP(x)
  }
  x
}

#' @rdname sirup-analyses
#' @export
listAnalyses <- function(x) {
  validateSIRUP(x)
  names(.sirup_analyses(x)) %||% character()
}

#' @rdname sirup-analyses
#' @export
removeAnalysis <- function(x, name, validate = TRUE) {
  validateSIRUP(x)
  name <- .sirup_analysis_name(name)
  analyses <- .sirup_analyses(x)
  analyses[name] <- NULL
  x$analyses <- analyses
  if (isTRUE(validate)) {
    validateSIRUP(x)
  }
  x
}

#' Update data stored in a SIRUP object
#'
#' `updateReceptors`, `updateSpatial`, and `updateAnalysis` apply a function to
#' one part of a SIRUP object. They return the updated SIRUP object.
#'
#' @details
#' `updateReceptors` needs a function that returns a data frame. The function
#' stores the result as a tibble and keeps the existing field mappings. If the
#' object has no clone field and the result has a `clone_id` column, the function
#' records `clone_id` as the clone field. Use `receptor_fields` to set or clear a
#' mapping yourself. If you give `sample`, the function applies `fun` to the rows
#' of that sample only. It then puts the result back in the full table.
#'
#' `updateSpatial` needs a function that returns a Seurat object. In a
#' multi-sample object, it applies `fun` to every Seurat object. If you give
#' `sample`, it updates only that object.
#'
#' `updateAnalysis` applies `fun` to the stored result called `name`. To store a
#' new result, use [setAnalysis()].
#'
#' @param x A SIRUP object.
#' @param fun Function, or the name of a function. The function receives the
#'   selected part of the SIRUP object as its first argument.
#' @param ... Other arguments for `fun`.
#' @param receptor_fields Named list of receptor field mappings to change. The
#'   names can be `spatial_id`, `image`, `sample`, `clone`, `locus`, or `count`.
#'   A `NULL` value clears an optional mapping.
#' @param sample Name of the spatial sample to update. Use it only for
#'   multi-sample objects.
#' @param name Name of the stored analysis. Give one non-empty character string.
#' @param validate Validate the updated SIRUP object?
#'
#' @return The updated SIRUP object.
#' @seealso [getReceptors()], [getSpatial()], [setAnalysis()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#'
#' # Remove the sequences in clone 3.
#' x2 <- updateReceptors(x, dplyr::filter, clone_id != "3")
#' nrow(getReceptors(x2))
#'
#' # Add a Seurat metadata column.
#' x3 <- updateSpatial(x, function(spatial) {
#'   spatial$gene_rich <- spatial$nFeature_Spatial >= 8
#'   spatial
#' })
#' table(getSpatial(x3)$gene_rich)
#'
#' # Change a stored result.
#' x4 <- setAnalysis(x, "n_clones", 5)
#' x4 <- updateAnalysis(x4, "n_clones", function(value) value + 1)
#' getAnalysis(x4, "n_clones")
#' @name sirup-transform
NULL

#' @rdname sirup-transform
#' @export
updateReceptors <- function(x,
                            fun,
                            ...,
                            sample = NULL,
                            receptor_fields = NULL,
                            validate = TRUE) {
  validateSIRUP(x)
  if (is.null(sample)) {
    receptors <- .sirup_apply(fun, x$receptors, ...)
  } else {
    selected <- .sirup_receptors_for_sample(x, sample)
    updated <- .sirup_apply(fun, selected, ...)
    if (!is.data.frame(updated)) {
      rlang::abort("`fun` must return a data frame/tibble when used with `updateReceptors()`.")
    }
    sample_col <- .sirup_receptor_fields(x)$sample
    updated <- tibble::as_tibble(updated)
    updated[[sample_col]] <- sample
    original <- x$receptors
    keep <- as.character(original[[sample_col]]) != sample | is.na(original[[sample_col]])
    receptors <- .sirup_bind_receptor_tables(original[keep, , drop = FALSE], updated)
  }
  if (!is.data.frame(receptors)) {
    rlang::abort("`fun` must return a data frame/tibble when used with `updateReceptors()`.")
  }
  x$receptors <- tibble::as_tibble(receptors)
  x$metadata$receptor_fields <- .sirup_update_receptor_fields(
    x,
    receptor_fields
  )
  if (isTRUE(validate)) {
    validateSIRUP(x)
  }
  x
}

#' @rdname sirup-transform
#' @export
updateSpatial <- function(x, fun, ..., sample = NULL, validate = TRUE) {
  validateSIRUP(x)
  if (is.null(sample)) {
    if (.sirup_is_spatial_list(x$spatial)) {
      spatial <- lapply(x$spatial, function(spatial_object) .sirup_apply(fun, spatial_object, ...))
      names(spatial) <- names(x$spatial)
    } else {
      spatial <- .sirup_apply(fun, x$spatial, ...)
    }
  } else {
    if (!.sirup_is_spatial_list(x$spatial)) {
      rlang::abort("`sample` can only be used when this SIRUP object contains a spatial list.")
    }
    .sirup_validate_sample_name(x, sample)
    spatial <- x$spatial
    spatial[[sample]] <- .sirup_apply(fun, spatial[[sample]], ...)
  }
  tryCatch(
    .sirup_validate_spatial_component(spatial),
    error = function(e) {
      rlang::abort("`fun` must return a Seurat object or a named list of Seurat objects when used with `updateSpatial()`.")
    }
  )
  x$spatial <- spatial
  x$metadata$sample_names <- .sirup_spatial_names(x$spatial)
  if (isTRUE(validate)) {
    validateSIRUP(x)
  }
  x
}

#' @rdname sirup-transform
#' @export
updateAnalysis <- function(x, name, fun, ..., sample = NULL, validate = TRUE) {
  value <- .sirup_apply(fun, getAnalysis(x, name, sample = sample), ...)
  setAnalysis(x, name, value, sample = sample, validate = validate)
}

#' @export
print.sirup <- function(x, ...) {
  validateSIRUP(x)
  info <- .sirup_join_diagnostics(x)
  receptor_fields <- .sirup_receptor_fields(x)
  polygon_count <- nrow(.sirup_polygons(x))
  spatial_count <- length(.sirup_spatial_names(x$spatial))
  image_count <- sum(vapply(.sirup_spatial_objects(x$spatial), function(spatial) {
    length(SeuratObject::Images(spatial))
  }, integer(1)))

  cat("<SIRUP object>\n")
  cat("Spatial data:    Seurat, ", format(info$n_spatial_ids, big.mark = ","), " ",
      .sirup_count_label(info$n_spatial_ids, "spatial ID"), ", ",
      format(.sirup_spatial_feature_count(x$spatial), big.mark = ","), " ",
      .sirup_count_label(.sirup_spatial_feature_count(x$spatial), "gene"), sep = "")
  if (spatial_count > 1L) {
    cat(", ", spatial_count, " ", .sirup_count_label(spatial_count, "spatial sample"), sep = "")
  }
  if (!is.na(image_count)) {
    cat(", ", image_count, " ", .sirup_count_label(image_count, "image"), sep = "")
  }
  cat("\n")
  cat("Receptor data:   tibble, ", format(nrow(x$receptors), big.mark = ","), " ",
      .sirup_count_label(nrow(x$receptors), "receptor sequence"), " x ",
      format(ncol(x$receptors), big.mark = ","), " ",
      .sirup_count_label(ncol(x$receptors), "column"), "\n", sep = "")
  cat("Spatial overlap: ", format(info$n_matched_ids, big.mark = ","), "/",
      format(info$n_receptor_ids, big.mark = ","), " ",
      .sirup_count_label(info$n_receptor_ids, "receptor ID"), " in spatial data (",
      formatC(100 * info$id_match_rate, format = "f", digits = 1), "%)\n", sep = "")
  if (spatial_count > 1L) {
    cat("Samples:         ", spatial_count, "\n", sep = "")
  }
  cat("Spatial ID field: ", receptor_fields$spatial_id, "\n", sep = "")
  if (spatial_count > 1L && !is.null(receptor_fields$sample)) {
    cat("Sample ID field: ", receptor_fields$sample, "\n", sep = "")
  }
  if (!is.null(receptor_fields$clone)) {
    cat("Clone ID field:  ", receptor_fields$clone, "\n", sep = "")
  }
  cat("Analyses:        ", length(.sirup_analyses(x)), "\n", sep = "")
  if (polygon_count > 0L) {
    cat("Polygons:        ", format(polygon_count, big.mark = ","), "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.sirup <- function(object, ...) {
  validateSIRUP(object)
  info <- .sirup_join_diagnostics(object)
  receptor_fields <- .sirup_receptor_fields(object)
  sample_summary <- .sirup_sample_summary(object)
  structure(
    c(
      info,
      list(
        n_receptor_records = nrow(object$receptors),
        n_receptor_columns = ncol(object$receptors),
        n_spatial_genes = .sirup_spatial_feature_count(object$spatial),
        spatial_id = receptor_fields$spatial_id,
        clone = receptor_fields$clone,
        image_scale = object$metadata$image_scale,
        swap_coords = object$metadata$swap_coords,
        loci = .sirup_field_values(object, "locus"),
        n_analyses = length(.sirup_analyses(object)),
        n_polygons = nrow(.sirup_polygons(object)),
        samples = sample_summary,
        analyses = listAnalyses(object)
      )
    ),
    class = c("summary.sirup", "list")
  )
}

#' @export
print.summary.sirup <- function(x, ...) {
  cat("SIRUP summary\n")
  labels <- c(
    "Spatial IDs:",
    "Receptor records:",
    "Receptor IDs:",
    "Matched IDs:",
    "Matched records:",
    "Image scale:",
    "Swap coords:",
    "Stored analyses:"
  )
  values <- c(
    format(x$n_spatial_ids, big.mark = ","),
    format(x$n_receptor_records, big.mark = ","),
    format(x$n_receptor_ids, big.mark = ","),
    paste0(format(x$n_matched_ids, big.mark = ","), " (",
           formatC(100 * x$id_match_rate, format = "f", digits = 1), "%)"),
    format(x$n_matched_records, big.mark = ","),
    .sirup_format_summary_option(x$image_scale),
    .sirup_format_summary_option(x$swap_coords),
    format(x$n_analyses, big.mark = ",")
  )
  if (!is.null(x$n_polygons) && x$n_polygons > 0L) {
    labels <- c(labels, "Polygons:")
    values <- c(values, format(x$n_polygons, big.mark = ","))
  }
  label_width <- max(nchar(labels))
  for (i in seq_along(labels)) {
    cat(sprintf("  %-*s %s\n", label_width, labels[[i]], values[[i]]))
  }
  samples <- x$samples
  if (is.list(samples) && length(samples) == 1L && is.data.frame(samples[[1L]])) {
    samples <- samples[[1L]]
  }
  if (is.data.frame(samples) && nrow(samples) > 0L) {
    cat("\nPer-sample summary\n")
    print(samples, row.names = FALSE)
  }
  invisible(x)
}

.sirup_format_summary_option <- function(x) {
  if (is.null(x)) {
    return("<unset>")
  }
  if (is.list(x) && !is.data.frame(x)) {
    values <- unlist(x, recursive = FALSE, use.names = TRUE)
  } else {
    values <- x
  }
  values <- as.vector(values, mode = "any")
  value_names <- names(unlist(x, recursive = FALSE, use.names = TRUE))
  formatted <- as.character(values)
  if (!is.null(value_names) && length(value_names) == length(formatted) &&
      any(!is.na(value_names) & nzchar(value_names))) {
    named <- !is.na(value_names) & nzchar(value_names)
    formatted[named] <- paste0(value_names[named], "=", formatted[named])
  }
  paste(formatted, collapse = ", ")
}

# Validate the optional named analysis list. Unique names make retrieval and
# replacement unambiguous, while values are deliberately left unrestricted.
.sirup_validate_analyses <- function(analyses) {
  if (!is.list(analyses)) {
    rlang::abort("`analyses` must be a named list.")
  }
  if (length(analyses) > 0L &&
      (is.null(names(analyses)) ||
       any(is.na(names(analyses)) | !nzchar(names(analyses))) ||
       anyDuplicated(names(analyses)))) {
    rlang::abort("`analyses` must have unique, non-empty names.")
  }
  invisible(analyses)
}

# Return the stored analysis list, treating pre-analysis SIRUP objects as
# having an empty list for backward compatibility.
.sirup_analyses <- function(x) {
  if (is.null(x$analyses)) {
    return(list())
  }
  x$analyses
}

# Bring an on-disk or legacy SIRUP object up to the current schema. Older
# objects may lack the `polygons`/`analyses` components (the "legacy" name
# sets tolerated by validateSIRUP()) or the spatial-metadata defaults. This
# rebuilds the canonical five-component object and backfills defaults without
# re-guessing the receptor field mappings. Version-stamp attributes added by
# saveSIRUP() are intentionally dropped here; readSIRUP() reads them first.
.sirup_upgrade <- function(x) {
  if (!isSIRUP(x)) {
    rlang::abort("`x` must inherit from class 'sirup'.")
  }
  metadata <- if (is.list(x$metadata)) x$metadata else list()
  if (is.null(metadata$polygons_include_boundary)) {
    metadata$polygons_include_boundary <- TRUE
  }
  if (is.null(metadata$sample_names)) {
    metadata$sample_names <- .sirup_spatial_names(x$spatial)
  }
  structure(
    list(
      spatial = x$spatial,
      receptors = x$receptors,
      metadata = metadata,
      polygons = x$polygons %||% .sirup_empty_polygons(),
      analyses = x$analyses %||% list()
    ),
    class = c("sirup", "list")
  )
}

.sirup_polygons <- function(x) {
  if (is.null(x$polygons)) {
    return(.sirup_empty_polygons())
  }
  x$polygons
}

.sirup_subset_polygons_component <- function(polygons,
                                             sample = NULL,
                                             relabel_single_sample = FALSE) {
  if (is.null(polygons) || nrow(polygons) == 0L || is.null(sample)) {
    return(polygons %||% .sirup_empty_polygons())
  }
  sample <- unique(as.character(sample))
  sample <- sample[!is.na(sample) & nzchar(sample)]
  if (length(sample) == 0L) {
    return(polygons[0, , drop = FALSE])
  }
  top_sample_col <- intersect(c("sample", "spatial_sample"), colnames(polygons))
  top_sample_col <- if (length(top_sample_col) > 0L) top_sample_col[[1L]] else NULL
  keep <- vapply(seq_len(nrow(polygons)), function(i) {
    if (!is.null(top_sample_col)) {
      value <- as.character(polygons[[top_sample_col]][[i]])
      if (!is.na(value) && nzchar(value)) {
        return(value %in% sample)
      }
    }
    vertices <- polygons$vertices[[i]]
    values <- .sirup_polygon_sample_values(vertices)
    length(values) == 0L || any(values %in% sample)
  }, logical(1))
  out <- polygons[keep, , drop = FALSE]
  if (!isTRUE(relabel_single_sample) || nrow(out) == 0L) {
    return(out)
  }
  labels <- vapply(seq_len(nrow(out)), function(i) {
    vertices <- tibble::as_tibble(out$vertices[[i]])
    if (".sirup_polygon_label" %in% colnames(vertices)) {
      values <- unique(stats::na.omit(as.character(vertices$.sirup_polygon_label)))
      values <- values[nzchar(values)]
      if (length(values) > 0L) {
        return(values[[1L]])
      }
    }
    parsed <- .sirup_parse_spatial_image_label(as.character(out$polygon_id[[i]]))
    if (!is.null(parsed)) {
      return(parsed$image)
    }
    as.character(out$polygon_id[[i]])
  }, character(1))
  labels[is.na(labels) | !nzchar(labels)] <- as.character(out$polygon_id[is.na(labels) | !nzchar(labels)])
  if (!anyDuplicated(labels)) {
    out$polygon_id <- labels
  }
  out
}

.sirup_empty_polygons <- function() {
  tibble::tibble(
    polygon_id = character(),
    vertices = list()
  )
}

.sirup_check_merge_receptor_fields <- function(objects) {
  fields <- lapply(objects, .sirup_receptor_fields)
  field_names <- names(fields[[1L]])
  out <- fields[[1L]]
  for (field in field_names) {
    if (identical(field, "sample")) {
      next
    }
    values <- vapply(fields, function(x) x[[field]] %||% NA_character_, character(1))
    values[is.na(values)] <- NA_character_
    unique_values <- unique(values)
    unique_values <- unique_values[!is.na(unique_values)]
    if (length(unique_values) > 1L || (identical(field, "spatial_id") && any(is.na(values)))) {
      rlang::abort(paste0(
        "SIRUP objects cannot be merged because receptor field `", field,
        "` is not consistent across inputs: ",
        paste(ifelse(is.na(values), "NULL", values), collapse = ", "), "."
      ))
    }
    out[[field]] <- if (length(unique_values) == 0L) NULL else unique_values[[1L]]
  }
  out <- out[c("spatial_id", "image", "sample", "clone", "locus", "count")]
  names(out) <- c("spatial_id", "image", "sample", "clone", "locus", "count")
  out
}

.sirup_merge_receptors_for_sample <- function(object, sample, sample_field) {
  receptors <- tibble::as_tibble(getReceptors(object))
  if (sample_field %in% colnames(receptors)) {
    values <- unique(stats::na.omit(as.character(receptors[[sample_field]])))
    values <- values[nzchar(values)]
    if (length(values) > 0L && (!identical(length(values), 1L) || !identical(values[[1L]], sample))) {
      rlang::abort(paste0(
        "`sample_field = \"", sample_field, "\"` already exists in one receptor table, ",
        "but its non-missing values do not match the requested sample `", sample,
        "`. Found: ", paste(values, collapse = ", "), "."
      ))
    }
  }
  receptors[[sample_field]] <- sample
  receptors
}

.sirup_polygons_for_merge <- function(object, sample) {
  polygons <- .sirup_polygons(object)
  if (nrow(polygons) == 0L) {
    return(polygons)
  }
  polygons <- tibble::as_tibble(polygons)
  labels <- vapply(seq_len(nrow(polygons)), function(i) {
    vertices <- tibble::as_tibble(polygons$vertices[[i]])
    values <- if (".sirup_polygon_label" %in% colnames(vertices)) {
      unique(stats::na.omit(as.character(vertices$.sirup_polygon_label)))
    } else {
      character()
    }
    values <- values[nzchar(values)]
    if (length(values) > 0L) values[[1L]] else as.character(polygons$polygon_id[[i]])
  }, character(1))
  polygons$polygon_id <- paste(sample, as.character(polygons$polygon_id), sep = "::")
  polygons$sample <- sample
  polygons$vertices <- lapply(seq_len(nrow(polygons)), function(i) {
    vertices <- tibble::as_tibble(polygons$vertices[[i]])
    vertices$sample <- sample
    vertices$.sirup_polygon_label <- labels[[i]]
    vertices
  })
  polygons
}

.sirup_bind_polygons_for_merge <- function(polygons) {
  polygons <- Filter(function(x) !is.null(x) && nrow(x) > 0L, polygons)
  if (length(polygons) == 0L) {
    return(.sirup_empty_polygons())
  }
  out <- dplyr::bind_rows(polygons)
  .sirup_validate_polygons(out)
  tibble::as_tibble(out)
}

.sirup_merge_analyses <- function(out, analyses, sample) {
  if (length(analyses) == 0L) {
    return(out)
  }
  for (name in names(analyses)) {
    current <- out[[name]] %||% list()
    current[[sample]] <- analyses[[name]]
    out[[name]] <- current
  }
  out
}

.sirup_merge_metadata <- function(objects, samples) {
  reserved <- c("receptor_fields", "image_scale", "swap_coords", "polygons_include_boundary", "polygons_path")
  metadata <- lapply(objects, function(object) object$metadata[setdiff(names(object$metadata), reserved)])
  all_names <- unique(unlist(lapply(metadata, names), use.names = FALSE))
  out <- list()
  differing <- list()
  for (name in all_names) {
    values <- lapply(metadata, `[[`, name)
    present <- vapply(values, Negate(is.null), logical(1))
    if (all(present) && all(vapply(values[-1], identical, logical(1), values[[1L]]))) {
      out[[name]] <- values[[1L]]
    } else {
      differing[[name]] <- stats::setNames(values, samples)
    }
  }
  if (length(differing) > 0L) {
    out$merged_from <- differing
  }
  out
}

.sirup_merge_sample_options <- function(objects, samples, field) {
  values <- lapply(seq_along(objects), function(i) {
    value <- objects[[i]]$metadata[[field]] %||% NULL
    if (is.null(value)) {
      return(NULL)
    }
    old_sample <- .sirup_spatial_names(objects[[i]]$spatial)[[1L]]
    .sirup_sample_option(value, sample = old_sample, argument = field)
  })
  if (all(vapply(values, is.null, logical(1)))) {
    return(NULL)
  }
  if (any(vapply(values, is.null, logical(1)))) {
    rlang::abort(paste0("Cannot merge SIRUP objects because some inputs are missing `metadata$", field, "`."))
  }
  out <- unlist(values, recursive = FALSE, use.names = FALSE)
  names(out) <- samples
  out
}

.sirup_read_polygons_path <- function(path, sample = NULL) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    rlang::abort("`polygons_path` must be one non-empty character path.")
  }
  if (!file.exists(path)) {
    rlang::abort("`polygons_path` does not exist: ", path)
  }
  out <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(sample)) {
    out$sample <- sample
  }
  out
}

.sirup_read_polygons_paths <- function(path, sample = NULL) {
  if (is.list(path) && !is.data.frame(path)) {
    path_names <- names(path)
    path <- unlist(path, recursive = FALSE, use.names = FALSE)
    names(path) <- path_names
  }
  if (!is.character(path) || length(path) < 1L || any(is.na(path) | !nzchar(path))) {
    rlang::abort("`polygons_path` must be one path, or a named vector/list of paths keyed by spatial sample.")
  }
  if (!is.null(sample)) {
    if (!is.character(sample) || length(sample) != length(path) || any(is.na(sample) | !nzchar(sample))) {
      rlang::abort("`sample` must be one non-empty sample name per `polygons_path`.")
    }
  } else if (!is.null(names(path)) && all(!is.na(names(path)) & nzchar(names(path)))) {
    sample <- names(path)
  }
  out <- lapply(seq_along(path), function(i) {
    .sirup_read_polygons_path(path[[i]], sample = if (is.null(sample)) NULL else sample[[i]])
  })
  dplyr::bind_rows(out)
}

.sirup_validate_multisample_polygons <- function(x, polygons) {
  .sirup_validate_multisample_polygon_input(x$spatial, polygons)
}

.sirup_validate_multisample_polygon_input <- function(spatial, polygons) {
  if (!.sirup_is_spatial_list(spatial) || is.null(polygons)) {
    return(invisible(polygons))
  }
  polygons <- tibble::as_tibble(polygons)
  if (all(c("polygon_id", "vertices") %in% names(polygons))) {
    top_level_sample <- any(c("sample", "spatial_sample") %in% colnames(polygons))
    has_sample <- top_level_sample | vapply(polygons$vertices, function(vertices) {
      any(c("sample", "spatial_sample") %in% colnames(vertices))
    }, logical(1L))
    if (!all(has_sample)) {
      rlang::abort("Polygons added to a multi-spatial SIRUP object must include `sample` or `spatial_sample` for every polygon.")
    }
    return(invisible(polygons))
  }
  if (!any(c("sample", "spatial_sample") %in% colnames(polygons))) {
    rlang::abort("Polygons added to a multi-spatial SIRUP object must include `sample` or `spatial_sample`, or use named `polygons_path` values keyed by sample.")
  }
  invisible(polygons)
}

.sirup_populate_spatial_sample_metadata <- function(spatial, sample_col) {
  if (!.sirup_is_spatial_list(spatial) || is.null(sample_col)) {
    return(spatial)
  }
  for (sample_name in names(spatial)) {
    metadata <- spatial[[sample_name]][[]]
    if (!sample_col %in% colnames(metadata)) {
      spatial[[sample_name]][[sample_col]] <- sample_name
    }
  }
  spatial
}

.sirup_prepare_polygons <- function(polygons) {
  if (is.null(polygons)) {
    return(.sirup_empty_polygons())
  }
  polygons <- tibble::as_tibble(polygons)
  if (all(c("polygon_id", "vertices") %in% names(polygons))) {
    out <- polygons
    if (!any(c("sample", "spatial_sample") %in% colnames(out))) {
      samples <- vapply(out$vertices, function(vertices) {
        values <- .sirup_polygon_sample_values(vertices)
        if (length(values) == 0L) NA_character_ else values[[1L]]
      }, character(1))
      if (any(!is.na(samples))) {
        out$sample <- samples
      }
    }
  } else {
    polygon_col <- if ("polygon_id" %in% names(polygons)) "polygon_id" else "polygon"
    if (!all(c(polygon_col, "x", "y") %in% names(polygons))) {
      rlang::abort("Flat polygons must contain `polygon_id` or `polygon`, plus `x` and `y`.")
    }
    sample_col <- intersect(c("sample", "spatial_sample"), names(polygons))
    sample_col <- if (length(sample_col) > 0L) sample_col[[1L]] else NULL
    polygon_values <- as.character(polygons[[polygon_col]])
    split_values <- polygon_values
    if (!is.null(sample_col)) {
      polygon_samples <- as.character(polygons[[sample_col]])
      polygon_sample_pairs <- unique(data.frame(
        polygon = polygon_values,
        sample = polygon_samples,
        stringsAsFactors = FALSE
      ))
      duplicated_across_samples <- any(duplicated(polygon_sample_pairs$polygon)) ||
        any(duplicated(polygon_sample_pairs$polygon, fromLast = TRUE))
      if (isTRUE(duplicated_across_samples)) {
        split_values <- paste(polygon_samples, polygon_values, sep = "::")
      }
    }
    split_rows <- split(polygons, split_values, drop = TRUE)
    out <- tibble::tibble(
      polygon_id = names(split_rows),
      sample = vapply(split_rows, function(rows) {
        if (is.null(sample_col)) {
          return(NA_character_)
        }
        values <- unique(stats::na.omit(as.character(rows[[sample_col]])))
        values <- values[nzchar(values)]
        if (length(values) == 0L) NA_character_ else values[[1L]]
      }, character(1)),
      vertices = lapply(split_rows, function(rows) {
        polygon_label <- unique(as.character(rows[[polygon_col]]))
        polygon_label <- polygon_label[!is.na(polygon_label) & nzchar(polygon_label)]
        polygon_label <- if (length(polygon_label) == 0L) NA_character_ else polygon_label[[1L]]
        keep <- intersect(c("spatial_id", "x", "y", "image", "sample", "spatial_sample"), names(rows))
        rows <- rows[, keep, drop = FALSE]
        if (!"spatial_id" %in% names(rows)) rows$spatial_id <- NA_character_
        rows$.sirup_polygon_label <- polygon_label
        rows[, c("spatial_id", "x", "y", setdiff(names(rows), c("spatial_id", "x", "y"))), drop = FALSE]
      })
    )
  }
  if ("sample" %in% names(out) && all(is.na(out$sample))) {
    out$sample <- NULL
  }
  out$vertices <- lapply(out$vertices, function(vertices) {
    vertices <- tibble::as_tibble(vertices)
    if (nrow(vertices) > 0L &&
        (vertices$x[[1L]] != vertices$x[[nrow(vertices)]] ||
         vertices$y[[1L]] != vertices$y[[nrow(vertices)]])) {
      vertices <- tibble::as_tibble(rbind(vertices, vertices[1, , drop = FALSE]))
    }
    vertices
  })
  .sirup_validate_polygons(out)
  out
}

# Add stored polygon IDs and per-polygon membership flags to both SIRUP data
# components. Receptors and Seurat metadata receive `polygon_id`, with `NA` for
# records outside all polygons.
.sirup_annotate_stored_polygons <- function(x) {
  if (nrow(.sirup_polygons(x)) == 0L) {
    return(x)
  }
  polygon <- .sirup_polygon_vertices_table(.sirup_polygons(x))
  x$spatial <- assignSpatialPolygons(
    object = x$spatial,
    polygon = polygon,
    image_scale = x$metadata$image_scale,
    swap_coords = x$metadata$swap_coords,
    polygon_id = "polygon",
    metadata = "polygon_id",
    include_boundary = isTRUE(x$metadata$polygons_include_boundary %||% TRUE)
  )
  x$spatial <- .sirup_ensure_spatial_metadata_column(x$spatial, "polygon_id", NA_character_)
  x$receptors <- assignReceptorPolygons(
    receptors = x$receptors,
    object = x,
    polygon = polygon,
    column = "polygon_id",
    spatial_id = .sirup_receptor_fields(x)$spatial_id,
    image_scale = x$metadata$image_scale,
    swap_coords = x$metadata$swap_coords,
    polygon_id = "polygon",
    include_boundary = isTRUE(x$metadata$polygons_include_boundary %||% TRUE)
  )
  x
}

.sirup_ensure_spatial_metadata_column <- function(spatial, column, value = NA_character_) {
  if (.sirup_is_spatial_list(spatial)) {
    for (sample in names(spatial)) {
      if (!column %in% colnames(spatial[[sample]][[]])) {
        spatial[[sample]][[column]] <- value
      }
    }
    return(spatial)
  }
  if (!column %in% colnames(spatial[[]])) {
    spatial[[column]] <- value
  }
  spatial
}

.sirup_validate_polygons <- function(polygons) {
  if (!tibble::is_tibble(polygons)) {
    rlang::abort("`polygons` must be a tibble.")
  }
  if (!all(c("polygon_id", "vertices") %in% names(polygons))) {
    rlang::abort("`polygons` must contain `polygon_id` and `vertices` columns.")
  }
  if (anyNA(polygons$polygon_id) || any(!nzchar(as.character(polygons$polygon_id))) ||
      anyDuplicated(as.character(polygons$polygon_id))) {
    rlang::abort("`polygons$polygon_id` must contain unique, non-empty IDs.")
  }
  if (!is.list(polygons$vertices) || length(polygons$vertices) != nrow(polygons)) {
    rlang::abort("`polygons$vertices` must be a list-column with one entry per polygon.")
  }
  for (i in seq_len(nrow(polygons))) {
    vertices <- tibble::as_tibble(polygons$vertices[[i]])
    if (!all(c("x", "y") %in% names(vertices))) {
      rlang::abort("Each polygon vertices table must contain `x` and `y`.")
    }
    if (nrow(vertices) < 3L || any(!is.finite(as.numeric(vertices$x))) ||
        any(!is.finite(as.numeric(vertices$y)))) {
      rlang::abort("Each polygon must contain at least three finite vertices.")
    }
  }
  invisible(polygons)
}

# Return TRUE for a single non-missing, non-empty character string.
.sirup_is_non_empty_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

# Validate and normalize one analysis name used for lookup or replacement.
.sirup_analysis_name <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    rlang::abort("`name` must be one non-empty character string.")
  }
  name
}

# Resolve a function supplied to the with/update APIs and apply it to the
# selected SIRUP component as its first argument.
.sirup_apply <- function(fun, sirup_data, ...) {
  fun <- match.fun(fun)
  fun(sirup_data, ...)
}

.sirup_resolve_image_scale <- function(image_scale) {
  if (is.null(image_scale)) {
    return(NULL)
  }
  if (is.list(image_scale) && !is.data.frame(image_scale)) {
    values <- unlist(image_scale, recursive = FALSE, use.names = FALSE)
    value_names <- names(image_scale)
    if (length(values) == 0L ||
        !all(vapply(values, function(value) is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value), logical(1)))) {
      rlang::abort("`image_scale` list values must be non-empty character strings.")
    }
    if (!is.null(value_names) && any(is.na(value_names))) {
      rlang::abort("`image_scale` names cannot be `NA`.")
    }
    image_scale <- unlist(image_scale, recursive = FALSE, use.names = TRUE)
    names(image_scale) <- value_names
    return(image_scale)
  }
  if (!is.character(image_scale) || length(image_scale) < 1L ||
      any(is.na(image_scale) | !nzchar(image_scale))) {
    rlang::abort("`image_scale` must be `NULL`, one non-empty character string, or a named character vector/list.")
  }
  if (length(image_scale) > 1L &&
      (is.null(names(image_scale)) || any(is.na(names(image_scale))) || any(!nzchar(names(image_scale))))) {
    rlang::abort("Multi-value `image_scale` must be named by spatial sample.")
  }
  image_scale
}

.sirup_resolve_swap_coords <- function(swap_coords) {
  if (is.null(swap_coords)) {
    return(NULL)
  }
  if (is.list(swap_coords) && !is.data.frame(swap_coords)) {
    values <- unlist(swap_coords, recursive = FALSE, use.names = FALSE)
    value_names <- names(swap_coords)
    if (length(values) == 0L ||
        !all(vapply(values, function(value) is.logical(value) && length(value) == 1L && !is.na(value), logical(1)))) {
      rlang::abort("`swap_coords` list values must be `TRUE` or `FALSE`.")
    }
    if (!is.null(value_names) && any(is.na(value_names))) {
      rlang::abort("`swap_coords` names cannot be `NA`.")
    }
    swap_coords <- unlist(swap_coords, recursive = FALSE, use.names = TRUE)
    names(swap_coords) <- value_names
    return(swap_coords)
  }
  if (!is.logical(swap_coords) || length(swap_coords) < 1L || any(is.na(swap_coords))) {
    rlang::abort("`swap_coords` must be `NULL`, `TRUE`, `FALSE`, or a named logical vector/list.")
  }
  if (length(swap_coords) > 1L &&
      (is.null(names(swap_coords)) || any(is.na(names(swap_coords))) || any(!nzchar(names(swap_coords))))) {
    rlang::abort("Multi-value `swap_coords` must be named by spatial sample.")
  }
  swap_coords
}

.sirup_validate_spatial_metadata <- function(metadata) {
  .sirup_resolve_image_scale(metadata$image_scale %||% NULL)
  .sirup_resolve_swap_coords(metadata$swap_coords %||% NULL)
  invisible(metadata)
}

.sirup_validate_metadata_sample_names <- function(x) {
  sample_names <- x$metadata$sample_names %||% NULL
  if (is.null(sample_names)) {
    return(invisible(x))
  }
  expected <- .sirup_spatial_names(x$spatial)
  if (!is.character(sample_names) ||
      length(sample_names) != length(expected) ||
      any(is.na(sample_names) | !nzchar(sample_names)) ||
      !identical(as.character(sample_names), expected)) {
    rlang::abort("`metadata$sample_names` must match the SIRUP object's spatial sample names.")
  }
  invisible(x)
}

.sirup_sample_option <- function(value, sample = NULL, argument = "option") {
  if (is.null(value) || length(value) == 0L) {
    return(NULL)
  }
  value_names <- names(value)
  if (!is.null(sample) && !is.null(value_names) && sample %in% value_names) {
    return(unname(value[[match(sample, value_names)]]))
  }
  unnamed <- if (is.null(value_names)) seq_along(value) else which(is.na(value_names) | !nzchar(value_names))
  if (length(unnamed) > 0L) {
    return(unname(value[[unnamed[[1L]]]]))
  }
  if (length(value) == 1L && (is.null(value_names) || is.na(value_names[[1L]]) || !nzchar(value_names[[1L]]))) {
    return(unname(value[[1L]]))
  }
  if (!is.null(sample)) {
    rlang::abort(paste0("`", argument, "` does not contain a value for spatial sample \"", sample, "\" and has no global fallback."))
  }
  if (length(value) == 1L && !is.null(value_names)) {
    return(unname(value[[1L]]))
  }
  value
}

.sirup_require_spatial_options <- function(image_scale, swap_coords, sample = NULL) {
  image_scale <- .sirup_resolve_image_scale(image_scale)
  swap_coords <- .sirup_resolve_swap_coords(swap_coords)
  image_scale <- .sirup_sample_option(image_scale, sample = sample, argument = "image_scale")
  swap_coords <- .sirup_sample_option(swap_coords, sample = sample, argument = "swap_coords")
  if (is.null(image_scale)) {
    rlang::abort("`image_scale` must be supplied explicitly or stored in a SIRUP object.")
  }
  if (is.null(swap_coords)) {
    rlang::abort("`swap_coords` must be supplied explicitly or stored in a SIRUP object.")
  }
  list(
    image_scale = image_scale,
    swap_coords = swap_coords
  )
}

.sirup_spatial_options <- function(object,
                                   image_scale = NULL,
                                   swap_coords = NULL) {
  image_scale <- .sirup_resolve_image_scale(image_scale)
  swap_coords <- .sirup_resolve_swap_coords(swap_coords)
  if (isSIRUP(object)) {
    validateSIRUP(object)
    metadata <- object$metadata
    if (is.null(image_scale)) {
      image_scale <- metadata$image_scale %||% NULL
    }
    if (is.null(swap_coords)) {
      swap_coords <- metadata$swap_coords %||% NULL
    }
  }
  .sirup_require_spatial_options(image_scale, swap_coords)
}

.sirup_spatial_context <- function(object, image_scale = NULL, swap_coords = NULL) {
  options <- .sirup_spatial_options(
    object = object,
    image_scale = image_scale,
    swap_coords = swap_coords
  )
  list(
    object = .sirup_spatial_input(object),
    image_scale = options$image_scale,
    swap_coords = options$swap_coords
  )
}

# Preserve receptor mappings after an update, infer conventional clone/locus
# fields when newly added, and apply any explicit mapping changes last.
.sirup_update_receptor_fields <- function(x, receptor_fields = NULL) {
  fields <- .sirup_receptor_fields(x)
  conventional <- c(clone = "clone_id", locus = "locus")
  for (field in names(conventional)) {
    if (is.null(fields[[field]]) && conventional[[field]] %in% names(x$receptors)) {
      fields[[field]] <- conventional[[field]]
    }
  }

  if (is.null(receptor_fields)) {
    return(fields)
  }
  if (!is.list(receptor_fields) || is.null(names(receptor_fields)) ||
      any(is.na(names(receptor_fields)) | !nzchar(names(receptor_fields))) ||
      anyDuplicated(names(receptor_fields))) {
    rlang::abort("`receptor_fields` must be a named list.")
  }
  allowed <- c("spatial_id", "image", "sample", "clone", "locus", "count")
  unknown <- setdiff(names(receptor_fields), allowed)
  if (length(unknown) > 0L) {
    rlang::abort(paste0(
      "Unknown receptor field mapping",
      if (length(unknown) == 1L) " " else "s ",
      paste0("'", unknown, "'", collapse = ", "),
      "."
    ))
  }
  for (field in names(receptor_fields)) {
    fields[field] <- list(receptor_fields[[field]])
  }
  fields
}

# Validate optional receptor field mappings while allowing default clone/locus
# fields to quietly become NULL when their conventional columns are absent.
.sirup_resolve_optional_col <- function(receptors, column, argument, allow_missing = FALSE) {
  if (is.null(column)) {
    return(NULL)
  }
  if (!is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column)) {
    rlang::abort(paste0("`", argument, "` must be `NULL` or one column name."))
  }
  if (!column %in% colnames(receptors)) {
    if (isTRUE(allow_missing)) {
      return(NULL)
    }
    rlang::abort(paste0("`", argument, "` refers to missing column '", column, "'."))
  }
  column
}

# Pluralize short labels in print methods.
.sirup_count_label <- function(n, singular) {
  if (identical(as.numeric(n), 1)) singular else paste0(singular, "s")
}

# Validate the spatial component of a SIRUP object. A single Seurat object is
# the legacy shape; a named list of Seurat objects represents multiple spatial
# contexts keyed by the receptor sample field.
.sirup_validate_spatial_component <- function(spatial) {
  if (inherits(spatial, "Seurat")) {
    return(invisible(spatial))
  }
  if (!is.list(spatial) || length(spatial) == 0L) {
    rlang::abort("`x$spatial` must be a Seurat object or a named list of Seurat objects.")
  }
  if (is.null(names(spatial)) || any(is.na(names(spatial)) | !nzchar(names(spatial))) ||
      anyDuplicated(names(spatial))) {
    rlang::abort("A spatial list must have unique, non-empty names.")
  }
  bad <- !vapply(spatial, inherits, logical(1), what = "Seurat")
  if (any(bad)) {
    rlang::abort("Every element of a spatial list must be a Seurat object.")
  }
  invisible(spatial)
}

.sirup_is_spatial_list <- function(spatial) {
  is.list(spatial) &&
    !inherits(spatial, "Seurat") &&
    length(spatial) > 0L &&
    !is.null(names(spatial)) &&
    all(!is.na(names(spatial)) & nzchar(names(spatial)))
}

.sirup_spatial_objects <- function(spatial) {
  if (.sirup_is_spatial_list(spatial)) spatial else list(spatial = spatial)
}

.sirup_spatial_names <- function(spatial) {
  if (.sirup_is_spatial_list(spatial)) names(spatial) else "spatial"
}

.sirup_spatial_feature_count <- function(spatial) {
  features <- unique(unlist(lapply(.sirup_spatial_objects(spatial), rownames), use.names = FALSE))
  length(features)
}

.sirup_spatial_cell_count <- function(spatial) {
  sum(vapply(.sirup_spatial_objects(spatial), ncol, integer(1)))
}

.sirup_spatial_cells <- function(spatial) {
  if (!.sirup_is_spatial_list(spatial)) {
    return(data.frame(.sirup_sample = "spatial", spatial_id = as.character(colnames(spatial)), stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(names(spatial), function(sample) {
    data.frame(
      .sirup_sample = sample,
      spatial_id = as.character(colnames(spatial[[sample]])),
      stringsAsFactors = FALSE
    )
  }))
}

.sirup_validate_sample_name <- function(x, sample) {
  if (!.sirup_is_non_empty_string(sample)) {
    rlang::abort("`sample` must be one non-empty spatial sample name.")
  }
  if (!.sirup_is_spatial_list(x$spatial)) {
    rlang::abort("`sample` can only be used when this SIRUP object contains a spatial list.")
  }
  if (!sample %in% names(x$spatial)) {
    rlang::abort(paste0("Unknown spatial sample: ", sample))
  }
  invisible(sample)
}

.sirup_receptors_for_sample <- function(x, sample) {
  .sirup_validate_sample_name(x, sample)
  sample_col <- .sirup_receptor_fields(x)$sample
  if (is.null(sample_col) || !sample_col %in% colnames(x$receptors)) {
    rlang::abort("This SIRUP object does not have a receptor sample field.")
  }
  x$receptors[!is.na(x$receptors[[sample_col]]) & as.character(x$receptors[[sample_col]]) == sample, , drop = FALSE]
}

.sirup_subset_spatial_component <- function(spatial,
                                            sample = NULL,
                                            spatial_id = NULL,
                                            drop_empty_spatial = TRUE) {
  if (.sirup_is_spatial_list(spatial)) {
    if (!is.null(sample)) {
      spatial <- spatial[sample]
    }
    if (!is.null(spatial_id)) {
      spatial <- lapply(spatial, function(spatial_object) {
        keep_cells <- intersect(colnames(spatial_object), spatial_id)
        if (length(keep_cells) == 0L) {
          return(NULL)
        }
        .sirup_subset_seurat_cells(spatial_object, keep_cells)
      })
      if (isTRUE(drop_empty_spatial)) {
        spatial <- Filter(Negate(is.null), spatial)
      }
    }
    if (length(spatial) == 0L) {
      rlang::abort("No spatial samples/cells remain after subsetting.")
    }
    if (length(spatial) == 1L) {
      return(spatial)
    }
    return(spatial)
  }
  if (!is.null(spatial_id)) {
    keep_cells <- intersect(colnames(spatial), spatial_id)
    if (length(keep_cells) == 0L) {
      rlang::abort("No spatial cells remain after subsetting by `spatial_id`.")
    }
    spatial <- .sirup_subset_seurat_cells(spatial, keep_cells)
  }
  spatial
}

.sirup_subset_seurat_cells <- function(spatial, cells) {
  cells <- intersect(colnames(spatial), as.character(cells))
  if (length(cells) == ncol(spatial)) {
    return(spatial)
  }
  subset(spatial, cells = cells)
}

.sirup_bind_receptor_tables <- function(x, y) {
  x <- tibble::as_tibble(x)
  y <- tibble::as_tibble(y)
  cols <- union(colnames(x), colnames(y))
  add_missing <- function(dat) {
    missing_cols <- setdiff(cols, colnames(dat))
    for (col in missing_cols) {
      dat[[col]] <- NA
    }
    dat[, cols, drop = FALSE]
  }
  dplyr::bind_rows(add_missing(x), add_missing(y))
}

.sirup_sample_summary <- function(x) {
  receptor_fields <- .sirup_receptor_fields(x)
  spatial_cells <- .sirup_spatial_cells(x$spatial)
  receptors <- x$receptors
  spatial_id_col <- receptor_fields$spatial_id
  sample_col <- receptor_fields$sample
  samples <- .sirup_spatial_names(x$spatial)
  if (!.sirup_is_spatial_list(x$spatial)) {
    receptor_values <- as.character(receptors[[spatial_id_col]])
    spatial_ids <- unique(spatial_cells$spatial_id)
    receptor_ids <- unique(receptor_values[!is.na(receptor_values) & nzchar(receptor_values)])
    matched_records <- receptor_values %in% spatial_ids
    matched_ids <- receptor_ids[receptor_ids %in% spatial_ids]
    return(tibble::tibble(
      sample = "spatial",
      n_spatial_ids = nrow(spatial_cells),
      n_receptor_records = nrow(receptors),
      n_receptor_ids = length(receptor_ids),
      n_matched_ids = length(matched_ids),
      n_matched_records = sum(matched_records, na.rm = TRUE),
      id_match_rate = if (length(receptor_ids) == 0L) NA_real_ else length(matched_ids) / length(receptor_ids)
    ))
  }
  dplyr::bind_rows(lapply(samples, function(sample) {
    sample_cells <- spatial_cells[spatial_cells$.sirup_sample == sample, , drop = FALSE]
    sample_receptors <- receptors[
      !is.na(receptors[[sample_col]]) & as.character(receptors[[sample_col]]) == sample,
      ,
      drop = FALSE
    ]
    receptor_values <- as.character(sample_receptors[[spatial_id_col]])
    spatial_ids <- unique(sample_cells$spatial_id)
    receptor_ids <- unique(receptor_values[!is.na(receptor_values) & nzchar(receptor_values)])
    matched_records <- receptor_values %in% spatial_ids
    matched_ids <- receptor_ids[receptor_ids %in% spatial_ids]
    tibble::tibble(
      sample = sample,
      n_spatial_ids = nrow(sample_cells),
      n_receptor_records = nrow(sample_receptors),
      n_receptor_ids = length(receptor_ids),
      n_matched_ids = length(matched_ids),
      n_matched_records = sum(matched_records, na.rm = TRUE),
      id_match_rate = if (length(receptor_ids) == 0L) NA_real_ else length(matched_ids) / length(receptor_ids)
    )
  }))
}

# Count spatial_id overlap between the spatial object and receptor table for
# validation, print, and summary methods.
.sirup_join_diagnostics <- function(x) {
  receptor_fields <- .sirup_receptor_fields(x)
  spatial_cells <- .sirup_spatial_cells(x$spatial)
  receptor_values <- as.character(x$receptors[[receptor_fields$spatial_id]])
  if (.sirup_is_spatial_list(x$spatial)) {
    sample_values <- as.character(x$receptors[[receptor_fields$sample]])
    spatial_key <- paste(spatial_cells$.sirup_sample, spatial_cells$spatial_id, sep = "\r")
    receptor_key <- paste(sample_values, receptor_values, sep = "\r")
    valid_receptor <- !is.na(sample_values) & nzchar(sample_values) &
      !is.na(receptor_values) & nzchar(receptor_values)
    receptor_ids <- unique(receptor_key[valid_receptor])
    spatial_receptor <- !is.na(sample_values) & sample_values %in% names(x$spatial) &
      !is.na(receptor_values) & nzchar(receptor_values)
    matched_records <- spatial_receptor & receptor_key %in% spatial_key
    matched_ids <- unique(receptor_key[matched_records])
  } else {
    receptor_ids <- unique(receptor_values[!is.na(receptor_values) & nzchar(receptor_values)])
    spatial_ids <- unique(spatial_cells$spatial_id)
    matched_records <- receptor_values %in% spatial_ids
    matched_ids <- receptor_ids[receptor_ids %in% spatial_ids]
  }
  list(
    n_spatial_ids = nrow(spatial_cells),
    n_receptor_ids = length(receptor_ids),
    n_matched_ids = length(matched_ids),
    n_matched_records = sum(matched_records, na.rm = TRUE),
    id_match_rate = length(matched_ids) / length(receptor_ids)
  )
}

# Retrieve the canonical receptor-field mapping list from a SIRUP object.
.sirup_receptor_fields <- function(x) {
  if (is.list(x$metadata) && is.list(x$metadata$receptor_fields)) {
    return(x$metadata$receptor_fields)
  }
  rlang::abort("SIRUP receptor field mappings are missing from `metadata$receptor_fields`.")
}

# Return sorted unique values for one mapped receptor field.
.sirup_field_values <- function(x, field) {
  column <- .sirup_receptor_fields(x)[[field]]
  if (is.null(column)) {
    return(character())
  }
  values <- unique(as.character(x$receptors[[column]]))
  sort(values[!is.na(values) & nzchar(values)])
}

# Accept either a Seurat object or SIRUP object where a spatial object is needed.
.sirup_spatial_input <- function(x) {
  if (!isSIRUP(x)) {
    return(x)
  }
  validateSIRUP(x)
  x$spatial
}

# Accept either a receptor table or SIRUP object where receptor data are needed.
.sirup_receptor_input <- function(x, spatial_id_col = NULL) {
  if (!isSIRUP(x)) {
    return(list(
      receptors = x,
      spatial_id_col = spatial_id_col,
      receptor_fields = NULL
    ))
  }
  validateSIRUP(x)
  receptor_fields <- .sirup_receptor_fields(x)
  list(
    receptors = x$receptors,
    spatial_id_col = spatial_id_col %||% receptor_fields$spatial_id,
    receptor_fields = receptor_fields
  )
}

# Resolve paired spatial/receptor inputs for functions that accept either a
# SIRUP object or separate Seurat-plus-receptor arguments.
.sirup_resolve_inputs <- function(object, receptors = NULL, spatial_id_col = NULL) {
  if (!isSIRUP(object)) {
    return(list(
      object = object,
      receptors = receptors,
      spatial_id_col = spatial_id_col,
      receptor_fields = NULL
    ))
  }

  validateSIRUP(object)
  container <- object
  receptor_fields <- .sirup_receptor_fields(container)
  using_container_receptors <- is.null(receptors)
  receptors <- receptors %||% container$receptors
  if (is.null(spatial_id_col)) {
    spatial_id_col <- if (using_container_receptors) {
      receptor_fields$spatial_id
    } else {
      .guess_spatial_id_col(receptors)
    }
  }
  list(
    object = container$spatial,
    receptors = receptors,
    spatial_id_col = spatial_id_col,
    receptor_fields = receptor_fields
  )
}

.sirup_receptor_sample_col <- function(receptors, receptor_fields = NULL, spatial = NULL) {
  mapped <- if (is.list(receptor_fields)) receptor_fields$sample else NULL
  if (!is.null(mapped) && mapped %in% colnames(receptors)) {
    return(mapped)
  }
  if (.sirup_is_spatial_list(spatial) && "sample" %in% colnames(receptors)) {
    return("sample")
  }
  NULL
}

.sirup_add_receptor_spatial_keys <- function(receptors, spatial_id_col, sample_col = NULL) {
  out <- .as_data_frame(receptors)
  out$.sirup_spatial_id <- as.character(out[[spatial_id_col]])
  if (!is.null(sample_col) && sample_col %in% colnames(out)) {
    out$.sirup_sample <- as.character(out[[sample_col]])
  }
  out
}

.sirup_add_coord_spatial_keys <- function(coords) {
  out <- .as_data_frame(coords)
  out$.sirup_spatial_id <- as.character(out$cell)
  if ("sample" %in% colnames(out)) {
    out$.sirup_sample <- as.character(out$sample)
  }
  out
}

.sirup_spatial_key_cols <- function(x, y = NULL) {
  cols <- ".sirup_spatial_id"
  if (".sirup_sample" %in% colnames(x) && (is.null(y) || ".sirup_sample" %in% colnames(y))) {
    cols <- c(".sirup_sample", cols)
  }
  cols
}
