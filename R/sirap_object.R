#' Create a SIRAP object
#'
#' A SIRAP object keeps a spatial Seurat object and an AIRR-style receptor
#' table together with receptor field mappings needed to join and interpret
#' them. Field mappings are stored in `metadata$receptor_fields`.
#' The receptor table is stored as a [tibble::tibble()]. Unmatched receptor
#' rows are retained, which allows a table to contain bulk or other
#' nonspatial receptor records alongside spatial records.
#'
#' @param spatial A Seurat object containing spatial expression data.
#' @param receptors A data-frame-like AIRR-style receptor table.
#' @param barcode_col Column in `receptors` matching spatial cell or spot
#'   barcodes. When `NULL`, a barcode-like column is guessed once and recorded
#'   in the object metadata.
#' @param image_col Optional receptor column identifying the Seurat image.
#' @param sample_col Optional sample column present in both `receptors` and the
#'   Seurat metadata. It is recorded for sample-aware analyses and validation;
#'   Seurat cell or spot names must still be unique within `spatial`.
#' @param clone_col Optional clone-assignment column. The default is recorded
#'   when present and otherwise stored as `NULL`.
#' @param locus_col Optional receptor-locus column. The default is recorded
#'   when present and otherwise stored as `NULL`.
#' @param count_col Optional column containing receptor observation weights or
#'   counts.
#' @param metadata Optional named list of dataset-level metadata.
#' @param validate Whether to validate the completed object.
#'
#' @return A SIRAP object with S3 class `sirap`.
#' @export
createSIRAP <- function(spatial,
                        receptors,
                        barcode_col = NULL,
                        image_col = NULL,
                        sample_col = NULL,
                        clone_col = "clone_id",
                        locus_col = "locus",
                        count_col = NULL,
                        metadata = list(),
                        validate = TRUE) {
  clone_missing <- missing(clone_col)
  locus_missing <- missing(locus_col)
  receptors <- tibble::as_tibble(receptors)
  barcode_col <- .guess_barcode_col(receptors, barcode_col)

  clone_col <- .sirap_resolve_optional_col(
    receptors,
    clone_col,
    "clone_col",
    allow_missing = clone_missing
  )
  locus_col <- .sirap_resolve_optional_col(
    receptors,
    locus_col,
    "locus_col",
    allow_missing = locus_missing
  )
  image_col <- .sirap_resolve_optional_col(receptors, image_col, "image_col")
  sample_col <- .sirap_resolve_optional_col(receptors, sample_col, "sample_col")
  count_col <- .sirap_resolve_optional_col(receptors, count_col, "count_col")

  if (!is.list(metadata)) {
    rlang::abort("`metadata` must be a named list.")
  }
  if (length(metadata) > 0L && (is.null(names(metadata)) || any(!nzchar(names(metadata))))) {
    rlang::abort("`metadata` must be a named list.")
  }
  if ("receptor_fields" %in% names(metadata)) {
    rlang::abort("`metadata$receptor_fields` is reserved; use the field-column arguments to `createSIRAP()`.")
  }
  metadata$receptor_fields <- list(
    barcode = barcode_col,
    image = image_col,
    sample = sample_col,
    clone = clone_col,
    locus = locus_col,
    count = count_col
  )

  out <- structure(
    list(
      spatial = spatial,
      receptors = receptors,
      metadata = metadata
    ),
    class = c("sirap", "list")
  )

  if (isTRUE(validate)) {
    validateSIRAP(out)
  }
  out
}

#' Validate a SIRAP object
#'
#' @param x Object to validate.
#' @param verbose Whether to print a compact validation summary.
#'
#' @return `x`, invisibly, when valid. Invalid objects raise an error.
#' @export
validateSIRAP <- function(x, verbose = FALSE) {
  if (!isSIRAP(x)) {
    rlang::abort("`x` must inherit from class 'sirap'.")
  }
  required_components <- c("spatial", "receptors", "metadata")
  if (!identical(names(x), required_components)) {
    rlang::abort(paste0(
      "A SIRAP object must contain, in order: ",
      paste(required_components, collapse = ", "),
      "."
    ))
  }
  if (!inherits(x$spatial, "Seurat")) {
    rlang::abort("`x$spatial` must be a Seurat object.")
  }
  if (!tibble::is_tibble(x$receptors)) {
    rlang::abort("`x$receptors` must be a tibble.")
  }
  if (!is.list(x$metadata)) {
    rlang::abort("`x$metadata` must be a list.")
  }
  if (length(x$metadata) > 0L && (is.null(names(x$metadata)) || any(!nzchar(names(x$metadata))))) {
    rlang::abort("`x$metadata` must be a named list.")
  }
  receptor_fields <- .sirap_receptor_fields(x)
  field_names <- c("barcode", "image", "sample", "clone", "locus", "count")
  if (!identical(names(receptor_fields), field_names)) {
    rlang::abort(paste0(
      "SIRAP receptor fields must contain, in order: ",
      paste(field_names, collapse = ", "),
      "."
    ))
  }

  barcode_col <- receptor_fields$barcode
  if (!is.character(barcode_col) || length(barcode_col) != 1L || is.na(barcode_col) || !nzchar(barcode_col)) {
    rlang::abort("The SIRAP barcode field must name one receptor column.")
  }
  .check_required_cols(x$receptors, barcode_col, "x$receptors")
  for (field in setdiff(field_names, "barcode")) {
    column <- receptor_fields[[field]]
    if (!is.null(column)) {
      if (!is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column)) {
        rlang::abort(paste0("The SIRAP ", field, " field must be `NULL` or one column name."))
      }
      .check_required_cols(x$receptors, column, "x$receptors")
    }
  }

  if (!is.null(receptor_fields$sample)) {
    spatial_metadata <- x$spatial[[]]
    .check_required_cols(spatial_metadata, receptor_fields$sample, "x$spatial metadata")
  }

  diagnostics <- .sirap_join_diagnostics(x)
  if (diagnostics$n_receptor_barcodes == 0L) {
    rlang::abort("`x$receptors` contains no non-missing spatial barcodes.")
  }
  if (diagnostics$n_matched_barcodes == 0L) {
    rlang::abort("No receptor barcodes match the spatial Seurat object.")
  }
  if (isTRUE(verbose)) {
    print(x)
  }
  invisible(x)
}

#' Test whether an object is a SIRAP object
#'
#' @param x An object.
#'
#' @return A single logical value.
#' @export
isSIRAP <- function(x) {
  inherits(x, "sirap")
}

#' Extract data from a SIRAP object
#'
#' @param x A SIRAP object.
#'
#' @return `spatialObject()` returns the Seurat object. `receptorTable()`
#'   returns the receptor tibble.
#' @name sirap-accessors
NULL

#' @rdname sirap-accessors
#' @export
spatialObject <- function(x) {
  validateSIRAP(x)
  x$spatial
}

#' @rdname sirap-accessors
#' @export
receptorTable <- function(x) {
  validateSIRAP(x)
  x$receptors
}

#' @export
print.sirap <- function(x, ...) {
  validateSIRAP(x)
  info <- .sirap_join_diagnostics(x)
  receptor_fields <- .sirap_receptor_fields(x)
  image_count <- if (requireNamespace("SeuratObject", quietly = TRUE)) {
    length(SeuratObject::Images(x$spatial))
  } else {
    NA_integer_
  }

  cat("<SIRAP object>\n")
  cat("Spatial data:    Seurat, ", format(info$n_spatial_barcodes, big.mark = ","), " ",
      .sirap_count_label(info$n_spatial_barcodes, "barcode"), ", ",
      format(nrow(x$spatial), big.mark = ","), " ",
      .sirap_count_label(nrow(x$spatial), "gene"), sep = "")
  if (!is.na(image_count)) {
    cat(", ", image_count, " ", .sirap_count_label(image_count, "image"), sep = "")
  }
  cat("\n")
  cat("Receptor data:   tibble, ", format(nrow(x$receptors), big.mark = ","), " ",
      .sirap_count_label(nrow(x$receptors), "receptor sequence"), " x ",
      format(ncol(x$receptors), big.mark = ","), " ",
      .sirap_count_label(ncol(x$receptors), "column"), "\n", sep = "")
  cat("Spatial overlap: ", format(info$n_matched_barcodes, big.mark = ","), "/",
      format(info$n_receptor_barcodes, big.mark = ","), " ",
      .sirap_count_label(info$n_receptor_barcodes, "receptor barcode"), " in spatial data (",
      formatC(100 * info$barcode_match_rate, format = "f", digits = 1), "%)\n", sep = "")
  cat("Barcode field:   ", receptor_fields$barcode, "\n", sep = "")
  if (!is.null(receptor_fields$clone)) {
    cat("Clone field:     ", receptor_fields$clone, "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.sirap <- function(object, ...) {
  validateSIRAP(object)
  info <- .sirap_join_diagnostics(object)
  receptor_fields <- .sirap_receptor_fields(object)
  structure(
    c(
      info,
      list(
        n_receptor_records = nrow(object$receptors),
        n_receptor_columns = ncol(object$receptors),
        n_spatial_genes = nrow(object$spatial),
        barcode_col = receptor_fields$barcode,
        clone_col = receptor_fields$clone,
        loci = .sirap_field_values(object, "locus")
      )
    ),
    class = c("summary.sirap", "list")
  )
}

#' @export
print.summary.sirap <- function(x, ...) {
  cat("SIRAP summary\n")
  cat("  Spatial barcodes:  ", format(x$n_spatial_barcodes, big.mark = ","), "\n", sep = "")
  cat("  Receptor records:  ", format(x$n_receptor_records, big.mark = ","), "\n", sep = "")
  cat("  Receptor barcodes: ", format(x$n_receptor_barcodes, big.mark = ","), "\n", sep = "")
  cat("  Matched barcodes:  ", format(x$n_matched_barcodes, big.mark = ","), " (",
      formatC(100 * x$barcode_match_rate, format = "f", digits = 1), "%)\n", sep = "")
  cat("  Matched records:   ", format(x$n_matched_records, big.mark = ","), "\n", sep = "")
  invisible(x)
}

.sirap_resolve_optional_col <- function(receptors, column, argument, allow_missing = FALSE) {
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

.sirap_count_label <- function(n, singular) {
  if (identical(as.numeric(n), 1)) singular else paste0(singular, "s")
}

.sirap_join_diagnostics <- function(x) {
  spatial_barcodes <- unique(as.character(colnames(x$spatial)))
  receptor_fields <- .sirap_receptor_fields(x)
  receptor_values <- as.character(x$receptors[[receptor_fields$barcode]])
  receptor_barcodes <- unique(receptor_values[!is.na(receptor_values) & nzchar(receptor_values)])
  matched <- receptor_barcodes[receptor_barcodes %in% spatial_barcodes]
  list(
    n_spatial_barcodes = length(spatial_barcodes),
    n_receptor_barcodes = length(receptor_barcodes),
    n_matched_barcodes = length(matched),
    n_matched_records = sum(receptor_values %in% spatial_barcodes, na.rm = TRUE),
    barcode_match_rate = length(matched) / length(receptor_barcodes)
  )
}

.sirap_receptor_fields <- function(x) {
  if (is.list(x$metadata) && is.list(x$metadata$receptor_fields)) {
    return(x$metadata$receptor_fields)
  }
  rlang::abort("SIRAP receptor field mappings are missing from `metadata$receptor_fields`.")
}

.sirap_field_values <- function(x, field) {
  column <- .sirap_receptor_fields(x)[[field]]
  if (is.null(column)) {
    return(character())
  }
  values <- unique(as.character(x$receptors[[column]]))
  sort(values[!is.na(values) & nzchar(values)])
}

.sirap_spatial_input <- function(x) {
  if (!isSIRAP(x)) {
    return(x)
  }
  validateSIRAP(x)
  x$spatial
}

.sirap_receptor_input <- function(x, barcode_col = NULL) {
  if (!isSIRAP(x)) {
    return(list(
      receptors = x,
      barcode_col = barcode_col,
      receptor_fields = NULL
    ))
  }
  validateSIRAP(x)
  receptor_fields <- .sirap_receptor_fields(x)
  list(
    receptors = x$receptors,
    barcode_col = barcode_col %||% receptor_fields$barcode,
    receptor_fields = receptor_fields
  )
}

.sirap_resolve_inputs <- function(object, receptors = NULL, barcode_col = NULL) {
  if (!isSIRAP(object)) {
    return(list(
      object = object,
      receptors = receptors,
      barcode_col = barcode_col,
      receptor_fields = NULL
    ))
  }

  validateSIRAP(object)
  container <- object
  receptor_fields <- .sirap_receptor_fields(container)
  using_container_receptors <- is.null(receptors)
  receptors <- receptors %||% container$receptors
  if (is.null(barcode_col)) {
    barcode_col <- if (using_container_receptors) {
      receptor_fields$barcode
    } else {
      .guess_barcode_col(receptors)
    }
  }
  list(
    object = container$spatial,
    receptors = receptors,
    barcode_col = barcode_col,
    receptor_fields = receptor_fields
  )
}
