#' Read a Visium HD spatial ID mapping file
#'
#' Reads a spatial ID mapping from CSV/TSV or parquet. Parquet support requires
#' the optional `arrow` package.
#'
#' @param path Path to a spatial ID mapping file.
#' @return A data frame.
#' @export
readSpatialIdMapping <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "parquet") {
    .require_namespace("arrow")
    return(.as_data_frame(arrow::read_parquet(path)))
  }
  if (ext %in% c("tsv", "tab")) {
    return(.as_data_frame(data.table::fread(path, sep = "\t")))
  }
  .as_data_frame(data.table::fread(path))
}

#' Extract spatial IDs from sequence IDs
#'
#' @param airr_df AIRR-style data frame.
#' @param id Column containing read or consensus identifiers.
#' @param spatial_id Name of the output spatial ID column.
#' @param regex Regular expression used to extract the spatial ID. If the regex
#'   contains a capture group, the first capture group is used. Otherwise the
#'   full match is used.
#' @param prefix Optional text added before the extracted match.
#' @param suffix Optional text added after the extracted match.
#' @return Input data with `spatial_id` added.
#' @export
extractSpatialIds <- function(airr_df,
                                    id = "sequence_id",
                                    spatial_id = "spatial_id",
                                    regex,
                                    prefix = "",
                                    suffix = "") {
  .check_required_cols(airr_df, id, "airr_df")
  out <- .as_data_frame(airr_df)
  ids <- as.character(out[[id]])

  # str_match makes a matrix with the full match in the first column and capture groups in subsequent columns
  match_matrix <- stringr::str_match(ids, regex)
  if (ncol(match_matrix) > 1) {
    extracted_spatial_id <- match_matrix[, 2]
  } else {
    extracted_spatial_id <- match_matrix[, 1]
  }
  out[[spatial_id]] <- ifelse(is.na(extracted_spatial_id), NA_character_, paste0(prefix, extracted_spatial_id, suffix))
  out
}

#' Add tissue-position coordinates to airr.tsv style data frames
#'
#' @param airr_df AIRR-style data frame containing a spatial_id column.
#' @param positions Tissue positions data frame or CSV/parquet path.
#' @param spatial_id Column in `airr_df` to match.
#' @param position_spatial_id Spatial ID column in `positions`.
#' @param x Coordinate column to copy or rename to `x`.
#' @param y Coordinate column to copy or rename to `y`.
#' @param keep_position_cols Whether to retain all position columns.
#' @return Data frame with `x` and `y` columns.
#' @export
addTissuePositions <- function(airr_df,
                                 positions,
                                 spatial_id = "spatial_id_8um",
                                 position_spatial_id = "barcode",
                                 x = "pxl_col_in_fullres",
                                 y = "pxl_row_in_fullres",
                                 keep_position_cols = TRUE) {
  .check_required_cols(airr_df, spatial_id, "airr_df")
  if (is.character(positions) && length(positions) == 1) {
    positions <- readSpatialIdMapping(positions)
  }
  .check_required_cols(positions, c(position_spatial_id, x, y), "positions")

  out <- dplyr::left_join(
    .as_data_frame(airr_df),
    .as_data_frame(positions),
    by = stats::setNames(position_spatial_id, spatial_id)
  )
  out$x <- out[[x]]
  out$y <- out[[y]]

  if (!keep_position_cols) {
    keep <- union(colnames(airr_df), c("x", "y"))
    out <- out[, keep, drop = FALSE]
  }
  out
}
